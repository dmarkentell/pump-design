"""3D finite-element stress analysis of the closed impeller.

Geometry (built with gmsh/OpenCASCADE from the hydraulic design):
  * back shroud: flat disc, hub OD to D2
  * front shroud: cone from the eye to the outlet, following the meridional passage, with the
    front wear-ring lip at the eye
  * hub: the bore fits the shaft; the nose reaches into the eye
  * blades: the same camber lines as the plan view (blade angle linear in radius), constant
    thickness, running from the leading edge at the mean eye diameter to D2, fused to both shrouds

Loads at running speed, all gauge pressures relative to suction:
  1. centrifugal body force rho w^2 r
  2. passage pressure on every wetted surface inside the impeller, from the rothalpy balance
     p - p_s = rho/2 [(u^2 - u1^2) - (w^2 - w1^2)], with w varying linearly from inlet to outlet
  3. side-gap pressure on the outside of both shrouds, p2 minus the drop from fluid rotating at
     half the impeller speed, the same model the leakage and thrust use
  4. blade loading: the circulation each blade carries, 2 pi r2 cu2 / z, spread over its length,
     gives a mean pressure difference rho w dw across the blade, added to the pressure side and
     taken off the suction side

Boundary condition: the bore is clamped to the shaft. Solved with quadratic tetrahedra (P2) in
scikit-fem and a conjugate-gradient solver with an algebraic-multigrid preconditioner (pyamg).

Units are SI throughout.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .meanline import G, Geometry

# ASTM A743 CA6NM (13Cr-4Ni cast stainless), the usual API impeller material: minimum properties
CA6NM = dict(name="ASTM A743 CA6NM", e=200e9, nu=0.29, rho=7700.0, sy=550e6, sut=755e6)


@dataclass(frozen=True)
class ImpellerSolid:
    t_shroud: float = 0.005      # shroud thickness
    lip_thick: float = 0.006     # wear-ring lip, radial
    lip_len: float = 0.015       # wear-ring lip, axial
    h: float = 0.004             # target element size
    h_blade: float = 0.0025      # element size on the blades


def blade_camber(g: Geometry, n: int = 40, extend: float = 0.0):
    """(r, theta) along the blade, blade angle linear in radius: d(theta) = dr / (r tan(beta)).

    extend continues the blade past D2 at the outlet angle (it is trimmed back to D2 later).
    """
    r1, r2 = g.d1m / 2, g.d2 / 2
    rs = np.linspace(r1, r2 + extend, n)
    beta = np.radians(g.beta1b + (g.beta2b - g.beta1b) * np.clip((rs - r1) / (r2 - r1), 0, 1))
    f = 1 / (rs * np.tan(beta))
    th = np.concatenate([[0], np.cumsum((f[1:] + f[:-1]) / 2 * np.diff(rs))])
    return rs, -th            # negative: blades sweep back against +z (counter-clockwise) rotation


def front_profile(g: Geometry, s: ImpellerSolid):
    """Axial position of the front shroud's inner surface as a function of radius."""
    r1s, r2 = g.d1 / 2, g.d2 / 2
    b_in = max(g.b2, 0.6 * (g.d1 - g.dh) / 2)
    return lambda r: g.b2 + (b_in - g.b2) * (r2 - r) / (r2 - r1s)


def build_mesh(g: Geometry, d_bore: float, hub_len: float, s: ImpellerSolid = ImpellerSolid(),
               msh_path: str | None = None):
    """Build the solid with gmsh and return node coordinates (3 x N) and tetrahedra (4 x M)."""
    import gmsh
    zf = front_profile(g, s)
    t = s.t_shroud
    r_bore, r_hub, r1s, r2 = d_bore / 2, g.dh / 2, g.d1 / 2, g.d2 / 2
    rx = r2 + 0.003                      # oversize outlet, turned to D2 at the end
    rs, th = blade_camber(g, extend=rx - r2)
    r_le = rs[0]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("impeller")
    occ = gmsh.model.occ

    def revolve_polygon(pts):
        tags = [occ.addPoint(x, 0, z) for x, z in pts]
        lines = [occ.addLine(tags[i], tags[(i + 1) % len(tags)]) for i in range(len(tags))]
        surf = occ.addPlaneSurface([occ.addCurveLoop(lines)])
        out = occ.revolve([(2, surf)], 0, 0, 0, 0, 0, 1, 2 * math.pi)
        return [d for d in out if d[0] == 3]

    back = revolve_polygon([(r_bore, -t), (rx, -t), (rx, 0), (r_bore, 0)])
    hub = revolve_polygon([(r_bore, -t), (r_hub, -t), (r_hub, hub_len - t), (r_bore, hub_len - t)])
    front = revolve_polygon([(r1s, zf(r1s)), (rx, zf(rx)), (rx, zf(rx) + t), (r1s + s.lip_thick, zf(r1s) + t),
                             (r1s + s.lip_thick, zf(r1s) + s.lip_len), (r1s, zf(r1s) + s.lip_len)])
    # passage envelope used to trim the blades; it overlaps both shrouds by half their thickness
    env = revolve_polygon([(r_le, -t / 2), (rx, -t / 2), (rx, zf(rx) + t / 2), (r_le, zf(r_le) + t / 2)])

    # blade: camber line thickened in plan view, extruded through the passage, trimmed to it
    x, y = rs * np.cos(th), rs * np.sin(th)
    dx, dy = np.gradient(x), np.gradient(y)
    nrm = np.hypot(dx, dy)
    nx, ny = -dy / nrm, dx / nrm
    half = g.t / 2
    side_a = [occ.addPoint(x[i] + half * nx[i], y[i] + half * ny[i], -t) for i in range(len(rs))]
    side_b = [occ.addPoint(x[i] - half * nx[i], y[i] - half * ny[i], -t) for i in range(len(rs))]
    ca = occ.addSpline(side_a)
    cb = occ.addSpline(side_b)
    c_out = occ.addLine(side_a[-1], side_b[-1])
    c_in = occ.addLine(side_b[0], side_a[0])
    sb = occ.addPlaneSurface([occ.addCurveLoop([ca, c_out, -cb, c_in])])
    ext = occ.extrude([(2, sb)], 0, 0, zf(r_le) + 2 * t)
    blade_raw = [d for d in ext if d[0] == 3]
    env_copy = occ.copy(env)
    blade, _ = occ.intersect(blade_raw, env_copy)
    blades = [blade[0]]
    for i in range(1, g.z):
        cp = occ.copy([blade[0]])
        occ.rotate(cp, 0, 0, 0, 0, 0, 1, 2 * math.pi * i / g.z)
        blades.append(cp[0])
    occ.remove(env, recursive=True)
    solid, _ = occ.fuse(back, hub + front + blades)
    # everything was built oversize at the outlet; turn the OD to D2 in one cut, like a lathe,
    # so the shrouds and blade tips share a single outer face
    lathe = occ.addCylinder(0, 0, -2 * t, 0, 0, hub_len + zf(r1s) + 4 * t, r2)
    solid, _ = occ.intersect(solid, [(3, lathe)])
    occ.synchronize()
    gmsh.model.addPhysicalGroup(3, [d[1] for d in solid], 1)
    n_volumes = len(solid)

    # finer elements on the blades, coarser elsewhere
    gmsh.option.setNumber("Mesh.MeshSizeMax", s.h)
    gmsh.option.setNumber("Mesh.MeshSizeMin", s.h_blade * 0.8)
    fld = gmsh.model.mesh.field
    cyl = fld.add("Cylinder")             # coarse inside the eye, fine where the blades are
    fld.setNumber(cyl, "Radius", r_le - 0.002)
    fld.setNumber(cyl, "ZCenter", 0.0)
    fld.setNumber(cyl, "ZAxis", 1.0)
    fld.setNumber(cyl, "VIn", s.h)
    fld.setNumber(cyl, "VOut", s.h_blade)
    fld.setAsBackgroundMesh(cyl)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)           # Delaunay (HXT trips on boolean slivers)
    gmsh.model.mesh.generate(3)
    if msh_path:
        gmsh.write(msh_path)
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    etypes, _, enodes = gmsh.model.mesh.getElements(3)
    gmsh.finalize()

    tet = enodes[list(etypes).index(4)].reshape(-1, 4)
    p = coords.reshape(-1, 3)
    index = np.zeros(int(node_tags.max()) + 1, dtype=np.int64)
    index[node_tags.astype(np.int64)] = np.arange(len(node_tags))
    t4 = index[tet.astype(np.int64)]
    used = np.unique(t4)
    remap = -np.ones(len(p), dtype=np.int64)
    remap[used] = np.arange(len(used))
    return p[used].T.copy(), remap[t4].T.copy()


def _fix_orientation(p, t):
    a, b, c, d = (p[:, t[i]] for i in range(4))
    vol = np.einsum("ij,ij->j", np.cross(b - a, c - a, axis=0), d - a)
    bad = vol < 0
    t2 = t.copy()
    t2[1, bad], t2[2, bad] = t[2, bad], t[1, bad]
    return t2


def solve(p, t, body_force, traction, dirichlet, mat=CA6NM, tol=1e-8):
    """Linear elasticity with P2 tetrahedra.

    body_force(x, y, z) -> (fx, fy, fz) arrays; traction(x, y, z, nx, ny, nz) -> (tx, ty, tz) on
    the boundary; dirichlet(x, y, z) -> bool mask of clamped boundary facets.
    Returns (basis, u, mesh).
    """
    import pyamg
    from scipy.sparse.linalg import cg
    from skfem import (Basis, BilinearForm, ElementTetP2, ElementVector, FacetBasis, LinearForm,
                       MeshTet, asm, condense)
    from skfem.helpers import ddot, sym_grad, trace

    mesh = MeshTet(p, _fix_orientation(p, t))
    e = ElementVector(ElementTetP2())
    basis = Basis(mesh, e, intorder=2)
    lam = mat["e"] * mat["nu"] / ((1 + mat["nu"]) * (1 - 2 * mat["nu"]))
    mu = mat["e"] / (2 * (1 + mat["nu"]))

    @BilinearForm
    def stiffness(u, v, w):
        eu, ev = sym_grad(u), sym_grad(v)
        return 2 * mu * ddot(eu, ev) + lam * trace(eu) * trace(ev)

    @LinearForm
    def body(v, w):
        fx, fy, fz = body_force(*w.x)
        return fx * v[0] + fy * v[1] + fz * v[2]

    K = asm(stiffness, basis)
    f = asm(body, basis)

    facets = mesh.boundary_facets()
    fb = FacetBasis(mesh, e, facets=facets, intorder=3)

    @LinearForm
    def surf(v, w):
        tx, ty, tz = traction(*w.x, *w.n)
        return tx * v[0] + ty * v[1] + tz * v[2]

    f = f + asm(surf, fb)

    mid = mesh.p[:, mesh.facets[:, facets]].mean(axis=1)
    clamp = facets[dirichlet(*mid)]
    D = basis.get_dofs(clamp).all()
    Kc, fc, x0, I = condense(K, f, D=D)
    # rigid-body modes (3 translations, 3 rotations) as near-null space for the multigrid
    xyz = basis.doflocs
    B = np.zeros((basis.N, 6))
    for comp, idx in enumerate(basis.split_indices()):
        B[idx, comp] = 1.0
    ix, iy, iz = basis.split_indices()
    B[ix, 3], B[iy, 3] = -xyz[1, ix], xyz[0, iy]
    B[iy, 4], B[iz, 4] = -xyz[2, iy], xyz[1, iz]
    B[ix, 5], B[iz, 5] = xyz[2, ix], -xyz[0, iz]
    ml = pyamg.smoothed_aggregation_solver(Kc.tocsr(), B=B[I], max_coarse=2000,
                                           strength=("symmetric", {"theta": 0.0}))
    M = ml.aspreconditioner(cycle="V")
    sol, info = cg(Kc, fc, rtol=tol, maxiter=2000, M=M)
    if info != 0:
        raise RuntimeError(f"CG did not converge (info={info})")
    u = np.zeros(basis.N)
    u[I] = sol
    return basis, u, mesh


def von_mises(basis, u, mat=CA6NM):
    """Von Mises stress at each element's quadrature points: (elements x points), plus the points."""
    from skfem.helpers import sym_grad, trace, eye
    lam = mat["e"] * mat["nu"] / ((1 + mat["nu"]) * (1 - 2 * mat["nu"]))
    mu = mat["e"] / (2 * (1 + mat["nu"]))
    w = basis.interpolate(u)
    eps = sym_grad(w)
    s = 2 * mu * eps + lam * trace(eps) * eye(np.ones_like(eps[0, 0]), 3)
    sxx, syy, szz = s[0, 0], s[1, 1], s[2, 2]
    sxy, syz, sxz = s[0, 1], s[1, 2], s[0, 2]
    vm = np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
                 + 3 * (sxy ** 2 + syz ** 2 + sxz ** 2))
    xq = basis.mapping.F(basis.X)
    return vm, xq, s


def hoop_stress(s, xq):
    """Tangential (hoop) stress at the quadrature points."""
    x, y = xq[0], xq[1]
    r = np.hypot(x, y)
    c, sn = x / r, y / r
    return c * c * s[1, 1] + sn * sn * s[0, 0] - 2 * c * sn * s[0, 1]


def disk_mesh(a: float, b: float, thick: float, h: float):
    """Annular disc a < r < b, for verification against the closed-form rotating-disc solution."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    occ = gmsh.model.occ
    outer = occ.addCylinder(0, 0, 0, 0, 0, thick, b)
    inner = occ.addCylinder(0, 0, 0, 0, 0, thick, a)
    occ.cut([(3, outer)], [(3, inner)])
    occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.MeshSizeMin", h / 2)
    gmsh.model.mesh.generate(3)
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    etypes, _, enodes = gmsh.model.mesh.getElements(3)
    gmsh.finalize()
    tet = enodes[list(etypes).index(4)].reshape(-1, 4)
    index = np.zeros(int(node_tags.max()) + 1, dtype=np.int64)
    index[node_tags.astype(np.int64)] = np.arange(len(node_tags))
    t4 = index[tet.astype(np.int64)]
    p = coords.reshape(-1, 3)
    used = np.unique(t4)
    remap = -np.ones(len(p), dtype=np.int64)
    remap[used] = np.arange(len(used))
    return p[used].T.copy(), remap[t4].T.copy()


def disk_clamped_bore_exact(r, a, b, omega, mat=CA6NM):
    """Plane-stress rotating disc, u(a) = 0, sigma_r(b) = 0. Returns hoop stress at r."""
    E, nu, rho = mat["e"], mat["nu"], mat["rho"]
    k = (1 - nu ** 2) * rho * omega ** 2 / (8 * E)
    # u = C1 r + C2 / r - k r^3 ;  sigma_r = E/(1-nu^2) [(1+nu) C1 - (1-nu) C2/r^2 - (3+nu) k r^2]
    A = np.array([[a, 1 / a], [(1 + nu), -(1 - nu) / b ** 2]])
    rhs = np.array([k * a ** 3, (3 + nu) * k * b ** 2])
    c1, c2 = np.linalg.solve(A, rhs)
    return E / (1 - nu ** 2) * ((1 + nu) * c2 / r ** 2 * 0 + (c1 + c2 / r ** 2 - k * r ** 2)
                                + nu * (c1 - c2 / r ** 2 - 3 * k * r ** 2))
