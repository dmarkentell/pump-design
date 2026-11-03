"""3D stress analysis of the optimized impeller (see pump/impeller_fea.py for the model).

Load cases, both with the rated impeller at running speed:
  rated    rated flow
  shutoff  zero flow: the highest pressures and the heaviest blade loading
plus the rated case again on a coarser mesh, to show which results are mesh-converged.

    python run_impeller_fea.py        # about 15 minutes, writes outputs/fea_*.png and fea.json
"""
from __future__ import annotations

import json
import math
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from pump.impeller_fea import CA6NM, ImpellerSolid, blade_camber, build_mesh, front_profile, solve, von_mises
from pump.meanline import G, evaluate
from run_drive_end import build
from run_study import C, DUTY, OUT


def load_functions(g, k, q, solid: ImpellerSolid):
    """Body force and boundary tractions for one operating point."""
    rho_f, mat = 998.2, CA6NM
    omega = 2 * math.pi * DUTY.n_rpm / 60
    r = evaluate(g, q, DUTY.n_rpm, k=k)
    t = solid.t_shroud
    zf = front_profile(g, solid)
    rs, _ = blade_camber(g)
    r_le, r1s, r2 = rs[0], g.d1 / 2, g.d2 / 2
    r_ring = g.ring_dia / 2
    u1m = omega * r_le
    w1 = math.hypot(r["cm1"], u1m)
    w2 = r["w2"]
    p2 = rho_f * G * r["h_p"]                                  # static rise across the impeller

    def p_pass_raw(rr):
        w = w1 + (w2 - w1) * np.clip((rr - r_le) / (r2 - r_le), 0, 1)
        return rho_f / 2 * ((omega * rr) ** 2 - u1m ** 2 - (w ** 2 - w1 ** 2))

    scale = p2 / p_pass_raw(r2)          # inviscid rothalpy, scaled so it lands on the model's p2

    def p_pass(rr):
        return np.where(rr > r_le, scale * p_pass_raw(rr), 0.0)

    def p_side(rr):
        return p2 - rho_f * (0.5 * omega) ** 2 * (r2 ** 2 - rr ** 2) / 2

    l_blade = (r2 - r_le) / math.sin(math.radians((g.beta1b + g.beta2b) / 2))
    dw = 2 * math.pi * r2 * r["cu2"] / (g.z * l_blade)        # circulation per blade / length
    dp_blade = rho_f * (w1 + w2) / 2 * dw
    tol = 3e-4

    def body(x, y, z):
        c = mat["rho"] * omega ** 2
        return c * x, c * y, 0 * z

    def traction(x, y, z, nx, ny, nz):
        rr = np.hypot(x, y)
        p = np.zeros_like(x)
        back = z < -t + tol
        p = np.where(back & (rr >= r_ring), p_side(rr), p)
        p = np.where(back & (rr < r_ring), 0.1 * (p_side(r_ring) - 0.0), p)   # behind the back ring
        front_out = (z > zf(rr) + t - tol) & (rr > r1s + solid.lip_thick + tol) & (rr < r2 - tol)
        p = np.where(front_out, p_side(rr), p)
        wetted = (z > tol) & (z < zf(rr) - tol) & (rr > r_le - tol) & (rr < r2 - tol)
        p = np.where(wetted, p_pass(rr), p)
        blade_face = wetted & (np.abs(nz) < 0.5) & (rr > r_le + tol)
        e_th = (-y * nx + x * ny) / rr                       # n . e_theta
        p = np.where(blade_face, p + np.sign(e_th) * dp_blade / 2, p)
        return -p * nx, -p * ny, -p * nz

    info = dict(q=q, p2_kpa=p2 / 1e3, dp_blade_kpa=dp_blade / 1e3, omega=omega, w1=w1, w2=w2)
    return body, traction, info


def regions(xq, g, solid):
    """Classify quadrature points: blades, back shroud, front shroud (incl. lip), hub."""
    t = solid.t_shroud
    zf = front_profile(g, solid)
    rr = np.hypot(xq[0], xq[1])
    z = xq[2]
    rs, _ = blade_camber(g)
    hub = rr < g.dh / 2 + 1e-4
    back = (z < 0) & ~hub
    front = (z > zf(rr)) & ~hub & (rr > g.d1 / 2 - 1e-4)
    blade = ~hub & ~back & ~front
    return dict(blades=blade, back_shroud=back, front_shroud=front, hub=hub)


def run_case(name, g, k, q, d_bore, hub_len, solid):
    t0 = time.time()
    p, t = build_mesh(g, d_bore, hub_len, solid)
    body, traction, info = load_functions(g, k, q, solid)
    clamp = lambda x, y, z: np.hypot(x, y) < d_bore / 2 + 1e-4
    basis, u, mesh = solve(p, t, body, traction, clamp)
    vm, xq, s = von_mises(basis, u)
    reg = regions(xq, g, solid)
    out = dict(case=name, nodes=int(p.shape[1]), elements=int(t.shape[1]), dofs=int(basis.N),
               seconds=round(time.time() - t0, 1), **info)
    for key, m in reg.items():
        v = vm[m]
        out[key] = dict(max=float(v.max() / 1e6), p99=float(np.percentile(v, 99) / 1e6),
                        p95=float(np.percentile(v, 95) / 1e6))
    # displacement at the mesh vertices (first 3*Nv dofs of P2 are vertex dofs per component)
    ux, uy, uz = (u[idx] for idx in basis.split_indices())
    xyz = basis.doflocs[:, basis.split_indices()[0]]
    r = np.hypot(xyz[0], xyz[1])
    ur = (ux * xyz[0] + uy * xyz[1]) / np.maximum(r, 1e-9)
    zf = front_profile(g, solid)
    lip = (r > g.d1 / 2 + solid.lip_thick - 2e-4) & (xyz[2] > zf(g.d1 / 2) + solid.t_shroud)
    tip = r > g.d2 / 2 - 2e-4
    out["lip_radial_um"] = dict(mean=float(ur[lip].mean() * 1e6), max=float(ur[lip].max() * 1e6))
    out["tip_radial_um"] = float(ur[tip].mean() * 1e6)
    out["max_disp_um"] = float(np.sqrt(ux ** 2 + uy ** 2 + uz ** 2).max() * 1e6)
    print(f"{name}: {out['elements']} tets, {out['dofs']} dofs, {out['seconds']} s; "
          + ", ".join(f"{kk} p99 {out[kk]['p99']:.0f} max {out[kk]['max']:.0f} MPa" for kk in reg), flush=True)
    return out, (mesh, vm, xq, u, basis)


def fig_stress(mesh, vm, fname, title):
    """Von Mises on the surface, averaged to the vertices, two views."""
    nt = mesh.t.shape[1]
    vm_el = vm.mean(axis=1)
    acc = np.zeros(mesh.p.shape[1]); cnt = np.zeros_like(acc)
    for i in range(4):
        np.add.at(acc, mesh.t[i], vm_el); np.add.at(cnt, mesh.t[i], 1)
    v_node = acc / np.maximum(cnt, 1) / 1e6
    f = mesh.facets[:, mesh.boundary_facets()]
    tri = mesh.p[:, f].transpose(2, 1, 0) * 1e3
    val = v_node[f].mean(axis=0)
    vmax = np.percentile(val, 99.5)
    cmap = plt.get_cmap("turbo")
    fig = plt.figure(figsize=(11, 5.4))
    for j, (elev, azim, sub) in enumerate(((35, -60, "eye side"), (-35, -60, "back side"))):
        ax = fig.add_subplot(1, 2, j + 1, projection="3d")
        pc = Poly3DCollection(tri, facecolors=cmap(np.clip(val / vmax, 0, 1)), edgecolor="none")
        ax.add_collection3d(pc)
        lim = tri[..., :2].max()
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim * 0.6, lim * 0.6)
        ax.set_box_aspect((1, 1, 0.6)); ax.view_init(elev, azim); ax.set_axis_off()
        ax.set_title(sub, fontsize=10)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, vmax))
    fig.colorbar(sm, ax=fig.axes, shrink=0.7, label="von Mises stress (MPa)")
    fig.suptitle(title, fontweight="bold")
    fig.savefig(OUT / fname, bbox_inches="tight")
    plt.close(fig)


def main():
    g, k = build()
    sh = json.load(open(OUT / "drive_end.json"))["shaft"]
    fine, coarse = ImpellerSolid(), ImpellerSolid(h=0.006, h_blade=0.004)
    cases = {}
    rated, data = run_case("rated", g, k, DUTY.q, sh["d_imp"], sh["hub_len"], fine)
    cases["rated"] = rated
    fig_stress(data[0], data[1], "fea_01_rated.png", "Von Mises stress, rated flow, 3560 rpm")
    shut, data = run_case("shutoff", g, k, 0.0, sh["d_imp"], sh["hub_len"], fine)
    cases["shutoff"] = shut
    fig_stress(data[0], data[1], "fea_02_shutoff.png", "Von Mises stress, shutoff, 3560 rpm")
    cases["rated_coarse"], _ = run_case("rated_coarse", g, k, DUTY.q, sh["d_imp"], sh["hub_len"], coarse)
    (OUT / "fea.json").write_text(json.dumps(dict(material=CA6NM, cases=cases), indent=1))


if __name__ == "__main__":
    main()
