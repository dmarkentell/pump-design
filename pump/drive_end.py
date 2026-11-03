"""Drive end of an OH2 pump: loads, shaft, bearings, critical speed, key, coupling and motor.

Layout, measured along the shaft from the impeller centerline (x = 0) toward the motor:

    impeller | seal chamber | inboard bearing ======= span ======= outboard bearings | coupling
       x=0       x_seal          x = a                               x = a + L          x = a + L + c

The impeller overhangs the inboard bearing by a. The inboard bearing is a deep-groove ball bearing
that carries radial load only. The outboard position is a back-to-back pair of 40 degree
angular-contact ball bearings, which carries the axial thrust and its share of the radial load.
That is the usual OH2 arrangement.

Loads come from the hydraulic model in meanline.py:
  * Radial: Stepanoff's single-volute correlation F = K_R rho g H D2 B2, with
    K_R = 0.36 |1 - (Q/Q_BEP)^2|, floored at 0.05 because a real volute never balances perfectly.
  * Axial: closed impeller with front and back wear rings at the same diameter and balance holes
    through the back shroud. The shroud pressures then cancel, leaving the momentum of the inlet
    flow, suction pressure on the shaft cross-section, the pressure drop across the balance holes,
    and an allowance of +/-10% of one shroud's force for unequal fluid rotation in the side gaps.

Everything is SI.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

import numpy as np
from scipy.linalg import eigh

from .meanline import G, Fluid, Geometry, LossCoeffs, bep, evaluate
from .shaftrule import HP_W, motor_hp, motor_torque, shaft_at_impeller
from .sizing import Duty

P_ATM = 101325.0


@dataclass(frozen=True)
class Material:
    name: str
    sut: float
    sy: float
    e: float = 205e9
    rho: float = 7850.0


# Shigley, Mechanical Engineering Design, Table A-20/A-21
SHAFT_STEEL = Material("AISI 4140, Q&T at 650 C (1200 F)", sut=758e6, sy=655e6)
KEY_STEEL = Material("AISI 1018, cold drawn", sut=440e6, sy=370e6)
IMPELLER_DENSITY = 7700.0   # CA6NM 13% chromium cast steel, common for API impellers

# DIN 6885 / ISO 773 parallel keys: (shaft over, up to] -> (width b, height h, shaft depth t1), mm
KEYS = [(17, 22, 6, 6, 3.5), (22, 30, 8, 7, 4.0), (30, 38, 10, 8, 5.0), (38, 44, 12, 8, 5.0),
        (44, 50, 14, 9, 5.5), (50, 58, 16, 10, 6.0), (58, 65, 18, 11, 7.0), (65, 75, 20, 12, 7.5)]


@dataclass(frozen=True)
class Requirements:
    seal_deflection_um: float = 50.0      # API 610 limit at the primary seal faces
    l10_rated_h: float = 25_000.0         # API 610 bearing life, continuous at rated conditions
    l10_max_load_h: float = 16_000.0      # API 610 bearing life at maximum loads
    critical_margin: float = 1.20         # first dry critical >= 1.2 x running speed (design rule)
    fatigue_sf: float = 2.0               # own design target
    key_sf: float = 2.0                   # own design target
    mcsf_frac: float = 0.40               # assumed minimum continuous stable flow, x BEP
    aor_max_frac: float = 1.20            # end of allowable operating region, x BEP
    head_increase: float = 1.05           # API 610: room for +5% head with a larger impeller


@dataclass(frozen=True)
class Layout:
    """Axial layout. Defaults are a plausible arrangement for this size, not a real drawing."""
    seal_from_hub: float = 0.045     # impeller hub face to primary seal faces
    seal_chamber: float = 0.120      # hub face to inboard bearing center (seal chamber + adapter)
    span: float = 0.190              # inboard to outboard bearing centers
    coupling_overhang: float = 0.090
    bearing_width: float = 0.027
    shroud_thickness: float = 0.005
    coupling_mass: float = 4.0       # shaft-side coupling hub, kg
    bearing_stiffness: float = 3.0e8  # N/m per bearing position: preloaded ball bearings, stiff housing
    suction_gauge: float | None = None   # Pa; None = from NPSH available (cold water)


@dataclass
class Shaft:
    d_imp: float
    d_seal: float
    d_brg: float
    d_span: float
    d_cpl: float
    hub_len: float

    def stations(self, lay: Layout) -> list[tuple[float, float, float]]:
        """(x_start, x_end, diameter) segments from the impeller end."""
        a = self.hub_len / 2 + lay.seal_chamber
        x1 = self.hub_len / 2
        x2 = a - lay.bearing_width / 2
        x3 = a + lay.bearing_width / 2
        x4 = a + lay.span - lay.bearing_width
        x5 = a + lay.span + lay.bearing_width / 2
        x6 = a + lay.span + lay.coupling_overhang
        return [(0, x1, self.d_imp), (x1, x2, self.d_seal), (x2, x3, self.d_brg), (x3, x4, self.d_span),
                (x4, x5, self.d_brg), (x5, x6, self.d_cpl)]


def round5(d: float) -> float:
    return math.ceil(d * 1e3 / 5 - 1e-9) * 5 / 1e3


def make_shaft(d_imp: float, d_seal: float, impeller_d2: float) -> Shaft:
    d_brg = round5(d_seal + 0.005)
    return Shaft(d_imp=d_imp, d_seal=d_seal, d_brg=d_brg, d_span=d_brg + 0.010, d_cpl=d_imp,
                 hub_len=max(1.6 * d_imp, 0.22 * impeller_d2))


# ---------------------------------------------------------------------------------- loads
def radial_force(g: Geometry, q: float, q_bep: float, head: float, lay: Layout, fluid: Fluid = Fluid()) -> float:
    k_r = max(0.36 * abs(1 - (q / q_bep) ** 2), 0.05)
    b2_total = g.b2 + 2 * lay.shroud_thickness
    return k_r * fluid.rho * G * head * g.d2 * b2_total


def suction_gauge(duty: Duty, lay: Layout, fluid: Fluid = Fluid()) -> float:
    if lay.suction_gauge is not None:
        return lay.suction_gauge
    return fluid.rho * G * duty.npsha + fluid.p_vapor - P_ATM


def axial_force(g: Geometry, r: dict, n_rpm: float, p_s: float, d_seal: float,
                fluid: Fluid = Fluid(), k: LossCoeffs = LossCoeffs()) -> dict:
    """Axial thrust on the rotor; positive pulls the rotor toward the suction (away from motor)."""
    omega = 2 * math.pi * n_rpm / 60
    r2, r_ring, r_sh = g.d2 / 2, g.ring_dia / 2, d_seal / 2
    p2 = p_s + fluid.rho * G * r["h_p"]                    # static pressure at the impeller outlet
    c = fluid.rho * (k.k_rot * omega) ** 2 / 2
    # one shroud, ring to outlet, gauge to suction: integral of (p(r) - p_s) 2 pi r dr
    shroud = 2 * math.pi * ((p2 - p_s - c * r2 ** 2) * (r2 ** 2 - r_ring ** 2) / 2
                            + c * (r2 ** 4 - r_ring ** 4) / 4)
    p_ring = p2 - c * (r2 ** 2 - r_ring ** 2)
    dp_holes = 0.10 * (p_ring - p_s)                       # balance holes ~5x the ring gap area
    f_holes = dp_holes * math.pi * (r_ring ** 2 - r_sh ** 2)
    f_shaft = -p_s * math.pi * r_sh ** 2
    cm1 = r["q_imp"] / (math.pi / 4 * (g.d1 ** 2 - g.dh ** 2))
    f_mom = -fluid.rho * r["q_imp"] * cm1
    nominal = f_holes + f_shaft + f_mom
    band = 0.10 * shroud
    return dict(nominal=nominal, lo=nominal - band, hi=nominal + band, worst=abs(nominal) + band,
                shroud=shroud, holes=f_holes, shaft=f_shaft, momentum=f_mom)


def impeller_mass(g: Geometry, shaft: Shaft, lay: Layout) -> tuple[float, float]:
    """Mass (kg) and diametral moment of inertia (kg m^2) of the impeller, from its geometry."""
    t = lay.shroud_thickness
    back = math.pi / 4 * (g.d2 ** 2 - g.dh ** 2) * t
    front = math.pi / 4 * (g.d2 ** 2 - g.d1 ** 2) * t
    l_blade = (g.d2 - g.d1m) / (2 * math.sin(math.radians((g.beta1b + g.beta2b) / 2)))
    blades = g.z * g.t * l_blade * (g.b2 + (g.d1 - g.dh) / 2) / 2
    hub = math.pi / 4 * (g.dh ** 2 - shaft.d_imp ** 2) * shaft.hub_len
    m = IMPELLER_DENSITY * (back + front + blades + hub)
    i_d = m * (g.d2 / 2) ** 2 / 8      # thin disc about a diameter is m r^2/4; the mass sits inboard
    return m, i_d


# ---------------------------------------------------------------------------- beam model
def _mesh(shaft: Shaft, lay: Layout, h: float = 0.004):
    segs = shaft.stations(lay)
    xs, ds = [0.0], []
    for x0, x1, d in segs:
        if x1 - x0 < 1e-9:
            continue
        n = max(1, math.ceil((x1 - x0) / h))
        for i in range(1, n + 1):
            xs.append(x0 + (x1 - x0) * i / n)
            ds.append(d)
    return np.array(xs), np.array(ds)


def beam(shaft: Shaft, lay: Layout, mat: Material = SHAFT_STEEL, point_loads=(), lumped=(),
         k_brg: float | None = None):
    """Euler-Bernoulli beam finite elements, two DOF per node (deflection, slope).

    point_loads: (x, F) pairs. lumped: (x, mass, diametral inertia) for the eigenproblem.
    Returns the mesh, static deflection, and the first two natural frequencies with mode shapes.
    """
    xs, ds = _mesh(shaft, lay)
    a = shaft.hub_len / 2 + lay.seal_chamber
    k_brg = lay.bearing_stiffness if k_brg is None else k_brg
    n = len(xs)
    K = np.zeros((2 * n, 2 * n))
    M = np.zeros_like(K)
    for e in range(n - 1):
        le = xs[e + 1] - xs[e]
        I = math.pi * ds[e] ** 4 / 64
        A = math.pi * ds[e] ** 2 / 4
        ke = mat.e * I / le ** 3 * np.array([[12, 6 * le, -12, 6 * le], [6 * le, 4 * le ** 2, -6 * le, 2 * le ** 2],
                                              [-12, -6 * le, 12, -6 * le], [6 * le, 2 * le ** 2, -6 * le, 4 * le ** 2]])
        me = mat.rho * A * le / 420 * np.array([[156, 22 * le, 54, -13 * le], [22 * le, 4 * le ** 2, 13 * le, -3 * le ** 2],
                                                 [54, 13 * le, 156, -22 * le], [-13 * le, -3 * le ** 2, -22 * le, 4 * le ** 2]])
        idx = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3]
        K[np.ix_(idx, idx)] += ke
        M[np.ix_(idx, idx)] += me
    nodes = lambda x: int(np.argmin(abs(xs - x)))
    i_a, i_b = nodes(a), nodes(a + lay.span)
    K[2 * i_a, 2 * i_a] += k_brg
    K[2 * i_b, 2 * i_b] += k_brg
    f = np.zeros(2 * n)
    for x, F in point_loads:
        f[2 * nodes(x)] += F
    for x, m, i_d in lumped:
        j = nodes(x)
        M[2 * j, 2 * j] += m
        M[2 * j + 1, 2 * j + 1] += i_d
    u = np.linalg.solve(K, f)
    w2, vec = eigh(K, M, subset_by_index=[0, 1])
    return dict(x=xs, d=ds, v=u[0::2], slope=u[1::2], omega=np.sqrt(w2), modes=vec[0::2].T,
                i_a=i_a, i_b=i_b, a=a)


def reactions(F: float, a: float, span: float) -> tuple[float, float]:
    """Bearing reactions for a load F at the overhung end (statics). Inboard, outboard."""
    return F * (a + span) / span, -F * a / span


# ---------------------------------------------------------------------------- fatigue
def marin_se(sut: float, d: float, reliability: float = 0.99) -> float:
    ka = 4.51 * (sut / 1e6) ** -0.265                     # machined surface
    dm = d * 1e3
    kb = 1.24 * dm ** -0.107 if dm <= 51 else 1.51 * dm ** -0.157
    ke = {0.5: 1.0, 0.9: 0.897, 0.99: 0.814, 0.999: 0.753}[reliability]
    return ka * kb * ke * 0.5 * sut


def de_goodman(m_alt: float, t_mean: float, d: float, kf: float, kfs: float, mat: Material = SHAFT_STEEL) -> float:
    """Safety factor, rotating bending (fully reversed) with steady torque. Shigley eq. 7-8."""
    se = marin_se(mat.sut, d)
    inv = 16 / (math.pi * d ** 3) * (2 * kf * m_alt / se + math.sqrt(3) * kfs * t_mean / mat.sut)
    return 1 / inv


# ---------------------------------------------------------------------------- bearings
def angular_pair_equivalent(fr: float, fa: float) -> float:
    """Back-to-back pair of 40 deg angular-contact ball bearings (SKF catalogue factors)."""
    if fa / max(fr, 1e-9) <= 1.14:
        return fr + 0.55 * fa
    return 0.57 * fr + 0.93 * fa


def c_required(p: float, n_rpm: float, hours: float) -> float:
    """Basic dynamic load rating needed for an L10 life in hours (ball bearings, exponent 3)."""
    return p * (hours * 60 * n_rpm / 1e6) ** (1 / 3)


# ---------------------------------------------------------------------------- the design
@dataclass
class DriveEnd:
    shaft: Shaft
    layout: Layout
    results: dict = field(default_factory=dict)


def load_cases(g: Geometry, duty: Duty, k: LossCoeffs, lay: Layout, d_seal: float, req: Requirements,
               fluid: Fluid = Fluid()) -> dict:
    """Radial and axial loads across the operating range, rated and maximum impellers."""
    q_bep, _ = bep(g, duty.q, duty.n_rpm, fluid, k)
    g_max = replace(g, d2=g.d2 * math.sqrt(req.head_increase))
    q_bep_max, _ = bep(g_max, duty.q, duty.n_rpm, fluid, k)
    p_s = suction_gauge(duty, lay, fluid)
    out = {"q_bep": q_bep, "q_bep_max": q_bep_max, "p_suction_gauge": p_s, "g_max": g_max}
    for name, geom, qb in (("rated_imp", g, q_bep), ("max_imp", g_max, q_bep_max)):
        fr, fa, qs = [], [], np.linspace(0.0, 1.3, 53) * qb
        for q in qs:
            r = evaluate(geom, float(q), duty.n_rpm, fluid, k)
            fr.append(radial_force(geom, q, qb, r["head"], lay, fluid))
            fa.append(axial_force(geom, r, duty.n_rpm, p_s, d_seal, fluid, k))
        out[name] = dict(q=qs, q_bep=qb, fr=np.array(fr), fa=fa)
    # rated point (rated impeller, rated flow) and the worst point in the allowable range
    r_rated = evaluate(g, duty.q, duty.n_rpm, fluid, k)
    out["rated"] = dict(fr=radial_force(g, duty.q, q_bep, r_rated["head"], lay, fluid),
                        fa=axial_force(g, r_rated, duty.n_rpm, p_s, d_seal, fluid, k))
    worst_fr, worst_fa = 0.0, 0.0
    for q in np.linspace(req.mcsf_frac, req.aor_max_frac, 41) * q_bep_max:
        r = evaluate(g_max, float(q), duty.n_rpm, fluid, k)
        worst_fr = max(worst_fr, radial_force(g_max, q, q_bep_max, r["head"], lay, fluid))
        worst_fa = max(worst_fa, axial_force(g_max, r, duty.n_rpm, p_s, d_seal, fluid, k)["worst"])
    out["worst"] = dict(fr=worst_fr, fa=worst_fa)
    return out


def design_drive_end(g: Geometry, duty: Duty, k: LossCoeffs, lay: Layout = Layout(),
                     req: Requirements = Requirements(), fluid: Fluid = Fluid()) -> DriveEnd:
    omega = 2 * math.pi * duty.n_rpm / 60
    hp = motor_hp(duty.q, duty.head, duty.n_rpm)
    torque = motor_torque(duty.q, duty.head, duty.n_rpm)
    d_imp = shaft_at_impeller(duty.q, duty.head, duty.n_rpm)

    # Step the seal section up until the seal faces meet the deflection limit (with 10% margin)
    # and the inboard bearing shoulder, where the bending moment peaks, meets the fatigue target.
    d_seal = d_imp + 0.010
    for _ in range(40):
        shaft = make_shaft(d_imp, d_seal, g.d2)
        loads = load_cases(g, duty, k, lay, d_seal, req, fluid)
        x_seal = shaft.hub_len / 2 + lay.seal_from_hub
        bm = beam(shaft, lay, point_loads=[(0.0, loads["worst"]["fr"])])
        dfl = abs(np.interp(x_seal, bm["x"], bm["v"]))
        f_tot = loads["worst"]["fr"] + impeller_mass(g, shaft, lay)[0] * G
        x_sh = shaft.hub_len / 2 + lay.seal_chamber - lay.bearing_width / 2
        sf_sh = de_goodman(f_tot * x_sh, torque, d_seal, 2.7, 2.2)
        if dfl <= 0.9 * req.seal_deflection_um * 1e-6 and sf_sh >= req.fatigue_sf:
            break
        d_seal += 0.005
    else:
        raise ValueError("seal deflection limit not reachable; shorten the overhang")

    a = shaft.hub_len / 2 + lay.seal_chamber
    m_imp, id_imp = impeller_mass(g, shaft, lay)
    w_imp = m_imp * G
    x_seal = shaft.hub_len / 2 + lay.seal_from_hub
    fr_worst = loads["worst"]["fr"]

    # deflection: hydraulic load alone (the API "dynamic" case) and with impeller weight added
    bm_h = beam(shaft, lay, point_loads=[(0.0, fr_worst)])
    bm_w = beam(shaft, lay, point_loads=[(0.0, fr_worst + w_imp)])
    seal_um = abs(np.interp(x_seal, bm_h["x"], bm_h["v"])) * 1e6
    seal_um_w = abs(np.interp(x_seal, bm_w["x"], bm_w["v"])) * 1e6

    # critical speeds, dry (no wear-ring Lomakin stiffness, which only raises them)
    lumped = [(0.0, m_imp, id_imp), (a + lay.span + lay.coupling_overhang, lay.coupling_mass, 0.0)]
    md = beam(shaft, lay, lumped=lumped)
    md_rigid = beam(shaft, lay, lumped=lumped, k_brg=1e12)
    crit_rpm = md["omega"][0] * 60 / (2 * math.pi)
    crit_rigid_rpm = md_rigid["omega"][0] * 60 / (2 * math.pi)

    # bearings: rated case and worst case (radial + weight, worst thrust)
    def bearing_loads(fr, fa):
        ra, rb = reactions(fr, a, lay.span)
        p_in = abs(ra)
        p_out = angular_pair_equivalent(abs(rb), fa)
        return abs(ra), abs(rb), p_in, p_out

    rated = bearing_loads(loads["rated"]["fr"] + w_imp, abs(loads["rated"]["fa"]["nominal"]))
    worst = bearing_loads(fr_worst + w_imp, loads["worst"]["fa"])
    c_in = max(c_required(rated[2], duty.n_rpm, req.l10_rated_h), c_required(worst[2], duty.n_rpm, req.l10_max_load_h))
    c_out = max(c_required(rated[3], duty.n_rpm, req.l10_rated_h), c_required(worst[3], duty.n_rpm, req.l10_max_load_h))

    # shaft fatigue at the three critical sections (rotating bending, steady motor torque)
    F = fr_worst + w_imp
    x_key = shaft.hub_len / 2
    x_shoulder = a - lay.bearing_width / 2
    sections = [
        dict(name="Impeller fit, keyseat", x=x_key, d=d_imp, kf=2.14, kfs=3.0, m=F * x_key),
        dict(name="Impeller fit to seal step, 3 mm fillet", x=x_key, d=d_imp, kf=1.7, kfs=1.5, m=F * x_key),
        dict(name="Inboard bearing shoulder, sharp fillet", x=x_shoulder, d=d_seal, kf=2.7, kfs=2.2, m=F * x_shoulder),
    ]
    for s in sections:
        s["sf"] = de_goodman(s["m"], torque, s["d"], s["kf"], s["kfs"])
        s["sigma_a"] = 32 * s["kf"] * s["m"] / (math.pi * s["d"] ** 3)
        s["tau_m"] = 16 * s["kfs"] * torque / (math.pi * s["d"] ** 3)
    # first-cycle yield at the worst section (von Mises on peak stresses)
    for s in sections:
        s["yield_sf"] = SHAFT_STEEL.sy / math.sqrt(s["sigma_a"] ** 2 + 3 * s["tau_m"] ** 2)

    # key at the impeller
    dmm = d_imp * 1e3
    b, h, t1 = next((b, h, t1) for lo, hi, b, h, t1 in KEYS if lo < dmm <= hi)
    key_len = shaft.hub_len - 0.006
    tau_key = 2 * torque / (d_imp * b * 1e-3 * key_len)
    sig_key = 2 * torque / (d_imp * (h - t1) * 1e-3 * key_len)
    key_sf = min(0.577 * KEY_STEEL.sy / tau_key, KEY_STEEL.sy / sig_key)

    # motor: API factor at rated, and end-of-curve power with the maximum impeller
    p_rated = evaluate(g, duty.q, duty.n_rpm, fluid, k)["p_shaft"]
    p_eoc = evaluate(loads["g_max"], req.aor_max_frac * loads["q_bep_max"], duty.n_rpm, fluid, k)["p_shaft"]

    res = dict(
        motor=dict(hp=hp, kw=hp * HP_W / 1e3, torque=torque, p_rated_kw=p_rated / 1e3, p_eoc_kw=p_eoc / 1e3,
                   factor_needed=None),
        loads=dict(rated_fr=loads["rated"]["fr"], rated_fa=loads["rated"]["fa"]["nominal"],
                   worst_fr=fr_worst, worst_fa=loads["worst"]["fa"], impeller_weight=w_imp,
                   impeller_mass=m_imp, p_suction_gauge=loads["p_suction_gauge"],
                   d2_max=loads["g_max"].d2),
        seal=dict(x=x_seal, deflection_um=seal_um, deflection_with_weight_um=seal_um_w,
                  tip_um=abs(bm_h["v"][0]) * 1e6),
        critical=dict(rpm=crit_rpm, rigid_rpm=crit_rigid_rpm, ratio=crit_rpm / duty.n_rpm,
                      second_rpm=md["omega"][1] * 60 / (2 * math.pi)),
        bearings=dict(inboard=dict(type="deep-groove ball", bore_mm=shaft.d_brg * 1e3,
                                   rated_load=rated[0], worst_load=worst[0], c_req=c_in),
                      outboard=dict(type="paired 40 deg angular-contact ball, back to back",
                                    bore_mm=shaft.d_brg * 1e3, rated_radial=rated[1], worst_radial=worst[1],
                                    rated_p=rated[3], worst_p=worst[3], c_req_pair=c_out,
                                    c_req_single=c_out / 1.62)),
        fatigue=sections,
        key=dict(size=f"{b} x {h} mm (DIN 6885)", length_mm=key_len * 1e3, tau=tau_key, sigma=sig_key, sf=key_sf),
        coupling=dict(torque_nameplate=torque, torque_rating_min=1.5 * torque),
        omega=omega,
    )
    res["motor"]["factor_needed"] = api_factor = res["motor"]["kw"] / res["motor"]["p_rated_kw"]
    de = DriveEnd(shaft=shaft, layout=lay, results=res)
    de.results["checks"] = checks(de, duty, req)
    de.results["beam_hyd"] = bm_h
    de.results["beam_modes"] = md
    de.results["load_curves"] = loads
    return de


def checks(de: DriveEnd, duty: Duty, req: Requirements) -> list[tuple[str, bool, str]]:
    r = de.results
    from .shaftrule import api_driver_factor
    need = api_driver_factor(r["motor"]["p_rated_kw"])
    fat = min(s["sf"] for s in r["fatigue"])
    return [
        ("Seal-face deflection <= 50 um (API 610)", r["seal"]["deflection_um"] <= req.seal_deflection_um,
         f"{r['seal']['deflection_um']:.0f} um"),
        ("First dry critical >= 1.2 x running speed", r["critical"]["ratio"] >= req.critical_margin,
         f"{r['critical']['rpm']:.0f} rpm ({r['critical']['ratio']:.2f} x)"),
        ("Shaft fatigue safety factor >= 2 (DE-Goodman)", fat >= req.fatigue_sf, f"{fat:.2f}"),
        ("Key safety factor >= 2", r["key"]["sf"] >= req.key_sf, f"{r['key']['sf']:.1f}"),
        (f"Motor >= {need:.0%} of rated power (API 610)", r["motor"]["factor_needed"] >= need - 1e-9,
         f"{r['motor']['hp']:.0f} hp = {r['motor']['factor_needed']:.0%}"),
        ("Motor covers end of curve, max impeller", r["motor"]["kw"] >= r["motor"]["p_eoc_kw"],
         f"{r['motor']['p_eoc_kw'] * 1e3 / HP_W:.1f} hp needed"),
    ]
