"""1D meanline performance model for a single-stage overhung centrifugal pump.

The model targets the kind of pump described by API 610 as type OH2: a single-stage,
end-suction, closed impeller hung off the end of a shaft that runs in its own bearing housing.

What it does, at each flow rate:
  1. Velocity triangles at the impeller inlet (eye) and outlet.
  2. Ideal (Euler) head with slip from the Wiesner correlation.
  3. Subtracts hydraulic losses: inlet incidence, impeller channel friction, blade loading, volute mismatch,
     volute friction, and the discharge diffuser cone.
  4. Leakage through the wear rings, solved iteratively, since it depends on the pressure
     the impeller makes.
  5. Disk friction on the impeller shrouds (Daily-Nece) and a mechanical loss fraction.
  6. NPSH required from the classic inlet-velocity correlation.

Everything is SI internally: m, s, kg, Pa, m^3/s, W. Angles are degrees, measured from the
tangential direction, which is the usual pump convention.

This is an engineering estimate, not CFD. Expect roughly 10-20% error on absolute head and
efficiency against a real pump. Trends between design variations are more trustworthy than
the absolute numbers.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace

import numpy as np

G = 9.80665


@dataclass(frozen=True)
class Fluid:
    rho: float = 998.2        # kg/m^3, water at 20 C
    nu: float = 1.004e-6      # m^2/s kinematic viscosity
    p_vapor: float = 2339.0   # Pa


@dataclass(frozen=True)
class Geometry:
    d2: float                 # impeller outlet diameter
    b2: float                 # impeller outlet width
    d1: float                 # eye diameter at the shroud
    dh: float                 # hub diameter at the eye
    beta1b: float             # blade inlet angle at the mean inlet diameter, deg
    beta2b: float             # blade outlet angle, deg
    z: int                    # number of blades
    t: float = 0.004          # blade thickness
    a_throat: float = 0.0     # volute throat area, m^2 (set by the sizing routine)
    d_discharge: float = 0.08 # discharge flange bore
    d_ring: float = 0.0       # wear ring diameter, 0 means d1 + 12 mm
    ring_clearance: float = 0.40e-3  # diametral running clearance
    ring_length: float = 0.015
    n_leak_paths: int = 2     # front wear ring + back wear ring with balance holes
    side_gap_ratio: float = 0.03     # axial gap / r2, for disk friction

    @property
    def d1m(self) -> float:
        return math.sqrt((self.d1 ** 2 + self.dh ** 2) / 2.0)

    @property
    def ring_dia(self) -> float:
        return self.d_ring if self.d_ring > 0 else self.d1 + 0.012


@dataclass(frozen=True)
class LossCoeffs:
    """Loss coefficients. Defaults are typical textbook values, not fitted to a real pump."""
    k_inc_high: float = 0.7   # incidence loss, flow above shock-free
    k_inc_low: float = 0.3    # incidence loss, flow below shock-free (recirculation softens it)
    cf: float = 0.006         # channel skin-friction coefficient (includes roughness)
    k_vol_fric: float = 0.10  # volute friction, fraction of throat velocity head
    k_vol_mis: float = 0.8    # volute mismatch: impeller swirl vs throat velocity
    k_cone: float = 0.25      # discharge cone, fraction of the diffused velocity head
    k_rot: float = 0.5        # fluid rotation factor in the shroud side room
    mech_fraction: float = 0.01   # bearings + mechanical seal, fraction of shaft power
    lam_c: float = 1.1        # NPSHr inlet coefficients
    lam_w: float = 0.25
    ring_friction: float = 0.04
    k_load: float = 0.0       # blade-loading loss (Coppage), x Df^2 u2^2. Off by default: it
                              # overpredicts low-flow losses (3% hump near shutoff on a normal
                              # 6-blade design), so the diffusion factor is reported, not charged.
    hyd_scale: float = 1.0    # multiplies the friction losses; set by calibrate.py


def _slip_wiesner(beta2b_deg: float, z: int, d1m_over_d2: float) -> float:
    s = math.sin(math.radians(beta2b_deg))
    sigma = 1.0 - math.sqrt(s) / z ** 0.7
    eps_lim = math.exp(-8.16 * s / z)
    if d1m_over_d2 > eps_lim:
        sigma *= 1.0 - ((d1m_over_d2 - eps_lim) / (1.0 - eps_lim)) ** 3
    return sigma


def evaluate(geom: Geometry, q: float, n_rpm: float,
             fluid: Fluid = Fluid(), k: LossCoeffs = LossCoeffs()) -> dict:
    """Performance at one operating point. q is the delivered flow in m^3/s."""
    omega = 2 * math.pi * n_rpm / 60.0
    r2 = geom.d2 / 2
    u2 = omega * r2
    u1m = omega * geom.d1m / 2
    u1s = omega * geom.d1 / 2
    u_ring = omega * geom.ring_dia / 2
    b1r = math.radians(geom.beta1b)
    b2r = math.radians(geom.beta2b)

    a1 = math.pi / 4 * (geom.d1 ** 2 - geom.dh ** 2)
    tau1 = 1.0 / (1.0 - geom.z * geom.t / (math.pi * geom.d1m * math.sin(b1r)))
    tau2 = 1.0 / (1.0 - geom.z * geom.t / (math.pi * geom.d2 * math.sin(b2r)))
    sigma = _slip_wiesner(geom.beta2b, geom.z, geom.d1m / geom.d2)

    # impeller channel, for friction
    b1 = (geom.d1 - geom.dh) / 2
    b_avg = (b1 + geom.b2) / 2
    d_avg = (geom.d1m + geom.d2) / 2
    beta_avg = math.radians((geom.beta1b + geom.beta2b) / 2)
    a_avg = max(math.pi * d_avg * math.sin(beta_avg) / geom.z - geom.t, 1e-4)
    dhyd = 2 * a_avg * b_avg / (a_avg + b_avg)
    l_chan = (geom.d2 - geom.d1m) / (2 * math.sin(beta_avg))

    s_rad = geom.ring_clearance / 2
    a_gap = math.pi * geom.ring_dia * s_rad
    mu_gap = 1.0 / math.sqrt(1.5 + k.ring_friction * geom.ring_length / (2 * s_rad))

    q_leak = 0.04 * max(q, 1e-3)
    for _ in range(6):
        qi = q + q_leak
        cm1 = qi / a1 * tau1
        cm2 = qi / (math.pi * geom.d2 * geom.b2)
        cu2 = u2 * (sigma - cm2 * tau2 / (u2 * math.tan(b2r)))
        h_th = u2 * cu2 / G

        dwu = u1m - cm1 / math.tan(b1r)          # >0 below shock-free flow
        k_inc = k.k_inc_low if dwu > 0 else k.k_inc_high
        h_inc = k_inc * dwu ** 2 / (2 * G)

        w1m = math.hypot(cm1, u1m)
        w1s = math.hypot(cm1, u1s)
        w2 = math.hypot(cm2 * tau2, u2 - cu2)
        w_avg = (w1m + w2) / 2
        h_fric = k.hyd_scale * 4 * k.cf * l_chan / dhyd * w_avg ** 2 / (2 * G)

        # blade loading (Coppage diffusion factor, from centrifugal compressor practice):
        # fewer blades carry more load each, so the relative flow decelerates harder
        e = geom.d1 / geom.d2
        df = (1 - w2 / w1s + 0.75 * u2 * cu2 / u2 ** 2
              / ((w1s / w2) * ((geom.z / math.pi) * (1 - e) + 2 * e)))
        h_load = k.hyd_scale * k.k_load * df ** 2 * u2 ** 2 / G

        c2sq = cm2 ** 2 + cu2 ** 2
        h_p = h_th - h_inc - h_fric - h_load - (c2sq - cm1 ** 2) / (2 * G)   # static rise in impeller
        dh_ring = max(h_p - k.k_rot ** 2 * (u2 ** 2 - u_ring ** 2) / (2 * G), 0.0)
        q_leak = geom.n_leak_paths * mu_gap * a_gap * math.sqrt(2 * G * dh_ring)

    # volute + discharge cone (only the delivered flow goes through)
    r3 = 0.5 * geom.d2 * 1.05 + math.sqrt(geom.a_throat / math.pi)
    cu3 = cu2 * r2 / r3
    c_th = q / geom.a_throat
    c_dis = q / (math.pi / 4 * geom.d_discharge ** 2)
    h_vol_mis = k.k_vol_mis * (cu3 - c_th) ** 2 / (2 * G)
    h_vol_fric = k.hyd_scale * k.k_vol_fric * c_th ** 2 / (2 * G)
    h_cone = k.hyd_scale * k.k_cone * max(c_th - c_dis, 0.0) ** 2 / (2 * G)

    head = h_th - h_inc - h_fric - h_load - h_vol_mis - h_vol_fric - h_cone

    # power
    re_disk = omega * r2 ** 2 / fluid.nu
    cm_disk = 0.0255 * geom.side_gap_ratio ** 0.1 / re_disk ** 0.2
    p_disk = cm_disk * fluid.rho * omega ** 3 * r2 ** 5       # both shrouds
    p_euler = fluid.rho * G * qi * h_th
    p_shaft = (p_euler + p_disk) / (1.0 - k.mech_fraction)
    p_hyd = fluid.rho * G * q * max(head, 0.0)
    eta = p_hyd / p_shaft if p_shaft > 0 else 0.0

    npshr = k.lam_c * cm1 ** 2 / (2 * G) + k.lam_w * w1s ** 2 / (2 * G)

    return dict(
        q=q, head=head, eta=eta, eta_h=head / h_th if h_th > 0 else 0.0, p_shaft=p_shaft, npshr=npshr,
        h_th=h_th, h_p=h_p, q_leak=q_leak, q_imp=qi, eta_vol=q / qi, sigma=sigma,
        u2=u2, cu2=cu2, cm2=cm2, cm1=cm1, w1s=w1s, w2=w2, de_haller=w2 / w1s,
        beta1_flow=math.degrees(math.atan2(cm1, u1m)),
        diffusion_factor=df,
        losses=dict(incidence=h_inc, impeller_friction=h_fric, blade_loading=h_load, volute_mismatch=h_vol_mis,
                    volute_friction=h_vol_fric, discharge_cone=h_cone),
        p_disk=p_disk, c_throat=c_th,
    )


def curve(geom: Geometry, q_max: float, n_rpm: float, fluid: Fluid = Fluid(),
          k: LossCoeffs = LossCoeffs(), npts: int = 61) -> dict:
    """Performance curve from shutoff to q_max. Returns arrays."""
    qs = np.linspace(0.0, q_max, npts)
    rows = [evaluate(geom, float(q), n_rpm, fluid, k) for q in qs]
    out = {key: np.array([r[key] for r in rows]) for key in
           ("q", "head", "eta", "p_shaft", "npshr", "h_th", "q_leak")}
    ibep = int(np.argmax(out["eta"]))
    out["q_bep"] = float(qs[ibep])
    out["eta_bep"] = float(out["eta"][ibep])
    out["head_bep"] = float(out["head"][ibep])
    return out


def bep(geom: Geometry, q_guess: float, n_rpm: float, fluid: Fluid = Fluid(),
        k: LossCoeffs = LossCoeffs()) -> tuple[float, float]:
    """Best-efficiency flow and efficiency, by golden-section search around q_guess."""
    lo, hi = 0.3 * q_guess, 2.0 * q_guess
    gr = (math.sqrt(5) - 1) / 2
    f = lambda q: -evaluate(geom, q, n_rpm, fluid, k)["eta"]
    c, d = hi - gr * (hi - lo), lo + gr * (hi - lo)
    for _ in range(40):
        if f(c) < f(d):
            hi = d
        else:
            lo = c
        c, d = hi - gr * (hi - lo), lo + gr * (hi - lo)
    q = (lo + hi) / 2
    return q, evaluate(geom, q, n_rpm, fluid, k)["eta"]


def with_changes(geom: Geometry, **kw) -> Geometry:
    return replace(geom, **kw)
