"""Size an impeller and volute for a duty point, then trim the impeller to hit the rated head.

Starting values come from standard statistical correlations for single-stage end-suction pumps
(head coefficient and outlet width vs specific speed, as tabulated in Gulich, "Centrifugal
Pumps"). The outlet diameter is then solved with the meanline model so the predicted head at the
rated flow equals the target, and the volute throat is matched to the impeller at that flow.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace

from scipy.optimize import brentq

from .meanline import G, Fluid, Geometry, LossCoeffs, bep, evaluate


@dataclass(frozen=True)
class Duty:
    q: float          # rated flow, m^3/s
    head: float       # rated head, m
    n_rpm: float
    npsha: float      # NPSH available at the suction flange, m

    @property
    def nq(self) -> float:
        """Specific speed, metric (rpm, m^3/s, m)."""
        return self.n_rpm * math.sqrt(self.q) / self.head ** 0.75

    @property
    def ns_us(self) -> float:
        """Specific speed, US units (rpm, gpm, ft)."""
        return self.n_rpm * math.sqrt(self.q * 15850.32) / (self.head * 3.28084) ** 0.75


@dataclass(frozen=True)
class DesignChoices:
    z: int = 6
    beta2b: float = 22.5        # deg
    b2_factor: float = 1.0      # x correlation value
    d1_factor: float = 1.0      # x eye diameter for minimum relative velocity
    incidence: float = 3.0      # deg added to the flow angle at rated flow
    throat_factor: float = 1.0  # throat velocity / impeller swirl at rated flow


def _b2_over_d2(nq: float) -> float:
    x = nq / 100.0
    return 0.017 + 0.262 * x - 0.08 * x ** 2 + 0.0093 * x ** 3


def _psi(nq: float) -> float:
    return 1.21 * math.exp(-0.77 * nq / 100.0)


def design(duty: Duty, ch: DesignChoices = DesignChoices(), fluid: Fluid = Fluid(),
           k: LossCoeffs = LossCoeffs()) -> Geometry:
    nq = duty.nq
    q_imp = duty.q * 1.05                       # allow ~5% leakage when sizing the eye
    d1 = ch.d1_factor * 4.5 * (q_imp / duty.n_rpm) ** (1 / 3)
    dh = max(0.030, 0.30 * d1)                  # hub sized by the overhung shaft
    d1m = math.sqrt((d1 ** 2 + dh ** 2) / 2)
    u1m = math.pi * d1m * duty.n_rpm / 60
    cm1 = q_imp / (math.pi / 4 * (d1 ** 2 - dh ** 2))
    beta1b = math.degrees(math.atan(cm1 / u1m)) + ch.incidence

    u2_0 = math.sqrt(2 * G * duty.head / _psi(nq))
    d2_0 = 60 * u2_0 / (math.pi * duty.n_rpm)
    b2_ratio = ch.b2_factor * _b2_over_d2(nq)
    d_dis = math.sqrt(4 * duty.q / (math.pi * 7.0))   # ~7 m/s at the discharge flange

    def build(d2: float, beta1b: float) -> Geometry:
        g = Geometry(d2=d2, b2=b2_ratio * d2, d1=d1, dh=dh, beta1b=beta1b,
                     beta2b=ch.beta2b, z=ch.z, a_throat=1.0, d_discharge=d_dis)
        r = evaluate(g, duty.q, duty.n_rpm, fluid, k)
        # match the throat to the impeller swirl at rated flow (free-vortex volute)
        r3_est = 0.5 * d2 * 1.05
        if r["cu2"] <= 0:
            return replace(g, a_throat=1.0)
        for _ in range(5):
            cu3 = r["cu2"] * (d2 / 2) / r3_est
            a_th = duty.q / (ch.throat_factor * cu3)
            r3_est = 0.5 * d2 * 1.05 + math.sqrt(a_th / math.pi)
        return replace(g, a_throat=a_th)

    def solve_d2(beta1b: float) -> float:
        def resid(d2: float) -> float:
            g = build(d2, beta1b)
            if g.a_throat == 1.0:
                return -duty.head
            return evaluate(g, duty.q, duty.n_rpm, fluid, k)["head"] - duty.head

        lo, hi = 0.75 * d2_0, 1.3 * d2_0
        while resid(hi) < 0:
            hi *= 1.15
            if hi > 3 * d2_0:
                raise ValueError("no impeller diameter reaches the rated head")
        while resid(lo) > 0:
            lo *= 0.9
        return brentq(resid, lo, hi, xtol=1e-6)

    # The real inlet flow angle includes blade blockage and the wear-ring leakage that recirculates
    # through the eye, so set the blade angle from the computed flow angle and repeat.
    for _ in range(4):
        d2 = solve_d2(beta1b)
        flow = evaluate(build(d2, beta1b), duty.q, duty.n_rpm, fluid, k)["beta1_flow"]
        new = flow + ch.incidence
        if abs(new - beta1b) < 0.01:
            break
        beta1b = new
    d2 = solve_d2(beta1b)
    return build(d2, beta1b)


def summarize(geom: Geometry, duty: Duty, fluid: Fluid = Fluid(), k: LossCoeffs = LossCoeffs()) -> dict:
    """Numbers used to judge a design. API 610 checks here are screening checks only."""
    r = evaluate(geom, duty.q, duty.n_rpm, fluid, k)
    q_bep, eta_bep = bep(geom, duty.q, duty.n_rpm, fluid, k)
    h0 = evaluate(geom, 0.0, duty.n_rpm, fluid, k)["head"]
    # Stability: head must keep rising as flow drops, down to 10% of rated. Below that the 1D
    # model cannot resolve low-flow recirculation, so the last 10% is not judged.
    h10 = evaluate(geom, 0.1 * duty.q, duty.n_rpm, fluid, k)["head"]
    h_low = max(evaluate(geom, f * duty.q, duty.n_rpm, fluid, k)["head"] for f in (0.15, 0.2, 0.3, 0.4))
    return dict(
        d2_mm=geom.d2 * 1e3, b2_mm=geom.b2 * 1e3, d1_mm=geom.d1 * 1e3, z=geom.z,
        beta1b=geom.beta1b, beta2b=geom.beta2b,
        head=r["head"], eta=r["eta"], p_shaft_kw=r["p_shaft"] / 1e3, npshr=r["npshr"],
        npsh_margin=duty.npsha - r["npshr"],
        q_bep=q_bep, eta_bep=eta_bep, rated_pct_bep=100 * duty.q / q_bep,
        head_rise_pct=100 * (h0 - r["head"]) / r["head"],
        droop_pct=100 * max(h_low - h10, 0.0) / h10,
        incidence_rated=geom.beta1b - evaluate(geom, duty.q, duty.n_rpm, fluid, k)["beta1_flow"],
        de_haller=r["de_haller"], eta_vol=r["eta_vol"],
        u2=r["u2"],
    )


def design_bep_matched(duty: Duty, ch: DesignChoices = DesignChoices(), fluid: Fluid = Fluid(),
                       k: LossCoeffs = LossCoeffs(), target_pct: float = 100.0) -> tuple[Geometry, DesignChoices]:
    """Size the volute throat so the rated flow sits at target_pct of the BEP flow.

    A bigger throat moves the best efficiency point to higher flow. Designers normally place
    the rated point at or a little left of BEP; API 610 asks for rated flow within 80-110% of BEP.
    """
    def resid(tf: float) -> float:
        c = replace(ch, throat_factor=tf)
        g = design(duty, c, fluid, k)
        q_bep, _ = bep(g, duty.q, duty.n_rpm, fluid, k)
        return 100 * duty.q / q_bep - target_pct

    lo, hi = 0.65, 1.1   # realistic throat sizing range (throat velocity vs swirl)
    if resid(lo) * resid(hi) > 0:
        return design(duty, ch, fluid, k), ch
    tf = brentq(resid, lo, hi, xtol=1e-3)
    c = replace(ch, throat_factor=tf)
    return design(duty, c, fluid, k), c
