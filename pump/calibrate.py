"""Calibrate the single hydraulic-loss multiplier against the statistical correlation.

The loss coefficients in LossCoeffs are generic textbook values. Rather than trusting them,
we scale the friction and blade-loading losses (impeller channel, blades, volute, discharge cone) by one factor so that the baseline design's hydraulic efficiency
at BEP matches what real pumps of the same size and specific speed achieve. The factor is then
held FIXED for every design variation, so comparisons between designs come from the physics,
not from re-tuning.

Only friction losses are scaled because they set the efficiency at the best efficiency point.
The off-design losses (incidence, volute mismatch) are zero or small at the design point, so the
correlation says nothing about them; they keep their textbook coefficients.
"""
from dataclasses import replace

from scipy.optimize import brentq

from .correlations import eta_hydraulic
from .meanline import LossCoeffs, bep, evaluate
from .sizing import DesignChoices, Duty, design_bep_matched


def calibrate(duty: Duty, ch: DesignChoices = DesignChoices(), k: LossCoeffs = LossCoeffs()) -> LossCoeffs:
    target = eta_hydraulic(duty.q, duty.nq)

    def resid(scale: float) -> float:
        kk = replace(k, hyd_scale=scale)
        g, _ = design_bep_matched(duty, ch, k=kk)
        q_bep, _ = bep(g, duty.q, duty.n_rpm, k=kk)
        return evaluate(g, q_bep, duty.n_rpm, k=kk)["eta_h"] - target

    scale = brentq(resid, 0.5, 4.0, xtol=1e-4)
    return replace(k, hyd_scale=scale)
