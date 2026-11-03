"""Statistical efficiency correlations for single-stage volute pumps.

These are the curve fits to large sets of measured pumps given by J. F. Gulich,
"Centrifugal Pumps" (Springer), chapter 3, for single-stage, single-entry pumps with
nq between roughly 10 and 100. They give what a well-designed pump of a given size and
specific speed typically achieves at its best efficiency point. Use them as a reality check,
not as a prediction for one particular pump.
"""
import math

Q_REF = 1.0  # m^3/s


def _m(q: float, nq: float, a: float = 1.0) -> float:
    return 0.08 * a * (Q_REF / q) ** 0.15 * (45.0 / nq) ** 0.06


def eta_hydraulic(q: float, nq: float) -> float:
    """Typical hydraulic efficiency at BEP. q in m^3/s (flow per impeller eye)."""
    return (1 - 0.055 * (Q_REF / q) ** _m(q, nq)
            - 0.2 * (0.26 - math.log10(nq / 25.0)) ** 2 * (Q_REF / q) ** 0.1)


def eta_overall(q: float, nq: float) -> float:
    """Typical overall (pump) efficiency at BEP for a single-stage volute pump."""
    return (1 - 0.095 * (Q_REF / q) ** _m(q, nq)
            - 0.3 * (0.35 - math.log10(nq / 23.0)) ** 2 * (Q_REF / q) ** 0.05)
