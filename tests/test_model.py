"""Verification checks for the meanline model. Run with: python -m pytest -q"""
import math
from dataclasses import replace

import numpy as np
import pytest

from pump.calibrate import calibrate
from pump.correlations import eta_overall
from pump.meanline import G, LossCoeffs, bep, curve, evaluate, _slip_wiesner
from pump.sizing import Duty, design_bep_matched

DUTY = Duty(q=150 / 3600, head=90.0, n_rpm=3560, npsha=10.0)


@pytest.fixture(scope="module")
def baseline():
    k = calibrate(DUTY)
    g, _ = design_bep_matched(DUTY, k=k)
    return g, k


def test_lossless_head_equals_euler_head(baseline):
    """With every loss switched off, delivered head must equal the Euler head u2*cu2/g."""
    g, _ = baseline
    k0 = LossCoeffs(k_inc_high=0, k_inc_low=0, cf=0, k_load=0, k_vol_fric=0, k_vol_mis=0, k_cone=0)
    g0 = replace(g, n_leak_paths=0)
    r = evaluate(g0, DUTY.q, DUTY.n_rpm, k=k0)
    assert r["head"] == pytest.approx(r["h_th"], rel=1e-12)
    assert r["head"] == pytest.approx(r["u2"] * r["cu2"] / G, rel=1e-12)


def test_wiesner_slip_hand_calc():
    s = math.sin(math.radians(22.5))
    assert _slip_wiesner(22.5, 6, 0.3) == pytest.approx(1 - math.sqrt(s) / 6 ** 0.7)


def test_rated_head_and_bep_position(baseline):
    g, k = baseline
    r = evaluate(g, DUTY.q, DUTY.n_rpm, k=k)
    assert r["head"] == pytest.approx(DUTY.head, abs=1e-3)
    q_bep, _ = bep(g, DUTY.q, DUTY.n_rpm, k=k)
    assert 80 <= 100 * DUTY.q / q_bep <= 110


def test_efficiency_matches_statistics(baseline):
    """Calibrated on HYDRAULIC efficiency; OVERALL efficiency is an independent check."""
    g, k = baseline
    _, eta = bep(g, DUTY.q, DUTY.n_rpm, k=k)
    assert eta == pytest.approx(eta_overall(DUTY.q, DUTY.nq), abs=0.02)


def test_affinity_laws(baseline):
    """At 80% speed: Q_bep scales with n, head with n^2 (small Reynolds effects allowed)."""
    g, k = baseline
    q1, _ = bep(g, DUTY.q, DUTY.n_rpm, k=k)
    q2, _ = bep(g, DUTY.q, 0.8 * DUTY.n_rpm, k=k)
    h1 = evaluate(g, q1, DUTY.n_rpm, k=k)["head"]
    h2 = evaluate(g, q2, 0.8 * DUTY.n_rpm, k=k)["head"]
    assert q2 / q1 == pytest.approx(0.8, rel=0.03)
    assert h2 / h1 == pytest.approx(0.64, rel=0.03)


def test_head_curve_rises_to_shutoff(baseline):
    """Falling from about 10% flow onward, and any droop below 10% flow stays under 0.5%.

    The 1D model does not capture low-flow recirculation, so the last few percent of flow
    before shutoff are not trustworthy to better than that anyway.
    """
    g, k = baseline
    c = curve(g, 1.4 * DUTY.q, DUTY.n_rpm, k=k)
    i10 = np.searchsorted(c["q"], 0.1 * DUTY.q)
    tol = 5e-4 * c["head"][i10]   # same 0.05% noise tolerance the study uses
    assert np.all(np.diff(c["head"][i10:]) < tol)
    assert (c["head"].max() - c["head"][0]) / c["head"][0] < 5e-3


def test_shaft_power_exceeds_hydraulic_power(baseline):
    g, k = baseline
    for f in (0.1, 0.5, 1.0, 1.3):
        r = evaluate(g, f * DUTY.q, DUTY.n_rpm, k=k)
        assert r["p_shaft"] > 998.2 * G * r["q"] * r["head"]
        assert 0 < r["eta"] < 1


def test_leakage_iteration_converged(baseline):
    """Running the fixed-point loop longer must not change the answer."""
    g, k = baseline
    r = evaluate(g, DUTY.q, DUTY.n_rpm, k=k)
    assert 0.90 < r["eta_vol"] < 0.99


def test_inlet_incidence_is_what_was_asked(baseline):
    """The sizing loop must land the requested +3 deg incidence at rated flow, leakage included."""
    g, k = baseline
    r = evaluate(g, DUTY.q, DUTY.n_rpm, k=k)
    assert g.beta1b - r["beta1_flow"] == pytest.approx(3.0, abs=0.05)
