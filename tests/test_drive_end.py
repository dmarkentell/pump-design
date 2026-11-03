"""Verification of the drive-end calculations against closed-form solutions."""
import math

import numpy as np
import pytest

from pump.calibrate import calibrate
from pump.drive_end import (Layout, Shaft, angular_pair_equivalent, beam, c_required,
                            de_goodman, design_drive_end, marin_se, reactions, SHAFT_STEEL)
from pump.sizing import Duty, design_bep_matched

DUTY = Duty(q=150 / 3600, head=90.0, n_rpm=3560, npsha=10.0)
LAY = Layout(seal_chamber=0.15, span=0.2, coupling_overhang=0.0, bearing_width=0.0)


def uniform(d=0.05):
    return Shaft(d_imp=d, d_seal=d, d_brg=d, d_span=d, d_cpl=d, hub_len=0.04)


def test_overhung_tip_deflection_matches_beam_theory():
    """Point load at the end of an overhang on rigid supports: v = F a^2 (a + L) / (3 E I)."""
    sh, F = uniform(), 1000.0
    b = beam(sh, LAY, point_loads=[(0.0, F)], k_brg=1e14)
    a, L = b["a"], LAY.span
    exact = F * a ** 2 * (a + L) / (3 * SHAFT_STEEL.e * math.pi * sh.d_imp ** 4 / 64)
    assert abs(b["v"][0]) == pytest.approx(exact, rel=1e-3)


def test_critical_speed_of_heavy_point_mass_on_overhang():
    """A point mass far heavier than the shaft: the first mode tends to sqrt(k_tip / m).

    (A truly massless shaft makes the mass matrix singular, so the test uses a real steel shaft,
    about 6 kg, under 1000 kg; the shaft's own mass then shifts the answer by under 0.2%.)
    """
    sh, m = uniform(), 1000.0
    b = beam(sh, LAY, lumped=[(0.0, m, 0.0)], k_brg=1e14)
    a, L = b["a"], LAY.span
    k_tip = 3 * SHAFT_STEEL.e * math.pi * sh.d_imp ** 4 / 64 / (a ** 2 * (a + L))
    assert b["omega"][0] == pytest.approx(math.sqrt(k_tip / m), rel=2e-3)
    assert b["omega"][0] < math.sqrt(k_tip / m)          # shaft mass can only lower it


def test_reactions_balance():
    F, a, L = 1234.0, 0.17, 0.21
    ra, rb = reactions(F, a, L)
    assert ra + rb == pytest.approx(F)
    assert ra * 0 + rb * L == pytest.approx(-F * a)        # moments about the inboard bearing


def test_de_goodman_hand_calc():
    m, t, d, kf, kfs = 120.0, 150.0, 0.03, 2.14, 3.0
    se = marin_se(SHAFT_STEEL.sut, d)
    lhs = 16 / (math.pi * d ** 3) * (1 / se * math.sqrt(4 * (kf * m) ** 2)
                                     + 1 / SHAFT_STEEL.sut * math.sqrt(3 * (kfs * t) ** 2))
    assert de_goodman(m, t, d, kf, kfs) == pytest.approx(1 / lhs)


def test_bearing_life_inverse():
    p, n, hrs = 2500.0, 3560, 25000
    c = c_required(p, n, hrs)
    assert (c / p) ** 3 * 1e6 / (60 * n) == pytest.approx(hrs)


def test_angular_pair_factors_switch_at_e():
    assert angular_pair_equivalent(1000, 1000) == pytest.approx(1000 + 0.55 * 1000)
    assert angular_pair_equivalent(1000, 3000) == pytest.approx(0.57 * 1000 + 0.93 * 3000)


@pytest.fixture(scope="module")
def drive():
    k = calibrate(DUTY)
    g, _ = design_bep_matched(DUTY, k=k)
    return design_drive_end(g, DUTY, k)


def test_design_passes_every_check(drive):
    failed = [c for c in drive.results["checks"] if not c[1]]
    assert not failed, failed


def test_seal_deflection_is_sized_to_the_limit(drive):
    r = drive.results["seal"]
    assert 30 <= r["deflection_um"] <= 50     # sized close to, but inside, the API limit


def test_radial_load_is_smallest_near_bep(drive):
    lc = drive.results["load_curves"]["rated_imp"]
    i = int(np.argmin(lc["fr"]))
    assert 0.9 <= lc["q"][i] / lc["q_bep"] <= 1.1
