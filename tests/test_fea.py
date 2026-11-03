"""Verification of the 3D solver against the closed-form rotating disc (clamped bore, free rim)."""
import math

import numpy as np
import pytest

pytest.importorskip("gmsh")
pytest.importorskip("pyamg")

from pump.impeller_fea import CA6NM, disk_clamped_bore_exact, disk_mesh, hoop_stress, solve, von_mises


def test_rotating_disc_hoop_stress_matches_closed_form():
    a, b, thick, om = 0.02, 0.11, 0.005, 2 * math.pi * 3560 / 60
    p, t = disk_mesh(a, b, thick, 0.006)
    c = CA6NM["rho"] * om ** 2
    basis, u, _ = solve(p, t, lambda x, y, z: (c * x, c * y, 0 * z),
                        lambda x, y, z, nx, ny, nz: (0 * x, 0 * x, 0 * x),
                        lambda x, y, z: np.hypot(x, y) < a + 1e-4)
    _, xq, s = von_mises(basis, u)
    hs, r = hoop_stress(s, xq), np.hypot(xq[0], xq[1])
    for rr in (0.05, 0.07, 0.09):
        m = abs(r - rr) < 0.003
        assert hs[m].mean() == pytest.approx(disk_clamped_bore_exact(rr, a, b, om), rel=0.02)
