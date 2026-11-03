"""The phone page runs a JavaScript port of the model. Check it gives the same numbers."""
import json
import shutil
import subprocess
from pathlib import Path

import pytest

from pump.calibrate import calibrate
from pump.sizing import DesignChoices, Duty, design_bep_matched, summarize

ROOT = Path(__file__).resolve().parents[1]
DUTY = Duty(q=150 / 3600, head=90.0, n_rpm=3560, npsha=10.0)
CASES = [dict(z=6, beta2b=22.5), dict(z=5, beta2b=18.6, b2_factor=1.38, d1_factor=0.88, incidence=5.8),
         dict(z=8, beta2b=30.0, b2_factor=0.8, d1_factor=1.1, incidence=1.0)]


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_js_matches_python():
    script = f"""
    const P = require({json.dumps(str(ROOT / 'web' / 'pump.js'))});
    const duty = {{q: {DUTY.q}, head: {DUTY.head}, n: {DUTY.n_rpm}, npsha: {DUTY.npsha}}};
    const k = P.calibrate(duty);
    const out = {json.dumps(CASES)}.map(c => P.summarize(P.designBepMatched(duty, c, k).g, duty, k));
    console.log(JSON.stringify({{scale: k.hyd_scale, out}}));
    """
    js = json.loads(subprocess.run(["node", "-e", script], capture_output=True, text=True, check=True).stdout)
    k = calibrate(DUTY)
    assert js["scale"] == pytest.approx(k.hyd_scale, rel=1e-3)
    for case, got in zip(CASES, js["out"]):
        g, _ = design_bep_matched(DUTY, DesignChoices(**case), k=k)
        want = summarize(g, DUTY, k=k)
        for key in ("d2_mm", "eta", "npshr", "rated_pct_bep", "head_rise_pct", "p_shaft_kw"):
            assert got[key] == pytest.approx(want[key], rel=2e-3, abs=2e-3), key
