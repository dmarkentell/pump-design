// JavaScript port of pump/meanline.py, pump/sizing.py and pump/calibrate.py.
// Kept line-for-line close to the Python so the two can be checked against each other
// (tests/test_js_port.py runs this file under Node and compares the numbers).
(function (root) {
  "use strict";
  const G = 9.80665;
  const FLUID = { rho: 998.2, nu: 1.004e-6 };
  const LOSS = {
    k_inc_high: 0.7, k_inc_low: 0.3, cf: 0.006, k_vol_fric: 0.10, k_vol_mis: 0.8, k_cone: 0.25,
    k_rot: 0.5, mech_fraction: 0.01, lam_c: 1.1, lam_w: 0.25, ring_friction: 0.04, k_load: 0.0,
    hyd_scale: 1.0,
  };
  const rad = (d) => (d * Math.PI) / 180;
  const deg = (r) => (r * 180) / Math.PI;

  function geometry(o) {
    const g = Object.assign({
      t: 0.004, a_throat: 0, d_discharge: 0.08, d_ring: 0, ring_clearance: 0.40e-3,
      ring_length: 0.015, n_leak_paths: 2, side_gap_ratio: 0.03,
    }, o);
    g.d1m = Math.sqrt((g.d1 * g.d1 + g.dh * g.dh) / 2);
    g.ring_dia = g.d_ring > 0 ? g.d_ring : g.d1 + 0.012;
    return g;
  }

  function slipWiesner(beta2b, z, ratio) {
    const s = Math.sin(rad(beta2b));
    let sigma = 1 - Math.sqrt(s) / Math.pow(z, 0.7);
    const epsLim = Math.exp(-8.16 * s / z);
    if (ratio > epsLim) sigma *= 1 - Math.pow((ratio - epsLim) / (1 - epsLim), 3);
    return sigma;
  }

  function evaluate(g, q, n, k) {
    k = k || LOSS;
    const omega = 2 * Math.PI * n / 60;
    const r2 = g.d2 / 2, u2 = omega * r2;
    const u1m = omega * g.d1m / 2, u1s = omega * g.d1 / 2, uRing = omega * g.ring_dia / 2;
    const b1r = rad(g.beta1b), b2r = rad(g.beta2b);
    const a1 = Math.PI / 4 * (g.d1 * g.d1 - g.dh * g.dh);
    const tau1 = 1 / (1 - g.z * g.t / (Math.PI * g.d1m * Math.sin(b1r)));
    const tau2 = 1 / (1 - g.z * g.t / (Math.PI * g.d2 * Math.sin(b2r)));
    const sigma = slipWiesner(g.beta2b, g.z, g.d1m / g.d2);

    const b1 = (g.d1 - g.dh) / 2, bAvg = (b1 + g.b2) / 2, dAvg = (g.d1m + g.d2) / 2;
    const betaAvg = rad((g.beta1b + g.beta2b) / 2);
    const aAvg = Math.max(Math.PI * dAvg * Math.sin(betaAvg) / g.z - g.t, 1e-4);
    const dhyd = 2 * aAvg * bAvg / (aAvg + bAvg);
    const lChan = (g.d2 - g.d1m) / (2 * Math.sin(betaAvg));

    const sRad = g.ring_clearance / 2;
    const aGap = Math.PI * g.ring_dia * sRad;
    const muGap = 1 / Math.sqrt(1.5 + k.ring_friction * g.ring_length / (2 * sRad));

    let qLeak = 0.04 * Math.max(q, 1e-3);
    let qi, cm1, cm2, cu2, hTh, hInc, hFric, hLoad, w1s, w2, df;
    for (let i = 0; i < 6; i++) {
      qi = q + qLeak;
      cm1 = qi / a1 * tau1;
      cm2 = qi / (Math.PI * g.d2 * g.b2);
      cu2 = u2 * (sigma - cm2 * tau2 / (u2 * Math.tan(b2r)));
      hTh = u2 * cu2 / G;
      const dwu = u1m - cm1 / Math.tan(b1r);
      const kInc = dwu > 0 ? k.k_inc_low : k.k_inc_high;
      hInc = kInc * dwu * dwu / (2 * G);
      const w1m = Math.hypot(cm1, u1m);
      w1s = Math.hypot(cm1, u1s);
      w2 = Math.hypot(cm2 * tau2, u2 - cu2);
      const wAvg = (w1m + w2) / 2;
      hFric = k.hyd_scale * 4 * k.cf * lChan / dhyd * wAvg * wAvg / (2 * G);
      const e = g.d1 / g.d2;
      df = 1 - w2 / w1s + 0.75 * u2 * cu2 / (u2 * u2) / ((w1s / w2) * ((g.z / Math.PI) * (1 - e) + 2 * e));
      hLoad = k.hyd_scale * k.k_load * df * df * u2 * u2 / G;
      const c2sq = cm2 * cm2 + cu2 * cu2;
      const hP = hTh - hInc - hFric - hLoad - (c2sq - cm1 * cm1) / (2 * G);
      const dhRing = Math.max(hP - k.k_rot * k.k_rot * (u2 * u2 - uRing * uRing) / (2 * G), 0);
      qLeak = g.n_leak_paths * muGap * aGap * Math.sqrt(2 * G * dhRing);
    }
    const r3 = 0.5 * g.d2 * 1.05 + Math.sqrt(g.a_throat / Math.PI);
    const cu3 = cu2 * r2 / r3;
    const cTh = q / g.a_throat;
    const cDis = q / (Math.PI / 4 * g.d_discharge * g.d_discharge);
    const hVolMis = k.k_vol_mis * (cu3 - cTh) ** 2 / (2 * G);
    const hVolFric = k.hyd_scale * k.k_vol_fric * cTh * cTh / (2 * G);
    const hCone = k.hyd_scale * k.k_cone * Math.max(cTh - cDis, 0) ** 2 / (2 * G);
    const head = hTh - hInc - hFric - hLoad - hVolMis - hVolFric - hCone;

    const reDisk = omega * r2 * r2 / FLUID.nu;
    const cmDisk = 0.0255 * Math.pow(g.side_gap_ratio, 0.1) / Math.pow(reDisk, 0.2);
    const pDisk = cmDisk * FLUID.rho * omega ** 3 * r2 ** 5;
    const pEuler = FLUID.rho * G * qi * hTh;
    const pShaft = (pEuler + pDisk) / (1 - k.mech_fraction);
    const pHyd = FLUID.rho * G * q * Math.max(head, 0);
    const eta = pShaft > 0 ? pHyd / pShaft : 0;
    const npshr = k.lam_c * cm1 * cm1 / (2 * G) + k.lam_w * w1s * w1s / (2 * G);
    return {
      q, head, eta, eta_h: hTh > 0 ? head / hTh : 0, p_shaft: pShaft, npshr, h_th: hTh, q_leak: qLeak,
      eta_vol: q / qi, sigma, u2, cu2, cm2, cm1, w1s, w2, de_haller: w2 / w1s,
      beta1_flow: deg(Math.atan2(cm1, u1m)), diffusion_factor: df, p_disk: pDisk, c_throat: cTh,
      losses: { incidence: hInc, impeller_friction: hFric, volute_mismatch: hVolMis,
                volute_friction: hVolFric, discharge_cone: hCone },
    };
  }

  function brentq(f, a, b, xtol) {
    // Brent's method, as in scipy.optimize.brentq
    let fa = f(a), fb = f(b);
    if (fa * fb > 0) throw new Error("root not bracketed");
    let c = a, fc = fa, d = b - a, e = d;
    for (let it = 0; it < 100; it++) {
      if (fb * fc > 0) { c = a; fc = fa; d = b - a; e = d; }
      if (Math.abs(fc) < Math.abs(fb)) { a = b; b = c; c = a; fa = fb; fb = fc; fc = fa; }
      const tol = 2 * 4.4e-16 * Math.abs(b) + 0.5 * xtol;
      const m = 0.5 * (c - b);
      if (Math.abs(m) <= tol || fb === 0) return b;
      if (Math.abs(e) >= tol && Math.abs(fa) > Math.abs(fb)) {
        let p, qq, r;
        const s = fb / fa;
        if (a === c) { p = 2 * m * s; qq = 1 - s; }
        else { qq = fa / fc; r = fb / fc; p = s * (2 * m * qq * (qq - r) - (b - a) * (r - 1)); qq = (qq - 1) * (r - 1) * (s - 1); }
        if (p > 0) qq = -qq; else p = -p;
        if (2 * p < Math.min(3 * m * qq - Math.abs(tol * qq), Math.abs(e * qq))) { e = d; d = p / qq; }
        else { d = m; e = m; }
      } else { d = m; e = m; }
      a = b; fa = fb;
      b += Math.abs(d) > tol ? d : (m > 0 ? tol : -tol);
      fb = f(b);
    }
    return b;
  }

  function bep(g, qGuess, n, k) {
    let lo = 0.3 * qGuess, hi = 2.0 * qGuess;
    const gr = (Math.sqrt(5) - 1) / 2;
    const f = (q) => -evaluate(g, q, n, k).eta;
    let c = hi - gr * (hi - lo), d = lo + gr * (hi - lo);
    for (let i = 0; i < 40; i++) {
      if (f(c) < f(d)) hi = d; else lo = c;
      c = hi - gr * (hi - lo); d = lo + gr * (hi - lo);
    }
    const q = (lo + hi) / 2;
    return [q, evaluate(g, q, n, k).eta];
  }

  const nqOf = (duty) => duty.n * Math.sqrt(duty.q) / Math.pow(duty.head, 0.75);
  const b2OverD2 = (nq) => { const x = nq / 100; return 0.017 + 0.262 * x - 0.08 * x * x + 0.0093 * x ** 3; };
  const psiOf = (nq) => 1.21 * Math.exp(-0.77 * nq / 100);

  function design(duty, ch, k) {
    ch = Object.assign({ z: 6, beta2b: 22.5, b2_factor: 1, d1_factor: 1, incidence: 3, throat_factor: 1 }, ch);
    const nq = nqOf(duty);
    const qImp = duty.q * 1.05;
    const d1 = ch.d1_factor * 4.5 * Math.cbrt(qImp / duty.n);
    const dh = Math.max(0.030, 0.30 * d1);
    const d1m = Math.sqrt((d1 * d1 + dh * dh) / 2);
    const u1m = Math.PI * d1m * duty.n / 60;
    const cm1 = qImp / (Math.PI / 4 * (d1 * d1 - dh * dh));
    let beta1b = deg(Math.atan(cm1 / u1m)) + ch.incidence;
    const u20 = Math.sqrt(2 * G * duty.head / psiOf(nq));
    const d20 = 60 * u20 / (Math.PI * duty.n);
    const b2r = ch.b2_factor * b2OverD2(nq);
    const dDis = Math.sqrt(4 * duty.q / (Math.PI * 7.0));

    function build(d2, b1b) {
      const g = geometry({ d2, b2: b2r * d2, d1, dh, beta1b: b1b, beta2b: ch.beta2b, z: ch.z, a_throat: 1.0, d_discharge: dDis });
      const r = evaluate(g, duty.q, duty.n, k);
      if (r.cu2 <= 0) return g;
      let r3 = 0.5 * d2 * 1.05, aTh = 1;
      for (let i = 0; i < 5; i++) {
        const cu3 = r.cu2 * (d2 / 2) / r3;
        aTh = duty.q / (ch.throat_factor * cu3);
        r3 = 0.5 * d2 * 1.05 + Math.sqrt(aTh / Math.PI);
      }
      g.a_throat = aTh;
      return g;
    }
    function solveD2(b1b) {
      const resid = (d2) => {
        const g = build(d2, b1b);
        if (g.a_throat === 1.0) return -duty.head;
        return evaluate(g, duty.q, duty.n, k).head - duty.head;
      };
      let lo = 0.75 * d20, hi = 1.3 * d20;
      while (resid(hi) < 0) { hi *= 1.15; if (hi > 3 * d20) throw new Error("no impeller diameter reaches the rated head"); }
      while (resid(lo) > 0) lo *= 0.9;
      return brentq(resid, lo, hi, 1e-6);
    }
    let d2;
    for (let i = 0; i < 4; i++) {
      d2 = solveD2(beta1b);
      const flow = evaluate(build(d2, beta1b), duty.q, duty.n, k).beta1_flow;
      const nb = flow + ch.incidence;
      if (Math.abs(nb - beta1b) < 0.01) break;
      beta1b = nb;
    }
    d2 = solveD2(beta1b);
    return build(d2, beta1b);
  }

  function designBepMatched(duty, ch, k, target) {
    target = target || 100;
    const resid = (tf) => {
      const g = design(duty, Object.assign({}, ch, { throat_factor: tf }), k);
      return 100 * duty.q / bep(g, duty.q, duty.n, k)[0] - target;
    };
    const lo = 0.65, hi = 1.1;
    if (resid(lo) * resid(hi) > 0) return { g: design(duty, ch, k), throat_factor: ch.throat_factor || 1 };
    const tf = brentq(resid, lo, hi, 1e-3);
    return { g: design(duty, Object.assign({}, ch, { throat_factor: tf }), k), throat_factor: tf };
  }

  function etaHydraulicCorr(q, nq) {
    const m = 0.08 * Math.pow(1 / q, 0.15) * Math.pow(45 / nq, 0.06);
    return 1 - 0.055 * Math.pow(1 / q, m) - 0.2 * (0.26 - Math.log10(nq / 25)) ** 2 * Math.pow(1 / q, 0.1);
  }
  function etaOverallCorr(q, nq) {
    const m = 0.08 * Math.pow(1 / q, 0.15) * Math.pow(45 / nq, 0.06);
    return 1 - 0.095 * Math.pow(1 / q, m) - 0.3 * (0.35 - Math.log10(nq / 23)) ** 2 * Math.pow(1 / q, 0.05);
  }

  function calibrate(duty) {
    const target = etaHydraulicCorr(duty.q, nqOf(duty));
    const resid = (scale) => {
      const k = Object.assign({}, LOSS, { hyd_scale: scale });
      const g = designBepMatched(duty, {}, k).g;
      const qb = bep(g, duty.q, duty.n, k)[0];
      return evaluate(g, qb, duty.n, k).eta_h - target;
    };
    return Object.assign({}, LOSS, { hyd_scale: brentq(resid, 0.5, 4.0, 1e-4) });
  }

  function summarize(g, duty, k) {
    const r = evaluate(g, duty.q, duty.n, k);
    const [qBep, etaBep] = bep(g, duty.q, duty.n, k);
    const h0 = evaluate(g, 0, duty.n, k).head;
    const h10 = evaluate(g, 0.1 * duty.q, duty.n, k).head;
    const hLow = Math.max(...[0.15, 0.2, 0.3, 0.4].map((f) => evaluate(g, f * duty.q, duty.n, k).head));
    return {
      d2_mm: g.d2 * 1e3, b2_mm: g.b2 * 1e3, d1_mm: g.d1 * 1e3, z: g.z, beta1b: g.beta1b, beta2b: g.beta2b,
      head: r.head, eta: r.eta, p_shaft_kw: r.p_shaft / 1e3, npshr: r.npshr, npsh_margin: duty.npsha - r.npshr,
      q_bep: qBep, eta_bep: etaBep, rated_pct_bep: 100 * duty.q / qBep,
      head_rise_pct: 100 * (h0 - r.head) / r.head, droop_pct: 100 * Math.max(hLow - h10, 0) / h10,
      incidence_rated: g.beta1b - r.beta1_flow, de_haller: r.de_haller, eta_vol: r.eta_vol, u2: r.u2,
    };
  }

  function curve(g, qMax, n, k, npts) {
    npts = npts || 61;
    const out = [];
    for (let i = 0; i < npts; i++) out.push(evaluate(g, (qMax * i) / (npts - 1), n, k));
    return out;
  }

  const api = { G, LOSS, geometry, evaluate, bep, design, designBepMatched, calibrate, summarize, curve,
                nqOf, etaHydraulicCorr, etaOverallCorr, slipWiesner };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.Pump = api;
})(typeof window !== "undefined" ? window : globalThis);
