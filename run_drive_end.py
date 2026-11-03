"""Drive-end design for the optimized impeller from run_study.py.

Reads outputs/results.json (run run_study.py first), rebuilds the optimized impeller, and sizes
the shaft, bearings, key, coupling and motor. Writes figures and outputs/drive_end.json.

    python run_drive_end.py
"""
from __future__ import annotations

import json
import math
from dataclasses import asdict, replace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pump.calibrate import calibrate
from pump.drive_end import Layout, Requirements, beam, design_drive_end
from pump.meanline import G
from pump.sizing import DesignChoices, design_bep_matched
from run_study import C, DUTY, FT, GPM, HP, OUT

LBF = 0.224809


def build():
    res = json.loads((OUT / "results.json").read_text())
    k = calibrate(DUTY)
    ch = DesignChoices(**{key: v for key, v in res["optimized"]["choices"].items() if key != "throat_factor"})
    g, _ = design_bep_matched(DUTY, ch, k=k)
    return g, k


def fig_rotor(de, fname):
    r, sh, lay = de.results, de.shaft, de.layout
    bm = r["beam_hyd"]
    a = bm["a"]
    F = r["loads"]["worst_fr"] + r["loads"]["impeller_weight"]
    x = bm["x"]
    moment = np.where(x <= a, F * x, F * a * np.clip((a + lay.span - x) / lay.span, 0, None))

    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(10, 8.6), sharex=True,
                                        gridspec_kw=dict(height_ratios=[1.6, 1, 1]))
    mm = 1e3
    for x0, x1, d in sh.stations(lay):
        ax0.add_patch(plt.Rectangle((x0 * mm, -d / 2 * mm), (x1 - x0) * mm, d * mm, facecolor="#cbd5e1",
                                    edgecolor=C["ink"], lw=1))
        if x1 - x0 > 0.012:
            ax0.text((x0 + x1) / 2 * mm, 0, f"Ø{d * mm:.0f}", ha="center", va="center", fontsize=8)
    g_d2 = r["loads"]["d2_max"] / math.sqrt(1.05)
    ax0.add_patch(plt.Rectangle((-0.5 * sh.hub_len * mm, -g_d2 / 2 * mm), sh.hub_len * mm, g_d2 * mm,
                                facecolor=C["base"], alpha=0.25, edgecolor=C["base"]))
    ax0.text(0, g_d2 / 2 * mm + 6, "impeller", ha="center", color=C["base"], fontsize=8.5)
    xs = r["seal"]["x"] * mm
    ax0.plot([xs, xs], [-sh.d_seal * mm, sh.d_seal * mm], color=C["opt"], lw=2)
    ax0.text(xs, sh.d_seal * mm + 4, "seal faces", ha="center", color=C["opt"], fontsize=8.5)
    for xb, label in ((a, "inboard\nball bearing"), (a + lay.span, "outboard\nangular-contact pair")):
        ax0.add_patch(plt.Rectangle(((xb - lay.bearing_width / 2) * mm, sh.d_brg / 2 * mm), lay.bearing_width * mm,
                                    18, facecolor=C["amber"], alpha=0.7))
        ax0.add_patch(plt.Rectangle(((xb - lay.bearing_width / 2) * mm, -sh.d_brg / 2 * mm - 18),
                                    lay.bearing_width * mm, 18, facecolor=C["amber"], alpha=0.7))
        ax0.text(xb * mm, sh.d_brg / 2 * mm + 22, label, ha="center", fontsize=8)
    xc = a + lay.span + lay.coupling_overhang
    ax0.text(xc * mm, sh.d_cpl * mm / 2 + 8, "coupling", ha="right", fontsize=8.5)
    ax0.annotate("", xy=(0, -g_d2 / 2 * mm - 4), xytext=(0, -g_d2 / 2 * mm - 40),
                 arrowprops=dict(arrowstyle="-|>", color=C["opt"], lw=2))
    ax0.text(6, -g_d2 / 2 * mm - 32, f"{F:.0f} N ({F * LBF:.0f} lbf)\nworst radial + weight", fontsize=8,
             color=C["opt"])
    ax0.set_ylim(-g_d2 / 2 * mm - 48, g_d2 / 2 * mm + 30); ax0.grid(False)
    ax0.set_yticks([]); ax0.set_title("Rotor layout (lengths to scale, diameters labeled, vertical scale exaggerated)")

    ax1.fill_between(x * mm, moment, color=C["base"], alpha=0.25)
    ax1.plot(x * mm, moment, color=C["base"], lw=2)
    ax1.set_ylabel("Bending moment (N m)")
    ax1.set_title(f"Bending moment, peak {F * a:.0f} N m at the inboard bearing")

    v = bm["v"] * 1e6
    ax2.plot(x * mm, v, color=C["base"], lw=2, label="hydraulic load, max impeller, worst flow")
    ax2.axhline(0, color=C["muted"], lw=0.8)
    for sgn in (1, -1):
        ax2.axhline(sgn * 50, color=C["opt"], ls="--", lw=1)
    ax2.plot(xs, np.interp(r["seal"]["x"], x, bm["v"]) * 1e6, "o", color=C["opt"], ms=7,
             label=f"seal faces: {r['seal']['deflection_um']:.0f} µm (limit 50 µm)")
    ax2.set_ylabel("Deflection (µm)"); ax2.set_xlabel("Distance from impeller centerline (mm)")
    ax2.set_title("Shaft deflection"); ax2.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT / fname)
    plt.close(fig)


def fig_loads(de, fname):
    lc = de.results["load_curves"]
    req = Requirements()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4))
    for key, color, label in (("rated_imp", C["base"], "rated impeller"), ("max_imp", C["opt"], "max impeller (+5% head)")):
        d = lc[key]
        qg = d["q"] * GPM
        ax1.plot(qg, d["fr"], color=color, lw=2, label=label)
        nom = np.array([f["nominal"] for f in d["fa"]])
        lo = np.array([f["lo"] for f in d["fa"]])
        hi = np.array([f["hi"] for f in d["fa"]])
        ax2.fill_between(qg, lo, hi, color=color, alpha=0.15)
        ax2.plot(qg, nom, color=color, lw=2, label=label)
    for ax in (ax1, ax2):
        qb = lc["max_imp"]["q_bep"] * GPM
        ax.axvspan(req.mcsf_frac * qb, req.aor_max_frac * qb, color=C["por"], alpha=0.5, lw=0)
        ax.axvline(DUTY.q * GPM, color=C["ink"], ls="--", lw=0.9)
        ax.set_xlabel("Flow (gpm)")
    ax1.set_ylabel("Radial force on impeller (N)")
    ax1.set_title("Radial load (single volute, Stepanoff)")
    ax1.legend(fontsize=8)
    ax2.axhline(0, color=C["muted"], lw=0.8)
    ax2.set_ylabel("Axial thrust (N, + toward suction)")
    ax2.set_title("Axial thrust, with ±10% shroud allowance")
    ax2.text(0.02, 0.04, "shaded x-band: assumed operating range 40-120% BEP", transform=ax2.transAxes,
             fontsize=7.5, color=C["muted"])
    fig.tight_layout()
    fig.savefig(OUT / fname)
    plt.close(fig)


def fig_modes(de, fname):
    md = de.results["beam_modes"]
    lay = de.layout
    fig, ax = plt.subplots(figsize=(10, 3.4))
    for i, (color, ls) in enumerate(((C["base"], "-"), (C["alt"], "--"))):
        shape = md["modes"][i] / np.max(np.abs(md["modes"][i]))
        rpm = md["omega"][i] * 60 / (2 * math.pi)
        ax.plot(md["x"] * 1e3, shape, color=color, ls=ls, lw=2, label=f"mode {i + 1}: {rpm:,.0f} rpm")
    for xb in (md["a"], md["a"] + lay.span):
        ax.axvline(xb * 1e3, color=C["amber"], lw=6, alpha=0.4)
    ax.axhline(0, color=C["muted"], lw=0.8)
    ax.set_xlabel("Distance from impeller centerline (mm)"); ax.set_ylabel("Relative deflection")
    ax.set_title(f"Dry lateral modes (running speed {DUTY.n_rpm:.0f} rpm; bearings shaded)")
    ax.legend(fontsize=8.5)
    fig.tight_layout()
    fig.savefig(OUT / fname)
    plt.close(fig)


def stiffness_study(g, k, fname):
    ks = [5e7, 1e8, 1.5e8, 2e8, 3e8, 5e8, 1e9, 1e12]
    rows = []
    for kb in ks:
        try:
            de = design_drive_end(g, DUTY, k, Layout(bearing_stiffness=kb))
        except ValueError:      # the bearings alone deflect more than the limit allows
            rows.append(dict(k=kb, d_seal=None, d_brg=None, seal_um=None, crit=None))
            continue
        rows.append(dict(k=kb, d_seal=de.shaft.d_seal * 1e3, d_brg=de.shaft.d_brg * 1e3,
                         seal_um=de.results["seal"]["deflection_um"], crit=de.results["critical"]["rpm"]))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.6))
    lab = [f"{r['k']:.0e}".replace("e+0", "e").replace("e+", "e") if r["k"] < 1e11 else "rigid" for r in rows]
    xs = np.arange(len(rows))
    ax1.bar(xs, [r["d_seal"] or 0 for r in rows], color=[C["opt"] if r["k"] == 3e8 else C["base"] for r in rows])
    for i, r in enumerate(rows):
        if r["d_seal"] is None:
            ax1.text(i, 2, "not\nreachable", ha="center", fontsize=7.5, color=C["opt"])
        else:
            ax1.text(i, r["d_seal"] + 1, f"{r['d_seal']:.0f}", ha="center", fontsize=8)
    ax1.set_xticks(xs, lab, fontsize=8); ax1.set_xlabel("Bearing radial stiffness (N/m)")
    ax1.set_ylabel("Shaft dia. at seal (mm)"); ax1.set_title("Shaft needed to hold 50 µm at the seal")
    ax2.bar(xs, [(r["crit"] or 0) / 1e3 for r in rows], color=[C["opt"] if r["k"] == 3e8 else C["base"] for r in rows])
    ax2.axhline(1.2 * DUTY.n_rpm / 1e3, color=C["amber"], ls="--", lw=1.2, label="1.2 x running speed")
    ax2.set_xticks(xs, lab, fontsize=8); ax2.set_xlabel("Bearing radial stiffness (N/m)")
    ax2.set_ylabel("First dry critical (1000 rpm)"); ax2.set_title("First critical speed"); ax2.legend(fontsize=8)
    fig.suptitle("How much the bearing stiffness assumption matters (red = value used)", fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / fname)
    plt.close(fig)
    return rows


def plain(o):
    if isinstance(o, dict):
        return {k: plain(v) for k, v in o.items() if k not in ("beam_hyd", "beam_modes", "load_curves")}
    if isinstance(o, (list, tuple)):
        return [plain(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def main():
    g, k = build()
    de = design_drive_end(g, DUTY, k)
    fig_rotor(de, "drive_01_rotor.png")
    fig_loads(de, "drive_02_loads.png")
    fig_modes(de, "drive_03_modes.png")
    rows = stiffness_study(g, k, "drive_04_bearing_stiffness.png")
    out = dict(shaft=asdict(de.shaft), layout=asdict(de.layout), requirements=asdict(Requirements()),
               results=plain(de.results), stiffness_study=rows)
    (OUT / "drive_end.json").write_text(json.dumps(out, indent=1, default=float))
    r = de.results
    print(f"shaft: impeller {de.shaft.d_imp * 1e3:.0f} / seal {de.shaft.d_seal * 1e3:.0f} / "
          f"bearings {de.shaft.d_brg * 1e3:.0f} / span {de.shaft.d_span * 1e3:.0f} mm")
    for name, ok, val in r["checks"]:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}: {val}")


if __name__ == "__main__":
    main()
