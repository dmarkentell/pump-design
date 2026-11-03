"""Impeller design study for an API 610 OH2-style pump.

Runs the whole step-1 study and writes every figure and number to outputs/:
  1. Calibrate the loss model and size a baseline impeller for the duty point.
  2. Performance curves, loss breakdown, velocity triangles, blade plan view.
  3. One-at-a-time sweeps of blade count, outlet angle, outlet width and eye diameter.
  4. Full blade-count x outlet-angle map.
  5. Random design-space sample: efficiency vs NPSHr trade-off.
  6. Constrained optimization (differential evolution) for best efficiency at the rated point.

    python run_study.py
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
from scipy.optimize import differential_evolution
from scipy.stats import qmc

from pump.calibrate import calibrate
from pump.correlations import eta_hydraulic, eta_overall
from pump.meanline import G, bep, curve, evaluate
from pump.sizing import DesignChoices, Duty, design_bep_matched, summarize

OUT = Path(__file__).parent / "outputs"
OUT.mkdir(exist_ok=True)

# ---- duty point: an illustrative refinery/process service, not a real pump ------------
DUTY = Duty(q=150 / 3600, head=90.0, n_rpm=3560, npsha=10.0)

# ---- screening requirements -------------------------------------------------------------
REQ = dict(
    npsh_margin_m=1.0,      # NPSHa - NPSHr at rated flow
    rated_pct_bep=(80, 110),  # API 610 rated-flow window
    head_rise_pct=10.0,     # rise from rated to shutoff (commonly specified)
    droop_pct=0.05,         # head keeps rising as flow drops (judged down to 10% of rated);
                            # 0.05% tolerance absorbs numerical noise on a flat curve
    de_haller=0.70,         # w2/w1 relative-velocity ratio, a diffusion limit
    b2_min_mm=8.0,          # narrowest outlet a sand casting comfortably makes
    z_range=(5, 9),         # usual practice at this specific speed; 1D losses cannot see the
                            # flow-guidance and pressure-pulsation penalties of fewer blades
)

# ---- units --------------------------------------------------------------------------------
GPM = 15850.32
FT = 3.28084
HP = 1 / 745.7

C = dict(ink="#1f2a37", muted="#6b7280", grid="#e5e7eb", base="#2563eb", opt="#dc2626",
         alt="#059669", amber="#d97706", por="#dbeafe")
plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 160, "font.size": 9.5, "axes.titlesize": 10.5,
    "axes.titleweight": "bold", "axes.edgecolor": C["muted"], "axes.labelcolor": C["ink"],
    "axes.grid": True, "grid.color": C["grid"], "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "xtick.color": C["muted"],
    "ytick.color": C["muted"], "legend.frameon": False,
})


def violations(s: dict) -> float:
    v = 0.0
    v += max(0.0, REQ["npsh_margin_m"] - s["npsh_margin"])
    lo, hi = REQ["rated_pct_bep"]
    v += max(0.0, lo - s["rated_pct_bep"]) / 10 + max(0.0, s["rated_pct_bep"] - hi) / 10
    v += max(0.0, REQ["head_rise_pct"] - s["head_rise_pct"]) / 5
    v += max(0.0, s["droop_pct"] - REQ["droop_pct"])
    v += max(0.0, REQ["de_haller"] - s["de_haller"]) * 10
    v += max(0.0, REQ["b2_min_mm"] - s["b2_mm"]) / 2
    v += 1.0 if not (REQ["z_range"][0] <= s["z"] <= REQ["z_range"][1]) else 0.0
    return v


def evaluate_choices(ch: DesignChoices, k) -> dict | None:
    try:
        g, ch2 = design_bep_matched(DUTY, ch, k=k)
        s = summarize(g, DUTY, k=k)
    except (ValueError, ZeroDivisionError):
        return None
    s["viol"] = violations(s)
    s["throat_factor"] = ch2.throat_factor
    s["geom"] = g
    s["choices"] = ch2
    return s


# ---------------------------------------------------------------------------------------
def fig_performance(designs: list[tuple[str, object, str]], k, fname: str, title: str):
    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    qmax = 1.5 * DUTY.q
    for label, g, color in designs:
        c = curve(g, qmax, DUTY.n_rpm, k=k, npts=91)
        q_bep, _ = bep(g, DUTY.q, DUTY.n_rpm, k=k)
        qg = c["q"] * GPM
        axs[0, 0].plot(qg, c["head"] * FT, color=color, lw=2, label=label)
        axs[0, 1].plot(qg, 100 * c["eta"], color=color, lw=2, label=label)
        axs[1, 0].plot(qg, c["p_shaft"] * HP, color=color, lw=2, label=label)
        axs[1, 1].plot(qg, c["npshr"] * FT, color=color, lw=2, label=label)
        if designs[0][1] is g:
            for ax in axs.flat:
                ax.axvspan(0.7 * q_bep * GPM, 1.2 * q_bep * GPM, color=C["por"], alpha=0.55, lw=0,
                           label="Preferred operating region (70-120% BEP)")
    axs[1, 1].axhline(DUTY.npsha * FT, color=C["amber"], ls="--", lw=1.4, label="NPSH available")
    axs[1, 1].axhline((DUTY.npsha - REQ["npsh_margin_m"]) * FT, color=C["amber"], ls=":", lw=1.2,
                      label="NPSHa - 1 m margin")
    for ax in axs.flat:
        ax.axvline(DUTY.q * GPM, color=C["ink"], lw=0.9, ls="--")
        ax.set_xlabel("Flow (gpm)")
        ax.set_xlim(0, qmax * GPM)
    axs[0, 0].plot(DUTY.q * GPM, DUTY.head * FT, "o", color=C["ink"], ms=7, zorder=5, label="Rated point")
    axs[0, 0].set_ylabel("Total head (ft)"); axs[0, 0].set_title("Head")
    axs[0, 1].set_ylabel("Pump efficiency (%)"); axs[0, 1].set_title("Efficiency"); axs[0, 1].set_ylim(0, 90)
    axs[1, 0].set_ylabel("Shaft power (hp)"); axs[1, 0].set_title("Power")
    axs[1, 1].set_ylabel("NPSH required (ft)"); axs[1, 1].set_title("NPSH required (3% head drop, est.)")
    axs[0, 0].legend(loc="lower left", fontsize=8)
    axs[1, 1].legend(loc="upper left", fontsize=8)
    fig.suptitle(title, fontweight="bold", fontsize=12, color=C["ink"])
    fig.tight_layout()
    fig.savefig(OUT / fname)
    plt.close(fig)


def fig_losses(g, k, fname: str):
    qs = np.linspace(0.02, 1.5, 80) * DUTY.q
    rows = [evaluate(g, q, DUTY.n_rpm, k=k) for q in qs]
    names = ["incidence", "impeller_friction", "volute_mismatch", "volute_friction", "discharge_cone"]
    labels = ["Inlet incidence", "Impeller friction", "Volute mismatch", "Volute friction", "Discharge cone"]
    colors = ["#f59e0b", "#2563eb", "#dc2626", "#7c3aed", "#059669"]
    stack = np.array([[r["losses"][n] for r in rows] for n in names]) * FT
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8))
    ax1.stackplot(qs * GPM, stack, labels=labels, colors=colors, alpha=0.85)
    ax1.axvline(DUTY.q * GPM, color=C["ink"], lw=0.9, ls="--")
    ax1.set_xlabel("Flow (gpm)"); ax1.set_ylabel("Head lost (ft)")
    ax1.set_title("Where the head goes"); ax1.legend(fontsize=8, loc="upper center", ncol=2)

    r = evaluate(g, DUTY.q, DUTY.n_rpm, k=k)
    p_shaft = r["p_shaft"]
    rho = 998.2
    budget = {
        "Delivered to fluid": rho * G * DUTY.q * r["head"],
        "Hydraulic losses": rho * G * (DUTY.q + r["q_leak"]) * (r["h_th"] - r["head"]),
        "Wear-ring leakage": rho * G * r["q_leak"] * r["head"],
        "Disk friction": r["p_disk"],
        "Bearings + seal": p_shaft * 0.01,
    }
    tot = sum(budget.values())
    vals = [100 * v / tot for v in budget.values()]
    bcol = [C["base"], "#93c5fd", "#fca5a5", "#fcd34d", "#d1d5db"]
    left = 0
    for (name, v), col in zip(zip(budget, vals), bcol):
        ax2.barh(0, v, left=left, color=col, edgecolor="white")
        if v > 3:
            ax2.text(left + v / 2, 0, f"{v:.1f}%", ha="center", va="center", fontsize=8.5,
                     color="white" if col == C["base"] else C["ink"])
        left += v
    ax2.set_yticks([]); ax2.set_xlim(0, 100); ax2.set_xlabel("% of shaft power at rated flow")
    ax2.set_title(f"Power budget at rated flow ({p_shaft * HP:.0f} hp)")
    ax2.legend([plt.Rectangle((0, 0), 1, 1, color=c) for c in bcol], budget.keys(), fontsize=8,
               loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=3)
    ax2.grid(False)
    fig.tight_layout()
    fig.savefig(OUT / fname)
    plt.close(fig)


def fig_triangles(g, k, fname: str):
    r = evaluate(g, DUTY.q, DUTY.n_rpm, k=k)
    u1m = math.pi * g.d1m * DUTY.n_rpm / 60
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))

    def tri(ax, u, cu, cm, title, beta_b):
        ax.annotate("", xy=(u, 0), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=C["ink"], lw=2))
        ax.annotate("", xy=(cu, cm), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=C["base"], lw=2))
        ax.annotate("", xy=(cu, cm), xytext=(u, 0), arrowprops=dict(arrowstyle="-|>", color=C["opt"], lw=2))
        bl = 0.45 * u
        ax.plot([u, u - bl * math.cos(math.radians(beta_b))], [0, bl * math.sin(math.radians(beta_b))],
                color=C["muted"], ls=":", lw=1.5)
        ax.text(u / 2, -0.07 * u, f"u = {u:.1f} m/s (blade speed)", ha="center", va="top", fontsize=8.5)
        ax.text(cu / 2 - 0.02 * u, cm / 2 + 0.03 * u, f"c = {math.hypot(cu, cm):.1f}", color=C["base"],
                ha="right", fontsize=8.5)
        ax.text((cu + u) / 2 + 0.02 * u, cm / 2 + 0.03 * u, f"w = {math.hypot(u - cu, cm):.1f}",
                color=C["opt"], fontsize=8.5)
        ax.set_aspect("equal"); ax.set_xlim(-0.05 * u, 1.08 * u); ax.set_ylim(-0.18 * u, max(cm, 0.2 * u) * 1.5)
        ax.set_title(title); ax.grid(False); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)

    tri(axs[0], u1m, 0.0, r["cm1"], f"Inlet (mean eye dia.) - flow angle {r['beta1_flow']:.1f}°, "
        f"blade {g.beta1b:.1f}°", g.beta1b)
    tri(axs[1], r["u2"], r["cu2"], r["cm2"], f"Outlet - blade {g.beta2b:.1f}°, slip factor {r['sigma']:.3f}",
        g.beta2b)
    fig.legend([plt.Line2D([], [], color=C["ink"], lw=2), plt.Line2D([], [], color=C["base"], lw=2),
                plt.Line2D([], [], color=C["opt"], lw=2), plt.Line2D([], [], color=C["muted"], ls=":")],
               ["Blade speed u", "Absolute velocity c", "Relative velocity w", "Blade direction"],
               loc="lower center", ncol=4, fontsize=8.5)
    fig.suptitle("Velocity triangles at rated flow (m/s)", fontweight="bold")
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(OUT / fname)
    plt.close(fig)


def blade_curve(g, npts=60):
    """Blade camber line in plan view with blade angle varying linearly with radius.

    d(theta) = dr / (r tan(beta)): the classic way to lay out a 2D (cylindrical) blade.
    """
    r1, r2 = g.d1m / 2, g.d2 / 2
    rs = np.linspace(r1, r2, npts)
    betas = np.radians(g.beta1b + (g.beta2b - g.beta1b) * (rs - r1) / (r2 - r1))
    dth = np.gradient(rs) / (rs * np.tan(betas))
    th = np.concatenate([[0], np.cumsum((dth[1:] + dth[:-1]) / 2 * 1)])
    return rs, th


def fig_planview(designs, fname):
    fig, axs = plt.subplots(1, len(designs), figsize=(4.6 * len(designs), 4.6))
    axs = np.atleast_1d(axs)
    for ax, (label, g, color) in zip(axs, designs):
        rs, th = blade_curve(g)
        for i in range(g.z):
            a = 2 * math.pi * i / g.z
            ax.plot(rs * 1e3 * np.cos(-th + a), rs * 1e3 * np.sin(-th + a), color=color, lw=2.4)
        t = np.linspace(0, 2 * math.pi, 200)
        for d, ls in ((g.d2, "-"), (g.d1, "--"), (g.dh, "-")):
            ax.plot(d / 2 * 1e3 * np.cos(t), d / 2 * 1e3 * np.sin(t), color=C["muted"], lw=1, ls=ls)
        ax.fill(g.dh / 2 * 1e3 * np.cos(t), g.dh / 2 * 1e3 * np.sin(t), color=C["grid"])
        ax.set_aspect("equal"); ax.grid(False)
        ax.set_title(f"{label}\nD2 {g.d2 * 1e3:.0f} mm, {g.z} blades, β2 {g.beta2b:.1f}°", color=color)
        ax.set_xlabel("mm")
        ra = g.d2 / 2 * 1e3 * 1.12
        t2 = np.linspace(math.radians(20), math.radians(70), 30)
        ax.plot(ra * np.cos(t2), ra * np.sin(t2), color=C["ink"], lw=1.2)
        ax.annotate("", xy=(ra * math.cos(t2[-1] + 0.05), ra * math.sin(t2[-1] + 0.05)),
                    xytext=(ra * math.cos(t2[-1]), ra * math.sin(t2[-1])),
                    arrowprops=dict(arrowstyle="-|>", color=C["ink"], lw=1.2))
        ax.text(ra * 0.98 * math.cos(math.radians(45)) + 6, ra * 0.98 * math.sin(math.radians(45)) + 6,
                "rotation", fontsize=8, color=C["ink"])
        lim = g.d2 / 2 * 1e3 * 1.3
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    fig.suptitle("Impeller plan view (blade camber lines, cylindrical blades)", fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / fname)
    plt.close(fig)


def fig_sweeps(base_ch: DesignChoices, k, fname: str) -> dict:
    sweeps = {
        "z": ("Blade count", [4, 5, 6, 7, 8, 9]),
        "beta2b": ("Outlet blade angle β2 (deg)", list(np.arange(14, 37, 2.0))),
        "b2_factor": ("Outlet width (x correlation)", list(np.linspace(0.6, 1.4, 11))),
        "d1_factor": ("Eye diameter (x min-NPSHr value)", list(np.linspace(0.85, 1.25, 11))),
    }
    metrics = [("eta", "Efficiency at rated (%)", 100), ("d2_mm", "Impeller dia. D2 (mm)", 1),
               ("head_rise_pct", "Head rise to shutoff (%)", 1), ("npshr", "NPSHr at rated (ft)", FT)]
    fig, axs = plt.subplots(len(metrics), len(sweeps), figsize=(12, 9.5), sharey="row")
    data = {}
    for j, (var, (label, vals)) in enumerate(sweeps.items()):
        rows = []
        for v in vals:
            s = evaluate_choices(replace(base_ch, **{var: type(getattr(base_ch, var))(v)}), k)
            rows.append(s)
        data[var] = [(float(v), {m: (s[m] if s else None) for m, _, _ in metrics} | {"viol": s["viol"] if s else None})
                     for v, s in zip(vals, rows)]
        for i, (m, mlabel, scale) in enumerate(metrics):
            ax = axs[i, j]
            xs = [v for v, s in zip(vals, rows) if s]
            ys = [s[m] * scale for s in rows if s]
            ok = [s["viol"] < 1e-9 for s in rows if s]
            ax.plot(xs, ys, color=C["base"], lw=1.8)
            ax.scatter([x for x, o in zip(xs, ok) if o], [y for y, o in zip(ys, ok) if o], color=C["base"], s=18, zorder=3)
            ax.scatter([x for x, o in zip(xs, ok) if not o], [y for y, o in zip(ys, ok) if not o],
                       facecolor="white", edgecolor=C["opt"], s=22, zorder=3)
            ax.axvline(getattr(base_ch, var), color=C["muted"], ls=":", lw=1)
            if j == 0:
                ax.set_ylabel(mlabel)
            if i == len(metrics) - 1:
                ax.set_xlabel(label)
        axs[2, 0].axhline(REQ["head_rise_pct"], color=C["amber"], ls="--", lw=1)
        axs[3, 0].axhline((DUTY.npsha - REQ["npsh_margin_m"]) * FT, color=C["amber"], ls="--", lw=1)
    for ax in axs[2]:
        ax.axhline(REQ["head_rise_pct"], color=C["amber"], ls="--", lw=1)
    for ax in axs[3]:
        ax.axhline((DUTY.npsha - REQ["npsh_margin_m"]) * FT, color=C["amber"], ls="--", lw=1)
    fig.suptitle("One change at a time from the baseline (every point re-sized to hit the rated head; "
                 "open red = fails a requirement; dotted = baseline)", fontweight="bold", fontsize=10.5)
    fig.tight_layout()
    fig.savefig(OUT / fname)
    plt.close(fig)
    return data


def fig_map(base_ch, k, fname):
    zs = [4, 5, 6, 7, 8, 9]
    betas = np.arange(14, 37, 2.0)
    eta = np.full((len(zs), len(betas)), np.nan)
    rise = np.full_like(eta, np.nan)
    feas = np.zeros_like(eta, dtype=bool)
    for i, z in enumerate(zs):
        for j, b in enumerate(betas):
            s = evaluate_choices(replace(base_ch, z=z, beta2b=float(b)), k)
            if s:
                eta[i, j] = 100 * s["eta"]; rise[i, j] = s["head_rise_pct"]; feas[i, j] = s["viol"] < 1e-9
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    for ax, arr, title, cmap in ((axs[0], eta, "Efficiency at rated flow (%)", "viridis"),
                                 (axs[1], rise, "Head rise to shutoff (%)", "magma")):
        im = ax.imshow(arr, origin="lower", aspect="auto", cmap=cmap,
                       extent=(betas[0] - 1, betas[-1] + 1, zs[0] - 0.5, zs[-1] + 0.5))
        fig.colorbar(im, ax=ax, shrink=0.9)
        norm = im.norm
        cm = plt.get_cmap(cmap)
        for i, z in enumerate(zs):
            for j, b in enumerate(betas):
                if np.isnan(arr[i, j]):
                    continue
                if not feas[i, j]:   # fade cells that fail a check
                    ax.add_patch(plt.Rectangle((b - 1, z - 0.5), 2, 1, facecolor="white", alpha=0.6,
                                               edgecolor="none"))
                    color = "#6b7280"
                else:
                    r_, g_, b_, _ = cm(norm(arr[i, j]))
                    color = "white" if 0.299 * r_ + 0.587 * g_ + 0.114 * b_ < 0.55 else "#111827"
                ax.text(b, z, f"{arr[i, j]:.1f}" if arr is eta else f"{arr[i, j]:.0f}", ha="center",
                        va="center", fontsize=7, color=color,
                        fontweight="bold" if feas[i, j] else "normal")
        ax.set_xlabel("Outlet blade angle β2 (deg)"); ax.set_ylabel("Blade count"); ax.set_title(title)
        ax.grid(False)
    fig.suptitle("Blade count x outlet angle (bold = meets every requirement, faded = fails one)",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / fname)
    plt.close(fig)
    return dict(z=zs, beta2b=betas.tolist(), eta=np.nan_to_num(eta).tolist(), rise=np.nan_to_num(rise).tolist(),
                feasible=feas.tolist())


BOUNDS = dict(beta2b=(14.0, 36.0), b2_factor=(0.6, 1.4), d1_factor=(0.85, 1.25), incidence=(0.0, 6.0))


def sample_space(k, n=1000, seed=1):
    sampler = qmc.LatinHypercube(d=5, seed=seed)
    u = sampler.random(n)
    pts = []
    for row in u:
        z = int(4 + min(5, int(row[0] * 6)))   # 4-9, so the plot shows what 4 blades would do
        vals = {key: lo + (hi - lo) * x for (key, (lo, hi)), x in zip(BOUNDS.items(), row[1:])}
        s = evaluate_choices(DesignChoices(z=z, **vals), k)
        if s:
            pts.append(s)
    return pts


def fig_tradeoff(pts, base, opt, fname):
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ok = [p for p in pts if p["viol"] < 1e-9]
    bad = [p for p in pts if p["viol"] >= 1e-9]
    ax.scatter([p["npshr"] * FT for p in bad], [100 * p["eta"] for p in bad], s=8, color="#d1d5db",
               label=f"Fails a requirement ({len(bad)})")
    sc = ax.scatter([p["npshr"] * FT for p in ok], [100 * p["eta"] for p in ok], s=10,
                    c=[p["z"] for p in ok], cmap="viridis", label=f"Meets all ({len(ok)})")
    fig.colorbar(sc, ax=ax, label="Blade count")
    ax.scatter(base["npshr"] * FT, 100 * base["eta"], marker="s", s=70, color=C["base"], edgecolor="white",
               zorder=5, label="Baseline")
    ax.scatter(opt["npshr"] * FT, 100 * opt["eta"], marker="*", s=200, color=C["opt"], edgecolor="white",
               zorder=6, label="Optimized")
    ax.axvline((DUTY.npsha - REQ["npsh_margin_m"]) * FT, color=C["amber"], ls="--", lw=1.2,
               label="NPSHr limit (NPSHa - 1 m)")
    ax.set_xlabel("NPSH required at rated flow (ft)"); ax.set_ylabel("Efficiency at rated flow (%)")
    ax.set_title(f"{len(pts)} random impeller designs, each sized to the same duty")
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT / fname)
    plt.close(fig)


def optimize(k, seed=7):
    best = None
    per_z = {}
    keys = list(BOUNDS)
    for z in range(REQ["z_range"][0], REQ["z_range"][1] + 1):
        def obj(x):
            s = evaluate_choices(DesignChoices(z=z, **dict(zip(keys, x))), k)
            if s is None:
                return 10.0
            return -s["eta"] + 2.0 * s["viol"]
        res = differential_evolution(obj, list(BOUNDS.values()), seed=seed, popsize=10, maxiter=40,
                                     tol=1e-7, polish=False)
        s = evaluate_choices(DesignChoices(z=z, **dict(zip(keys, res.x))), k)
        per_z[z] = dict(eta=s["eta"], viol=s["viol"], x=dict(zip(keys, map(float, res.x))), nfev=int(res.nfev))
        if s["viol"] < 1e-9 and (best is None or s["eta"] > best["eta"]):
            best = s
    return best, per_z


def clean(s: dict) -> dict:
    out = {k: v for k, v in s.items() if k not in ("geom", "choices")}
    out["geometry"] = asdict(s["geom"])
    out["choices"] = asdict(s["choices"])
    return out


def main():
    k = calibrate(DUTY)
    base_ch = DesignChoices()
    base = evaluate_choices(base_ch, k)
    base_ch = base["choices"]
    print(f"calibration: friction-loss scale {k.hyd_scale:.3f}")
    print(f"baseline: D2 {base['d2_mm']:.1f} mm, eta {100 * base['eta']:.1f}%, NPSHr {base['npshr']:.2f} m")

    sweeps = fig_sweeps(base_ch, k, "04_sweeps.png")
    grid = fig_map(base_ch, k, "05_blade_map.png")
    pts = sample_space(k)
    opt, per_z = optimize(k)
    print(f"optimized: z {opt['z']}, beta2 {opt['beta2b']:.1f}, D2 {opt['d2_mm']:.1f} mm, eta {100 * opt['eta']:.2f}%")

    designs = [("Baseline", base["geom"], C["base"]), ("Optimized", opt["geom"], C["opt"])]
    fig_performance(designs[:1], k, "01_baseline_curves.png",
                    f"Baseline impeller: {DUTY.q * GPM:.0f} gpm at {DUTY.head * FT:.0f} ft, {DUTY.n_rpm:.0f} rpm")
    fig_losses(base["geom"], k, "02_losses.png")
    fig_triangles(base["geom"], k, "03_velocity_triangles.png")
    fig_tradeoff(pts, base, opt, "06_tradeoff.png")
    fig_performance(designs, k, "07_baseline_vs_optimized.png", "Baseline vs optimized impeller")
    fig_planview(designs, "08_plan_view.png")

    # curves for the interactive page and the write-up
    curves = {}
    for name, s in (("baseline", base), ("optimized", opt)):
        c = curve(s["geom"], 1.5 * DUTY.q, DUTY.n_rpm, k=k, npts=46)
        curves[name] = {key: np.round(c[key], 5).tolist() for key in ("q", "head", "eta", "p_shaft", "npshr")}

    feas = [p for p in pts if p["viol"] < 1e-9]
    results = dict(
        duty=dict(q_m3h=DUTY.q * 3600, q_gpm=DUTY.q * GPM, head_m=DUTY.head, head_ft=DUTY.head * FT,
                  n_rpm=DUTY.n_rpm, npsha_m=DUTY.npsha, nq=DUTY.nq, ns_us=DUTY.ns_us),
        requirements=REQ,
        calibration=dict(friction_scale=k.hyd_scale, eta_h_target=eta_hydraulic(DUTY.q, DUTY.nq),
                         eta_overall_correlation=eta_overall(DUTY.q, DUTY.nq)),
        loss_coeffs=asdict(k),
        baseline=clean(base), optimized=clean(opt), optimizer_per_z=per_z,
        sample=dict(n=len(pts), feasible=len(feas),
                    eta_best_feasible=max(p["eta"] for p in feas),
                    eta_median_feasible=float(np.median([p["eta"] for p in feas]))),
        sweeps=sweeps, blade_map=grid, curves=curves,
    )
    (OUT / "results.json").write_text(json.dumps(results, indent=1, default=float))
    print(f"sample: {len(feas)}/{len(pts)} feasible")


if __name__ == "__main__":
    main()
