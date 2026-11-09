"""Coverage map + model-vs-experiment figure for the glass forward model.

    python scripts/plot_coverage_validation.py [forward_results]
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

R = sys.argv[1] if len(sys.argv) > 1 else "forward_results"
G = json.load(open(f"{R}/grid.json"))
V = json.load(open(f"{R}/validation.json"))
I = json.load(open(f"{R}/inputs.json"))
HW = {f: d["hw_max"] for f, d in I["families"].items()}

IONS = ["Pr3+", "Nd3+", "Yb3+", "Er3+", "Tm3+", "Ho3+", "Dy3+"]
COL = dict(zip(IONS, ["#7b3294", "#e66101", "#1b7837", "#d6604d", "#2166ac", "#b8860b", "#8c510a"]))
fams = sorted(HW, key=lambda f: (HW[f], f))           # low phonon energy first
LMIN, LMAX = 450, 4800
TICKS = [500, 600, 800, 1000, 1300, 1550, 2000, 2700, 3500, 4500]


def xaxis(ax, label=True):
    ax.set_xscale("log"); ax.set_xlim(LMIN, LMAX)
    ax.set_xticks(TICKS); ax.set_xticklabels([f"{t/1000:g}" for t in TICKS])
    ax.xaxis.set_minor_formatter(NullFormatter())
    if label:
        ax.set_xlabel("wavelength (µm)")
    ax.grid(axis="x", alpha=0.25)


fig = plt.figure(figsize=(10, 19))
gs = fig.add_gridspec(4, 2, height_ratios=[6.2, 1.15, 2.6, 2.6], hspace=0.45, wspace=0.32)

# ---------------------------------------------------------------- A: coverage map
ax = fig.add_subplot(gs[0, :])
lane = 0.095
for row, fam in enumerate(fams):
    y0 = row
    ax.axhspan(y0 - 0.46, y0 + 0.46, color="0.96" if row % 2 else "white", zorder=0)
    for k, ion in enumerate(IONS):
        y = y0 - 0.33 + k * lane
        g = next((x for x in G if x["family"] == fam and x["ion"] == ion), None)
        for w in (g or {}).get("windows", []):
            if w["hi_nm"] - w["lo_nm"] >= 5:
                ax.plot([w["lo_nm"], w["hi_nm"]], [y, y], color=COL[ion], lw=4.2,
                        solid_capstyle="butt", zorder=2)
    # measured data for this family, on a lane above the ions
    yv = y0 + 0.36
    for r in V:
        if r.get("family") != fam:
            continue
        b = r.get("band_exp_nm")
        if b:
            ax.plot(b, [yv, yv], color="k", lw=2.2, zorder=4)
            ax.plot(b, [yv, yv], "|", color="k", ms=7, zorder=4)
        pk = r.get("exp_lambda_em_peak_nm") or r.get("exp_lambda_abs_peak_nm")
        if pk:
            ax.plot(pk, yv, "D", color="k", ms=5, zorder=5)
        mp = r.get("model_em_peak_nm")
        if mp and (b or pk):
            ax.plot(mp, yv, "x", color="red", ms=7, mew=2, zorder=6)
ax.set_yticks(range(len(fams)))
ax.set_yticklabels([f"{f}\n{HW[f]:.0f} cm⁻¹" for f in fams], fontsize=8.5)
ax.set_ylim(-0.6, len(fams) - 0.4)
ax.set_ylabel("glass family  (sorted by max phonon energy)")
xaxis(ax)
ax.set_title("A. Where the model predicts gain (≥0.1 dB/m above loss), and where measured fibres sit",
             fontsize=10.5, loc="left")
handles = [Patch(color=COL[i], label=i.replace("3+", "")) for i in IONS]
handles += [Line2D([], [], color="k", lw=2.2, label="measured gain band"),
            Line2D([], [], color="k", marker="D", ls="", ms=5, label="measured peak"),
            Line2D([], [], color="red", marker="x", ls="", ms=7, mew=2, label="model peak")]
ax.legend(handles=handles, ncol=5, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.07),
          frameon=False)

# ---------------------------------------------------------------- B: coverage strip
axb = fig.add_subplot(gs[1, :])
grid_l = np.geomspace(LMIN, LMAX, 3000)


def covered(families):
    c = np.zeros_like(grid_l, bool)
    for g in G:
        if g["family"] in families:
            for w in g.get("windows", []):
                if w["hi_nm"] - w["lo_nm"] >= 5:
                    c |= (grid_l >= w["lo_nm"]) & (grid_l <= w["hi_nm"])
    return c


strips = [("any glass", list(HW)),
          ("silica-based\n(MCVD today)", ["silica", "aluminosilicate", "germanosilicate", "phosphosilicate"])]
for k, (lab, fl) in enumerate(strips):
    c = covered(fl)
    y = 1 - k
    for i in range(len(grid_l) - 1):
        axb.fill_between(grid_l[i:i + 2], y - 0.35, y + 0.35,
                         color="#4daf4a" if c[i] else "#f4a3a3", lw=0)
axb.set_yticks([1, 0]); axb.set_yticklabels([s[0] for s in strips], fontsize=8.5)
axb.set_ylim(-0.6, 1.6)
xaxis(axb)
axb.set_title("B. Union coverage: green = some ion gives gain, red = gap", fontsize=10.5, loc="left")

# ---------------------------------------------------------------- C: parity plots
def parity(ax, xs, ys, cols, labels, lo, hi, log, title, unit):
    ax.plot([lo, hi], [lo, hi], color="0.5", lw=1)
    if log:
        ax.fill_between([lo, hi], [lo / 2, hi / 2], [lo * 2, hi * 2], color="0.92", lw=0)
        ax.set_xscale("log"); ax.set_yscale("log")
    for x, y, c, lab in zip(xs, ys, cols, labels):
        ax.plot(x, y, "o", color=c, ms=6, mec="k", mew=0.5)
        ax.annotate(lab, (x, y), fontsize=6.5, xytext=(3, 3), textcoords="offset points")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    from matplotlib.ticker import FuncFormatter
    fmt = FuncFormatter(lambda v, _: f"{v:g}")
    ax.xaxis.set_major_formatter(fmt); ax.yaxis.set_major_formatter(fmt)
    ax.xaxis.set_minor_formatter(NullFormatter()); ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(f"measured ({unit})", fontsize=8.5); ax.set_ylabel(f"model ({unit})", fontsize=8.5)
    ax.set_title(title, fontsize=9.5, loc="left")
    ax.tick_params(labelsize=7.5)


def short(r):
    return r["entry"].replace(" band", "").replace("aluminosilicate", "Al-sil").replace("phosphosilicate", "P-sil")


# C1 wavelength offsets (nm): peak, short edge, long edge
axc = fig.add_subplot(gs[2, :])
rows = []
for r in V:
    pk = r.get("exp_lambda_em_peak_nm") or r.get("exp_lambda_abs_peak_nm")
    mp = r.get("model_em_peak_nm") if r.get("exp_lambda_em_peak_nm") else r.get("model_abs_peak_nm")
    b, w = r.get("band_exp_nm"), r.get("window_model_nm") or []
    d = {}
    if pk and mp:
        d["peak"] = mp - pk
    if b and w:
        lo = min(x[0] for x in w); hi = max(x[1] for x in w)
        d["short edge"], d["long edge"] = lo - b[0], hi - b[1]
    if d:
        rows.append((short(r), r["ion"], d))
mk = {"peak": ("D", "k"), "short edge": ("<", "#2166ac"), "long edge": (">", "#d6604d")}
for i, (lab, ion, d) in enumerate(rows):
    for k, v in d.items():
        axc.plot(v, i, mk[k][0], color=mk[k][1], ms=7, mec="k", mew=0.4)
axc.axvline(0, color="0.4", lw=1)
axc.axvspan(-50, 50, color="0.92", zorder=0)
axc.set_yticks(range(len(rows))); axc.set_yticklabels([r[0] for r in rows], fontsize=8)
axc.invert_yaxis(); axc.set_xlim(-300, 350)
axc.set_xlabel("model − measured (nm)", fontsize=9)
axc.set_title("C1. How far off: peak position and gain-band edges (grey = within ±50 nm)", fontsize=10, loc="left")
axc.legend(handles=[Line2D([], [], marker=m, color=c, ls="", ms=7, label=k) for k, (m, c) in mk.items()],
           fontsize=8, loc="lower right")
axc.grid(axis="x", alpha=0.3)

# C2 emission cross section
xs, ys, cs, ls = [], [], [], []
for r in V:
    e, m = r.get("exp_sigma_em"), r.get("model_sigma_em")
    if isinstance(e, (list, tuple)):
        e = float(np.mean([v for v in e if v]))
    if e and m:
        xs.append(e * 1e21); ys.append(m * 1e21); cs.append(COL.get(r["ion"], "k")); ls.append(short(r))
parity(fig.add_subplot(gs[3, 0]), xs, ys, cs, ls, 0.2, 50, True,
       "C2. Peak emission cross section", "×10⁻²¹ cm²")

# C3 upper-level lifetime
xs, ys, cs, ls = [], [], [], []
for r in V:
    e, m = r.get("exp_tau_ms"), r.get("model_tau_ms")
    if isinstance(e, (list, tuple)):
        e = float(np.mean([v for v in e if v]))
    if e and m and "Yb" not in r["ion"]:            # Yb lifetime is an input
        xs.append(e); ys.append(m); cs.append(COL.get(r["ion"], "k")); ls.append(short(r))
parity(fig.add_subplot(gs[3, 1]), xs, ys, cs, ls, 0.05, 40, True, "C3. Upper-level lifetime", "ms")

fig.text(0.01, 0.005, "Grey band in C2/C3 = within 2× of measured. Yb lifetime is an input (not shown). "
         "Measured values: gainmat/refs/fiber_validation.py; model: forward_results/.", fontsize=7.5)
out = f"{R}/coverage_vs_experiment.png"
fig.savefig(out, dpi=150, bbox_inches="tight"); print(out)
