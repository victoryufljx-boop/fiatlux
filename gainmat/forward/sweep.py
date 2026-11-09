"""Sweeps of the glass forward model.

  grid_sweep      every ion x glass family -> gain windows (the "what wavelengths can
                  rare-earth-doped fibre support" map)
  phonon_sweep    a continuous host knob: max phonon energy (with the IR edge and the
                  energy-gap law following it), everything else fixed -> which windows open
  mixture_sweep   linear interpolation of host descriptors between two families
                  (n, hw_max, site geometry, Judd-Ofelt Omega) as a stand-in for a
                  composition series, e.g. silica -> tellurite, ZBLAN -> fluoroindate
"""
from __future__ import annotations

import json
import os
import time

import numpy as np

from .glass import run, FAMILIES, omegas_for, PUMPS
from .ions import IONS


def _summ(r):
    return dict(ion=r["ion"], family=r["family"], omega_src=r["omega_src"], omegas=r["omegas"],
                n=r["n"], hw_max=r["hw_max"],
                windows=r["windows"],
                pumps={str(p): {k: v for k, v in d.items() if k in ("pump_level", "windows", "error")}
                       for p, d in r["pumps"].items()},
                tau_ms={p: d.get("tau_ms") for p, d in r["pumps"].items()}.popitem()[1] if r["pumps"] else None,
                eta={p: d.get("eta") for p, d in r["pumps"].items()}.popitem()[1] if r["pumps"] else None)


def grid_sweep(ions=None, families=None, n_sites=30, out=None, verbose=True):
    ions = ions or IONS
    families = families or list(FAMILIES)
    res = []
    for ion in ions:
        for fam in families:
            t = time.time()
            try:
                r = run(ion, fam, n_sites=n_sites)
                s = _summ(r); s["t_s"] = round(time.time() - t, 1)
            except Exception as e:
                s = dict(ion=ion, family=fam, error=repr(e))
            res.append(s)
            if verbose:
                print(ion, fam, s.get("t_s"), [(w["lo_nm"], w["hi_nm"]) for w in s.get("windows", [])], flush=True)
            if out:
                json.dump(res, open(out, "w"), indent=1, default=float)
    return res


def phonon_sweep(ion, base="ZBLAN", hws=(350, 450, 550, 650, 750, 850, 950, 1100, 1200),
                 n_sites=25, verbose=True):
    res = []
    for hw in hws:
        r = run(ion, base, n_sites=n_sites, overrides=dict(hw_max=float(hw), gap_family=None))
        res.append(dict(hw_max=hw, windows=r["windows"], eta=list(r["pumps"].values())[0].get("eta")))
        if verbose:
            print(ion, hw, [(w["lo_nm"], w["hi_nm"]) for w in r["windows"]], flush=True)
    return res


def mixture_sweep(ion, fam_a, fam_b, xs=np.linspace(0, 1, 6), n_sites=25, verbose=True):
    A, B = FAMILIES[fam_a], FAMILIES[fam_b]
    oa, _ = omegas_for(ion, fam_a); ob, _ = omegas_for(ion, fam_b)
    res = []
    for x in xs:
        ov = {k: (1 - x) * A[k] + x * B[k] for k in ("n", "hw_max", "CN", "R", "R_std", "ang_std", "loss_dB_m")}
        ov["omegas"] = tuple((1 - x) * np.array(oa) + x * np.array(ob))
        ov["gap_family"] = None                          # use the phonon-order law throughout
        if A["ligand"] != B["ligand"]:
            ov["ligand"] = A["ligand"] if x < 0.5 else B["ligand"]
        r = run(ion, fam_a, n_sites=n_sites, overrides=ov)
        res.append(dict(x=float(x), hw_max=ov["hw_max"], n=ov["n"], windows=r["windows"]))
        if verbose:
            print(ion, fam_a, "->", fam_b, round(x, 2), [(w["lo_nm"], w["hi_nm"]) for w in r["windows"]], flush=True)
    return res


def coverage_plot(grid, path, threshold_label="gain >= 0.1 dB/m above loss"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = [g for g in grid if g.get("windows")]
    ions = sorted({g["ion"] for g in rows}, key=lambda s: IONS.index(s) if s in IONS else 99)
    cmap = plt.get_cmap("tab10")
    fig, ax = plt.subplots(figsize=(11, 0.28 * len(rows) + 1.5))
    labels = []
    for y, g in enumerate(rows):
        c = cmap(ions.index(g["ion"]) % 10)
        for w in g["windows"]:
            ax.plot([w["lo_nm"], w["hi_nm"]], [y, y], lw=6, color=c, solid_capstyle="butt")
        labels.append(f"{g['ion'].replace('3+', '')} : {g['family']}")
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(labels, fontsize=7)
    ax.set_xscale("log"); ax.set_xlim(400, 5000)
    ticks = [400, 500, 600, 800, 1000, 1300, 1550, 2000, 2700, 3500, 5000]
    ax.set_xticks(ticks); ax.set_xticklabels([str(t) for t in ticks])
    from matplotlib.ticker import NullFormatter
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("wavelength (nm)")
    ax.set_title(f"Predicted fibre gain windows ({threshold_label}; site-ensemble lineshapes + rate equations)", fontsize=9)
    ax.grid(axis="x", alpha=0.3); ax.invert_yaxis()
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def phonon_plot(sw, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    fig, axs = plt.subplots(1, len(sw), figsize=(3.2 * len(sw), 4), sharey=True)
    for ax, (ion, rows) in zip(axs, sw.items()):
        for r in rows:
            for w in r["windows"]:
                ax.plot([w["lo_nm"], w["hi_nm"]], [r["hw_max"]] * 2, lw=5, color="C0", solid_capstyle="butt")
        ax.set_xscale("log"); ax.set_xlim(450, 5000); ax.set_title(ion.replace("3+", ""))
        ax.set_xticks([500, 1000, 2000, 3000, 5000]); ax.set_xticklabels(["0.5", "1", "2", "3", "5"])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel("wavelength (um)"); ax.grid(axis="x", alpha=0.3)
    axs[0].set_ylabel("host max phonon energy (cm$^{-1}$)")
    fig.suptitle("Predicted gain windows vs host phonon energy (fluoride-like site, everything else fixed)", fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)
