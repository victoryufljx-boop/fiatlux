"""Run the glass forward model end to end: validation, ion x glass grid, continuous sweeps.

    python scripts/run_forward.py [--sites 30] [--out results/forward]
Writes JSON results, coverage / sweep plots and a markdown summary.
"""
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np

from gainmat.forward.glass import (load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ,
                                   FAMILIES, OMEGA, SPM_DEFAULT)
from gainmat.forward.sweep import grid_sweep, phonon_sweep, mixture_sweep, coverage_plot

ap = argparse.ArgumentParser()
ap.add_argument("--sites", type=int, default=30)
ap.add_argument("--out", default="results/forward")
ap.add_argument("--skip", default="")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
skip = set(a.skip.split(","))
t0 = time.time()
loaded = load_refs()
site_used = apply_site_refs()
cal = calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ)
meta = dict(loaded_refs=loaded, site_refs_used=site_used, spm_calibration=cal,
            families=FAMILIES, spm=SPM_DEFAULT, spm_scale=__import__("gainmat.forward.glass", fromlist=["SPM_SCALE"]).SPM_SCALE,
            omega={f"{k[0]}:{k[1]}": v for k, v in OMEGA.items()})
json.dump(meta, open(f"{a.out}/inputs.json", "w"), indent=1, default=str)
print("refs loaded:", loaded, flush=True)

if "validate" not in skip:
    from gainmat.forward.validate import validate
    try:
        rows = validate(n_sites=a.sites)
    except ImportError as e:
        rows = [dict(error=f"no validation refs: {e}")]
    json.dump(rows, open(f"{a.out}/validation.json", "w"), indent=1, default=float)

if "grid" not in skip:
    g = grid_sweep(n_sites=a.sites, out=f"{a.out}/grid.json")
    coverage_plot(g, f"{a.out}/coverage.png")

if "sweeps" not in skip:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sw = {}
    for ion in ("Er3+", "Ho3+", "Pr3+", "Tm3+", "Dy3+"):
        sw[ion] = phonon_sweep(ion, base="ZBLAN", n_sites=max(15, a.sites // 2))
    json.dump(sw, open(f"{a.out}/phonon_sweep.json", "w"), indent=1, default=float)
    fig, axs = plt.subplots(1, len(sw), figsize=(3.2 * len(sw), 4), sharey=True)
    for ax, (ion, rows) in zip(axs, sw.items()):
        for r in rows:
            for w in r["windows"]:
                ax.plot([w["lo_nm"], w["hi_nm"]], [r["hw_max"]] * 2, lw=5, solid_capstyle="butt")
        ax.set_xscale("log"); ax.set_xlim(450, 5000); ax.set_title(ion.replace("3+", ""))
        ax.set_xticks([500, 1000, 2000, 3000, 5000]); ax.set_xticklabels(["0.5", "1", "2", "3", "5"])
        ax.set_xlabel("wavelength (um)")
    axs[0].set_ylabel("host max phonon energy (cm$^{-1}$)")
    fig.suptitle("Gain windows vs host phonon energy (ZBLAN-like site, all else fixed)", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{a.out}/phonon_sweep.png", dpi=140); plt.close(fig)

    mix = {}
    for ion, fa, fb in (("Er3+", "silica", "tellurite"), ("Er3+", "ZBLAN", "fluoroindate"),
                        ("Tm3+", "silica", "germanate"), ("Pr3+", "ZBLAN", "chalcogenide")):
        mix[f"{ion}:{fa}->{fb}"] = mixture_sweep(ion, fa, fb, n_sites=max(15, a.sites // 2))
    json.dump(mix, open(f"{a.out}/mixture_sweep.json", "w"), indent=1, default=float)

print("total %.0f s" % (time.time() - t0), flush=True)
print("DONE")
