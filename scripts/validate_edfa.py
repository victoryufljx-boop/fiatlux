"""Erbium-doped fibre: measured absorption alpha(lambda) and gain coefficient g*(lambda)
(Giles & Desurvire, JLT 9, 271 (1991), Fig. 2) vs the model.

alpha = Gamma N sigma_a and g* = Gamma N sigma_e, so g*/alpha = sigma_e/sigma_a needs no
concentration or overlap: a direct test of the model's Stark structure / McCumber ratio.
Shapes are compared peak-normalised.  The Giles gain-vs-pump curves (Figs. 4-8) are
calculations, not measurements, and are not used.

    python scripts/validate_edfa.py  ->  forward_results/edfa_validation.json, edfa_validation.png
"""
import json
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np

from gainmat.forward.glass import load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ, run
from gainmat.refs.edfa_refs import GILES1991 as G

load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ, verbose=False)
sp2 = G["spectra_fig2"]
lam = np.array(sp2["wavelength_nm"], float)
rows, plots = [], []
for key, fam, lab in (("a_Ge_silicate", "germanosilicate", "Ge-silicate (Fig. 2a)"),
                      ("b_AlGe_silicate", "aluminosilicate", "Al:Ge-silicate (Fig. 2b)")):
    if key not in sp2:
        key = [k for k in sp2 if k.startswith(key[0] + "_")][0]
    a = np.array(sp2[key]["alpha_dB_per_m"], float); g = np.array(sp2[key]["gstar_dB_per_m"], float)
    s = run("Er3+", fam, n_sites=30, pumps=[])["_spec"]
    x = 1e7 / lam
    sa = np.interp(x, s["nu"], s["sigma_abs"][("4I13/2", "4I15/2")])
    se = np.interp(x, s["nu"], s["sigma_em"][("4I13/2", "4I15/2")])
    ratio_m, ratio_e = se / sa, g / a
    sh_a = float(np.sqrt(np.mean((sa / sa.max() - a / a.max()) ** 2)))
    sh_g = float(np.sqrt(np.mean((se / se.max() - g / g.max()) ** 2)))
    lr = float(np.sqrt(np.mean(np.log10(ratio_m / ratio_e) ** 2)))
    i1550 = np.argmin(abs(lam - 1550)); i1530 = np.argmin(abs(lam - 1531))
    rows += [
        dict(case=f"{lab}: alpha shape rms (peak-normalised)", model=round(sh_a, 3), measured=0.0),
        dict(case=f"{lab}: g* shape rms (peak-normalised)", model=round(sh_g, 3), measured=0.0),
        dict(case=f"{lab}: sigma_e/sigma_a at 1531 nm", model=round(float(ratio_m[i1530]), 2),
             measured=round(float(ratio_e[i1530]), 2)),
        dict(case=f"{lab}: sigma_e/sigma_a at 1550 nm", model=round(float(ratio_m[i1550]), 2),
             measured=round(float(ratio_e[i1550]), 2)),
        dict(case=f"{lab}: sigma_e/sigma_a rms log10 over 1471-1569 nm", model=round(lr, 3), measured=0.0),
        dict(case=f"{lab}: absorption peak nm", model=float(lam[sa.argmax()]), measured=float(lam[a.argmax()])),
    ]
    plots.append((lab, a / a.max(), g / g.max(), sa / sa.max(), se / se.max(), ratio_e, ratio_m))
for r in rows:
    print(f"{r['case']:60s} model {r['model']}   measured {r['measured']}")
json.dump(rows, open("forward_results/edfa_validation.json", "w"), indent=1, default=float)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, axs = plt.subplots(2, 2, figsize=(10, 6.5))
for j, (lab, a, g, sa, se, re, rm) in enumerate(plots):
    ax = axs[0, j]
    ax.plot(lam, a, "k-", label="alpha (meas.)"); ax.plot(lam, g, "k--", label="g* (meas.)")
    ax.plot(lam, sa, color="tab:red", label="sigma_a (model)"); ax.plot(lam, se, "--", color="tab:red", label="sigma_e (model)")
    ax.set_title(f"Er {lab}, peak-normalised", fontsize=9); ax.legend(fontsize=7)
    ax = axs[1, j]
    ax.semilogy(lam, re, "k-", label="g*/alpha (meas.)"); ax.semilogy(lam, rm, color="tab:red", label="sigma_e/sigma_a (model)")
    ax.set_xlabel("wavelength (nm)"); ax.legend(fontsize=7); ax.set_title("emission / absorption ratio", fontsize=9)
fig.suptitle("Er-doped silica fibre vs Giles & Desurvire (1991) Fig. 2 (digitised)", fontsize=10)
fig.tight_layout(); fig.savefig("forward_results/edfa_validation.png", dpi=130)
print("wrote forward_results/edfa_validation.json/.png")
