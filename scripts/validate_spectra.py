"""Full cross-section spectra: model vs digitised measured curves (shape and scale).

Metrics per curve: scale = model/measured integrated cross section over the measured range;
shape = rms difference of peak-normalised curves; peak shift (nm); width ratio (FWHM).

    python scripts/validate_spectra.py  ->  forward_results/spectra_validation.json, .png
"""
import json
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np

from gainmat.forward.glass import load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ, run
from gainmat.refs.tm_silica_refs import AGGER2006
from gainmat.refs.midir_ho_refs import LEGOUET2020, JACKSON2012
from gainmat.refs.zblan_yb_refs import PASK1995

load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ, verbose=False)

curves = []   # (label, ion, family, (upper, lower), 'abs'|'em', lam_nm array, sigma cm^2 array, src)
ag = AGGER2006["cross_sections"]["spectra_digitized_1e-25_m2"]
for k, kind in (("Tm1_sigma_a", "abs"), ("Tm1_sigma_e", "em")):
    lam = np.array(ag["lambda_nm"], float); s = np.array([np.nan if v is None else v for v in ag[k]], float)
    curves.append((f"Tm:Al-silica {kind} (Agger 2006)", "Tm3+", "aluminosilicate", ("3F4", "3H6"), kind,
                   lam, s * 1e-21, "Agger & Povlsen, OE 14, 50 (2006) Fig. 4"))
lg = np.array(LEGOUET2020["XS_CURVE_DIGITIZED_nm_abs_em_1e-25m2"], float)
curves.append(("Ho:Al-silica abs (Le Gouet 2020)", "Ho3+", "aluminosilicate", ("5I7", "5I8"), "abs",
               lg[:, 0], lg[:, 1] * 1e-21, "Le Gouet et al., arXiv:2003.14188"))
curves.append(("Ho:Al-silica em (Le Gouet 2020)", "Ho3+", "aluminosilicate", ("5I7", "5I8"), "em",
               lg[:, 0], lg[:, 2] * 1e-21, "Le Gouet et al., arXiv:2003.14188"))
pk = PASK1995["cross_sections"]["curve_DIGITIZED"]
lam = np.array(pk["wavelength_nm"], float)
for k, kind in (("sigma_abs_1e-24m2", "abs"), ("sigma_em_1e-24m2", "em")):
    s = np.array([np.nan if v is None else v for v in pk[k]], float)
    curves.append((f"Yb:Ge-silica {kind} (Pask 1995)", "Yb3+", "germanosilicate", ("2F5/2", "2F7/2"), kind,
                   lam, s * 1e-20, "Pask et al., JSTQE 1, 2 (1995) Fig. 1(b)"))
J = JACKSON2012["FIG4_EMISSION_XS_CURVES_um_1e-25m2"]
jm = {"Ho_5I7-5I8_aluminosilicate": ("Ho3+", "aluminosilicate", ("5I7", "5I8")),
      "Tm_3F4-3H6_aluminosilicate": ("Tm3+", "aluminosilicate", ("3F4", "3H6")),
      "Er_4I13/2-4I15/2_aluminosilicate": ("Er3+", "aluminosilicate", ("4I13/2", "4I15/2")),
      "Er_4I11/2-4I13/2_ZBLAN": ("Er3+", "ZBLAN", ("4I11/2", "4I13/2")),
      "Ho_5I6-5I7_ZBLAN": ("Ho3+", "ZBLAN", ("5I6", "5I7")),
      "Dy_6H13/2-6H15/2_ZBLAN": ("Dy3+", "ZBLAN", ("6H13/2", "6H15/2"))}
for key, (ion, fam, pair) in jm.items():
    d = np.array(J[key], float)
    curves.append((f"{ion[:-2]}:{fam[:5]} {pair[0]} em (Jackson 2012)", ion, fam, pair, "em",
                   d[:, 0] * 1e3, d[:, 1] * 1e-21, "Jackson, Nat. Photon. 6, 423 (2012) Fig. 4"))

cache = {}
rows = []
for label, ion, fam, pair, kind, lam, s_meas, src in curves:
    if (ion, fam) not in cache:
        cache[(ion, fam)] = run(ion, fam, n_sites=30, pumps=[])["_spec"]
    sp = cache[(ion, fam)]
    s_mod = np.interp(1e7 / lam, sp["nu"], sp["sigma_abs" if kind == "abs" else "sigma_em"][pair])
    ok = np.isfinite(s_meas)
    l, a, b = lam[ok], s_meas[ok], s_mod[ok]
    order = np.argsort(l); l, a, b = l[order], a[order], b[order]
    scale = np.trapezoid(b, l) / np.trapezoid(a, l)
    shape = float(np.sqrt(np.mean((b / b.max() - a / a.max()) ** 2)))

    def fwhm(y):
        h = y.max() / 2; idx = np.where(y >= h)[0]
        return l[idx.max()] - l[idx.min()]
    rows.append(dict(label=label, ion=ion, family=fam, pair=pair, kind=kind, src=src,
                     scale_integrated=float(scale), peak_ratio=float(b.max() / a.max()),
                     shape_rms=shape, peak_meas_nm=float(l[a.argmax()]), peak_model_nm=float(l[b.argmax()]),
                     fwhm_meas_nm=float(fwhm(a)), fwhm_model_nm=float(fwhm(b)),
                     lam_nm=l.tolist(), meas=a.tolist(), model=b.tolist()))
    print(f"{label:42s} scale {scale:5.2f}  peak x{b.max()/a.max():4.2f}  shape rms {shape:4.2f}  "
          f"peak {l[b.argmax()]:6.0f}/{l[a.argmax()]:6.0f} nm  FWHM {fwhm(b):5.0f}/{fwhm(a):5.0f} nm")
json.dump(rows, open("forward_results/spectra_validation.json", "w"), indent=1, default=float)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
n = len(rows); nc = 4; nr = int(np.ceil(n / nc))
fig, axs = plt.subplots(nr, nc, figsize=(4 * nc, 2.8 * nr))
for ax, r in zip(axs.flat, rows):
    ax.plot(r["lam_nm"], np.array(r["meas"]) * 1e21, "k.-", lw=1, ms=3, label="measured")
    ax.plot(r["lam_nm"], np.array(r["model"]) * 1e21, color="tab:red", lw=1.5, label="model")
    ax.set_title(r["label"], fontsize=8); ax.tick_params(labelsize=7)
    ax.set_xlabel("nm", fontsize=7); ax.set_ylabel("$\\sigma$ (1e-21 cm$^2$)", fontsize=7)
for ax in list(axs.flat)[n:]:
    ax.axis("off")
axs.flat[0].legend(fontsize=7)
fig.suptitle("Cross-section spectra: model (red) vs digitised measurements (black)", fontsize=10)
fig.tight_layout(); fig.savefig("forward_results/spectra_validation.png", dpi=130)
print("wrote forward_results/spectra_validation.json/.png")
