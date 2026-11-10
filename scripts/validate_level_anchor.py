"""Measured-level mode: Stark patterns anchored to measured levels from ONE source, tested on
spectra and devices from OTHER sources that were never used.

  Yb3+  levels: Newell et al., Opt. Commun. 273, 256 (2007) (77 K, silica fibre, composition not
        stated): 2F7/2 0, 425, 530, 865; 2F5/2 gap 630; zero-phonon line 977 nm.
        tests: Pask et al. 1995 germanosilicate emission/absorption curves and two Pask lasers.
  Er3+  levels: Desurvire & Simpson, OL 15, 547 (1990) / Desurvire FLN (4.2 K silicates): 4I15/2
        0, 30, 55, 129, 234 ... total 420; 4I13/2 total 330.
        tests: Giles & Desurvire sigma_e/sigma_a at 1500-1564 nm (Al- and Ge-silica) and the Er
        absorption peak / FWHM.
  Tm3+  levels: Lincoln et al., J. Lumin. 50, 297 (1991) (9 K FLN, Ge-silica): 3H6 ~0, ~250, ~430.
        tests: Simpson thesis 3F4 absorption peak (Ge-silica 1578 nm), 3H4 FWHM, Agger 2006 curves.

    python scripts/validate_level_anchor.py [--sideband S] -> forward_results/level_anchor_validation.json
"""
import argparse, json, os, sys, time
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np
from gainmat.forward import glass as g
from gainmat.forward.spectra import fwhm_nm
from gainmat.forward.fiber import local_model_for, slope_and_threshold
from gainmat.refs.zblan_yb_refs import PASK1995
from gainmat.refs.stark_anchor_refs import ANCHORS

ap = argparse.ArgumentParser(); ap.add_argument("--sideband", type=float, default=0.0)
ap.add_argument("--out", default="forward_results/level_anchor_validation.json"); a = ap.parse_args()
g.load_refs(); g.apply_site_refs(); g.calibrate_spm(370.0, g.CARNALL_ER_LAF3_BKQ, verbose=False)
res = dict(sideband_S=a.sideband, yb=[], er=[], tm=[], devices=[])


def ov_for(ion, anchored):
    if not anchored:
        return dict(level_anchor=None)          # anchors are on by default in glass.run; switch them off
    o = dict(level_anchor={k: v for k, v in ANCHORS[ion]["levels"].items()})
    if a.sideband:
        o["sideband_S"] = a.sideband
    return o


# ---------------- Yb vs Pask curves
pk = PASK1995["cross_sections"]["curve_DIGITIZED"]; lam_m = np.array(pk["wavelength_nm"], float)
cur = {k: np.array([np.nan if v is None else v for v in pk[f"sigma_{k}_1e-24m2"]]) * 1e-20 for k in ("abs", "em")}
for anchored in (False, True):
    sp = g.run("Yb3+", "germanosilicate", pumps=[], n_sites=60, overrides=ov_for("Yb3+", anchored))["_spec"]
    p = ("2F5/2", "2F7/2"); r = dict(variant="anchored" if anchored else "model")
    for k in ("abs", "em"):
        y = sp["sigma_" + k][p]; mod = np.interp(1e7 / lam_m, sp["nu"], y)
        ok = np.isfinite(cur[k]) & (cur[k] >= 0.05e-20)
        r[f"{k}_peak_nm"] = float(1e7 / sp["nu"][np.argmax(y)]); r[f"{k}_peak"] = float(y.max())
        r[f"{k}_logerr_dex_median"] = float(np.median(np.abs(np.log10(np.maximum(mod[ok], 1e-24) / cur[k][ok]))))
        r[f"{k}_rel"] = {int(l): [round(float(np.interp(1e7 / l, sp["nu"], y) / y.max()), 3),
                                  round(float(np.interp(l, lam_m, cur[k]) / np.nanmax(cur[k])), 3)]
                         for l in ((912, 960, 1000, 1035) if k == "abs" else (1000, 1035, 1064, 1100, 1150))}
    res["yb"].append(r); print("Yb", json.dumps(r), flush=True)

# ---------------- Er vs Giles ratio
GILES = {"germanosilicate": [-0.399, -0.037, 0.079, 0.222], "aluminosilicate": [-0.372, -0.039, 0.042, 0.159]}
LG = [1500, 1531, 1550, 1564]
for fam in ("aluminosilicate", "germanosilicate"):
    for anchored in (False, True):
        sp = g.run("Er3+", fam, pumps=[], n_sites=60, overrides=ov_for("Er3+", anchored))["_spec"]
        p = ("4I13/2", "4I15/2"); sa, se = sp["sigma_abs"][p], sp["sigma_em"][p]
        lr = [float(np.log10(np.interp(1e7 / l, sp["nu"], se) / np.interp(1e7 / l, sp["nu"], sa))) for l in LG]
        r = dict(family=fam, variant="anchored" if anchored else "model", log_ratio=lr, giles=GILES[fam],
                 mean_err_dex=float(np.mean(np.subtract(lr, GILES[fam]))), rms_err_dex=float(np.sqrt(np.mean(np.subtract(lr, GILES[fam]) ** 2))),
                 abs_peak_nm=float(1e7 / sp["nu"][np.argmax(sa)]), abs_fwhm_nm=float(fwhm_nm(sp["nu"], sa)))
        res["er"].append(r); print("Er", json.dumps(r), flush=True)

# ---------------- Tm
for fam, meas in (("germanosilicate", 1578.5), ("aluminosilicate", 1632.5)):
    for anchored in (False, True):
        sp = g.run("Tm3+", fam, pumps=[], n_sites=60, overrides=ov_for("Tm3+", anchored))["_spec"]
        lam = 1e7 / sp["nu"]
        def band(pair, which, lo, hi):
            y = sp["sigma_abs" if which == "abs" else "sigma_em"][pair]; m = (lam > lo) & (lam < hi)
            return float(lam[m][np.argmax(y[m])]), float(fwhm_nm(sp["nu"], np.where(m, y, 0)))
        r = dict(family=fam, variant="anchored" if anchored else "model", abs_3F4=band(("3F4", "3H6"), "abs", 1450, 1950),
                 em_3F4=band(("3F4", "3H6"), "em", 1600, 2150), abs_3H4=band(("3H4", "3H6"), "abs", 740, 850),
                 abs_3H5=band(("3H5", "3H6"), "abs", 1100, 1300), measured_3F4_abs=meas)
        res["tm"].append(r); print("Tm", json.dumps(r), flush=True)

# ---------------- Pask lasers (never used to set anything)
N = 550e-6 * 2.2e22
for anchored in (False, True):
    lm = local_model_for("Yb3+", "germanosilicate", n_sites=30, overrides=ov_for("Yb3+", anchored))
    t0 = time.time()
    s, th, cv = slope_and_threshold(lm, 90, 3.75, 0.17, N, 1047, [0.05, 0.1, 0.2, 0.3], 1102,
                                    R1=0.98, R2=0.04, alpha_bg_dB_m=0.005, nz=60)
    res["devices"].append(dict(case="Yb 1102 nm laser, 1047 pump (Pask 1995)", variant="anchored" if anchored else "model",
                               slope=s, threshold_W=th, measured_slope=0.90, measured_threshold_W=0.030, seconds=time.time() - t0))
    print("dev", json.dumps(res["devices"][-1], default=float), flush=True)
    s, th, cv = slope_and_threshold(lm, 5, 3.75, 0.17, N, 840, [0.08, 0.15, 0.25, 0.35], 1040,
                                    R1=0.04, R2=0.04, alpha_bg_dB_m=0.005, nz=40)
    res["devices"].append(dict(case="Yb free-running 5 m, 840 pump (Pask 1995)", variant="anchored" if anchored else "model",
                               slope=s, threshold_W=th, measured_threshold_W=0.045))
    print("dev", json.dumps(res["devices"][-1], default=float), flush=True)
json.dump(res, open(a.out, "w"), indent=1, default=float)
