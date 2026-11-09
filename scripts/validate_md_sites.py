"""Score dopant-site geometries from MD against measured Stark structure and spectra.

    python scripts/validate_md_sites.py --sites md_sites/Er-SiO2.json:silica [more path:family ...]
        [--max-sites 120] [--out forward_results/md_sites_validation.json]

One crystal-field strength is set ONCE: the oxygen scale that makes Er 4I15/2 in the pure-silica
MD sites match the measured total splitting (Zemon 1991 silicate 400 cm^-1, Desurvire 395-500).
Everything else is a prediction, for every glass given, compared with the same quantities from the
random-polyhedron sites the forward model uses by default:
  - Yb3+ Stark pattern: 2F7/2 levels (Newell 2007: 0, 425, 530, 865) and 2F5/2 first gap (630)
  - Yb3+ zero-phonon line (975-977 nm), emission tail at 1064/1100 nm relative to peak
  - Er3+ 4I15/2 pattern (Desurvire FLN 0, 30, 55, 129, 234...) and 4I13/2 total splitting
  - Tm3+ 3H6 -> 3F4 absorption peak and FWHM vs composition (Simpson: Ge-silica 1578 nm;
    Al2O3 3.4 % 1632; 9 % 1660), 3H4 band FWHM (18-19 nm), 3H5 peak (1209 nm)
  - Er sigma_e/sigma_a against Giles (log10 ratio error at the Giles wavelengths)
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np

from gainmat.forward.glass import (load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ, run,
                                   FAMILIES, spm_for, load_md_sites)
from gainmat.forward.spectra import random_sites, spm_bkq, fwhm_nm
from gainmat.forward.ions import get_ion

ap = argparse.ArgumentParser()
ap.add_argument("--sites", nargs="+", required=True, help="path.json:family")
ap.add_argument("--max-sites", type=int, default=120)
ap.add_argument("--calib", default="", help="path.json used to set the O scale (default: first given)")
ap.add_argument("--out", default="forward_results/md_sites_validation.json")
a = ap.parse_args()
load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ, verbose=False)

MEAS = dict(
    yb_7_2=[0, 425, 530, 865], yb_5_2_gap=630, yb_zpl_nm=(975, 977),
    er_4I15_low=[0, 30, 55, 129, 234], er_4I15_total=(395, 500), er_4I13_total=(305, 357),
    tm_3F4_abs={"germanosilicate": 1578.5, "Al3.4": 1632.5, "Al9": 1660.5, "aluminosilicate": 1641.5},
    tm_3F4_fwhm=(160, 182), tm_3H4_fwhm=(18, 19), tm_3H5_peak=1209,
)


def subsample(sites, n, seed=0):
    if len(sites) <= n:
        return sites
    idx = np.random.default_rng(seed).choice(len(sites), n, replace=False)
    return [sites[i] for i in sorted(idx)]


def pattern(ion, label, sites, s, lig="O"):
    m = get_ion(ion)
    p = dict(spm_for(ion, lig)); [p.__setitem__(f"B{k}", p[f"B{k}"] * s) for k in (2, 4, 6)]
    out = []
    for st in sites:
        E, W, par = m.stark(spm_bkq(st, {lig: p}))
        e = np.sort(E[par == m.mult_index(label)]); out.append((e - e[0])[::2] if ion != "Tm3+" else (e - e[0]))
    return np.mean(out, axis=0)


def peak(r, pair, which, lo, hi):
    s = r["_spec"]; lam = 1e7 / s["nu"]; y = s["sigma_abs" if which == "abs" else "sigma_em"][pair]
    m = (lam > lo) & (lam < hi)
    return float(lam[m][np.argmax(y[m])]), float(fwhm_nm(s["nu"], np.where(m, y, 0.0)))


def at_rel(r, pair, which, l):
    s = r["_spec"]; y = s["sigma_abs" if which == "abs" else "sigma_em"][pair]
    return float(np.interp(1e7 / l, s["nu"], y) / y.max())


def evaluate(tag, family, sites, s):
    ov = dict(cf_scale=s) if sites is None else dict(sites=sites, cf_scale=s)
    res = dict(tag=tag, family=family, n_sites=None if sites is None else len(sites), o_scale=s)
    geo = sites if sites is not None else random_sites(60, "O", FAMILIES[family]["CN"], FAMILIES[family]["R"],
                                                       FAMILIES[family]["R_std"], FAMILIES[family]["ang_std"], seed=0)
    res["yb_2F7_2"] = pattern("Yb3+", "2F7/2", geo, s).round(0).tolist()
    res["yb_2F5_2"] = pattern("Yb3+", "2F5/2", geo, s).round(0).tolist()
    res["er_4I15_2"] = pattern("Er3+", "4I15/2", geo, s).round(0).tolist()
    res["er_4I13_2"] = pattern("Er3+", "4I13/2", geo, s).round(0).tolist()
    y = run("Yb3+", family, n_sites=40, pumps=[], overrides=ov); p = ("2F5/2", "2F7/2")
    res["yb_zpl_abs_nm"], _ = peak(y, p, "abs", 940, 1000)
    res["yb_em_peak_nm"], _ = peak(y, p, "em", 940, 1150)
    res["yb_em_1064_rel"] = at_rel(y, p, "em", 1064); res["yb_em_1100_rel"] = at_rel(y, p, "em", 1100)
    t = run("Tm3+", family, n_sites=40, pumps=[], overrides=ov)
    res["tm_3F4_abs_peak_nm"], res["tm_3F4_abs_fwhm_nm"] = peak(t, ("3F4", "3H6"), "abs", 1450, 1950)
    res["tm_3F4_em_peak_nm"], _ = peak(t, ("3F4", "3H6"), "em", 1600, 2150)
    res["tm_3H4_abs_peak_nm"], res["tm_3H4_abs_fwhm_nm"] = peak(t, ("3H4", "3H6"), "abs", 740, 850)
    res["tm_3H5_abs_peak_nm"], _ = peak(t, ("3H5", "3H6"), "abs", 1120, 1300)
    e = run("Er3+", family, n_sites=40, pumps=[], overrides=ov); pe = ("4I13/2", "4I15/2")
    res["er_abs_peak_nm"], res["er_abs_fwhm_nm"] = peak(e, pe, "abs", 1450, 1620)
    lam_g = [1500, 1531, 1550, 1564]
    res["er_log_se_sa"] = {l: float(np.log10(at_rel(e, pe, "em", l) * e["_spec"]["sigma_em"][pe].max()
                                              / (at_rel(e, pe, "abs", l) * e["_spec"]["sigma_abs"][pe].max())))
                           for l in lam_g}
    print(json.dumps({k: v for k, v in res.items()}, default=float), flush=True)
    return res


specs = [x.rsplit(":", 1) for x in a.sites]
loaded = {p: subsample(load_md_sites(p), a.max_sites) for p, _ in specs}
# one O scale, set on Er 4I15/2 total splitting in the calibration glass (target 400 cm^-1)
cal = a.calib or specs[0][0]
tot = pattern("Er3+", "4I15/2", loaded[cal], 1.0)[-1]
s_md = 400.0 / tot
print(f"O scale for MD sites: {s_md:.3f} (Er 4I15/2 {tot:.0f} -> 400 cm^-1 in {cal})")
rows = []
for p, fam in specs:
    rows.append(evaluate(f"MD {os.path.basename(p)}", fam, loaded[p], s_md))
fams = sorted({fam for _, fam in specs})
for fam in fams:
    rows.append(evaluate(f"random polyhedra ({fam})", fam, None, 1.0))
GILES = {"germanosilicate": [-0.399, -0.037, 0.079, 0.222], "aluminosilicate": [-0.372, -0.039, 0.042, 0.159]}


def score(r, glass_key):
    """Errors against the measurements that apply to this glass."""
    sc = {}
    yb = np.array(r["yb_2F7_2"]); sc["yb_2F7_2_rms_cm1"] = float(np.sqrt(np.mean((yb - MEAS["yb_7_2"]) ** 2)))
    sc["yb_2F5_2_gap_err_cm1"] = r["yb_2F5_2"][1] - MEAS["yb_5_2_gap"]
    sc["yb_zpl_err_nm"] = r["yb_zpl_abs_nm"] - 976.0
    er = np.array(r["er_4I15_2"][:5]); sc["er_4I15_low5_rms_cm1"] = float(np.sqrt(np.mean((er - MEAS["er_4I15_low"]) ** 2)))
    sc["tm_3H4_fwhm_err_nm"] = r["tm_3H4_abs_fwhm_nm"] - 18.7
    sc["tm_3H5_peak_err_nm"] = r["tm_3H5_abs_peak_nm"] - 1209
    if glass_key in MEAS["tm_3F4_abs"]:
        sc["tm_3F4_abs_peak_err_nm"] = r["tm_3F4_abs_peak_nm"] - MEAS["tm_3F4_abs"][glass_key]
    g = GILES.get(r["family"])
    if g:
        sc["er_giles_ratio_err_dex"] = float(np.mean([r["er_log_se_sa"][l] - v for l, v in zip([1500, 1531, 1550, 1564], g)]))
    return sc


GKEY = {"Er-Al3-SiO2": "Al3.4", "Er-Al9-SiO2": "Al9", "Er-Ge10-SiO2": "germanosilicate"}
for r in rows:
    key = next((v for k, v in GKEY.items() if k in r["tag"]), r["family"])
    r["score"] = score(r, key)
    print(r["tag"], json.dumps({k: round(v, 1) if isinstance(v, float) else v for k, v in r["score"].items()}))
json.dump(dict(measured=MEAS, o_scale_md=s_md, calibration=cal, rows=rows), open(a.out, "w"), indent=1, default=float)
print("wrote", a.out)
