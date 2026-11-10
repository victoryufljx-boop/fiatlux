"""Crystal-field intrinsic parameters on MD glass sites, fitted jointly across Er, Yb and Tm with
leave-one-ion-out tests.

Sites: Pedone-potential melt-quench MD of Er-doped silica glasses (branches md2/*, Er sites used
for every ion).  Parameters (oxygen ligands, superposition model): overall scale s and rank
multipliers r2, r6 on the handbook B2, B6 (r4 = 1).  Targets:
  Er (pure-silica sites): 4I15/2 levels 0, 30, 55, 129, 234 (Desurvire FLN) and total 420;
                          4I13/2 total 330 (silicates 305-357)
  Yb (pure-silica sites): 2F7/2 levels 0, 425, 530, 865 and 2F5/2 first gap 630 (Newell 2007)
  Tm (3.4 % Al2O3 sites): 3H6 -> 3F4 absorption peak 1632.5 nm and FWHM 181.5 nm (Simpson),
                          3F4 emission peak 1825 nm (Agger 2006, Al-silica)
Held-out checks for every fold: the held-out ion's targets, plus (never fitted) the Tm absorption
peak on the Ge-silica sites (Simpson 1578 nm) and the Er sigma_e/sigma_a against Giles.

    python -m gainmat.ml.md_cf_fit --sites-dir DIR [--holdout Yb3+] --out ml_results/md_cf
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

T_ER = dict(low=[0, 30, 55, 129, 234], total=420, t13=330)
T_YB = dict(levels=[0, 425, 530, 865], gap52=630)
T_TM = dict(abs=1632.5, fwhm=181.5, em=1825.0)


def pattern(ion, label, sites, s, rk, kramers=True):
    from ..forward.glass import spm_for
    from ..forward.spectra import spm_bkq
    from ..forward.ions import get_ion
    m = get_ion(ion); p = dict(spm_for(ion, "O"))
    for i, k in enumerate((2, 4, 6)):
        p[f"B{k}"] = p[f"B{k}"] * s * rk[i]
    out = []
    for st in sites:
        E, W, par = m.stark(spm_bkq(st, {"O": p}))
        e = np.sort(E[par == m.mult_index(label)]); e = e - e[0]
        out.append(e[::2] if kramers else e)
    return np.mean(out, axis=0)


def tm_bands(sites, family, s, rk):
    from ..forward import glass as g
    from ..forward.spectra import fwhm_nm
    old = g.RANK_SCALE.get("O")
    g.RANK_SCALE["O"] = tuple(rk)
    try:
        r = g.run("Tm3+", family, pumps=[], overrides=dict(sites=sites, cf_scale=s))
    finally:
        g.RANK_SCALE["O"] = old if old else (1.0, 1.0, 1.0)
    sp = r["_spec"]; lam = 1e7 / sp["nu"]
    out = {}
    for which, lo, hi in (("abs", 1450, 1950), ("em", 1600, 2150)):
        y = sp["sigma_abs" if which == "abs" else "sigma_em"][("3F4", "3H6")]; m = (lam > lo) & (lam < hi)
        out[which] = float(lam[m][np.argmax(y[m])])
        if which == "abs":
            out["fwhm"] = float(fwhm_nm(sp["nu"], np.where(m, y, 0.0)))
    return out


def residuals(params, S, ions):
    s, r2, r6 = params
    rk = (r2, 1.0, r6)
    res = {}
    if "Er3+" in ions:
        e15 = pattern("Er3+", "4I15/2", S["silica"], s, rk); e13 = pattern("Er3+", "4I13/2", S["silica"], s, rk)
        res["Er3+"] = list((e15[:5] - T_ER["low"]) / 30.0) + [(e15[-1] - T_ER["total"]) / 30.0, (e13[-1] - T_ER["t13"]) / 30.0]
    if "Yb3+" in ions:
        y7 = pattern("Yb3+", "2F7/2", S["silica"], s, rk); y5 = pattern("Yb3+", "2F5/2", S["silica"], s, rk)
        res["Yb3+"] = list((y7 - T_YB["levels"]) / 40.0) + [(y5[1] - T_YB["gap52"]) / 50.0]
    if "Tm3+" in ions:
        t = tm_bands(S["al3"], "aluminosilicate", s, rk)
        res["Tm3+"] = [(t["abs"] - T_TM["abs"]) / 10.0, (t["fwhm"] - T_TM["fwhm"]) / 15.0, (t["em"] - T_TM["em"]) / 15.0]
    return res


def loss(params, S, ions):
    r = residuals(params, S, ions)
    return float(sum(np.sum(np.square(v)) for v in r.values()))


def fit(S, ions, maxiter=70, log=print):
    from scipy.optimize import minimize
    hist = []

    def f(x):
        p = np.exp(x)
        L = loss(p, S, ions) + 0.5 * float(np.sum((x[1:] / 1.0) ** 2))
        hist.append((L, p)); log(f"  {len(hist):3d} loss {L:8.2f} s {p[0]:.3f} r2 {p[1]:.3f} r6 {p[2]:.3f}")
        return L
    x0 = np.log([1.15, 1.0, 1.0])
    minimize(f, x0, method="Nelder-Mead", options=dict(maxfev=maxiter, xatol=0.01, fatol=0.05,
             initial_simplex=np.array([x0, x0 + [0.3, 0, 0], x0 + [0, 0.6, 0], x0 + [0, 0, 0.6]])))
    return min(hist, key=lambda h: h[0])[1]


def report(params, S):
    s, r2, r6 = params; rk = (r2, 1.0, r6)
    out = dict(params=dict(s=s, r2=r2, r6=r6))
    out["Er_4I15_2"] = pattern("Er3+", "4I15/2", S["silica"], s, rk).round(0).tolist()
    out["Er_4I13_2"] = pattern("Er3+", "4I13/2", S["silica"], s, rk).round(0).tolist()
    out["Yb_2F7_2"] = pattern("Yb3+", "2F7/2", S["silica"], s, rk).round(0).tolist()
    out["Yb_2F5_2"] = pattern("Yb3+", "2F5/2", S["silica"], s, rk).round(0).tolist()
    out["Tm_Al3"] = tm_bands(S["al3"], "aluminosilicate", s, rk)
    if "ge10" in S:
        out["Tm_Ge10"] = tm_bands(S["ge10"], "germanosilicate", s, rk)
    if "al9" in S:
        out["Tm_Al9"] = tm_bands(S["al9"], "aluminosilicate", s, rk)
    out["residuals"] = {k: [round(float(x), 2) for x in v] for k, v in residuals(params, S, ("Er3+", "Yb3+", "Tm3+")).items()}
    out["chi2_by_ion"] = {k: float(np.sum(np.square(v))) for k, v in residuals(params, S, ("Er3+", "Yb3+", "Tm3+")).items()}
    return out


def main():
    from ..forward.glass import load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ, load_md_sites
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites-dir", required=True)
    ap.add_argument("--holdout", default="")
    ap.add_argument("--n", type=int, default=60)
    ap.add_argument("--maxiter", type=int, default=70)
    ap.add_argument("--offset", type=int, default=0, help="start index of the strided site subsample")
    ap.add_argument("--out", default="ml_results/md_cf")
    a = ap.parse_args()
    load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ, verbose=False)
    S = {}
    for key, name in (("silica", "Er-SiO2"), ("al3", "Er-Al3-SiO2"), ("ge10", "Er-Ge10-SiO2"), ("al9", "Er-Al9-SiO2")):
        p = os.path.join(a.sites_dir, name + ".json")
        if os.path.exists(p):
            st = load_md_sites(p)
            S[key] = st[a.offset % max(1, len(st) // a.n)::max(1, len(st) // a.n)][:a.n]
    ions = tuple(i for i in ("Er3+", "Yb3+", "Tm3+") if i != a.holdout)
    os.makedirs(a.out, exist_ok=True)
    tag = ("all" if not a.holdout else "minus_" + a.holdout.replace("+", "")) + (f"_off{a.offset}" if a.offset else "")
    logf = open(os.path.join(a.out, tag + ".log"), "w")
    log = lambda s: (print(s, flush=True), logf.write(s + "\n"), logf.flush())
    t0 = time.time()
    base = report((1.154, 1.0, 1.0), S)
    best = fit(S, ions, a.maxiter, log)
    rep = report(tuple(best), S)
    out = dict(holdout=a.holdout or None, fitted_on=list(ions), handbook=base, fitted=rep, seconds=time.time() - t0)
    json.dump(out, open(os.path.join(a.out, tag + ".json"), "w"), indent=1, default=float)
    log(f"{tag}: params {rep['params']} chi2 by ion: handbook {base['chi2_by_ion']} -> fitted {rep['chi2_by_ion']}")
    log(f"   Tm abs: Al3 {rep['Tm_Al3']}  Ge10 {rep.get('Tm_Ge10')}  Al9 {rep.get('Tm_Al9')}")


if __name__ == "__main__":
    main()
