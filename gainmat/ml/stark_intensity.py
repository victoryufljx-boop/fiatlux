"""Stark-resolved intensities (odd-rank crystal field) learned across ions in silica glasses.

The isotropic Judd-Ofelt form gives every Stark component the same Omega-weighted |U|^2 share;
the explicit form (spectra.aniso_S) computes each component's strength from the odd-rank
crystal field of the same site, renormalised to the measured Omega (band integrals unchanged).
Its free parameters are the ratios of the two intrinsic intensity amplitudes per lambda,
a(lambda, lambda+1) / a(lambda, lambda-1), for lambda = 2, 4, 6: properties of the ligand, so
shared by every ion.  Fit them on the line shapes of Er, Yb and Ho; predict Tm (never used):
line shapes (Agger 2006, Jackson 2012) and peak positions (Jackson & King 1999).

    python -m gainmat.ml.stark_intensity [--holdout Tm3+] [--maxiter 60] --out ml_results/stark
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

FAMS = ("silica", "aluminosilicate", "germanosilicate")
KINDS = ("curve", "peak_nm", "ratio_log10")


def observations():
    from .observations import build
    from ..refs.tm_silica_refs import JACKSON1999 as J
    obs = [o for o in build() if o["family"] in FAMS and o["kind"] in KINDS]
    d = next((v for v in J.values() if isinstance(v, dict) and "sigma_a_3H6_3F4_peak_nm" in v), J)
    src = "Jackson & King, JLT 17, 948 (1999) (digitised peaks)"
    for pair, which, key, band, sig in [
            (("3F4", "3H6"), "abs", "sigma_a_3H6_3F4_peak_nm", (1450, 1900), 20.0),
            (("3H5", "3H6"), "abs", "sigma_a_3H6_3H5_peak_nm", (1100, 1300), 10.0),
            (("3H4", "3H6"), "abs", "sigma_a_3H6_3H4_peak_nm", (740, 850), 6.0),
            (("3F4", "3H6"), "em", "sigma_e_3F4_3H6_peak_nm", (1600, 2150), 20.0)]:
        if key in d:
            obs.append(dict(family="silica", ion="Tm3+", kind="peak_nm", pair=pair, which=which, band_nm=band,
                            value=float(d[key]), sigma_abs=sig, src=src))
    return obs


def predict(o, r):
    from .host_calib import predict as p0, _band
    if o["kind"] != "peak_nm":
        return p0(o, r)
    lam, s = _band(r["_spec"], tuple(o["pair"]), o["which"])
    lo, hi = o["band_nm"]; m = (lam >= lo) & (lam <= hi)
    v = float(lam[m][np.argmax(s[m])])
    return np.array([(v - o["value"]) / o["sigma_abs"]]), dict(model=v)


def evaluate(obs, intensity, n_sites=16, seed=0, detail=False):
    from ..forward.glass import run
    res, info = [], []
    for fam, ion in sorted({(o["family"], o["ion"]) for o in obs}):
        ov = {} if intensity is None else dict(intensity=intensity)
        r = run(ion, fam, n_sites=n_sites, pumps=[], overrides=ov, seed=seed, dnu=4.0)
        r["_spec"]["ion"] = ion
        for o in obs:
            if (o["family"], o["ion"]) != (fam, ion):
                continue
            v, d = predict(o, r)
            res.append(v)
            if detail:
                info.append(dict(family=fam, ion=ion, kind=o["kind"], pair=o.get("pair"), which=o.get("which"),
                                 resid=[float(x) for x in v], **d, src=o["src"]))
    r = np.concatenate(res)
    return (float(np.sum(r ** 2)), info) if detail else float(np.sum(r ** 2))


def components(info):
    acc = {}
    for i in info:
        names = ["shape", "width", "peak_shift"] if i["kind"] == "curve" else [i["kind"]]
        for n, v in zip(names, i["resid"]):
            acc[n] = acc.get(n, 0.0) + v * v
    return {k: round(v, 1) for k, v in acc.items()}


def main():
    from scipy.optimize import minimize
    from ..forward.glass import load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ
    ap = argparse.ArgumentParser()
    ap.add_argument("--holdout", default="Tm3+")
    ap.add_argument("--sites", type=int, default=16)
    ap.add_argument("--maxiter", type=int, default=60)
    ap.add_argument("--out", default="ml_results/stark")
    a = ap.parse_args()
    load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ, verbose=False)
    os.makedirs(a.out, exist_ok=True)
    obs = observations()
    tr = [o for o in obs if o["ion"] != a.holdout]; te = [o for o in obs if o["ion"] == a.holdout]
    tag = "all" if not a.holdout else "minus_" + a.holdout.replace("+", "")
    logf = open(os.path.join(a.out, tag + ".log"), "w")
    log = lambda s: (print(s, flush=True), logf.write(s + "\n"), logf.flush())
    t0 = time.time()
    base = {}
    for name, it in (("isotropic JO", None), ("odd-rank, equal amplitudes", dict(ratio={2: 1.0, 4: 1.0, 6: 1.0}))):
        Ltr, itr = evaluate(tr, it, a.sites, detail=True)
        base[name] = dict(train=Ltr, train_components=components(itr))
        if te:
            Lte, ite = evaluate(te, it, a.sites, detail=True)
            base[name].update(test=Lte, test_components=components(ite), test_detail=ite)
        log(f"{name}: train {Ltr:.1f} {components(itr)}" + (f" | held-out {a.holdout} {Lte:.1f} {components(ite)}" if te else ""))
    hist = []

    def f(x):
        it = dict(ratio={2: float(np.exp(x[0])), 4: float(np.exp(x[1])), 6: float(np.exp(x[2]))})
        L = evaluate(tr, it, a.sites) + 0.5 * float(np.sum((x / 1.5) ** 2))
        hist.append((L, it)); log(f"  eval {len(hist):3d} loss {L:8.2f} {json.dumps(it['ratio'])}")
        return L
    minimize(f, np.zeros(3), method="Nelder-Mead", options=dict(maxfev=a.maxiter, xatol=0.05, fatol=0.2,
             initial_simplex=np.array([[0, 0, 0], [1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]])))
    L, best = min(hist, key=lambda h: h[0])
    Ltr, itr = evaluate(tr, best, a.sites, detail=True)
    out = dict(holdout=a.holdout or None, baselines=base, learned=best, train=Ltr, train_components=components(itr),
               evals=len(hist), seconds=time.time() - t0)
    if te:
        Lte, ite = evaluate(te, best, a.sites, detail=True)
        out.update(test=Lte, test_components=components(ite), test_detail=ite)
        log(f"learned {best['ratio']}: train {Ltr:.1f} | held-out {a.holdout} {Lte:.1f} {components(ite)}")
        for i in ite:
            log(f"    {i['kind']:9s} {i['pair']} {i.get('which')}: {i.get('model') or i.get('peak_shift_nm')}")
    json.dump(out, open(os.path.join(a.out, tag + ".json"), "w"), indent=1, default=float)
    log(f"done {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
