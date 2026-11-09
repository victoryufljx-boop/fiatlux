"""Uncertainty propagation: gain windows as probabilities instead of yes/no.

Samples the forward model's learned/uncertain inputs and reruns it:
  Omega_2,4,6  log-normal around the value used.  Spread: if measured for that glass, the
               literature reproducibility between publications (0.25 / 0.10 / 0.13 dex, from
               LOMS replicates); if predicted, the out-of-sample error of the L1 model
               (0.37 / 0.30 / 0.27 dex, publication-grouped CV).
  host theta   (cf_scale, disorder, neph_std, gamma_h) drawn from the spread of the L2 fits
               for that family when available (ml_results/l2/<family>_*.json), else fixed.
Returns, on a wavelength grid, the probability that gain - loss >= threshold (dB/m) under
the family's standard pumps.

    python -m gainmat.ml.uq --ion Er3+ --family ZBLAN --n 16
"""
from __future__ import annotations

import argparse
import glob
import json
import os

import numpy as np

SIG_MEASURED = np.array([0.25, 0.10, 0.13])
SIG_PREDICTED = np.array([0.37, 0.30, 0.27])


def theta_samples(family, l2_dir="ml_results/l2"):
    out = []
    for p in glob.glob(os.path.join(l2_dir, f"{family}_*.json")):
        try:
            d = json.load(open(p))
            out.append(d["theta_fit"])
        except Exception:
            pass
    return out


def window_probability(ion, family, n=16, n_sites=20, seed=0, threshold_dB_m=0.1, use_theta=True,
                       lam=(400.0, 5000.0), dnu=4.0):
    from ..forward.glass import run, omegas_for, FAMILIES
    om, src = omegas_for(ion, family)
    sig = SIG_MEASURED if src == "measured" else SIG_PREDICTED
    ths = theta_samples(family) if use_theta else []
    f = FAMILIES[family]
    rng = np.random.default_rng(seed)
    hits, nu = None, None
    for s in range(n):
        o = tuple(float(x) for x in np.array(om) * 10 ** (rng.normal(0, sig)))
        ov = dict(omegas=o)
        if ths:
            th = ths[rng.integers(len(ths))]
            ov.update(cf_scale=th["cf_scale"], R_std=f["R_std"] * th["disorder"], ang_std=f["ang_std"] * th["disorder"],
                      neph_std=th["neph_std"], gamma_h=th["gamma_h"])
        r = run(ion, family, n_sites=n_sites, overrides=ov, seed=seed + s, lam=lam, dnu=dnu)
        if r["_env"] is None:
            continue
        ok = (r["_env"] - r["_loss"]) >= threshold_dB_m
        hits = ok.astype(float) if hits is None else hits + ok
        nu = r["_spec"]["nu"]
    if hits is None:
        return None
    P = hits / n
    lamg = 1e7 / nu
    return dict(ion=ion, family=family, omega=om, omega_src=src, sigma_log10=sig.tolist(), n=n,
                n_theta=len(ths), lam_nm=lamg[::-1].tolist(), P=P[::-1].tolist())


def bands(res, p_lo=0.2, p_hi=0.8, min_width_nm=2.0):
    """Robust (P >= p_hi) and possible (P >= p_lo) wavelength ranges."""
    lam = np.array(res["lam_nm"]); P = np.array(res["P"])
    out = {}
    for name, thr in (("robust", p_hi), ("possible", p_lo)):
        ok = P >= thr; rngs = []; i = 0
        while i < len(ok):
            if ok[i]:
                j = i
                while j + 1 < len(ok) and ok[j + 1]:
                    j += 1
                if lam[j] - lam[i] >= min_width_nm:
                    rngs.append((round(float(lam[i]), 1), round(float(lam[j]), 1)))
                i = j + 1
            else:
                i += 1
        out[name] = rngs
    return out


def main():
    from ..forward.glass import load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", default="")
    ap.add_argument("--family", default="")
    ap.add_argument("--n", type=int, default=16)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--families", default="", help="comma list restricting --all")
    ap.add_argument("--no-ml", action="store_true", help="explicit model only (no learned theta / L3 quenching)")
    ap.add_argument("--out", default="ml_results/uq")
    a = ap.parse_args()
    load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ, verbose=False)
    from ..forward import glass as _g
    _g.USE_ML = not a.no_ml
    os.makedirs(a.out, exist_ok=True)
    from ..forward.glass import FAMILIES, EMAX
    fams = [f for f in FAMILIES if not a.families or f in a.families.split(",")]
    combos = [(i, f) for i in EMAX for f in fams] if a.all else [(a.ion, a.family)]
    for ion, fam in combos:
        tag = f"{ion.replace('+', '')}_{fam}"
        if os.path.exists(os.path.join(a.out, tag + ".json")):
            continue
        try:
            r = window_probability(ion, fam, n=a.n)
        except Exception as e:
            print(tag, "error", e); continue
        if r is None:
            continue
        r["bands"] = bands(r)
        json.dump(r, open(os.path.join(a.out, tag + ".json"), "w"))
        print(tag, r["omega_src"], "robust", r["bands"]["robust"], "possible", r["bands"]["possible"], flush=True)


if __name__ == "__main__":
    main()
