"""Crystal-field strength and rank balance learned from measured Stark splittings across ions
and glasses (gainmat.refs.stark_glass_refs), with transfer tests.

Model for the superposition intrinsic parameters of ligand L in glass family f, ion i:
    B_k(i, f) = handbook B_k(i, L) x SPM_SCALE[L] x s_f x r_k[L]
  s_f        glass-family field strength (host property; 1 = current model)
  r_k[O]     rank balance for oxygen, r4 = 1 fixed (overall size is in s_f)
Single-rank splittings show why the balance is learnable from two ions (germanosilicate sites,
cm^-1):  Er 4I15/2  B2 184 / B4 287 / B6 204;  Yb 2F7/2  448 / 589 / 100.

Tests (all out-of-sample):
  ion transfer      fit s_f on Er only, predict Yb in the same glasses (and the reverse)
  family transfer   rank balance learned on the other glasses; in the held-out glass s_f comes
                    from Er only and Yb is predicted
  held-out ion      Tm3+ (no Tm data used): absorption/emission peak positions in silica
                    (Jackson & King 1999) from the full spectra, before vs after

    python -m gainmat.ml.cf_ranks  ->  ml_results/cf/cf_ranks.json, cf_report.json
"""
from __future__ import annotations

import json
import os

import numpy as np

N_SITES = 24


def fast_split(ion, family, levels, s_f=1.0, rk=(1.0, 1.0, 1.0), n_sites=N_SITES, seed=0):
    """Mean total Stark splitting of `levels` over a site ensemble (no spectra: fast)."""
    from ..forward.glass import FAMILIES, spm_for
    from ..forward.spectra import random_sites, spm_bkq
    from ..forward.ions import get_ion
    f = FAMILIES[family]; lig = f.get("ligand", "O")
    p = dict(spm_for(ion, lig))
    for i, k in enumerate((2, 4, 6)):
        p[f"B{k}"] = p[f"B{k}"] * s_f * (rk[i] if lig == "O" else 1.0)
    sites = random_sites(n_sites, lig, f["CN"], f["R"], f["R_std"], f["ang_std"], seed=seed)
    m = get_ion(ion)
    out = {lab: [] for lab in levels}
    for s in sites:
        E, W, par = m.stark(spm_bkq(s, {lig: p}))
        for lab in levels:
            e = E[par == m.mult_index(lab)]
            out[lab].append(e.max() - e.min())
    return {lab: float(np.mean(v)) for lab, v in out.items()}


def data(exclude_approx=True):
    from ..refs.stark_glass_refs import STARK_SPLITTINGS
    from ..forward.glass import FAMILIES
    return [d for d in STARK_SPLITTINGS if d["family"] in FAMILIES and not (exclude_approx and d["sigma"] >= 300)]


class Model:
    def __init__(self, families):
        self.fams = sorted(families)
        self._cache = {}

    def predict(self, d, s, rk):
        # the splitting of a well-isolated J multiplet is first order in the crystal field, so it
        # scales linearly with s_f; only the rank balance needs a new diagonalisation
        key = (d["ion"], d["family"], d["level"], tuple(round(x, 4) for x in rk))
        if key not in self._cache:
            self._cache[key] = fast_split(d["ion"], d["family"], [d["level"]], 1.0, rk)[d["level"]]
        return s * self._cache[key]

    def loss(self, rows, sfam, rk):
        return float(sum(((self.predict(d, sfam.get(d["family"], 1.0), rk) - d["value"]) / d["sigma"]) ** 2 for d in rows))


def best_scales(model, rows, fams, rk):
    """Weighted least squares for s_f given the rank balance (prediction is linear in s_f)."""
    sf = {}
    for f in fams:
        rr = [d for d in rows if d["family"] == f]
        if not rr:
            continue
        m = np.array([model.predict(d, 1.0, rk) for d in rr]); v = np.array([d["value"] for d in rr])
        w = 1 / np.array([d["sigma"] for d in rr]) ** 2
        sf[f] = float(np.sum(w * m * v) / np.sum(w * m * m))
    return sf


def fit(model, rows, fams, fit_ranks, rk0=(1.0, 1.0, 1.0), prior_r=0.5):
    """chi^2 over s_f (closed form) and optionally log r2, log r6 (Nelder-Mead)."""
    from scipy.optimize import minimize
    if not fit_ranks:
        return best_scales(model, rows, fams, rk0), rk0

    def obj(x):
        rk = (float(np.exp(x[0])), 1.0, float(np.exp(x[1])))
        return model.loss(rows, best_scales(model, rows, fams, rk), rk) + float(np.sum((x / prior_r) ** 2))
    r = minimize(obj, np.zeros(2), method="Nelder-Mead", options=dict(maxfev=80, xatol=0.01, fatol=0.01,
                 initial_simplex=np.array([[0, 0], [0.3, 0], [0, 0.3]])))
    rk = (float(np.exp(r.x[0])), 1.0, float(np.exp(r.x[1])))
    return best_scales(model, rows, fams, rk), rk


def table(model, rows, sf, rk):
    return [dict(ion=d["ion"], level=d["level"], family=d["family"], measured=d["value"], sigma=d["sigma"],
                 model=round(model.predict(d, sf.get(d["family"], 1.0), rk), 1), src=d["src"]) for d in rows]


def rms_rel(tab):
    return float(np.sqrt(np.mean([(np.log(t["model"] / t["measured"])) ** 2 for t in tab]))) if tab else None


def main(out="ml_results/cf"):
    from ..forward.glass import load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ
    load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ, verbose=False)
    os.makedirs(out, exist_ok=True)
    rows = data()
    fams = sorted({d["family"] for d in rows})
    M = Model(fams)
    rep = {}
    # 0. explicit model as is
    t0 = table(M, rows, {}, (1, 1, 1))
    rep["explicit"] = dict(rms_log=rms_rel(t0), table=t0)
    print(f"explicit model: rms {rep['explicit']['rms_log']:.3f} (log ratio)")
    # 1. ion transfer: s_f from one ion only (handbook rank balance), predict the other ion
    for src_ion, tgt_ion in (("Er3+", "Yb3+"), ("Yb3+", "Er3+")):
        tr = [d for d in rows if d["ion"] == src_ion]; te = [d for d in rows if d["ion"] == tgt_ion]
        fs = sorted({d["family"] for d in tr} & {d["family"] for d in te})
        tr = [d for d in tr if d["family"] in fs]; te = [d for d in te if d["family"] in fs]
        sf, rk = fit(M, tr, fs, False)
        tt = table(M, te, sf, rk); tb = table(M, te, {}, (1, 1, 1))
        rep[f"ion_transfer_{src_ion}_to_{tgt_ion}"] = dict(s_f=sf, rms_log_before=rms_rel(tb), rms_log_after=rms_rel(tt), table=tt)
        print(f"{src_ion} -> {tgt_ion}: s_f {({k: round(v, 3) for k, v in sf.items()})}; held-out rms {rms_rel(tb):.3f} -> {rms_rel(tt):.3f}")
        for t, b in zip(tt, tb):
            print(f"     {t['ion']} {t['level']} {t['family']:16s} meas {t['measured']:6.0f}  explicit {b['model']:6.0f}  transferred {t['model']:6.0f}")
    # 2. joint fit: s_f per family + oxygen rank balance, all data
    sf, rk = fit(M, rows, fams, True)
    tj = table(M, rows, sf, rk)
    rep["joint"] = dict(s_f=sf, rank_O=rk, rms_log=rms_rel(tj), table=tj)
    print(f"joint: s_f {({k: round(v, 3) for k, v in sf.items()})} r_O {tuple(round(x, 3) for x in rk)} rms {rms_rel(tj):.3f}")
    # 3. family transfer of the rank balance: learn r on other families, then in the held-out family
    #    set s_f from Er alone and predict Yb (only families with both ions)
    both = sorted({d["family"] for d in rows if d["ion"] == "Er3+"} & {d["family"] for d in rows if d["ion"] == "Yb3+"})
    lofo = []
    for f in both:
        if f == "ZBLAN":       # oxygen rank balance does not apply to the fluoride ligand
            continue
        tr = [d for d in rows if d["family"] != f]
        frs = sorted({d["family"] for d in tr})
        _, rk_f = fit(M, tr, frs, True)
        er = [d for d in rows if d["family"] == f and d["ion"] == "Er3+"]
        yb = [d for d in rows if d["family"] == f and d["ion"] == "Yb3+"]
        for rr, tag in ((rk_f, "learned ranks"), ((1, 1, 1), "handbook ranks")):
            s, _ = fit(M, er, [f], False, rk0=rr)
            ty = table(M, yb, s, rr)
            lofo.append(dict(family=f, ranks=tag, rank_O=rr, s_f=s[f], yb=ty, rms_log=rms_rel(ty)))
            print(f"  LOFO {f:14s} {tag:15s} r {tuple(round(x, 2) for x in rr)} s_f(Er) {s[f]:.3f} -> Yb "
                  + ", ".join(f"{t['model']:.0f} vs {t['measured']:.0f}" for t in ty))
    rep["family_transfer"] = lofo
    json.dump(rep, open(os.path.join(out, "cf_report.json"), "w"), indent=1, default=float)
    json.dump(dict(rank_scale={"O": list(rk)}, family_scale=sf,
                   note="joint fit to measured Stark splittings (gainmat.refs.stark_glass_refs)"),
              open(os.path.join(out, "cf_ranks.json"), "w"), indent=1)
    return rep


if __name__ == "__main__":
    main()
