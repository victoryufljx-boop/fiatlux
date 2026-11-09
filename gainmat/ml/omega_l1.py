"""L1: Judd-Ofelt Omega_2,4,6 from glass composition, built to generalise across glass families.

Structure (log Omega_lambda, lambda = 2, 4, 6):
    log Omega = a_ion + w . phi(host) + f_GBM(phi, ion)
  a_ion     per-ion offset (the 4f wavefunction part)
  w . phi   a shared linear host effect on physical descriptors (optical basicity, field
            strength, anion type, polarisability, former fractions) -- the same for every ion,
            so a host measured with one ion informs every other ion, and it extrapolates
  f_GBM     shallow boosted trees on the residual (interactions, saturation), kept small
Uncertainty: bootstrap ensemble (resampling publications).

Evaluation: publication-grouped CV; leave-one-glass-family-out; basicity extrapolation;
cross-ion transfer within a host.  Baselines: per-ion median, per-ion-and-family median,
gradient boosting on element fractions (the earlier model), gradient boosting on descriptors.

The LOMS data (Hrabovsky et al., Sci. Rep. 2025; www.LOMS.cz) are not redistributed; set
LOMS_CSV.  Outputs: ml_results/l1/*.json (+ plots), model pickle in private_data.
"""
from __future__ import annotations

import json
import os
import pickle

import numpy as np

from .descriptors import describe, NAMES
from ..forward.omega_model import parse_components, featurize, IONS_FIT

LOMS = os.environ.get("LOMS_CSV", "/home/claude/private_data/loms_jo.csv")


def load(path=LOMS):
    import pandas as pd
    d = pd.read_csv(path)
    d = d[d.host_type.str.startswith("Glass")]
    d = d[d.re_ion.isin(IONS_FIT)].dropna(subset=["omega2", "omega4", "omega6", "components"])
    d = d[(d.omega2 > 0) & (d.omega4 > 0) & (d.omega6 > 0) & (d.omega2 < 40) & (d.omega4 < 40) & (d.omega6 < 40)]
    rows, phis, fams, elx = [], [], [], []
    for i, r in d.iterrows():
        comp = parse_components(r.components)
        v, fam = describe(comp, dopant=r.re_ion)
        if v is None or not np.all(np.isfinite(v)):
            continue
        rows.append(i); phis.append(v); fams.append(fam); elx.append(featurize(comp, r.re_ion))
    d = d.loc[rows].copy()
    d["family"] = fams
    d["pub"] = d.pub_doi.fillna(d.pub_title).values
    d["comp_key"] = d.components.astype(str)
    Y = np.log10(d[["omega2", "omega4", "omega6"]].values)
    return d, np.array(phis), np.array(elx), Y


def ion_onehot(ions):
    return np.array([[1.0 if i == k else 0.0 for k in IONS_FIT] for i in ions], float).reshape(len(ions), len(IONS_FIT))


class Hybrid:
    """Shared linear host effect + per-ion offsets + small residual GBM, per Omega_lambda."""

    def __init__(self, alpha=30.0, gbm=True, seed=0):
        self.alpha, self.gbm, self.seed = alpha, gbm, seed

    def fit(self, phi, ions, y):
        from sklearn.linear_model import Ridge
        from sklearn.ensemble import GradientBoostingRegressor
        self.mu, self.sd = phi.mean(0), phi.std(0) + 1e-9
        z = (phi - self.mu) / self.sd
        self.zlo, self.zhi = z.min(0) - 0.5, z.max(0) + 0.5     # limit linear extrapolation
        X = np.hstack([ion_onehot(ions), z])
        self.lin = Ridge(alpha=self.alpha, fit_intercept=False).fit(X, y)
        if self.gbm:
            r = y - self.lin.predict(X)
            self.res = GradientBoostingRegressor(n_estimators=150, max_depth=2, learning_rate=0.05,
                                                 subsample=0.8, min_samples_leaf=8,
                                                 random_state=self.seed).fit(np.hstack([ion_onehot(ions), z]), r)
        return self

    def predict(self, phi, ions):
        z = np.clip((phi - self.mu) / self.sd, self.zlo, self.zhi)
        X = np.hstack([ion_onehot(ions), z])
        p = self.lin.predict(X)
        if self.gbm:
            p = p + self.res.predict(X)
        return p


class GBM:
    def __init__(self, seed=0):
        self.seed = seed

    def fit(self, X, ions, y):
        from sklearn.ensemble import GradientBoostingRegressor
        self.m = GradientBoostingRegressor(n_estimators=300, max_depth=3, learning_rate=0.05, subsample=0.8,
                                           random_state=self.seed).fit(np.hstack([X, ion_onehot(ions)]), y)
        return self

    def predict(self, X, ions):
        return self.m.predict(np.hstack([X, ion_onehot(ions)]))


def _median_by(keys_tr, y_tr, keys_te, fallback):
    out = np.empty(len(keys_te))
    for i, k in enumerate(keys_te):
        m = [y for kk, y in zip(keys_tr, y_tr) if kk == k]
        out[i] = np.median(m) if m else fallback[i]
    return out


def evaluate_split(d, phi, elx, Y, tr, te, n_boot=0):
    """Predictions of every model for one train/test split; returns dict name -> (n_te, 3)."""
    ions = d.re_ion.values; fam = d.family.values
    P = {k: np.zeros((len(te), 3)) for k in ("ion_median", "ion_family_median", "gbm_elements",
                                             "gbm_physics", "linear_shared", "hybrid")}
    S = np.zeros((len(te), 3))
    for j in range(3):
        y = Y[:, j]
        ion_med = _median_by(ions[tr], y[tr], ions[te], np.full(len(te), np.median(y[tr])))
        P["ion_median"][:, j] = ion_med
        P["ion_family_median"][:, j] = _median_by(list(zip(ions[tr], fam[tr])), y[tr],
                                                  list(zip(ions[te], fam[te])), ion_med)
        P["gbm_elements"][:, j] = GBM().fit(elx[tr], ions[tr], y[tr]).predict(elx[te], ions[te])
        P["gbm_physics"][:, j] = GBM().fit(phi[tr], ions[tr], y[tr]).predict(phi[te], ions[te])
        P["linear_shared"][:, j] = Hybrid(gbm=False).fit(phi[tr], ions[tr], y[tr]).predict(phi[te], ions[te])
        P["hybrid"][:, j] = Hybrid().fit(phi[tr], ions[tr], y[tr]).predict(phi[te], ions[te])
        if n_boot:
            pubs = np.unique(d.pub.values[tr]); rng = np.random.default_rng(j)
            bs = []
            for b in range(n_boot):
                pick = set(rng.choice(pubs, len(pubs), replace=True))
                idx = tr[np.isin(d.pub.values[tr], list(pick))]
                bs.append(Hybrid(seed=b).fit(phi[idx], ions[idx], y[idx]).predict(phi[te], ions[te]))
            S[:, j] = np.std(bs, axis=0)
    return P, S


def metrics(P, Yte, S=None):
    out = {}
    for k, p in P.items():
        e = p - Yte
        out[k] = dict(rmse_log10=np.sqrt(np.mean(e ** 2, axis=0)).round(3).tolist(),
                      median_abs_rel=np.median(np.abs(10 ** e - 1), axis=0).round(3).tolist())
    if S is not None:
        # coverage: |error| <= 1.64 sqrt(S^2 + s0^2), s0 = in-family scatter of the hybrid
        e = P["hybrid"] - Yte
        s0 = 0.12
        cov = np.mean(np.abs(e) <= 1.64 * np.sqrt(S ** 2 + s0 ** 2), axis=0)
        out["hybrid_90pc_interval_coverage"] = cov.round(3).tolist()
        out["hybrid_mean_sigma_log10"] = np.sqrt(np.mean(S ** 2 + s0 ** 2, axis=0)).round(3).tolist()
    return out


def run_all(out_dir="ml_results/l1", n_boot=12):
    from sklearn.model_selection import GroupKFold
    os.makedirs(out_dir, exist_ok=True)
    d, phi, elx, Y = load()
    n = len(Y); idx = np.arange(n)
    report = dict(n=n, n_publications=int(d.pub.nunique()), families=d.family.value_counts().to_dict(),
                  ions=d.re_ion.value_counts().to_dict(), descriptors=NAMES, splits={})
    print(f"{n} glass records, {d.pub.nunique()} publications")
    # 1. publication-grouped CV
    preds = {}
    allP = {k: np.zeros((n, 3)) for k in ("ion_median", "ion_family_median", "gbm_elements", "gbm_physics",
                                         "linear_shared", "hybrid")}
    allS = np.zeros((n, 3))
    for tr, te in GroupKFold(n_splits=5).split(phi, Y, d.pub.values):
        P, S = evaluate_split(d, phi, elx, Y, tr, te, n_boot=n_boot)
        for k in P:
            allP[k][te] = P[k]
        allS[te] = S
    report["splits"]["publication_5fold"] = metrics(allP, Y, allS)
    # 2. leave-one-family-out
    lofo = {}
    for fam in sorted(d.family.unique()):
        te = idx[d.family.values == fam]; tr = idx[d.family.values != fam]
        if len(te) < 8:
            continue
        P, S = evaluate_split(d, phi, elx, Y, tr, te, n_boot=n_boot)
        lofo[fam] = dict(n=len(te), **metrics(P, Y[te], S))
        print(f"LOFO {fam:18s} n={len(te):3d}  " + "  ".join(
            f"{k} {np.mean(v['rmse_log10']):.3f}" for k, v in lofo[fam].items() if isinstance(v, dict) and "rmse_log10" in v))
    report["splits"]["leave_one_family_out"] = lofo
    # 3. extrapolation in optical basicity (oxides only)
    b = phi[:, NAMES.index("basicity")]; ox = phi[:, NAMES.index("x_O")] > 0.5
    for name, mask_te in (("basicity_top20pc", ox & (b > np.quantile(b[ox], 0.8))),
                          ("basicity_bottom20pc", ox & (b < np.quantile(b[ox], 0.2)))):
        te = idx[mask_te]; tr = idx[~mask_te]
        P, S = evaluate_split(d, phi, elx, Y, tr, te, n_boot=n_boot)
        report["splits"][name] = dict(n=len(te), **metrics(P, Y[te], S))
    # 4. cross-ion transfer within a host: for compositions measured with >= 2 ions, hold out one
    #    ion's record(s); the other ions' records of that same glass stay in training
    keyc = d.comp_key.values; ions = d.re_ion.values
    multi = [k for k in set(keyc) if len(set(ions[keyc == k])) >= 2]
    te_list = []
    rng = np.random.default_rng(0)
    for k in multi:
        ii = idx[keyc == k]
        ion_out = rng.choice(sorted(set(ions[ii])))
        te_list += list(ii[ions[ii] == ion_out])
    te = np.array(sorted(set(te_list)))
    if len(te) >= 10:
        # also remove other records of the held-out (glass, ion) pairs; keep the other ions
        tr = np.setdiff1d(idx, te)
        P, S = evaluate_split(d, phi, elx, Y, tr, te, n_boot=n_boot)
        report["splits"]["cross_ion_same_glass"] = dict(n=len(te), n_glasses=len(multi), **metrics(P, Y[te], S))
    json.dump(report, open(os.path.join(out_dir, "omega_l1_report.json"), "w"), indent=1, default=float)
    return report, d, phi, Y


def fit_final(path_out="/home/claude/private_data/omega_l1.pkl", n_boot=20):
    d, phi, elx, Y = load()
    ions = d.re_ion.values
    models = []
    pubs = np.unique(d.pub.values); rng = np.random.default_rng(1)
    for b in range(n_boot):
        pick = set(rng.choice(pubs, len(pubs), replace=True))
        idx = np.where(np.isin(d.pub.values, list(pick)))[0]
        models.append([Hybrid(seed=b).fit(phi[idx], ions[idx], Y[idx, j]) for j in range(3)])
    full = [Hybrid().fit(phi, ions, Y[:, j]) for j in range(3)]
    pickle.dump(dict(full=full, boot=models, names=NAMES), open(path_out, "wb"))
    return full, models


_M = None


def predict(components: dict, ion: str, path="/home/claude/private_data/omega_l1.pkl"):
    """(Omega2, Omega4, Omega6) in 1e-20 cm^2 and their 1-sigma log10 spread."""
    global _M
    if _M is None:
        _M = pickle.load(open(path, "rb"))
    el = ion.replace("3+", "")
    v, fam = describe(components, dopant=el)
    phi = v[None, :]
    mu = np.array([m.predict(phi, [el])[0] for m in _M["full"]])
    bs = np.array([[m.predict(phi, [el])[0] for m in ms] for ms in _M["boot"]])
    return tuple(10 ** mu), tuple(bs.std(0)), fam
