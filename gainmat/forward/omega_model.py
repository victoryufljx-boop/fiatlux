"""Judd-Ofelt Omega predictor: glass composition (mol%) + ion -> (Omega2, Omega4, Omega6).

Trained on the LOMS.cz Judd-Ofelt dataset (Hrabovsky, Varak, Krystufek, Sci. Rep. 2025,
doi:10.1038/s41598-025-13620-0; data: figshare 10.6084/m9.figshare.c.7904465, www.LOMS.cz).
The raw dataset is NOT redistributed in this repository (its licence is personal / educational /
internal use); point LOMS_CSV at a local export to train.  Only derived quantities are committed.

Features: element fractions computed from the mol% components (cations normalised among
cations, anions among anions), plus the ion.  Model: gradient boosting on log(Omega + 0.1),
one model per Omega_lambda, all ions pooled with an ion one-hot.  Validation is grouped by
publication, so a paper's samples never sit in both train and test.
"""
from __future__ import annotations

import json
import os
import re

import numpy as np

ANIONS = ["O", "F", "S", "Se", "Cl", "Br", "I", "N", "Te_anion"]
CATIONS = ["Si", "Al", "B", "P", "Ge", "Te", "Bi", "Pb", "Ga", "Nb", "Ti", "Zn", "Zr", "Hf", "In",
           "Sb", "As", "W", "Mo", "V", "Ta", "Li", "Na", "K", "Rb", "Cs", "Mg", "Ca", "Sr", "Ba",
           "Cd", "La", "Y", "Gd", "Lu", "Sc", "Ce", "Ag", "Tl", "Sn", "Th"]
RE = {"Pr", "Nd", "Pm", "Sm", "Eu", "Gd_dop", "Tb", "Dy", "Ho", "Er", "Tm", "Yb"}
IONS_FIT = ["Pr", "Nd", "Sm", "Eu", "Tb", "Dy", "Ho", "Er", "Tm"]
_SUB = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")


def parse_formula(f: str) -> dict:
    """Element counts of a chemical formula with digits, decimals and parentheses."""
    f = f.translate(_SUB).replace(" ", "")
    tokens = re.findall(r"([A-Z][a-z]?|\(|\)|\d+(?:\.\d+)?)", f)
    stack = [{}]
    i = 0
    while i < len(tokens):
        t = tokens[i]
        if t == "(":
            stack.append({})
        elif t == ")":
            grp = stack.pop()
            mult = 1.0
            if i + 1 < len(tokens) and re.fullmatch(r"\d+(?:\.\d+)?", tokens[i + 1]):
                mult = float(tokens[i + 1]); i += 1
            for k, v in grp.items():
                stack[-1][k] = stack[-1].get(k, 0) + v * mult
        elif re.fullmatch(r"[A-Z][a-z]?", t):
            n = 1.0
            if i + 1 < len(tokens) and re.fullmatch(r"\d+(?:\.\d+)?", tokens[i + 1]):
                n = float(tokens[i + 1]); i += 1
            stack[-1][t] = stack[-1].get(t, 0) + n
        i += 1
    return stack[0]


def parse_components(s: str) -> dict:
    """'SiO2=60 mol% | Na2O=40 mol%' -> {'SiO2': 60, 'Na2O': 40}."""
    out = {}
    for part in str(s).split("|"):
        m = re.match(r"\s*(.+?)=\s*([\d.]+)\s*mol%", part)
        if m:
            out[m.group(1).strip()] = float(m.group(2))
    return out


def featurize(components: dict, ion: str) -> np.ndarray:
    """Element-fraction features for a composition {formula: mol%}; the dopant itself
    (as RE fluoride/oxide) is excluded from the host fractions."""
    el = {}
    for f, x in components.items():
        try:
            c = parse_formula(f)
        except Exception:
            continue
        if ion in c and len(c) <= 3:                 # the dopant compound itself
            continue
        for k, v in c.items():
            el[k] = el.get(k, 0) + v * x
    an = {k: el.get(k, 0.0) for k in ANIONS if k in el}
    ca = {k: v for k, v in el.items() if k not in ANIONS and k not in RE and k != ion}
    sa, sc = sum(an.values()) or 1.0, sum(ca.values()) or 1.0
    xa = [an.get(k, 0.0) / sa for k in ANIONS]
    xc = [ca.get(k, 0.0) / sc for k in CATIONS]
    other = 1.0 - sum(xc)
    onehot = [1.0 if ion == i else 0.0 for i in IONS_FIT]
    return np.array(xa + xc + [other] + onehot)


def load_dataset(path: str, glass_only=True):
    import pandas as pd
    d = pd.read_csv(path)
    if glass_only:
        d = d[d.host_type.str.startswith("Glass")]
    d = d[d.re_ion.isin(IONS_FIT)].dropna(subset=["omega2", "omega4", "omega6", "components"])
    d = d[(d.omega2 >= 0) & (d.omega4 > 0) & (d.omega6 > 0) & (d.omega2 < 40) & (d.omega4 < 40) & (d.omega6 < 40)]
    X = np.array([featurize(parse_components(c), ion) for c, ion in zip(d.components, d.re_ion)])
    Y = np.log(d[["omega2", "omega4", "omega6"]].values + 0.1)
    groups = d.pub_doi.fillna(d.pub_title).values
    return d, X, Y, groups


def fit(path=None, out=None, verbose=True):
    from sklearn.ensemble import GradientBoostingRegressor
    from sklearn.model_selection import GroupKFold
    path = path or os.environ.get("LOMS_CSV")
    d, X, Y, groups = load_dataset(path)
    gkf = GroupKFold(n_splits=5)
    pred = np.zeros_like(Y)
    base = np.zeros_like(Y)
    ions = d.re_ion.values
    for tr, te in gkf.split(X, Y, groups):
        for j in range(3):
            m = GradientBoostingRegressor(n_estimators=300, max_depth=3, learning_rate=0.05,
                                          subsample=0.8, random_state=0)
            m.fit(X[tr], Y[tr, j]); pred[te, j] = m.predict(X[te])
            for ion in set(ions[te]):          # baseline: per-ion median of the training folds
                mask_tr = ions[tr] == ion
                base[te[ions[te] == ion], j] = np.median(Y[tr][mask_tr, j]) if mask_tr.any() else np.median(Y[tr, j])
    err = lambda P: np.median(np.abs(np.exp(P) - np.exp(Y)) / np.exp(Y), axis=0)
    report = dict(n=len(Y), n_publications=len(set(groups)),
                  median_rel_err_model=dict(zip(("O2", "O4", "O6"), err(pred).round(3).tolist())),
                  median_rel_err_ion_median=dict(zip(("O2", "O4", "O6"), err(base).round(3).tolist())))
    models = []
    for j in range(3):
        m = GradientBoostingRegressor(n_estimators=300, max_depth=3, learning_rate=0.05,
                                      subsample=0.8, random_state=0)
        m.fit(X, Y[:, j]); models.append(m)
    if verbose:
        print(json.dumps(report, indent=1))
    if out:
        import pickle
        pickle.dump(dict(models=models, report=report), open(out, "wb"))
    return models, report


_MODELS = None


def predict(components: dict, ion: str, model_path=None):
    """Omega (1e-20 cm^2) for a host composition {formula: mol%} and ion ('Er3+' or 'Er')."""
    global _MODELS
    ion = ion.replace("3+", "")
    if _MODELS is None:
        import pickle
        p = model_path or os.environ.get("OMEGA_MODEL", "/home/claude/private_data/omega_model.pkl")
        _MODELS = pickle.load(open(p, "rb"))["models"]
    x = featurize(components, ion)[None, :]
    return tuple(float(np.exp(m.predict(x)[0]) - 0.1) for m in _MODELS)


# representative compositions (mol%) of the forward model's glass families
FAMILY_COMPOSITIONS = {
    "silica": {"SiO2": 100},
    "aluminosilicate": {"SiO2": 94, "Al2O3": 6},
    "germanosilicate": {"SiO2": 90, "GeO2": 10},
    "phosphosilicate": {"SiO2": 90, "P2O5": 10},
    "phosphate": {"P2O5": 60, "Al2O3": 10, "K2O": 15, "BaO": 15},
    "germanate": {"GeO2": 60, "PbO": 30, "Na2O": 10},
    "tellurite": {"TeO2": 75, "ZnO": 20, "Na2O": 5},
    "ZBLAN": {"ZrF4": 53, "BaF2": 20, "LaF3": 4, "AlF3": 3, "NaF": 20},
    "fluoroindate": {"InF3": 40, "ZnF2": 20, "SrF2": 20, "BaF2": 20},
    "chalcogenide": {"Ga2S3": 70, "La2S3": 30},
}
