"""L2: host microstructure, learned jointly across every ion measured in that host.

Parameters of a host (theta):
  cf_scale   multiplier on the crystal-field strength (superposition-model B_k)
  disorder   multiplier on the site disorder (radial R_std and angular ang_std together)
  neph_std   per-site fractional scatter of multiplet energies (nephelauxetic)
  gamma_h    homogeneous linewidth at room temperature (cm^-1)
These are properties of the glass, not of the ion: one theta must reproduce the widths,
splittings, shapes and sigma_e/sigma_a ratios of every ion doped into that glass.  The
leave-one-ion-out test fits theta on all ions but one and predicts the held-out ion.

    python -m gainmat.ml.host_calib --family ZBLAN [--holdout Er3+] [--sites 24] --out ml_results/l2
"""
from __future__ import annotations

import argparse
import copy
import json
import os
import time

import numpy as np

PARAMS = ("cf_scale", "disorder", "neph_std", "gamma_h")
BOUNDS = dict(cf_scale=(0.4, 2.5), disorder=(0.2, 3.0), neph_std=(0.0, 0.012), gamma_h=(5.0, 250.0))


def default_theta(family):
    from ..forward.glass import FAMILIES, NEPH_STD
    return dict(cf_scale=1.0, disorder=1.0, neph_std=NEPH_STD, gamma_h=FAMILIES[family].get("gamma_h", 49.0))


def _to_x(theta):
    return np.array([np.log(theta["cf_scale"]), np.log(theta["disorder"]),
                     np.sqrt(max(theta["neph_std"], 0.0) / 0.001), np.log(theta["gamma_h"])])


def _from_x(x):
    th = dict(cf_scale=float(np.exp(x[0])), disorder=float(np.exp(x[1])),
              neph_std=float(0.001 * x[2] ** 2), gamma_h=float(np.exp(x[3])))
    for k, (lo, hi) in BOUNDS.items():
        th[k] = float(np.clip(th[k], lo, hi))
    return th


DNU = 4.0          # spectral grid for calibration (cm^-1); widths of interest are >= 100 cm^-1


def simulate(family, ion, theta, n_sites=24, seed=0):
    from ..forward.glass import run, FAMILIES
    f = FAMILIES[family]
    ov = dict(cf_scale=theta["cf_scale"], R_std=f["R_std"] * theta["disorder"],
              ang_std=f["ang_std"] * theta["disorder"], neph_std=theta["neph_std"], gamma_h=theta["gamma_h"])
    r = run(ion, family, n_sites=n_sites, pumps=[], overrides=ov, seed=seed, dnu=DNU)
    r["_spec"]["ion"] = ion
    return r


def _band(spec, pair, which, T=295.0):
    """Cross section of a transition; an upper level written 'A,B' is a thermally coupled group
    (e.g. Ho 5S2,5F4): emission is the Boltzmann-weighted sum, (2J+1) exp(-E/kT) per multiplet."""
    nu = spec["nu"]
    key = "sigma_em" if which == "em" else "sigma_abs"
    ups = [u.strip() for u in pair[0].split(",")]
    if len(ups) == 1:
        return 1e7 / nu, spec[key][(ups[0], pair[1])]
    from ..forward.ions import get_ion
    m = get_ion(spec.get("ion", "")) if spec.get("ion") else None
    w, tot = [], np.zeros_like(nu)
    for u in ups:
        if (u, pair[1]) not in spec[key]:
            continue
        J = None
        if m is not None:
            J = m.multiplets[m.mult_index(u)]["J"]
        g = 2 * J + 1 if J is not None else 1.0
        w.append((g * np.exp(-spec["bary"][u] / (0.695 * T)), spec[key][(u, pair[1])]))
    z = sum(a for a, _ in w)
    for a, s in w:
        tot = tot + a / z * s
    return 1e7 / nu, tot


def predict(obs, r):
    """Model value(s) for one observation from a run() result; returns residual vector."""
    spec = r["_spec"]
    k = obs["kind"]
    if k == "eff_width_nm":
        lam, s = _band(spec, tuple(obs["pair"]), obs.get("which", "em"))
        lo, hi = obs["band_nm"]
        m = (lam >= lo) & (lam <= hi)
        if m.sum() < 3 or s[m].max() <= 0:
            return np.array([3.0]), dict(model=None)
        l, y = lam[m][::-1], s[m][::-1]
        w = np.trapezoid(y, l) / y.max()
        return np.array([np.log(w / obs["value"]) / obs["sigma_rel"]]), dict(model=float(w))
    if k == "curve":
        lam, s = _band(spec, tuple(obs["pair"]), obs["which"])
        x = np.asarray(obs["lam_nm"], float); y = np.asarray(obs["sigma"], float)
        ok = np.isfinite(y) & (y >= 0)
        x, y = x[ok], y[ok]
        b = np.interp(x, lam[::-1], s[::-1])
        if b.max() <= 0:
            return np.array([3.0, 3.0, 3.0]), dict(model=None)
        o = np.argsort(x); x, y, b = x[o], y[o], b[o]
        shape = np.sqrt(np.mean((b / b.max() - y / y.max()) ** 2))
        wm = np.trapezoid(b, x) / b.max(); we = np.trapezoid(y, x) / y.max()
        shift = x[b.argmax()] - x[y.argmax()]
        res = [shape / 0.05, np.log(wm / we) / obs["sigma_rel"], shift / 15.0]
        return np.array(res), dict(shape_rms=float(shape), width_model=float(wm), width_meas=float(we),
                                   peak_shift_nm=float(shift))
    if k == "split_cm1":
        v = r["split"].get(obs["label"])
        return np.array([(v - obs["value"]) / obs["sigma_abs"]]), dict(model=float(v))
    if k == "zpl_fwhm_cm1":
        se = spec["site_edges"]
        z = np.array([s[obs["upper"]][0] - s[obs["lower"]][0] for s in se])
        v = 2.355 * z.std()
        return np.array([(v - obs["value"]) / obs["sigma_abs"]]), dict(model=float(v))
    if k == "ratio_log10":
        nu = spec["nu"]; x = 1e7 / obs["lam_point"]
        se = np.interp(x, nu, spec["sigma_em"][tuple(obs["pair"])])
        sa = np.interp(x, nu, spec["sigma_abs"][tuple(obs["pair"])])
        v = np.log10(max(se, 1e-40) / max(sa, 1e-40))
        return np.array([(v - obs["value"]) / obs["sigma_abs"]]), dict(model=float(v))
    raise ValueError(k)


def evaluate(family, theta, obs, n_sites=24, seed=0, detail=False):
    ions = sorted({o["ion"] for o in obs})
    res, info = [], []
    for ion in ions:
        r = simulate(family, ion, theta, n_sites, seed)
        for o in obs:
            if o["ion"] != ion:
                continue
            v, d = predict(o, r)
            res.append(v)
            if detail:
                info.append(dict(ion=ion, kind=o["kind"], pair=o.get("pair") or o.get("label"),
                                 resid=v.tolist(), **d, src=o["src"]))
    r = np.concatenate(res) if res else np.zeros(1)
    return (float(np.sum(r ** 2)), info) if detail else float(np.sum(r ** 2))


def fit(family, obs, n_sites=24, seed=0, theta0=None, maxiter=120, log=None):
    from scipy.optimize import minimize
    th0 = theta0 or default_theta(family)
    x0 = _to_x(th0)
    hist = []

    def f(x):
        th = _from_x(x)
        L = evaluate(family, th, obs, n_sites, seed)
        prior = 0.5 * ((x[0] / 0.5) ** 2 + (x[1] / 1.0) ** 2)   # weak pull to the explicit model
        hist.append((L + prior, th))
        if log:
            log(f"  eval {len(hist):3d} loss {L:9.2f} prior {prior:5.2f} {json.dumps({k: round(v, 4) for k, v in th.items()})}")
        return L + prior
    simplex = np.array([x0] + [x0 + d for d in np.diag([0.35, 0.5, 0.8, 0.5])])
    out = minimize(f, x0, method="Nelder-Mead",
                   options=dict(maxfev=maxiter, xatol=0.02, fatol=0.05, initial_simplex=simplex))
    best = min(hist, key=lambda h: h[0])
    return best[1], best[0], len(hist)


def main():
    from ..forward.glass import load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ
    from .observations import build
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", required=True)
    ap.add_argument("--holdout", default="")
    ap.add_argument("--sites", type=int, default=24)
    ap.add_argument("--maxiter", type=int, default=120)
    ap.add_argument("--out", default="ml_results/l2")
    a = ap.parse_args()
    load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ, verbose=False)
    os.makedirs(a.out, exist_ok=True)
    allobs = [o for o in build() if o["family"] == a.family]
    train = [o for o in allobs if o["ion"] != a.holdout]
    test = [o for o in allobs if o["ion"] == a.holdout]
    tag = f"{a.family}_{'all' if not a.holdout else 'minus_' + a.holdout.replace('+', '')}"
    logf = open(os.path.join(a.out, tag + ".log"), "w")
    log = lambda s: (print(s, flush=True), logf.write(s + "\n"), logf.flush())
    t0 = time.time()
    th0 = default_theta(a.family)
    L0_train, info0_train = evaluate(a.family, th0, train, a.sites, detail=True)
    log(f"{tag}: {len(train)} train obs, {len(test)} test obs; default loss {L0_train:.2f}")
    th, Lbest, nev = fit(a.family, train, a.sites, maxiter=a.maxiter, log=log)
    L_train, info_train = evaluate(a.family, th, train, a.sites, detail=True)
    out = dict(family=a.family, holdout=a.holdout or None, theta_default=th0, theta_fit=th,
               loss_train_default=L0_train, loss_train_fit=L_train, n_train=len(train), n_test=len(test),
               evals=nev, seconds=time.time() - t0, detail_train=info_train)
    if test:
        L0, i0 = evaluate(a.family, th0, test, a.sites, detail=True)
        L1, i1 = evaluate(a.family, th, test, a.sites, detail=True)
        out.update(loss_test_default=L0, loss_test_fit=L1, detail_test_default=i0, detail_test_fit=i1)
        log(f"held-out {a.holdout}: loss default {L0:.2f} -> learned {L1:.2f}")
    json.dump(out, open(os.path.join(a.out, tag + ".json"), "w"), indent=1, default=float)
    log(f"done in {time.time() - t0:.0f} s; theta {th}")


if __name__ == "__main__":
    main()
