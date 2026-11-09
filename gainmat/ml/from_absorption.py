"""'Give it one measured absorption spectrum' mode, and the test of whether it helps.

For each glass where the absorption AND emission of the same band were measured, predict the
emission four ways and score against the measured emission (never used):
  M0  composition-only forward model (no spectrum used)
  M1  forward model with its host parameters (cf_scale, disorder, gamma_h) refit to the measured
      absorption band; emission then comes from the refit Stark structure
  M2  standard practice: McCumber transform of the measured absorption,
      sigma_e = sigma_a exp((eps - h nu)/kT), with eps taken from the M0 model
  M3  as M2 with eps from the M1 (absorption-refit) model
Scores: emission shape (rms of peak-normalised curves), emission peak position (nm), and the
emission/absorption peak ratio (log10, the McCumber-sensitive quantity).

    python -m gainmat.ml.from_absorption  ->  ml_results/from_absorption/report.json
"""
from __future__ import annotations

import json
import os

import numpy as np

KT = 0.6950348 * 295.0


def cases():
    from .observations import build
    obs = [o for o in build() if o["kind"] == "curve" and o["family"] != "ZBLAN"]
    out = []
    for a in obs:
        if a["which"] != "abs":
            continue
        e = next((o for o in obs if o["which"] == "em" and o["src"] == a["src"] and o["ion"] == a["ion"]
                  and o["family"] == a["family"]), None)
        if e:
            out.append((a, e))
    return out


def _model(ion, family, theta, n_sites=20):
    from .host_calib import simulate
    r = simulate(family, ion, theta, n_sites)
    return r


def _curves(r, pair):
    from .host_calib import _band
    lam, sa = _band(r["_spec"], tuple(pair), "abs")
    _, se = _band(r["_spec"], tuple(pair), "em")
    return lam, sa, se


def _eps(lam, sa, se):
    """McCumber energy from a model's own sigma_e/sigma_a (exact in the model) at the abs peak."""
    i = int(np.argmax(sa))
    nu = 1e7 / lam[i]
    return nu + KT * np.log(se[i] / sa[i])


def scores(x, em_pred, em_meas, ab_meas):
    o = np.argsort(x); x, p, m, a = x[o], em_pred[o], em_meas[o], ab_meas[o]
    shape = float(np.sqrt(np.mean((p / p.max() - m / m.max()) ** 2)))
    shift = float(x[p.argmax()] - x[m.argmax()])
    ratio_err = float(np.log10((p.max() / a.max()) / (m.max() / a.max())))
    return dict(shape_rms=shape, peak_shift_nm=shift, peak_ratio_err_log10=ratio_err)


def refit_to_absorption(ion, family, a, n_sites=16, maxiter=60):
    from scipy.optimize import minimize
    from .host_calib import default_theta, predict
    th0 = default_theta(family)

    def unpack(x):
        return dict(th0, cf_scale=float(np.exp(x[0])), disorder=float(np.exp(x[1])), gamma_h=float(th0["gamma_h"] * np.exp(x[2])))

    def f(x):
        r = _model(ion, family, unpack(x), n_sites)
        v, _ = predict(a, r)
        return float(np.sum(v ** 2)) + 0.5 * float(np.sum((x / 1.0) ** 2))
    res = minimize(f, np.zeros(3), method="Nelder-Mead", options=dict(maxfev=maxiter, xatol=0.02, fatol=0.05,
                   initial_simplex=np.array([[0, 0, 0], [0.3, 0, 0], [0, 0.5, 0], [0, 0, 0.5]])))
    return unpack(res.x)


def main(out="ml_results/from_absorption"):
    from ..forward.glass import load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ
    from .host_calib import default_theta
    load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ, verbose=False)
    os.makedirs(out, exist_ok=True)
    rep = []
    for a, e in cases():
        ion, fam, pair = a["ion"], a["family"], tuple(a["pair"])
        xa = np.asarray(a["lam_nm"], float); ya = np.asarray(a["sigma"], float)
        oka = np.isfinite(ya) & (ya >= 0); xa, ya = xa[oka], ya[oka]
        xe = np.asarray(e["lam_nm"], float); ye = np.asarray(e["sigma"], float)
        ok = np.isfinite(ye) & (ye >= 0); xe, ye = xe[ok], ye[ok]
        # score where the measured emission is above 5 % of its peak, inside the measured absorption range
        # (McCumber is unreliable in the far tails; the same window applies to every method)
        ok = (ye >= 0.05 * ye.max()) & (xe >= xa.min()) & (xe <= xa.max()); xe, ye = xe[ok], ye[ok]
        ab_on_e = np.interp(xe, *zip(*sorted(zip(xa, ya))))
        row = dict(ion=ion, family=fam, src=a["src"])
        th0 = default_theta(fam)
        r0 = _model(ion, fam, th0); l0, sa0, se0 = _curves(r0, pair)
        m0 = np.interp(xe, l0[::-1], se0[::-1])
        a0 = np.interp(xe, l0[::-1], sa0[::-1])
        # M0 uses its own absorption for the ratio
        s = scores(xe, m0, ye, ab_on_e); s["peak_ratio_err_log10"] = float(np.log10((m0.max() / a0.max()) / (ye.max() / ab_on_e.max())))
        row["M0_model_only"] = s
        th1 = refit_to_absorption(ion, fam, a)
        r1 = _model(ion, fam, th1); l1, sa1, se1 = _curves(r1, pair)
        m1 = np.interp(xe, l1[::-1], se1[::-1]); a1 = np.interp(xe, l1[::-1], sa1[::-1])
        s = scores(xe, m1, ye, ab_on_e); s["peak_ratio_err_log10"] = float(np.log10((m1.max() / a1.max()) / (ye.max() / ab_on_e.max())))
        row["M1_model_refit_to_absorption"] = dict(s, theta=th1)
        for tag, (l, sa, se) in (("M2_mccumber_eps_from_model", (l0, sa0, se0)), ("M3_mccumber_eps_from_refit", (l1, sa1, se1))):
            eps = _eps(l, sa, se)
            pred = ab_on_e * np.exp((eps - 1e7 / xe) / KT)
            row[tag] = dict(scores(xe, pred, ye, ab_on_e), eps_cm1=float(eps))
        rep.append(row)
        print(f"{ion} {fam} ({a['src'][:28]})")
        for k in ("M0_model_only", "M1_model_refit_to_absorption", "M2_mccumber_eps_from_model", "M3_mccumber_eps_from_refit"):
            v = row[k]
            print(f"   {k:30s} shape {v['shape_rms']:.3f}  peak shift {v['peak_shift_nm']:+6.1f} nm  em/abs ratio err {v['peak_ratio_err_log10']:+.3f} dex", flush=True)
    summ = {}
    for k in ("M0_model_only", "M1_model_refit_to_absorption", "M2_mccumber_eps_from_model", "M3_mccumber_eps_from_refit"):
        summ[k] = dict(shape_rms_mean=float(np.mean([r[k]["shape_rms"] for r in rep])),
                       peak_shift_abs_mean_nm=float(np.mean([abs(r[k]["peak_shift_nm"]) for r in rep])),
                       ratio_err_rms_dex=float(np.sqrt(np.mean([r[k]["peak_ratio_err_log10"] ** 2 for r in rep]))))
        print(f"{k:30s} mean shape {summ[k]['shape_rms_mean']:.3f}  |peak shift| {summ[k]['peak_shift_abs_mean_nm']:.1f} nm  ratio rms {summ[k]['ratio_err_rms_dex']:.3f} dex")
    json.dump(dict(cases=rep, summary=summ), open(os.path.join(out, "report.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
