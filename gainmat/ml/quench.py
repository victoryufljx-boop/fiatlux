"""L3: extrinsic quenching of the ~2 um levels in silica-family fibres.

What the data say (ml_results/l3/extrinsic_quenching.json: measured lifetime vs the
single-ion model's intrinsic rate, W_extra = 1/tau_obs - 1/tau_model):
  * fluoride (ZBLAN) glass: W_extra ~ 0 for Tm 3F4, Er 4I13/2, Ho 5I7, from 1.6e19 to 3.2e20 cm^-3
  * silica-family fibres: W_extra = 1.2e3 - 2.6e3 s^-1 for Tm 3F4 and Ho 5I7 (5300-6000 cm^-1),
    rising only weakly with concentration; Er 4I13/2 (6700 cm^-1) only ~20 s^-1.
Dry fluorides vs OH-containing silica, a resonance window near the OH combination bands, and a
sub-linear concentration dependence (energy migration to a fixed set of quenching sites) are
the physical constraints.  The model therefore has the form
    W_q = 0                                    (fluoride hosts)
    W_q = A (N / 1e20 cm^-3)^p, p fixed at 0   (silica family, level energy 4800-6300 cm^-1; see P_FIXED)
    W_q = 0 otherwise (Er 4I13/2 class: tens of s^-1, below the model's own error)
with A fitted.  Tested by predicting each point from the others, Ho 5I7 from Tm data only, and an
independent low-concentration series (Simpson thesis), which is what fixed p at 0.
"""
from __future__ import annotations

import json
import os

import numpy as np

WINDOW = (4800.0, 6300.0)
# Concentration exponent.  The first version fitted p = 0.43 on fibres with N = 5e19 - 3e20.  An
# independent series (Simpson, PhD thesis, Victoria University: nine Tm aluminosilicate fibres,
# Tm2O3 100 - 1950 ppm, N = 7e17 - 1.3e19 cm^-3) has 3F4 lifetimes of 384 - 573 us with no
# concentration trend; the N^0.43 law predicts 1.0 - 2.2 ms there (rms log error 1.32) while an
# N-independent rate fitted to the same original points gives rms 0.18.  The extra decay is a
# property of the silica host, not concentration quenching, so p is fixed at 0.
P_FIXED = 0.0
FLUORIDE = {"ZBLAN", "fluoroindate"}
PARAMS = dict(A=None, p=None)


def _fit(pts, p_fixed="default"):
    N = np.array([p["N_cm3"] for p in pts]); W = np.array([p["W_extra_s1"] for p in pts])
    x = np.log(N / 1e20); y = np.log(W)
    pf = P_FIXED if p_fixed == "default" else p_fixed
    if pf is not None:
        return float(np.exp(np.mean(y - pf * x))), float(pf)
    if len(pts) >= 3 and np.ptp(x) > 0.3:
        p, lnA = np.polyfit(x, y, 1)
    else:
        p, lnA = 0.45, np.mean(y - 0.45 * x)
    return float(np.exp(lnA)), float(p)


def load_params(path=None):
    path = path or os.path.join(os.path.dirname(__file__), "..", "..", "ml_results", "l3", "quench_model.json")
    if PARAMS["A"] is None and os.path.exists(path):
        d = json.load(open(path)); PARAMS.update(A=d["A_s1"], p=d["p"])
    return PARAMS


def extra_rates(family, bary: dict, N_cm3):
    """{level: W_q} for every multiplet of an ion in this host (levels outside the window: none)."""
    load_params()
    if PARAMS["A"] is None:
        return {}
    out = {}
    for lab, E in bary.items():
        if E and WINDOW[0] <= E <= WINDOW[1]:
            w = rate(family, E, N_cm3)
            if w > 0:
                out[lab] = w
    return out


def rate(family, level_energy_cm1, N_cm3, A=None, p=None):
    if family in FLUORIDE:
        return 0.0
    if not (WINDOW[0] <= level_energy_cm1 <= WINDOW[1]):
        return 0.0
    A = A if A is not None else PARAMS["A"]; p = p if p is not None else PARAMS["p"]
    return A * (N_cm3 / 1e20) ** p


def load_points(path="ml_results/l3/extrinsic_quenching.json"):
    return [r for r in json.load(open(path)) if r["family"] not in FLUORIDE
            and WINDOW[0] <= r["E_cm1"] <= WINDOW[1]]


def evaluate(path="ml_results/l3/extrinsic_quenching.json", out="ml_results/l3/quench_model.json"):
    pts = load_points(path)
    A, p = _fit(pts)
    PARAMS.update(A=A, p=p)
    loo = []
    for i, q in enumerate(pts):
        a, pp = _fit(pts[:i] + pts[i + 1:])
        w = rate(q["family"], q["E_cm1"], q["N_cm3"], a, pp)
        tau = 1.0 / (w + 1e3 / q["tau_model_ms"])
        loo.append(dict(ion=q["ion"], level=q["level"], N_cm3=q["N_cm3"], W_meas=q["W_extra_s1"], W_pred=w,
                        tau_obs_ms=q["tau_obs_ms"], tau_pred_ms=tau * 1e3, tau_single_ion_ms=q["tau_model_ms"], src=q["src"]))
    tm = [q for q in pts if q["ion"] == "Tm3+"]; ho = [q for q in pts if q["ion"] == "Ho3+"]
    cross = []
    if tm and ho:
        a, pp = _fit(tm)
        for q in ho:
            w = rate(q["family"], q["E_cm1"], q["N_cm3"], a, pp)
            cross.append(dict(fit_on="Tm3+ only", predict=f"{q['ion']} {q['level']}", N_cm3=q["N_cm3"],
                              W_meas=q["W_extra_s1"], W_pred=w, tau_obs_ms=q["tau_obs_ms"],
                              tau_pred_ms=1e3 / (w + 1e3 / q["tau_model_ms"]), tau_single_ion_ms=q["tau_model_ms"]))
    A_free, p_free = _fit(pts, p_fixed=None)
    res = dict(A_s1=A, p=p, window_cm1=WINDOW, n_points=len(pts), leave_one_out=loo, cross_ion=cross,
               free_exponent_fit=dict(A_s1=A_free, p=p_free),
               external_test="ml_results/l3/simpson_test.json (Simpson thesis, 9 fibres, N 7e17-1.3e19): "
                             "rms log error 1.32 with p = 0.43 vs 0.18 with p = 0")
    json.dump(res, open(out, "w"), indent=1, default=float)
    print(f"W_q = {A:.0f} (N/1e20)^{p:.2f} s^-1 for silica-family levels at {WINDOW[0]:.0f}-{WINDOW[1]:.0f} cm^-1")
    for r in loo:
        print(f"  LOO {r['ion']} {r['level']} N {r['N_cm3']:.1e}: tau obs {r['tau_obs_ms']:.3f}  pred {r['tau_pred_ms']:.3f}  "
              f"(single-ion {r['tau_single_ion_ms']:.2f}) ms")
    for r in cross:
        print(f"  Tm-only -> {r['predict']}: tau obs {r['tau_obs_ms']:.3f}  pred {r['tau_pred_ms']:.3f}  "
              f"(single-ion {r['tau_single_ion_ms']:.2f}) ms")
    return res


if __name__ == "__main__":
    evaluate()
