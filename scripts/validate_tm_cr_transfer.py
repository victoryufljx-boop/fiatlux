"""Tm:silica lasers (Jackson & Mossman 2003) with the pump band taken from measurements, and the
3H4,3H6 -> 3F4,3F4 cross-relaxation strength learned on one fibre and tested on the others.

The 3H6 -> 3H4 band is built from two published numbers: peak 8.6e-21 cm2 at 790 nm (Jackson &
King 1999) and ~4.8e-21 cm2 at 805 nm (Jackson & Mossman 2003 core absorption / Tm density); the
model's own band is right at the peak but too narrow at the 805 nm pump.

    python scripts/validate_tm_cr_transfer.py -> forward_results/tm_cr_transfer.json
"""
import sys, json, numpy as np
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from gainmat.forward.glass import load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ
from gainmat.forward.fiber import local_model_for, slope_and_threshold
from gainmat.refs.tm_silica_refs import JACKSON2003 as J
from gainmat.ml.quench import extra_rates
load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ, verbose=False)
nu0 = 1e7 / 790; sd = (1e7 / 790 - 1e7 / 805) / np.sqrt(-2 * np.log(4.8 / 8.6))
nu = np.linspace(nu0 - 5 * sd, nu0 + 5 * sd, 300)
band = (1e7 / nu, 8.6e-21 * np.exp(-0.5 * ((nu - nu0) / sd) ** 2))
ov = dict(measured_absorption={("3H4", "3H6"): band})
fib = {f["fibre"]: f for f in J["fibres"]}; fig5 = J["fig5_output_vs_launched_pump_digitized"]
cases = {4: ("fibre4_1.6wt", 0.46, 5.24), 3: ("fibre3_1.3wt", 0.415, 5.99), 1: ("fibre1_0.7wt", 0.213, 7.41)}
lm0 = local_model_for("Tm3+", "aluminosilicate", n_sites=30, overrides=ov)
def sim(k, c):
    key, s_meas, th_meas = cases[k]; f = fib[k]; L = fig5[key]["L_m_est"]; N = f["N_Tm_cm3_from_alpha_core"]
    core = min(max(f["core_diam_um_est_from_area_ratio"], 20.0), 25.0)
    Gp = (core / 300.0) ** 2 * 4 / np.pi / (1 - 0.1); k3101 = c * 1.8e-16 * N / 4.5e20
    et = [("3H4", "3F4", "3H6", "3F4", k3101), ("3F4", "3H4", "3F4", "3H6", 0.084 * k3101)]
    q = extra_rates("aluminosilicate", lm0.spec["bary"], N)
    lm = local_model_for("Tm3+", "aluminosilicate", n_sites=30, quench=q, energy_transfer=et, overrides=ov)
    sl, th, cv = slope_and_threshold(lm, L, core, 0.19, N, 805, [6, 9, 12, 15, 18, 22], 1990, R1=0.99, R2=0.035,
                                     alpha_bg_dB_m=0.1, nz=40, pump_overlap=Gp)
    return sl, th, s_meas, th_meas
res = {}; out = []
for c in (1, 2, 4, 8):
    sl, th, sm, tm = sim(4, c); res[c] = (sl, th)
    err = (np.log(sl / sm)) ** 2 + (np.log(th / tm)) ** 2 if sl and th else 99
    out.append(dict(fit_on="fibre 4", cr_scale=c, slope=sl, threshold_W=th, measured_slope=sm, measured_threshold_W=tm))
    print(f"fit on fibre 4: CR x{c}: slope {sl:.3f} ({sm}) threshold {th:.2f} ({tm})  err {err:.3f}", flush=True)
best = min(res, key=lambda c: (np.log(res[c][0] / 0.46)) ** 2 + (np.log(res[c][1] / 5.24)) ** 2)
print("best CR scale from fibre 4:", best)
for k in (3, 1):
    for c in (1, best):
        sl, th, sm, tm = sim(k, c)
        out.append(dict(predict=f"fibre {k}", cr_scale=c, slope=sl, threshold_W=th, measured_slope=sm, measured_threshold_W=tm))
        print(f"  predict fibre {k} with CR x{c}: slope {sl:.3f} ({sm}) threshold {th:.2f} ({tm})", flush=True)
json.dump(dict(best_cr_scale_from_fibre4=best, rows=out, pump_band="Gaussian in energy: 8.6e-21 at 790 nm, 4.8e-21 at 805 nm"),
          open("forward_results/tm_cr_transfer.json", "w"), indent=1, default=float)
