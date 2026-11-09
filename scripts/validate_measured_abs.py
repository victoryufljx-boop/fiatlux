"""Measured-absorption mode on devices: lasers rebuilt with the band's measured absorption
(emission by McCumber with the model's eps) vs the composition-only model.

    python scripts/validate_measured_abs.py -> forward_results/measured_abs_devices.json
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np

from gainmat.forward.glass import load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ
from gainmat.forward.fiber import local_model_for, slope_and_threshold
from gainmat.ml.observations import build

load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ, verbose=False)
obs = build()


def curve(ion, fam, src_key, which="abs"):
    o = next(o for o in obs if o["kind"] == "curve" and o["which"] == which and o["ion"] == ion
             and o["family"] == fam and src_key in o["src"])
    return o["lam_nm"], o["sigma"]


rows = []
# Le Gouet 2020 Ho:silica laser (as in validate_devices.py)
ma = {("5I7", "5I8"): curve("Ho3+", "aluminosilicate", "Le Gouet")}
me = {("5I7", "5I8"): curve("Ho3+", "aluminosilicate", "Le Gouet") + (curve("Ho3+", "aluminosilicate", "Le Gouet", "em"),)}
N = 5.3e19; core = 20.0; Gp = np.pi * (core / 2) ** 2 / (0.8284 * 250.0 ** 2)
for var, ov in (("model only", None), ("measured absorption", dict(measured_absorption=ma)),
                ("measured absorption + emission", dict(measured_absorption=me))):
    t0 = time.time()
    lm = local_model_for("Ho3+", "aluminosilicate", n_sites=30, overrides=ov)
    s, th, cv = slope_and_threshold(lm, 16, core, 0.08, N, 1940, [50, 100, 150, 200], 2120,
                                    R1=0.99, R2=0.03, alpha_bg_dB_m=0.04, nz=40, pump_overlap=Gp)
    rows.append(dict(case="Ho:silica laser (Le Gouet 2020)", variant=var, model_slope=s, measured_slope=0.49,
                     model_threshold_W=th, measured_threshold_W=13, curve=cv, seconds=time.time() - t0))
    print(json.dumps(rows[-1], default=float), flush=True)
# Pask 1995 Yb:germanosilicate lasers
ma = {("2F5/2", "2F7/2"): curve("Yb3+", "germanosilicate", "Pask")}
me = {("2F5/2", "2F7/2"): curve("Yb3+", "germanosilicate", "Pask") + (curve("Yb3+", "germanosilicate", "Pask", "em"),)}
N = 550e-6 * 2.2e22
for var, ov in (("model only", None), ("measured absorption", dict(measured_absorption=ma)),
                ("measured absorption + emission", dict(measured_absorption=me))):
    lm = local_model_for("Yb3+", "germanosilicate", n_sites=30, overrides=ov)
    t0 = time.time()
    s, th, cv = slope_and_threshold(lm, 90, 3.75, 0.17, N, 1047, [0.05, 0.1, 0.2, 0.3], 1102,
                                    R1=0.98, R2=0.04, alpha_bg_dB_m=0.005, nz=60)
    rows.append(dict(case="Yb:germanosilicate 1102 nm laser (Pask 1995)", variant=var, model_slope=s, measured_slope=0.90,
                     model_threshold_W=th, measured_threshold_W=0.030, curve=cv, seconds=time.time() - t0))
    print(json.dumps(rows[-1], default=float), flush=True)
    t0 = time.time()
    s, th, cv = slope_and_threshold(lm, 5, 3.75, 0.17, N, 840, [0.08, 0.15, 0.25, 0.35], 1040,
                                    R1=0.04, R2=0.04, alpha_bg_dB_m=0.005, nz=40)
    rows.append(dict(case="Yb:germanosilicate free-running 5 m (Pask 1995)", variant=var, model_slope=s,
                     model_threshold_W=th, measured_threshold_W=0.045, curve=cv, seconds=time.time() - t0))
    print(json.dumps(rows[-1], default=float), flush=True)
json.dump(rows, open("forward_results/measured_abs_devices.json", "w"), indent=1, default=float)
