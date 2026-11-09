"""Inter-ion energy transfer in Er:ZBLAN and Er,Pr:ZBLAN (Golding et al., PRB 62, 856 (2000)).

1. Effective Er 4I13/2 and 4I11/2 lifetimes with Pr co-doping: model intrinsic rates +
   Golding's Er->Pr coefficients, vs Golding Fig. 6(a).  (The coefficients were derived from
   these lifetimes with the measured intrinsic lifetimes, so this checks the model's intrinsic
   rates and the rate-law bookkeeping, not the coefficients.)
2. What the energy transfer does to the 2.8 um Er line under 975 nm CW pumping:
   single-ion vs + Er-Er upconversion / cross-relaxation vs + Pr deactivation.

    python scripts/validate_codoping.py  ->  forward_results/codoping_validation.json
"""
import json
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np

from gainmat.forward.glass import load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ
from gainmat.forward.fiber import local_model_for, CoDopedModel
from gainmat.refs.er_pr_zblan_refs import GOLDING2000 as G, ENERGY_TRANSFER as ET

load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ, verbose=False)
er = local_model_for("Er3+", "ZBLAN", n_sites=30)
pr = local_model_for("Pr3+", "ZBLAN", n_sites=30)
k = {e["name"].split()[0] + ("13" if "4I13/2 ->" in e["name"] else "11" if "4I11/2 ->" in e["name"] else ""): e
     for e in ET}
et_pr13 = next(e for e in ET if e["name"].startswith("ET Er 4I13/2"))["k_cm3_s"]
et_pr11 = next(e for e in ET if e["name"].startswith("ET Er 4I11/2"))["k_cm3_s"]
meas = G["lifetimes_Er_Pr_codoped"]
rows = []
print("1. Er lifetimes with Pr (us): model / measured")
for key in meas["tau_eff_4I13_2_us"]:
    NEr, NPr = key
    m = CoDopedModel({"Er": (er, NEr), "Pr": (pr, NPr)},
                     et=[("Er", "4I13/2", "4I15/2", "Pr", "3H4", "3F3", et_pr13[key]),
                         ("Er", "4I11/2", "4I15/2", "Pr", "3H4", "1G4", et_pr11[key])])
    t13 = m.lifetimes("Er", "4I13/2") * 1e6; t11 = m.lifetimes("Er", "4I11/2") * 1e6
    r = dict(N_Er=NEr, N_Pr=NPr, tau13_model_us=t13, tau13_meas_us=meas["tau_eff_4I13_2_us"][key],
             tau11_model_us=t11, tau11_meas_us=meas["tau_eff_4I11_2_us"][key])
    rows.append(r)
    print(f"  Er {NEr:.1e} Pr {NPr:.1e}:  4I13/2 {t13:6.0f} / {r['tau13_meas_us']:5.0f}   "
          f"4I11/2 {t11:6.0f} / {r['tau11_meas_us']:5.0f}")

# 2. 2.8 um line under CW pumping
NEr, NPr = 8e20, 4e19
etu1 = next(e for e in ET if e["name"].startswith("ETU1"))["k_cm3_s"][8e20]
etu2 = next(e for e in ET if e["name"].startswith("ETU2"))["k_cm3_s"][8e20]
cr = next(e for e in ET if e["name"].startswith("CR") and e["k_cm3_s"])["k_cm3_s"][8e20]
er_er = [("Er", "4I13/2", "4I15/2", "Er", "4I13/2", "4I9/2", etu1),
         ("Er", "4I11/2", "4I15/2", "Er", "4I11/2", "4F7/2", etu2),
         ("Er", "4S3/2", "4I9/2", "Er", "4I15/2", "4I13/2", cr)]
er_pr = [("Er", "4I13/2", "4I15/2", "Pr", "3H4", "3F3", et_pr13[(NEr, NPr)]),
         ("Er", "4I11/2", "4I15/2", "Pr", "3H4", "1G4", et_pr11[(NEr, NPr)])]
variants = {"single-ion": (CoDopedModel({"Er": (er, NEr)}), None),
            "+ Er-Er ETU/CR": (CoDopedModel({"Er": (er, NEr)}, er_er), None),
            "+ Er-Er + Pr": (CoDopedModel({"Er": (er, NEr), "Pr": (pr, NPr)}, er_er + er_pr), None)}
print("2. Er:ZBLAN 8e20 cm^-3, 975 nm CW: gain at 2750 nm (dB/m) / n(4I13/2) / n(4I11/2)")
lab = er.labels
demo = {}
for name, (m, _) in variants.items():
    out = []
    for I in (1e4, 1e5, 3e5, 1e6):
        n = m.populations([(975.0, I)])
        g = m.gain_cm(2750.0, n) * 4.343 * 100
        out.append(dict(I_W_cm2=I, gain_dB_m=float(g), n13=float(n["Er"][lab.index("4I13/2")]),
                        n11=float(n["Er"][lab.index("4I11/2")])))
    demo[name] = out
    print(f"  {name:16s} " + "  ".join(f"{o['gain_dB_m']:7.2f}/{o['n13']:.2f}/{o['n11']:.2f}" for o in out))
json.dump(dict(lifetimes=rows, er_2p8um_cw=demo,
               src="Golding, Jackson, King, Pollnau, PRB 62, 856 (2000)"),
          open("forward_results/codoping_validation.json", "w"), indent=1, default=float)
print("wrote forward_results/codoping_validation.json")
