"""Check the energy-gap-law parameters (gainmat/refs/ln_refs.py GAP_LAW) against measured
multiphonon rates (gainmat/refs/multiphonon_refs.py MULTIPHONON_DATA, gaps as stated in the sources).
    python scripts/check_gap_law.py  ->  forward_results/gap_law_check.json"""
import json, os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np
from gainmat.refs.ln_refs import GAP_LAW
from gainmat.refs.multiphonon_refs import MULTIPHONON_DATA

MAP = {"silicate": "silicate", "phosphate": "phosphate", "germanate": "germanate",
       "tellurite": "tellurite", "silica": "silicate"}
rows = []
for d in MULTIPHONON_DATA:
    fam = MAP.get(d["glass"])
    if fam is None or abs(d.get("T_K", 295) - 295) > 5 or not d.get("W_s1"):
        continue
    p = GAP_LAW[fam]
    w = p["B"] * np.exp(-p["alpha"] * d["dE_cm1"])
    rows.append(dict(ion=d["ion"], glass=d["glass"], level=d["level"], dE=d["dE_cm1"],
                     W_meas=d["W_s1"], W_law=float(w), ratio=float(w / d["W_s1"]), src=d["src"]))
r = np.array([x["ratio"] for x in rows])
lay = np.array([x["ratio"] for x in rows if x["glass"] != "silica"])
out = dict(n=len(rows), rms_log10_all=float(np.sqrt(np.mean(np.log10(r) ** 2))),
           rms_log10_layne=float(np.sqrt(np.mean(np.log10(lay) ** 2))), rows=rows)
# temperature dependence (Layne Fig. 9 series): W(T)/W(~300 K), measured vs Eq. 11 factor
from gainmat.ln4f.nonrad import temperature_factor
tser = {}
for d in MULTIPHONON_DATA:
    if "Fig. 9" in d["src"]:
        tser.setdefault(d["glass"], []).append((d["T_K"], d["W_s1"], d["dE_cm1"]))
trows = []
for g, pts in tser.items():
    pts.sort(); T0, W0, dE = pts[0]
    hw = GAP_LAW[MAP[g]]["hw"]
    for T, W, _ in pts[1:]:
        trows.append(dict(glass=g, T=T, ratio_meas=W / W0,
                          ratio_model=temperature_factor(dE, hw, T, T0)))
out["temperature"] = trows
os.makedirs("forward_results", exist_ok=True)
json.dump(out, open("forward_results/gap_law_check.json", "w"), indent=1)
for x in rows:
    print(f'{x["ion"]:5s} {x["glass"]:9s} {x["level"]:7s} dE {x["dE"]:6.0f}  meas {x["W_meas"]:9.3g}  law {x["W_law"]:9.3g}  x{x["ratio"]:.2f}')
for t in trows:
    print(f'  T-dependence {t["glass"]:9s} {t["T"]:4.0f} K  meas x{t["ratio_meas"]:.2f}  Eq.11 x{t["ratio_model"]:.2f}')
print("rms log10 error: all %.2f, Layne glasses %.2f (n=%d)" % (out["rms_log10_all"], out["rms_log10_layne"], len(rows)))
