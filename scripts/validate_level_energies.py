"""Multiplet energies (site-averaged Stark barycentres) vs the levels measured in ZBLAN for six
ions (Wetenkamp, West & Toebben 1992, Fig. 1, digitised; +-50-100 cm^-1).

    python scripts/validate_level_energies.py -> forward_results/level_energy_validation.json
"""
import json
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np

from gainmat.forward.glass import load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ, run
from gainmat.refs.zblan_yb_refs import WETENKAMP1992 as W

load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ, verbose=False)
key = next(k for k, v in W.items() if isinstance(v, dict) and isinstance(v.get("Tm3+"), dict)
           and isinstance(v["Tm3+"].get("3F4"), (int, float)))
rows = []
for ion, levels in W[key].items():
    if not ion.endswith("3+"):
        continue
    b = run(ion, "ZBLAN", n_sites=16, pumps=[])["_spec"]["bary"]
    for lab, E in levels.items():
        mb = [b[x.strip()] for x in lab.split(",") if x.strip() in b]
        if mb:
            m = float(np.mean(mb))
            rows.append(dict(ion=ion, level=lab, measured_cm1=E, model_cm1=round(m, 1), diff_cm1=round(m - E, 1),
                             diff_pct=round((m - E) / E * 100, 2)))
d = np.array([r["diff_cm1"] for r in rows]); p = np.array([r["diff_pct"] for r in rows])
summ = dict(n=len(rows), median_abs_cm1=float(np.median(np.abs(d))), median_abs_pct=float(np.median(np.abs(p))),
            within_1pct=float(np.mean(np.abs(p) <= 1)), within_2pct=float(np.mean(np.abs(p) <= 2)),
            per_ion_mean_pct={i: float(np.mean([r["diff_pct"] for r in rows if r["ion"] == i])) for i in sorted({r["ion"] for r in rows})})
json.dump(dict(summary=summ, rows=rows, src="Wetenkamp, West & Toebben, JNCS 140, 35 (1992) Fig. 1 (digitised)"),
          open("forward_results/level_energy_validation.json", "w"), indent=1)
print(json.dumps(summ, indent=1))
