"""Concentration self-quenching predicted from the model's spectra (gainmat.forward.quenching) vs measured
lifetime-vs-concentration series collected from open-access papers (gainmat/refs/data/et_macro_records.jsonl).

Test case with a cross-relaxation channel in the model: Tm 3H4 in 75TeO2-20ZnO-5Na2O tellurite
(Taher, PhD thesis, Swansea 2011; 0.82-9.06e20 cm^-3).  Two variants:
  predictive : tau0 = the model's own low-concentration lifetime
  given tau0 : tau0 = the measured lifetime of the most dilute sample (tests the transfer physics alone)
Series without a self-quenching channel in the model (Er 4I13/2, Tm 3F4: their quenching runs through
migration to OH / defects, which needs an impurity concentration) are listed with the model's flat
prediction, to show what the model leaves out.

    python scripts/validate_quenching.py -> forward_results/quenching_validation.json
"""
import json, os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np
from gainmat.forward import glass as g
from gainmat.forward.quenching import SELF_QUENCH, lifetime_vs_N

g.load_refs(); g.apply_site_refs(); g.calibrate_spm(370.0, g.CARNALL_ER_LAF3_BKQ, verbose=False)
recs = [json.loads(l) for l in open("gainmat/refs/data/et_macro_records.jsonl")]
series = {}
for r in recs:
    if r["quantity"] != "lifetime_ms" or not r.get("conc_cm3") or r.get("proxy"):
        continue
    k = (r["ion_system"], r["family"], r["level"], r["src"].split(",")[0])
    series.setdefault(k, []).append((float(r["conc_cm3"]), float(r["value"])))
out = []
for (ion, fam, lev, src), pts in series.items():
    pts = sorted(pts)
    if len(pts) < 3:
        continue
    N = [p[0] for p in pts]; tm = np.array([p[1] for p in pts])
    row = dict(ion=ion, family=fam, level=lev, src=src, N_cm3=N, tau_meas_ms=tm.tolist(),
               self_quenching_channel=lev in SELF_QUENCH.get(ion, {}))
    if row["self_quenching_channel"]:
        for var, t0 in (("predictive", None), ("given tau0", tm[0] * 1e-3)):
            q = lifetime_vs_N(ion, fam, N, level=lev, tau0_s=t0)
            tp = np.array([p["tau_mean_ms"] for p in q["points"]])
            row[var] = dict(tau_model_ms=tp.round(4).tolist(), tau0_ms=q["tau0_ms"], C_CR=q["C_CR_cm6_s"], C_mig=q["C_mig_cm6_s"],
                            log10_ratio=np.log10(tp / tm).round(3).tolist(),
                            median_abs_dex=float(np.median(np.abs(np.log10(tp / tm)))))
            # quenching rate beyond the dilute sample: measured vs model
            row[var]["k_q_meas_s"] = (1e3 / tm - 1e3 / tm[0]).round(1).tolist()
            row[var]["k_q_model_s"] = (1e3 / tp - 1e3 / tp[0]).round(1).tolist()
    else:
        row["note"] = "no self-quenching channel in the model: predicted lifetime independent of N"
    out.append(row); print(json.dumps(row), flush=True)
json.dump(out, open("forward_results/quenching_validation.json", "w"), indent=1, default=float)
