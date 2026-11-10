"""Model Stark structure in glass vs measured levels collected by extraction agents from open-access
papers (gainmat/refs/data/stark_records.jsonl): lowest level of a multiplet above the lowest ground
level (E0) and total splitting, for the glass records only (crystal records are listed but skipped).
Model = crystal-field ensemble with anchors switched off (pure prediction).

    python scripts/validate_stark_agent.py -> forward_results/stark_agent_validation.json
"""
import json, os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np
from gainmat.forward import glass as g
from gainmat.refs.oa_agent_data import load_stark

g.load_refs(); g.apply_site_refs(); g.calibrate_spm(370.0, g.CARNALL_ER_LAF3_BKQ, verbose=False)
GROUND = {"Er3+": "4I15/2", "Yb3+": "2F7/2", "Tm3+": "3H6", "Ho3+": "5I8", "Nd3+": "4I9/2", "Pr3+": "3H4", "Dy3+": "6H15/2"}
_sp = {}


def spec(ion, fam):
    if (ion, fam) not in _sp:
        _sp[(ion, fam)] = g.run(ion, fam, n_sites=60, pumps=[], overrides=dict(level_anchor=None))["_spec"]
    return _sp[(ion, fam)]


rows, skipped = [], []
for r in load_stark():
    ion, fam, mult = r.get("ion"), r.get("family"), r.get("multiplet")
    if fam not in g.FAMILIES or ion not in GROUND:
        skipped.append(f"{ion} {mult} {fam}"); continue
    try:
        sp = spec(ion, fam)
        lo, hi = sp["edges"][mult]; g0 = sp["edges"][GROUND[ion]][0]
    except Exception as e:
        skipped.append(f"{ion} {mult} {fam}: {e!r}"[:100]); continue
    out = dict(ion=ion, family=fam, multiplet=mult, src=str(r.get("src"))[:120], T_K=r.get("T_K"))
    if r.get("E0_cm1") not in (None, 0, "0") and mult != GROUND[ion]:
        out.update(E0_meas=float(r["E0_cm1"]), E0_model=round(lo - g0, 1), E0_err=round(lo - g0 - float(r["E0_cm1"]), 1))
    ts = r.get("total_splitting_cm1")
    if ts is None and r.get("levels_cm1"):
        lv = [x for x in r["levels_cm1"] if x is not None]
        ts = max(lv) if len(lv) >= 2 else None
    if ts:
        out.update(split_meas=float(ts), split_model=round(hi - lo, 1), split_ratio=round((hi - lo) / float(ts), 2))
    if len(out) > 5:
        rows.append(out); print(json.dumps(out), flush=True)
e0 = [abs(r["E0_err"]) for r in rows if "E0_err" in r]
sr = [r["split_ratio"] for r in rows if "split_ratio" in r]
summ = dict(n_E0=len(e0), E0_median_abs_err_cm1=float(np.median(e0)) if e0 else None,
            n_split=len(sr), split_ratio_median=float(np.median(sr)) if sr else None,
            split_within_25pct=float(np.mean(np.abs(np.log(sr)) <= np.log(1.25))) if sr else None, skipped=skipped)
print(json.dumps(summ, indent=1))
json.dump(dict(summary=summ, rows=rows), open("forward_results/stark_agent_validation.json", "w"), indent=1)
