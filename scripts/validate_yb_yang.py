"""Yb3+ level patterns in 8 glasses (Yang et al., OL 39, 1772 (2014)) vs the crystal-field model, and the
Yang phosphate pattern as an anchor tested on Wang et al. 2017 phosphate cross sections at 975/1005 nm.

    python scripts/validate_yb_yang.py -> forward_results/yb_yang_validation.json
"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from gainmat.forward import glass as g
from gainmat.forward.spectra import random_sites, spm_bkq
from gainmat.forward.ions import get_ion
from gainmat.refs.yang2014_yb_refs import YANG2014
from gainmat.refs.stark_anchor_refs import ANCHORS
g.load_refs(); g.apply_site_refs(); g.calibrate_spm(370.0, g.CARNALL_ER_LAF3_BKQ, verbose=False)
m = get_ion("Yb3+"); out = dict(patterns=[], wang=[])
for ml in (False, True):
    for k, v in YANG2014.items():
        fam = v["family"]
        if not fam:
            continue
        F = g.FAMILIES[fam]; sites = random_sites(60, F["ligand"], F["CN"], F["R"], F["R_std"], F["ang_std"], seed=0)
        sc = (g.cf_family_scale(fam) or 1.0) if ml else 1.0
        p = dict(g.spm_for("Yb3+", F["ligand"])); p = dict(p, **{f"B{j}": p[f"B{j}"] * sc for j in (2, 4, 6)})
        a7, a5 = [], []
        for s in sites:
            E, W, par = m.stark(spm_bkq(s, {F["ligand"]: p}))
            e7 = np.sort(E[par == m.mult_index("2F7/2")])[::2]; e5 = np.sort(E[par == m.mult_index("2F5/2")])[::2]
            a7.append(e7 - e7[0]); a5.append(e5 - e5[0])
        a7, a5 = np.mean(a7, 0), np.mean(a5, 0)
        m7 = np.array(v["f72"], float); m5 = np.array(v["f52"], float) - v["f52"][0]
        out["patterns"].append(dict(glass=k, family=fam, proxy=v["proxy"], ml_family_scale=ml, scale=sc,
                                    model_2F72=a7.round(0).tolist(), meas_2F72=m7.tolist(), model_2F52=a5.round(0).tolist(), meas_2F52=m5.tolist(),
                                    shape_b_over_d=[float(a7[1] / a7[3]), float(m7[1] / m7[3])], shape_c_over_d=[float(a7[2] / a7[3]), float(m7[2] / m7[3])],
                                    rms_after_scale_cm1=float(np.sqrt(np.mean((a7 / a7[3] * m7[3] - m7) ** 2)))))
p = ("2F5/2", "2F7/2")
for tag, ov in (("model", dict(level_anchor=None)), ("ZPL anchor (default)", {}), ("Yang phosphate anchor", dict(level_anchor=ANCHORS["Yb3+:phosphate"]["levels"]))):
    sp = g.run("Yb3+", "phosphate", pumps=[], n_sites=60, overrides=ov)["_spec"]
    f = lambda k, l: float(np.interp(1e7 / l, sp["nu"], sp["sigma_" + k][p]))
    out["wang"].append(dict(variant=tag, sa975=f("abs", 975), se975=f("em", 975), sa1005=f("abs", 1005), se1005=f("em", 1005),
                            em_ratio=f("em", 1005) / f("em", 975), abs_ratio=f("abs", 1005) / f("abs", 975)))
out["wang_measured"] = dict(src="Wang et al., Materials 10, 241 (2017) Table 1, BPS9 (71 P2O5)", sa975=1.09e-20, se975=1.46e-20,
                            sa1005=0.15e-20, se1005=0.86e-20, em_ratio=0.86 / 1.46, abs_ratio=0.15 / 1.09)
json.dump(out, open("forward_results/yb_yang_validation.json", "w"), indent=1, default=float)
for r in out["patterns"]:
    print(r["glass"], r["family"], "ML" if r["ml_family_scale"] else "  ", "b/d", np.round(r["shape_b_over_d"], 2).tolist(), "c/d", np.round(r["shape_c_over_d"], 2).tolist(), "rms after scale", round(r["rms_after_scale_cm1"]))
for r in out["wang"]:
    print(r["variant"], "em ratio", round(r["em_ratio"], 2), "abs ratio", round(r["abs_ratio"], 3), "(meas 0.59 / 0.138)")
