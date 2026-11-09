"""Compare Er3+ lineshapes from MLIP-MD dopant sites against the stochastic site model.

    python scripts/compare_md_sites.py md_sites/silica.json silica [md_sites/ZBLAN.json ZBLAN ...]
"""
import json
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np

from gainmat.forward.glass import (load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ,
                                   FAMILIES, load_md_sites, run, key_transitions)

load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ, verbose=False)
args = sys.argv[1:]
out = {}
for path, fam in zip(args[0::2], args[1::2]):
    sites = load_md_sites(path)
    if not sites:
        print(path, "no sites"); continue
    cn = [len(s) for s in sites]
    rr = [np.linalg.norm(l[1:]) for s in sites for l in s]
    res = {}
    for label, ov in (("stochastic", {}), ("md", {"sites": sites})):
        r = run("Er3+", fam, n_sites=30, overrides=ov)
        kt = key_transitions(r, [("4I13/2", "4I15/2"), ("4I11/2", "4I13/2")])
        res[label] = dict(split_4I15=r["split"].get("4I15/2"), split_4I13=r["split"].get("4I13/2"),
                          windows=[(w["lo_nm"], w["hi_nm"]) for w in r["windows"]], transitions=kt)
    out[fam] = dict(n_sites=len(sites), CN_mean=float(np.mean(cn)), R_mean=float(np.mean(rr)),
                    R_std=float(np.std(rr)), **res)
    print(json.dumps({fam: out[fam]}, indent=1, default=float))
json.dump(out, open("forward_results/md_vs_stochastic.json", "w"), indent=1, default=float)
