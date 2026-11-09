"""Summarise the UQ grid (ml_results/uq/*.json): robust (P >= 0.8) and possible (P >= 0.2) gain
windows per (ion, glass), and the probability inside each experimentally known band.

    python -m gainmat.ml.uq_summary  ->  ml_results/uq/SUMMARY.md, summary.json
"""
import glob
import json
import os

import numpy as np


def main(d="ml_results/uq"):
    rows = []
    for p in sorted(glob.glob(os.path.join(d, "*_*.json"))):
        if os.path.basename(p) == "summary.json":
            continue
        r = json.load(open(p))
        lam = np.array(r["lam_nm"]); P = np.array(r["P"])
        rob = r["bands"]["robust"]; pos = r["bands"]["possible"]
        # 'contested' = possible but not robust: where the uncertainty decides
        cont_nm = float(np.sum((P >= 0.2) & (P < 0.8)) * np.median(np.abs(np.diff(lam))))
        rows.append(dict(ion=r["ion"], family=r["family"], omega_src=r["omega_src"].split(" ")[0], n_theta=r["n_theta"],
                         robust=rob, possible=pos, contested_nm=cont_nm))
    known = []
    V = json.load(open("forward_results/validation.json"))
    for v in V:
        if "band_exp_nm" not in v:
            continue
        p = os.path.join(d, f"{v['ion'].replace('+', '')}_{v['family']}.json")
        if not os.path.exists(p):
            continue
        r = json.load(open(p)); lam = np.array(r["lam_nm"]); P = np.array(r["P"])
        lo, hi = v["band_exp_nm"]; m = (lam >= lo) & (lam <= hi)
        known.append(dict(entry=v["entry"], max_P=float(P[m].max()), frac_P50=float(np.mean(P[m] >= 0.5))))
    json.dump(dict(grid=rows, known_bands=known), open(os.path.join(d, "summary.json"), "w"), indent=1)
    L = ["# Gain-window probabilities", "",
         "P = fraction of Monte Carlo draws (Omega spread x host-theta spread from the L2 folds, ML layer on) with "
         "gain - background loss >= 0.1 dB/m under the family's standard pumps.  robust: P >= 0.8; possible: P >= 0.2.", "",
         "## Known fibre media", "", "| medium | max P in measured band | fraction of band with P >= 0.5 |", "|---|---|---|"]
    L += [f"| {k['entry']} | {k['max_P']:.2f} | {k['frac_P50']:.2f} |" for k in known]
    L += ["", "## Grid", "", "| ion | glass | Omega | theta draws | robust windows (nm) | possible windows (nm) | width decided by uncertainty (nm) |",
          "|---|---|---|---|---|---|---|"]
    fmt = lambda b: ", ".join(f"{a:.0f}-{c:.0f}" for a, c in b) or "-"
    for r in rows:
        L.append(f"| {r['ion']} | {r['family']} | {r['omega_src']} | {r['n_theta']} | {fmt(r['robust'])} | {fmt(r['possible'])} | {r['contested_nm']:.0f} |")
    open(os.path.join(d, "SUMMARY.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
