"""Summarise bake-off JSONs into markdown tables.

    python -m gainmat.host.bakeoff_report results/bakeoff > table.md
"""
from __future__ import annotations

import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np

from gainmat.refs.phonon_refs import PHONON_REFS, match_modes

HOSTS = ["CaF2", "Al2O3", "SiO2-alpha-quartz", "LiYF4"]
SHORT = {"SiO2-alpha-quartz": "quartz"}


def load(d):
    rows = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(d, "*.json")))]
    return [r for r in rows if "modes" in r]


def score(r):
    ref = PHONON_REFS[r["host"]]
    exp = defaultdict(list)
    for k in ("raman", "ir_to"):
        for irr, fl in ref[k].items():
            exp[irr] += fl
    pairs = match_modes(r["modes"], exp)
    err = [c - e for _, e, c in pairs if c is not None]
    nmiss = sum(c is None for _, _, c in pairs)
    out = dict(mae=float(np.mean(np.abs(err))) if err else None,
               mse=float(np.mean(err)) if err else None, n=len(err), nmiss=nmiss,
               max_to=max(m["freq_cm1"] for m in r["modes"]),
               exp_max_to=max(max(v) for v in exp.values()), pairs=pairs)
    if "nac_gamma" in r and ref["ir_lo"]:
        lo_exp = max(max(v) for v in ref["ir_lo"].values())
        lo_calc = max(max(v) for v in r["nac_gamma"].values())
        out.update(lo_calc=lo_calc, lo_exp=lo_exp)
    return out


def main(d):
    rows = load(d)
    models = sorted({r["model"] for r in rows})
    print("| model | protocol | " + " | ".join(
        f"{SHORT.get(h, h)} dV% / MAE / maxTO" for h in HOSTS) + " | mean MAE | mean bias | t (s) |")
    print("|---|---|" + "---|" * (len(HOSTS) + 3))
    summary = {}
    for m in models:
        for prot in ("relaxed", "expcell"):
            cells, maes, biases, ts = [], [], [], []
            for h in HOSTS:
                r = next((x for x in rows if x["model"] == m and x["host"] == h
                          and x["protocol"] == prot), None)
                if r is None:
                    cells.append("–"); continue
                s = score(r)
                maes.append(s["mae"]); biases.append(s["mse"]); ts.append(r["t_s"])
                cells.append(f"{r['vol_pct']:+.1f} / {s['mae']:.1f} / {s['max_to']:.0f}")
            if maes:
                summary[(m, prot)] = (np.mean(maes), len(maes))
                print(f"| {m} | {prot} | " + " | ".join(cells) +
                      f" | **{np.mean(maes):.1f}** | {np.mean(biases):+.1f} | {sum(ts):.0f} |")
    exp_line = " | ".join(f"exp maxTO {max(max(v) for k in ('raman','ir_to') for v in PHONON_REFS[h][k].values()):.0f}" for h in HOSTS)
    print(f"| experiment | | {exp_line} | | | |")
    print()
    print("LO (NAC, Born charges + eps_inf from refs) — highest Gamma LO, cm^-1")
    print()
    print("| model | protocol | CaF2 LO (exp 463) | Al2O3 LO (exp 906.6) |")
    print("|---|---|---|---|")
    for m in models:
        for prot in ("relaxed", "expcell"):
            c = []
            for h in ("CaF2", "Al2O3"):
                r = next((x for x in rows if x["model"] == m and x["host"] == h
                          and x["protocol"] == prot), None)
                s = score(r) if r else {}
                c.append(f"{s['lo_calc']:.0f}" if "lo_calc" in s else "–")
            print(f"| {m} | {prot} | " + " | ".join(c) + " |")
    print()
    best = min(summary, key=lambda k: summary[k][0])
    print(f"Best (lowest mean TO+Raman MAE): {best[0]} / {best[1]} "
          f"({summary[best][0]:.1f} cm^-1 over {summary[best][1]} hosts)")
    return summary


def detail(d, model, prot):
    for r in load(d):
        if r["model"] == model and r["protocol"] == prot:
            s = score(r)
            print(r["host"], [(i, e, c) for i, e, c in s["pairs"]])


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "results/bakeoff")
