"""Explicit model vs explicit + ML layer, check by check.

    python scripts/compare_ml.py [--off forward_results] [--on forward_results_ml]
        -> forward_results_ml/COMPARISON.md, comparison.json

Both folders hold residuals.jsonl written by scripts/validation_scorecard.py (the ML-on folder
from a run with GAINMAT_ML=1: learned ZBLAN host theta + L3 extrinsic quenching).  Rows whose
measurement was used to *fit* an ML component are flagged in-sample and reported separately,
so the out-of-sample comparison is not flattered by training points.
"""
import argparse
import json
import os
from collections import defaultdict

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--off", default="forward_results")
ap.add_argument("--on", default="forward_results_ml")
a = ap.parse_args()


def load(d):
    out = {}
    for line in open(os.path.join(d, "residuals.jsonl")):
        r = json.loads(line)
        out[(r["quantity"], r["ion"], r["host"], r["item"], r.get("src", ""))] = r
    return out


# measurements that trained an ML component
l3 = json.load(open("ml_results/l3/extrinsic_quenching.json"))
L3_TRAIN = {(p["ion"], round(p["tau_obs_ms"], 3)) for p in l3 if p["family"] not in ("ZBLAN", "fluoroindate")}
# L2 (ZBLAN theta) was fitted to band widths / peak-normalised shapes from these sources; an
# absolute cross section from the same spectra is a different number but not an independent one
L2_SOURCES = ("Wetenkamp", "Jackson, Nat. Photon")


def cls(q):
    q = q.split("(")[0]
    return {"lifetime:tau": "lifetime", "tau_ms": "lifetime", "lifetime:tau_rad": "radiative lifetime",
            "cross_section:em": "peak cross section", "cross_section:abs": "peak cross section",
            "spectrum_integrated": "integrated cross section", "vendor_xs_integrated": "integrated cross section",
            "stark_cm1": "Stark splitting"}.get(q, q)


def subset(r):
    if cls(r["quantity"]) == "lifetime" and r.get("measured") is not None and r["host"] not in ("ZBLAN", "fluoroindate"):
        if (r["ion"], round(float(r["measured"]), 3)) in L3_TRAIN:
            return "in-sample (trained L3)"
    if r["host"] == "ZBLAN" and any(k in r.get("src", "") for k in L2_SOURCES) and "cross" in cls(r["quantity"]):
        return "same source as L2 widths"
    return "independent"


off, on = load(a.off), load(a.on)
rows = []
for k, r0 in off.items():
    r1 = on.get(k)
    if r1 is None or r0.get("log10_ratio") is None or r1.get("log10_ratio") is None:
        continue
    if not (np.isfinite(r0["log10_ratio"]) and np.isfinite(r1["log10_ratio"])):
        continue
    rows.append(dict(cls=cls(r0["quantity"]), ion=r0["ion"], host=r0["host"], item=r0["item"],
                     measured=r0.get("measured"), off=r0["log10_ratio"], on=r1["log10_ratio"],
                     changed=abs(r1["log10_ratio"] - r0["log10_ratio"]) > 0.005, subset=subset(r0)))


def stats(v):
    v = np.asarray(v)
    return dict(n=len(v), rms_dex=float(np.sqrt(np.mean(v ** 2))) if len(v) else None,
                within30=float(np.mean(np.abs(v) <= np.log10(1.3))) if len(v) else None,
                within2x=float(np.mean(np.abs(v) <= np.log10(2))) if len(v) else None)


summ = defaultdict(dict)
for c in sorted({r["cls"] for r in rows}):
    for tag in ("independent", "same source as L2 widths", "in-sample (trained L3)"):
        v = [r for r in rows if r["cls"] == c and r["subset"] == tag]
        if not v:
            continue
        summ[c][tag] = dict(off=stats([r["off"] for r in v]), on=stats([r["on"] for r in v]),
                            n_changed=sum(r["changed"] for r in v))

os.makedirs(a.on, exist_ok=True)
json.dump(dict(summary=summ, rows=rows), open(os.path.join(a.on, "comparison.json"), "w"), indent=1)
L = ["# Explicit model vs explicit + ML layer", "",
     "Residual = log10(model / measured).  ML layer = learned ZBLAN host microstructure (L2, adopted because it "
     "transfers across ions) + extrinsic quenching rule for silica-family 2 um levels (L3).  Rows whose "
     "measurement trained an ML component are shown separately.", "",
     "| check | subset | n | changed by ML | rms dex: explicit | rms dex: +ML | within 30 %: explicit | +ML | within 2x: explicit | +ML |",
     "|---|---|---|---|---|---|---|---|---|---|"]
for c, d in summ.items():
    for tag, s in d.items():
        o, n = s["off"], s["on"]
        L.append(f"| {c} | {tag} | {o['n']} | {s['n_changed']} | {o['rms_dex']:.3f} | {n['rms_dex']:.3f} | "
                 f"{o['within30']:.0%} | {n['within30']:.0%} | {o['within2x']:.0%} | {n['within2x']:.0%} |")
L += ["", "## Rows the ML layer changed", "",
      "| check | ion | host | item | measured | x explicit | x +ML | subset |", "|---|---|---|---|---|---|---|---|"]
for r in sorted(rows, key=lambda r: (r["cls"], r["ion"], r["host"])):
    if r["changed"]:
        L.append(f"| {r['cls']} | {r['ion']} | {r['host']} | {r['item']} | {r['measured']} | {10 ** r['off']:.2f} | "
                 f"{10 ** r['on']:.2f} | {r['subset']} |")
open(os.path.join(a.on, "COMPARISON.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L[:len(L)]))
