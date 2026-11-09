"""Summarise the L2 host-microstructure fits (ml_results/l2/*.json):

  - per family: the fitted theta (all ions) and its spread across leave-one-ion-out folds
    (a parameter that moves a lot when one ion is dropped is set by that ion, not the glass)
  - leave-one-ion-out: held-out loss with the default theta vs with theta learned from the
    other ions; per-observation detail
  - cross-host transfer: theta of family A applied to family B's observations
  - writes ml_results/l2/SUMMARY.json and l2_summary.png; learned_theta.json is what the
    forward model loads with use_ml

    python -m gainmat.ml.l2_summary [--cross]
"""
from __future__ import annotations

import argparse
import glob
import json
import os

import numpy as np


def collect(d="ml_results/l2"):
    out = {}
    for p in sorted(glob.glob(os.path.join(d, "*.json"))):
        if os.path.basename(p) in ("SUMMARY.json", "learned_theta.json", "cross_host.json"):
            continue
        r = json.load(open(p))
        out.setdefault(r["family"], {})[r["holdout"] or "all"] = r
    return out


def summarise(res):
    S = {}
    for fam, runs in res.items():
        s = dict(theta_all=runs.get("all", {}).get("theta_fit"), loio={})
        ths = [r["theta_fit"] for k, r in runs.items()]
        if ths:
            s["theta_spread"] = {k: [float(min(t[k] for t in ths)), float(max(t[k] for t in ths))] for k in ths[0]}
        if "all" in runs:
            s["loss_all_default"] = runs["all"]["loss_train_default"]
            s["loss_all_fit"] = runs["all"]["loss_train_fit"]
        for ion, r in runs.items():
            if ion == "all":
                continue
            s["loio"][ion] = dict(loss_default=r.get("loss_test_default"), loss_learned=r.get("loss_test_fit"),
                                  n_test=r.get("n_test"), theta=r["theta_fit"],
                                  improved=bool(r.get("loss_test_fit", 1e9) < r.get("loss_test_default", 0)))
        S[fam] = s
    return S


def components(res):
    """Held-out loss split by what the residual measures.  theta is a broadening/crystal-field
    description of the host, so it should move *widths* and *splittings*; line-shape detail,
    peak positions and the McCumber sigma_e/sigma_a ratio depend on the ion's Stark structure."""
    out = {}
    for fam, runs in res.items():
        acc = {}
        for ion, r in runs.items():
            if ion == "all" or "detail_test_default" not in r:
                continue
            for a, b in zip(r["detail_test_default"], r["detail_test_fit"]):
                if a["kind"] == "curve":
                    parts = [("shape", 0), ("width", 1), ("peak_shift", 2)]
                elif a["kind"] in ("eff_width_nm", "zpl_fwhm_cm1"):
                    parts = [("width", 0)]
                elif a["kind"] == "split_cm1":
                    parts = [("splitting", 0)]
                else:
                    parts = [("sigma_e/sigma_a ratio", 0)]
                for name, i in parts:
                    d = acc.setdefault(name, dict(default=0.0, learned=0.0, n=0))
                    d["default"] += a["resid"][i] ** 2; d["learned"] += b["resid"][i] ** 2; d["n"] += 1
        out[fam] = acc
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="ml_results/l2")
    ap.add_argument("--cross", action="store_true", help="evaluate cross-host transfer (runs the model)")
    a = ap.parse_args()
    res = collect(a.dir)
    S = summarise(res)
    # Adoption gate: a family's learned theta goes into the forward model only if it transfers,
    # i.e. summed over the leave-one-ion-out folds the held-out loss falls AND most folds improve.
    # A theta that only fits the ions it was trained on is a per-ion fudge, not a host property.
    learned = {}
    for fam, s in S.items():
        f = s["loio"]
        if not s.get("theta_all") or not f:
            s["adopted"] = False; continue
        L0 = sum(v["loss_default"] for v in f.values()); L1 = sum(v["loss_learned"] for v in f.values())
        nb = sum(v["improved"] for v in f.values())
        s["loio_total"] = dict(default=L0, learned=L1, folds_improved=nb, folds=len(f))
        s["adopted"] = bool(L1 < L0 and nb > len(f) / 2)
        if s["adopted"]:
            learned[fam] = s["theta_all"]
    json.dump(learned, open(os.path.join(a.dir, "learned_theta.json"), "w"), indent=1)
    for fam, s in S.items():
        print(f"== {fam}: theta(all) {s.get('theta_all')}")
        if "loss_all_default" in s:
            print(f"   all-ion loss {s['loss_all_default']:.1f} -> {s['loss_all_fit']:.1f}")
        if "loio_total" in s:
            t = s["loio_total"]
            print(f"   held-out total {t['default']:.1f} -> {t['learned']:.1f}, {t['folds_improved']}/{t['folds']} folds better; "
                  f"{'ADOPTED' if s['adopted'] else 'not adopted (does not transfer across ions)'}")
        for ion, v in s["loio"].items():
            print(f"   leave-out {ion:5s}: held-out loss {v['loss_default']:.1f} (explicit) -> {v['loss_learned']:.1f} (learned from other ions)"
                  f"  {'BETTER' if v['improved'] else 'worse'}")
    comp = components(res)
    for fam, acc in comp.items():
        print(f"   {fam} held-out by component: " + "; ".join(
            f"{k} {v['default']:.0f}->{v['learned']:.0f} (n={v['n']})" for k, v in acc.items()))
        S[fam]["loio_components"] = acc
    if a.cross:
        from ..forward.glass import load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ
        from .observations import build
        from .host_calib import evaluate, default_theta
        load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ, verbose=False)
        obs = build()
        cross = {}
        allth = {fam: s["theta_all"] for fam, s in S.items() if s.get("theta_all")}
        for B in allth:
            ob = [o for o in obs if o["family"] == B]
            row = {"default": evaluate(B, default_theta(B), ob, 20)}
            for A, th in allth.items():
                row[A] = evaluate(B, th, ob, 20)
            cross[B] = row
            print("cross-host", B, {k: round(v, 1) for k, v in row.items()})
        S["cross_host"] = cross
    json.dump(S, open(os.path.join(a.dir, "SUMMARY.json"), "w"), indent=1, default=float)
    # plot
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fams = [f for f in S if f != "cross_host" and S[f]["loio"]]
    if fams:
        fig, axs = plt.subplots(1, len(fams), figsize=(4.5 * len(fams), 3.8), squeeze=False)
        for ax, f in zip(axs[0], fams):
            ions = list(S[f]["loio"]); x = np.arange(len(ions))
            ax.bar(x - 0.2, [S[f]["loio"][i]["loss_default"] for i in ions], 0.4, label="explicit model")
            ax.bar(x + 0.2, [S[f]["loio"][i]["loss_learned"] for i in ions], 0.4, label="host learned from other ions")
            ax.set_xticks(x); ax.set_xticklabels(ions); ax.set_yscale("log"); ax.set_title(f"{f}: held-out ion loss", fontsize=9)
            ax.legend(fontsize=7)
        fig.tight_layout(); fig.savefig(os.path.join(a.dir, "l2_summary.png"), dpi=130)


if __name__ == "__main__":
    main()
