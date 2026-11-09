"""L4: a learned error model for the forward model's outputs.

Rather than a free-form residual regressor (84 heterogeneous points would be memorised), the
residuals log10(model/measured) are modelled per output class c and host class h as
    r ~ Normal(bias_{c,h}, s_{c,h}^2)
with partial pooling (bias_{c,h} shrunk toward bias_c, s from the pooled scatter).  The model
answers two questions honestly:
  1. Is there a *systematic* bias worth correcting?  Only if correcting it improves residuals
     on held-out sources (leave-one-source-out), not just in-sample.
  2. What error bar should a new prediction carry?  Checked by the coverage of 68 / 90 %
     intervals on held-out sources.

    python -m gainmat.ml.error_model  ->  ml_results/l4/error_model.json
"""
from __future__ import annotations

import json
import os
from collections import defaultdict

import numpy as np

FLUORIDE = {"ZBLAN", "fluoroindate"}


def _hclass(host):
    return "fluoride" if host in FLUORIDE else "oxide"


def _cclass(q):
    q = q.split("(")[0]
    return {"lifetime:tau": "lifetime", "tau_ms": "lifetime", "lifetime:tau_rad": "radiative_lifetime",
            "cross_section:em": "cross_section", "cross_section:abs": "cross_section",
            "spectrum_integrated": "cross_section_integrated", "vendor_xs_integrated": "cross_section_integrated"
            }.get(q, None)


def load(path="forward_results/residuals.jsonl"):
    rows = []
    for line in open(path):
        r = json.loads(line)
        c = _cclass(r["quantity"])
        if c is None or r.get("log10_ratio") is None or not np.isfinite(r["log10_ratio"]):
            continue
        r["cclass"], r["hclass"] = c, _hclass(r["host"])
        r["source"] = r["src"].split(",")[0].split("(")[0].strip()[:40]
        rows.append(r)
    return rows


def fit(rows, k=3.0):
    """Partial pooling: bias_{c,h} = (n mean + k bias_c) / (n + k); scatter pooled per class."""
    by_c, by_ch = defaultdict(list), defaultdict(list)
    for r in rows:
        by_c[r["cclass"]].append(r["log10_ratio"]); by_ch[(r["cclass"], r["hclass"])].append(r["log10_ratio"])
    par = {}
    for c, v in by_c.items():
        bc = float(np.mean(v)); sc = float(np.std(v, ddof=1)) if len(v) > 1 else 0.3
        par[c] = dict(bias=bc, scatter=max(sc, 0.05), n=len(v))
        for h in ("oxide", "fluoride"):
            vv = by_ch.get((c, h), [])
            n = len(vv)
            b = (np.sum(vv) + k * bc) / (n + k) if n else bc
            s = np.sqrt((np.sum((np.array(vv) - b) ** 2) + k * sc ** 2) / (n + k)) if n else sc
            par[f"{c}|{h}"] = dict(bias=float(b), scatter=float(max(s, 0.05)), n=n)
    return par


def predict(par, c, h):
    p = par.get(f"{c}|{h}") or par.get(c)
    return (p["bias"], p["scatter"]) if p else (0.0, 0.3)


def loso(rows):
    """Leave-one-source-out: fit on other sources, predict held-out residuals."""
    srcs = sorted({r["source"] for r in rows})
    out = []
    for s in srcs:
        tr = [r for r in rows if r["source"] != s]; te = [r for r in rows if r["source"] == s]
        par = fit(tr)
        for r in te:
            b, sd = predict(par, r["cclass"], r["hclass"])
            out.append(dict(source=s, cclass=r["cclass"], hclass=r["hclass"], resid=r["log10_ratio"],
                            corrected=r["log10_ratio"] - b, sigma=sd))
    return out


def main(out_dir="ml_results/l4"):
    os.makedirs(out_dir, exist_ok=True)
    rows = load()
    par = fit(rows)
    lo = loso(rows)
    summ = {}
    for c in sorted({r["cclass"] for r in lo}):
        v = [r for r in lo if r["cclass"] == c]
        raw = np.array([r["resid"] for r in v]); cor = np.array([r["corrected"] for r in v])
        sd = np.array([r["sigma"] for r in v])
        summ[c] = dict(n=len(v), rms_raw=float(np.sqrt(np.mean(raw ** 2))), rms_bias_corrected=float(np.sqrt(np.mean(cor ** 2))),
                       coverage68=float(np.mean(np.abs(cor) <= sd)), coverage90=float(np.mean(np.abs(cor) <= 1.645 * sd)),
                       correction_helps_out_of_source=bool(np.sqrt(np.mean(cor ** 2)) < np.sqrt(np.mean(raw ** 2))))
        print(f"{c:26s} n={len(v):3d}  rms raw {summ[c]['rms_raw']:.3f}  bias-corrected {summ[c]['rms_bias_corrected']:.3f}  "
              f"cov68 {summ[c]['coverage68']:.2f} cov90 {summ[c]['coverage90']:.2f}")
    json.dump(dict(parameters=par, leave_one_source_out=summ, n_rows=len(rows)),
              open(os.path.join(out_dir, "error_model.json"), "w"), indent=1)
    for k, v in par.items():
        print(f"  {k:36s} bias {v['bias']:+.3f} dex  scatter {v['scatter']:.3f} dex  n={v['n']}")
    return par, summ


if __name__ == "__main__":
    main()
