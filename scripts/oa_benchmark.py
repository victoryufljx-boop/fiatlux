"""External benchmark on open-access literature (gainmat/refs/oa_benchmark_refs.py).

Three predictors per record, all blind to that record:
  model         the forward model as shipped (its own Omega: measured family value or ion median)
  model+Omega   the same, given the paper's own Judd-Ofelt Omegas (tests the spectroscopy physics
                when the intensity parameters are known, as they are for any measured glass)
  lookup        what a practitioner would do without a model: the mean of the same quantity for the
                same ion and transition in OTHER papers on the same glass family (leave-one-paper-out);
                if there is none, other papers on any family.  None if no other paper reports it.
Errors: log10(pred/meas) for cross sections, lifetimes and Omegas; pred - meas (nm) for peaks/widths.

    python scripts/oa_benchmark.py -> forward_results/oa_benchmark.json, oa_benchmark.md
"""
import json, os, sys
from collections import defaultdict
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np
from gainmat.forward import glass as g
from gainmat.forward.fiber import LocalModel
from gainmat.refs.oa_benchmark_refs import RECORDS

g.load_refs(); g.apply_site_refs(); g.calibrate_spm(370.0, g.CARNALL_ER_LAF3_BKQ, verbose=False)
LOG = {"sigma_em_cm2", "sigma_abs_cm2", "tau_rad_ms", "tau_ms"}
_cache = {}


def model_vals(ion, fam, omegas=None):
    key = (ion, fam, tuple(omegas) if omegas else None)
    if key in _cache:
        return _cache[key]
    ov = dict(omegas=tuple(omegas)) if omegas else None
    r = g.run(ion, fam, n_sites=60, pumps=[], overrides=ov)
    F = g.FAMILIES[fam]
    lm = LocalModel(ion, r["_spec"], r["omegas"], F["n"], F["hw_max"], F.get("gap_family"), emax=g.EMAX[ion])
    _cache[key] = (r, lm)
    return r, lm


def fwhm_local(nu, y):
    i = int(np.argmax(y)); h = y[i] / 2; lo = hi = i
    while lo > 0 and y[lo - 1] >= h: lo -= 1
    while hi < len(y) - 1 and y[hi + 1] >= h: hi += 1
    return float(1e7 / nu[lo] - 1e7 / nu[hi])


def predict_q(ion, fam, pair, q, omegas=None):
    r, lm = model_vals(ion, fam, omegas)
    sp = r["_spec"]; nu = sp["nu"]
    if q in ("tau_rad_ms", "tau_ms"):
        i = lm.labels.index(pair[0]); A = float(lm.A[i].sum())
        return 1e3 / A if q == "tau_rad_ms" else 1e3 / (A + float(lm.W[i]))
    if q == "omegas":
        return tuple(map(float, r["omegas"]))
    if pair not in sp["sigma_abs"]:
        return None
    y = sp["sigma_em" if q.startswith(("em", "sigma_em")) else "sigma_abs"][pair]
    return {"em_peak_nm": lambda: float(1e7 / nu[np.argmax(y)]), "abs_peak_nm": lambda: float(1e7 / nu[np.argmax(y)]),
            "em_fwhm_nm": lambda: fwhm_local(nu, y), "sigma_em_cm2": lambda: float(y.max()),
            "sigma_abs_cm2": lambda: float(y.max())}[q]()


def paper_omegas(rec):
    """Omegas the same paper gives for the same sample (matched by src and the note's sample tag)."""
    tag = rec["note"].split(";")[0].split(",")[0].strip()
    for o in RECORDS:
        if o["q"] == "omegas" and o["src"] == rec["src"] and o["ion"] == rec["ion"] and o["family"] == rec["family"]:
            otag = o["note"].split(";")[0].split(",")[0].strip()
            if not tag or not otag or tag == otag or tag in o["note"] or otag in rec["note"]:
                return o["value"]
    return None


def err(q, pred, meas):
    if pred is None or meas is None:
        return None
    if q == "omegas":
        return float(np.mean([np.log10(max(p, 1e-3) / max(m, 1e-3)) for p, m in zip(pred, meas)]))
    return float(np.log10(pred / meas)) if q in LOG else float(pred - meas)


def lookup(rec):
    same = [o["value"] for o in RECORDS if o["ion"] == rec["ion"] and o["pair"] == rec["pair"] and o["q"] == rec["q"]
            and o["src"] != rec["src"] and o["family"] == rec["family"]]
    scope = "same family"
    if not same:
        same = [o["value"] for o in RECORDS if o["ion"] == rec["ion"] and o["pair"] == rec["pair"] and o["q"] == rec["q"]
                and o["src"] != rec["src"]]
        scope = "any family"
    if not same:
        return None, None
    if rec["q"] == "omegas":
        return tuple(np.exp(np.mean(np.log(np.array(same, float)), axis=0))), scope
    return (float(np.exp(np.mean(np.log(same)))) if rec["q"] in LOG else float(np.mean(same))), scope


rows = []
for rec in RECORDS:
    q = rec["q"]
    m = predict_q(rec["ion"], rec["family"], rec["pair"], q)
    om = paper_omegas(rec) if q != "omegas" else None
    mo = predict_q(rec["ion"], rec["family"], rec["pair"], q, om) if om else None
    lk, scope = lookup(rec)
    row = dict(ion=rec["ion"], family=rec["family"], proxy=rec["proxy"], pair=rec["pair"], q=q, measured=rec["value"],
               model=m, model_paper_omega=mo, lookup=lk, lookup_scope=scope, src=rec["src"], note=rec["note"],
               err_model=err(q, m, rec["value"]), err_model_paper_omega=err(q, mo, rec["value"]),
               err_lookup=err(q, lk, rec["value"]))
    rows.append(row)
    print(json.dumps({k: row[k] for k in ("ion", "family", "q", "measured", "model", "model_paper_omega", "lookup",
                                          "err_model", "err_model_paper_omega", "err_lookup")}, default=float), flush=True)

# summary per quantity class: rows where the lookup exists too (fair comparison) and all rows
summ = {}
for q in sorted({r["q"] for r in rows}):
    R = [r for r in rows if r["q"] == q]
    def stat(key, sub):
        v = [abs(r[key]) for r in sub if r[key] is not None]
        return dict(n=len(v), median_abs=float(np.median(v)) if v else None, rms=float(np.sqrt(np.mean(np.square(v)))) if v else None)
    both = [r for r in R if r["err_lookup"] is not None and r["err_model"] is not None]
    summ[q] = dict(unit="dex" if (q in LOG or q == "omegas") else "nm",
                   all=dict(model=stat("err_model", R), model_paper_omega=stat("err_model_paper_omega", R)),
                   paired_with_lookup=dict(model=stat("err_model", both), model_paper_omega=stat("err_model_paper_omega", both),
                                           lookup=stat("err_lookup", both)),
                   model_beats_lookup=sum(abs(r["err_model"]) < abs(r["err_lookup"]) for r in both), n_paired=len(both))
json.dump(dict(summary=summ, rows=rows), open("forward_results/oa_benchmark.json", "w"), indent=1, default=float)
md = ["# External benchmark: open-access literature", "",
      "`scripts/oa_benchmark.py`, data in `gainmat/refs/oa_benchmark_refs.py` (never used to build the model).", "",
      "| quantity | unit | n | model (own Omega) median abs err | model + paper's Omega | lookup (other papers) | model beats lookup |",
      "|---|---|---|---|---|---|---|"]
for q, s in summ.items():
    a = s["all"]; p = s["paired_with_lookup"]
    f = lambda d: "—" if d["median_abs"] is None else f"{d['median_abs']:.3g} (n={d['n']})"
    md.append(f"| {q} | {s['unit']} | {a['model']['n']} | {f(a['model'])} | {f(a['model_paper_omega'])} | {f(p['lookup'])} | "
              f"{s['model_beats_lookup']}/{s['n_paired']} |")
open("forward_results/oa_benchmark.md", "w").write("\n".join(md) + "\n")
print("\n".join(md))
