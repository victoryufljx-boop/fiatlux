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
from gainmat.refs.oa_benchmark_refs import RECORDS as _HAND
from gainmat.refs.oa_agent_data import load_oa
import re as _re
_AG, _REJ = load_oa([m for r in _HAND for m in _re.findall(r"doi:(10\.\S+?)(?:\)|,|$| )", r["src"])])
RECORDS = sorted(list(_HAND) + _AG, key=lambda r: (r["ion"], r["family"]))   # grouped: cache reuse
print(f"records: {len(_HAND)} hand-collected + {len(_AG)} agent-collected (rejected {_REJ})", flush=True)

g.load_refs(); g.apply_site_refs(); g.calibrate_spm(370.0, g.CARNALL_ER_LAF3_BKQ, verbose=False)
LOG = {"sigma_em_cm2", "sigma_abs_cm2", "tau_rad_ms", "tau_ms"}   # n_refr, peaks, widths: linear errors
_cache = {}


def model_vals(ion, fam, omegas=None):
    key = (ion, fam, tuple(omegas) if omegas else None)
    if key in _cache:
        return _cache[key]
    ov = dict(omegas=tuple(omegas)) if omegas else None
    r = g.run(ion, fam, n_sites=60, pumps=[], overrides=ov)
    F = g.FAMILIES[fam]
    lm = LocalModel(ion, r["_spec"], r["omegas"], F["n"], F["hw_max"], F.get("gap_family"), emax=g.EMAX[ion])
    # keep only what the benchmark reads (full run dicts for ~150 (ion, glass, Omega) keys exhaust memory)
    r = dict(_spec=r["_spec"], omegas=r["omegas"])
    if len(_cache) >= 4:
        _cache.pop(next(iter(_cache)))
    import gc; gc.collect()
    _cache[key] = (r, lm)
    return r, lm


def fwhm_local(nu, y):
    i = int(np.argmax(y)); h = y[i] / 2; lo = hi = i
    while lo > 0 and y[lo - 1] >= h: lo -= 1
    while hi < len(y) - 1 and y[hi + 1] >= h: hi += 1
    return float(1e7 / nu[lo] - 1e7 / nu[hi])


def predict_q(ion, fam, pair, q, omegas=None):
    if q == "n_refr":
        return float(g.FAMILIES[fam]["n"])
    r, lm = model_vals(ion, fam, omegas)
    sp = r["_spec"]; nu = sp["nu"]
    if q in ("tau_rad_ms", "tau_ms"):
        if pair[0] not in lm.labels:
            return None
        i = lm.labels.index(pair[0]); A = float(lm.A[i].sum())
        return 1e3 / A if q == "tau_rad_ms" else 1e3 / (A + float(lm.W[i]))
    if q == "omegas":
        return tuple(map(float, r["omegas"]))
    if q == "n_refr":
        return float(g.FAMILIES[fam]["n"])
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


def library(rec, q):
    """Mean of q (em_peak_nm / em_fwhm_nm) for the same ion, transition and family in OTHER papers."""
    v = [o["value"] for o in RECORDS if o["ion"] == rec["ion"] and o["pair"] == rec["pair"] and o["q"] == q
         and o["family"] == rec["family"] and o["src"] != rec["src"]]
    return float(np.mean(v)) if v else None


def hybrid_q(rec, omegas=None):
    """Model band moved and stretched (in energy, about its emission peak) so that its emission peak and
    FWHM equal the library values from other papers; the model's integrated strength is kept, so the
    peak cross sections are predictions.  Falls back to the model where the library is empty."""
    ion, fam, pair, q = rec["ion"], rec["family"], rec["pair"], rec["q"]
    if q not in ("em_peak_nm", "em_fwhm_nm", "sigma_em_cm2", "sigma_abs_cm2", "abs_peak_nm"):
        return None
    r, _ = model_vals(ion, fam, omegas)
    sp = r["_spec"]; nu = sp["nu"]
    if pair not in sp["sigma_em"]:
        return None
    se, sa = sp["sigma_em"][pair], sp["sigma_abs"][pair]
    lp, lw = library(rec, "em_peak_nm"), library(rec, "em_fwhm_nm")
    nu0 = nu[np.argmax(se)]
    shift = (1e7 / lp - nu0) if lp else 0.0
    s = 1.0
    if lw:
        wm = fwhm_local(nu, se); lam0 = 1e7 / nu0
        # FWHM in energy, model vs library (library width given in nm at the library peak)
        s = (lw / (lp or lam0) ** 2) / (wm / lam0 ** 2)
    def tf(y):                                       # y'(nu) = y(nu0 + (nu - nu0 - shift) / s) / s
        return np.interp(nu0 + (nu - nu0 - shift) / s, nu, y, left=0, right=0) / s
    se2, sa2 = tf(se), tf(sa)
    return {"em_peak_nm": lambda: float(1e7 / nu[np.argmax(se2)]), "abs_peak_nm": lambda: float(1e7 / nu[np.argmax(sa2)]),
            "em_fwhm_nm": lambda: fwhm_local(nu, se2), "sigma_em_cm2": lambda: float(se2.max()),
            "sigma_abs_cm2": lambda: float(sa2.max())}[q]()


rows = []
for rec in RECORDS:
    q = rec["q"]
    try:
        m = predict_q(rec["ion"], rec["family"], rec["pair"], q)
    except Exception as e:
        print("skip", rec["ion"], rec["family"], rec["pair"], q, repr(e)[:80]); continue
    om = paper_omegas(rec) if q != "omegas" else None
    try:
        mo = predict_q(rec["ion"], rec["family"], rec["pair"], q, om) if om else None
    except Exception:
        mo = None
    lk, scope = lookup(rec)
    try:
        hy = hybrid_q(rec)
        hyo = hybrid_q(rec, om) if om else None
    except Exception:
        hy = hyo = None
    row = dict(ion=rec["ion"], family=rec["family"], proxy=rec["proxy"], pair=rec["pair"], q=q, measured=rec["value"],
               model=m, model_paper_omega=mo, lookup=lk, lookup_scope=scope, src=rec["src"], note=rec["note"],
               err_model=err(q, m, rec["value"]), err_model_paper_omega=err(q, mo, rec["value"]),
               agent=bool(rec.get("agent")),
               err_lookup=err(q, lk, rec["value"]), hybrid=hy, hybrid_paper_omega=hyo,
               err_hybrid=err(q, hy, rec["value"]), err_hybrid_paper_omega=err(q, hyo, rec["value"]))
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
    hyb = [r for r in R if r["err_hybrid"] is not None and r["err_lookup"] is not None]
    summ[q] = dict(unit="dex" if (q in LOG or q == "omegas") else ("-" if q == "n_refr" else "nm"),
                   all=dict(model=stat("err_model", R), model_paper_omega=stat("err_model_paper_omega", R)),
                   paired_with_lookup=dict(model=stat("err_model", both), model_paper_omega=stat("err_model_paper_omega", both),
                                           lookup=stat("err_lookup", both)),
                   model_beats_lookup=sum(abs(r["err_model"]) < abs(r["err_lookup"]) for r in both), n_paired=len(both),
                   hybrid_vs_lookup=dict(hybrid=stat("err_hybrid", hyb), lookup=stat("err_lookup", hyb), model=stat("err_model", hyb),
                                         hybrid_paper_omega=stat("err_hybrid_paper_omega", hyb),
                                         hybrid_beats_lookup=sum(abs(r["err_hybrid"]) < abs(r["err_lookup"]) for r in hyb), n=len(hyb)))
def _med(sub, key):
    v = [abs(r[key]) for r in sub if r[key] is not None]
    return (float(np.median(v)), len(v)) if v else (None, 0)
by_ion = {}
for q in sorted({r["q"] for r in rows}):
    for ion in sorted({r["ion"] for r in rows if r["q"] == q}):
        sub = [r for r in rows if r["q"] == q and r["ion"] == ion]
        both = [r for r in sub if r["err_lookup"] is not None and r["err_model"] is not None]
        by_ion[f"{q}|{ion}"] = dict(model=_med(sub, "err_model"), model_paper_omega=_med(sub, "err_model_paper_omega"),
                                    lookup=_med(both, "err_lookup"), model_on_paired=_med(both, "err_model"),
                                    model_beats_lookup=sum(abs(r["err_model"]) < abs(r["err_lookup"]) for r in both), n_paired=len(both))
no_tm = {}
for q in sorted({r["q"] for r in rows}):
    sub = [r for r in rows if r["q"] == q and r["ion"] != "Tm3+"]
    both = [r for r in sub if r["err_lookup"] is not None and r["err_model"] is not None]
    no_tm[q] = dict(model=_med(sub, "err_model"), lookup=_med(both, "err_lookup"), model_on_paired=_med(both, "err_model"),
                    model_beats_lookup=sum(abs(r["err_model"]) < abs(r["err_lookup"]) for r in both), n_paired=len(both))
json.dump(dict(summary=summ, by_ion=by_ion, without_tm=no_tm, n_hand=len(_HAND), n_agent=len(_AG), rejected=_REJ, rows=rows),
          open("forward_results/oa_benchmark.json", "w"), indent=1, default=float)
md = ["# External benchmark: open-access literature", "",
      "`scripts/oa_benchmark.py`, data in `gainmat/refs/oa_benchmark_refs.py` (never used to build the model).", "",
      "| quantity | unit | n | model (own Omega) median abs err | model + paper's Omega | lookup (other papers) | model beats lookup |",
      "|---|---|---|---|---|---|---|"]
for q, s in summ.items():
    a = s["all"]; p = s["paired_with_lookup"]
    f = lambda d: "—" if d["median_abs"] is None else f"{d['median_abs']:.3g} (n={d['n']})"
    md.append(f"| {q} | {s['unit']} | {a['model']['n']} | {f(a['model'])} | {f(a['model_paper_omega'])} | {f(p['lookup'])} | "
              f"{s['model_beats_lookup']}/{s['n_paired']} |")
md += ["", "Hybrid = model band shifted/stretched to the emission peak and FWHM reported by OTHER papers on the same "
       "ion and family (library), integrated strength from the model.  Rows where both hybrid and lookup exist:", "",
       "| quantity | n | model | hybrid | hybrid + paper's Omega | lookup | hybrid beats lookup |", "|---|---|---|---|---|---|---|"]
for q, s in summ.items():
    h = s["hybrid_vs_lookup"]
    if h["n"]:
        f = lambda d: "—" if d["median_abs"] is None else f"{d['median_abs']:.3g}"
        md.append(f"| {q} | {h['n']} | {f(h['model'])} | {f(h['hybrid'])} | {f(h['hybrid_paper_omega'])} | {f(h['lookup'])} | {h['hybrid_beats_lookup']}/{h['n']} |")
md += ["", "Without Tm (median abs error; paired = rows where a lookup exists):", "",
       "| quantity | model (all) | model (paired) | lookup (paired) | model beats lookup |", "|---|---|---|---|---|"]
for q, v in no_tm.items():
    f = lambda t: "—" if t[0] is None else f"{t[0]:.3g} (n={t[1]})"
    md.append(f"| {q} | {f(v['model'])} | {f(v['model_on_paired'])} | {f(v['lookup'])} | {v['model_beats_lookup']}/{v['n_paired']} |")
md += ["", f"Records: {len(_HAND)} hand-collected + {len(_AG)} collected by extraction agents (rejected by checks: {_REJ})."]
open("forward_results/oa_benchmark.md", "w").write("\n".join(md) + "\n")
print("\n".join(md))
