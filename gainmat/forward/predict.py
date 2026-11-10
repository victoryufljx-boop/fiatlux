"""One call for a prediction: the main optical properties of an ion in a glass, each with an
uncertainty taken from the validation residuals and a statement of what it rests on.

    from gainmat.forward.predict import predict
    p = predict("Er3+", "aluminosilicate")          # dict, JSON-serialisable
    python -m gainmat.forward.predict Er3+ aluminosilicate [--ml] [--N 1e20]

Uncertainty is empirical, not a model property: for each quantity class the rms log10 error (or
nm error) of the model against independent measurements in forward_results/ (residuals.jsonl,
spectra_validation.json, level_validation.json), per ion where there are at least MIN_N checks,
otherwise over all ions.  Known systematic biases found in validation (e.g. Tm 3F4 band position
in silicates) are attached as warnings rather than silently corrected.
"""
from __future__ import annotations

import json
import os

import numpy as np

from . import glass as g

R = os.path.join(os.path.dirname(__file__), "..", "..", "forward_results")
MIN_N = 3
MAIN = {
    "Er3+": [("4I13/2", "4I15/2"), ("4I11/2", "4I15/2"), ("4I11/2", "4I13/2")],
    "Yb3+": [("2F5/2", "2F7/2")],
    "Tm3+": [("3F4", "3H6"), ("3H4", "3H6"), ("3H4", "3F4")],
    "Ho3+": [("5I7", "5I8"), ("5I6", "5I7")],
    "Nd3+": [("4F3/2", "4I11/2"), ("4F3/2", "4I9/2"), ("4F3/2", "4I13/2")],
    "Pr3+": [("1G4", "3H5")],
    "Dy3+": [("6H13/2", "6H15/2")],
}
KNOWN_BIAS = [
    dict(ion="*", families=("*",), pair=None,
         text="external open-access benchmark (forward_results/oa_benchmark.md): emission peaks typically "
              "2-8 nm off for Er/Nd/Ho/Yb (Tm 3F4 up to ~90 nm); emission widths ~20 nm; peak cross "
              "sections ~0.19 dex (0.14 dex given the glass's measured Omegas); radiative lifetimes 0.12 dex "
              "(0.015 dex given measured Omegas).  Supplying measured Omegas is the single largest improvement."),
    dict(ion="Tm3+", families=("silica", "aluminosilicate", "germanosilicate", "phosphosilicate"),
         pair=("3F4", "3H6"),
         text="3H6->3F4 absorption peak is ~60-150 nm too long and its composition trend is not "
              "reproduced (measured 1578-1670 nm in silicate fibres); emission peak within ~25 nm. "
              "Use measured_absorption= if a measured band is available (ml_results/RESULTS.md, round 3-4)."),
    dict(ion="Yb3+", families=("silica", "aluminosilicate", "germanosilicate"), pair=("2F5/2", "2F7/2"),
         text="peak cross section ~0.8x Pask 1995 (radiative lifetime 0.92x); emission beyond "
              "~1080 nm underestimated (no vibronic sideband)."),
    dict(ion="Er3+", families=("silica", "aluminosilicate", "germanosilicate"), pair=("4I13/2", "4I15/2"),
         text="sigma_e/sigma_a about +0.1 dex high across 1500-1565 nm (Giles); absorption peak ~3-5 nm short."),
]


def _load(name):
    p = os.path.join(R, name)
    return json.load(open(p)) if os.path.exists(p) else None


def _uncertainty_table():
    """rms errors by (class, ion) from the validation files; '*' = all ions."""
    rows = []
    p = os.path.join(R, "residuals.jsonl")
    if os.path.exists(p):
        for line in open(p):
            r = json.loads(line)
            q = r["quantity"]
            cls = ("tau_rad" if q == "lifetime:tau_rad" else "tau" if q in ("lifetime:tau", "tau_ms(validation)")
                   else "sigma_peak" if q.startswith("cross_section") else "sigma_integrated"
                   if q in ("spectrum_integrated", "vendor_xs_integrated") else None)
            if cls and r.get("log10_ratio") is not None and np.isfinite(r["log10_ratio"]):
                rows.append((cls, r["ion"], float(r["log10_ratio"])))
    sp = _load("spectra_validation.json") or []
    for r in sp:
        for k_mod, k_meas, cls in (("peak_model_nm", "peak_meas_nm", "peak_nm"),):
            if r.get(k_mod) is not None and r.get(k_meas) is not None:
                rows.append((cls, r.get("ion", "*"), float(r[k_mod]) - float(r[k_meas])))
    lv = _load("level_validation.json") or []
    for r in lv:
        if r.get("group") == "cross_section" and r.get("peak_model_nm") and r.get("peak_meas_nm"):
            rows.append(("peak_nm", r["ion"], float(r["peak_model_nm"]) - float(r["peak_meas_nm"])))
    tab = {}
    for cls in {c for c, _, _ in rows}:
        allv = np.array([v for c, _, v in rows if c == cls])
        # fallback for ions with few checks: robust spread over all ions (1.4826 x median |err|), so
        # one ion's known systematic error (e.g. Tm band positions) does not inflate everyone's
        tab[(cls, "*")] = (float(1.4826 * np.median(np.abs(allv))), len(allv))
        for ion in {i for c, i, _ in rows if c == cls}:
            v = [x for c, i, x in rows if c == cls and i == ion]
            if len(v) >= MIN_N:
                tab[(cls, ion)] = (float(np.sqrt(np.mean(np.square(v)))), len(v))
    return tab


_UT = None


def _u(cls, ion):
    global _UT
    if _UT is None:
        _UT = _uncertainty_table()
    v = _UT.get((cls, ion)) or _UT.get((cls, "*"))
    return None if v is None else dict(rms=round(v[0], 3), n_checks=v[1],
                                        scope=ion if (cls, ion) in _UT else "all ions (robust)")


_READY = False


def _setup():
    global _READY
    if not _READY:
        g.load_refs(); g.apply_site_refs(); g.calibrate_spm(370.0, g.CARNALL_ER_LAF3_BKQ, verbose=False)
        _READY = True


def _fwhm_local(nu, y):
    """FWHM (nm) of the contiguous region around the maximum that stays above half maximum."""
    i = int(np.argmax(y)); h = y[i] / 2
    lo = i
    while lo > 0 and y[lo - 1] >= h:
        lo -= 1
    hi = i
    while hi < len(y) - 1 and y[hi + 1] >= h:
        hi += 1
    return float(1e7 / nu[lo] - 1e7 / nu[hi])


def predict(ion, family=None, N_cm3=None, ml=None, n_sites=60, measured_absorption=None, pairs=None,
            composition=None, omegas=None):
    """family: one of glass.FAMILIES, or give composition= (dict of mol%, or a formula such as
    '70TeO2-15K2O-15Nb2O5'); the family is then chosen by gainmat.forward.composition.family_for."""
    _setup()
    fam_rule = None
    if family is None:
        from .composition import family_for
        family, fam_rule = family_for(composition)
        if family is None:
            raise ValueError(f"no model glass family for this composition ({fam_rule})")
    from .fiber import LocalModel
    from .spectra import fwhm_nm
    old_ml = g.USE_ML
    if ml is not None:
        g.USE_ML = bool(ml)
    try:
        ov = {}
        if N_cm3:
            ov["N_cm3"] = N_cm3
        if measured_absorption:
            ov["measured_absorption"] = measured_absorption
        if omegas is not None:
            ov["omegas"] = tuple(float(x) for x in omegas)
        r = g.run(ion, family, n_sites=n_sites, overrides=ov or None)
    finally:
        g.USE_ML = old_ml
    spec = r["_spec"]; nu = spec["nu"]
    fam = g.FAMILIES[family]
    lm = LocalModel(ion, spec, r["omegas"], fam["n"], fam["hw_max"], fam.get("gap_family"), emax=g.EMAX[ion])
    from ..refs.stark_anchor_refs import ANCHORS
    anc = next((v for k, v in ANCHORS.items() if (k == ion or v.get("ion") == ion)
                and v.get("auto", True) and family in v["families"]), None)
    anchored = bool(g.LEVEL_ANCHORS and anc)
    prov = dict(
        stark=("measured levels: " + anc["src"]) if anchored else
        "crystal-field model (handbook superposition parameters x Er-calibrated O/F scale, random sites)",
        omegas=dict(values=list(map(float, r["omegas"])), source="supplied by user" if omegas is not None else r["omega_src"]),
        ml_layer=bool(g.USE_ML if ml is None else ml),
        measured_absorption=bool(measured_absorption),
        refractive_index=fam["n"], max_phonon_cm1=fam["hw_max"])
    bands = []
    for p in (pairs or MAIN.get(ion, [])):
        if p not in spec["sigma_abs"]:
            continue
        sa, se = spec["sigma_abs"][p], spec["sigma_em"][p]
        ia, ie = int(np.argmax(sa)), int(np.argmax(se))
        b = dict(transition=f"{p[0]}->{p[1]}", abs_peak_nm=round(float(1e7 / nu[ia]), 1), sigma_abs_peak_cm2=float(sa[ia]),
                 em_peak_nm=round(float(1e7 / nu[ie]), 1), sigma_em_peak_cm2=float(se[ie]),
                 em_fwhm_nm=round(_fwhm_local(nu, se), 1),
                 uncertainty=dict(sigma_peak_log10=_u("sigma_peak", ion), peak_nm=_u("peak_nm", ion)))
        warn = [k["text"] for k in KNOWN_BIAS if k["ion"] == ion and family in k["families"] and k["pair"] == p]
        if warn:
            b["known_bias"] = warn
        bands.append(b)
    levels = []
    for i, lab in enumerate(lm.labels):
        A = float(lm.A[i].sum())
        if A <= 0:
            continue
        tot = A + float(lm.W[i])
        levels.append(dict(level=lab, tau_rad_ms=1e3 / A, tau_ms=1e3 / tot, quantum_efficiency=A / tot,
                           uncertainty=dict(tau_rad_log10=_u("tau_rad", ion), tau_log10=_u("tau", ion))))
    if fam_rule:
        prov["family_from_composition"] = fam_rule
    general = [k["text"] for k in KNOWN_BIAS if k["ion"] == "*"]
    from .transfer import standard_parameters
    et = standard_parameters(ion, family, spec=spec, overrides=dict(n=fam["n"]))
    tv = _load("transfer_validation.json")
    et_unc = (dict(rms_log10=round(tv["summary"]["overlap / model defaults"]["rms_dex"], 3),
                   n_checks=tv["summary"]["overlap / model defaults"]["n"],
                   note="vs published Forster-Dexter overlaps of measured spectra (all 7 within 2x); a decay-fitted "
                        "Tm cross-relaxation value is 4x higher than predicted")
              if tv else None)
    return dict(ion=ion, family=family, N_cm3=r["N_cm3"], provenance=prov, bands=bands, levels=levels[:8],
                accuracy=general, energy_transfer=dict(parameters=et, uncertainty=et_unc),
                gain_windows=[{k: float(v) for k, v in w.items()} for w in r["windows"]],
                note="uncertainties are rms errors of this model against independent measurements "
                     "(forward_results/); 'known_bias' lists systematic errors found in validation")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("ion"); ap.add_argument("family")
    ap.add_argument("--ml", action="store_true"); ap.add_argument("--N", type=float, default=None)
    a = ap.parse_args()
    print(json.dumps(predict(a.ion, a.family, N_cm3=a.N, ml=a.ml or None), indent=1, default=float))
