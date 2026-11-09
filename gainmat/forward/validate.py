"""Validate the glass forward model against known fibre gain media.

For each entry of gainmat.refs.fiber_validation.VALIDATION, run the forward model for
(ion, host family) and compare:
  * gain window overlapping the measured band (edges, nm)
  * peak absorption / emission cross sections and their wavelengths for the transition
  * emission FWHM
  * upper-level lifetime (1 / (A_rad + W_nr)) against the measured lifetime
Nothing is fitted to these targets.
"""
from __future__ import annotations

import json
import re

import numpy as np

from .glass import run, key_transitions, FAMILIES

FAMILY_ALIASES = {
    "silica": "silica", "pure silica": "silica", "aluminosilicate": "aluminosilicate",
    "al-silica": "aluminosilicate", "alumino-silicate": "aluminosilicate",
    "germanosilicate": "germanosilicate", "phosphosilicate": "phosphosilicate", "phosphate": "phosphate", "germanate": "germanate",
    "tellurite": "tellurite", "zblan": "ZBLAN", "fluoride": "ZBLAN", "fluoroindate": "fluoroindate",
    "inf3": "fluoroindate", "chalcogenide": "chalcogenide", "sulfide": "chalcogenide",
}


def _val(x):
    return x.get("value") if isinstance(x, dict) else x


def _family(s):
    s = (s or "").lower()
    for k, v in FAMILY_ALIASES.items():
        if k in s:
            return v
    return None


def _ion(s):
    s = s.strip()
    return s if s.endswith("3+") else s + "3+"


def _pair(tr):
    m = re.split(r"\s*(?:->|→|-)\s*(?=\d)", tr.strip()) if tr else None
    if not m or len(m) != 2:
        return None
    return (m[0].strip(), m[1].strip())


def validate(entries=None, n_sites=40, verbose=True):
    if entries is None:
        from ..refs.fiber_validation import VALIDATION as entries
    rows, cache = [], {}
    for e in entries:
        ion = _ion(_val(e.get("ion")) or "")
        fam = _family(_val(e.get("host")) or _val(e.get("host_family")) or _val(e.get("family")))
        pair = _pair(_val(e.get("transition")) or "")
        if fam is None or fam not in FAMILIES:
            rows.append(dict(entry=e.get("name") or f"{ion}:{e.get('host')}", skipped="unknown host"))
            continue
        key = (ion, fam)
        try:
            if key not in cache:
                cache[key] = run(ion, fam, n_sites=n_sites)
            r = cache[key]
        except Exception as ex:                       # keep going
            rows.append(dict(entry=f"{ion}:{fam}", error=repr(ex)))
            continue
        row = dict(entry=e.get("name") or f"{ion}:{fam} {pair}", ion=ion, family=fam,
                   transition=pair, omega_src=r["omega_src"])
        band = _val(e.get("gain_band_nm")) or _val(e.get("band_nm"))
        if isinstance(band, (list, tuple)) and len(band) == 2:
            lo, hi = float(band[0]), float(band[1])
            ov = [w for w in r["windows"] if w["hi_nm"] >= lo and w["lo_nm"] <= hi]
            row["band_exp_nm"] = (lo, hi)
            row["window_model_nm"] = [(w["lo_nm"], w["hi_nm"]) for w in ov]
            row["band_reproduced"] = bool(ov)
        if pair:
            kt = key_transitions(r, [pair]).get(f"{pair[0]}->{pair[1]}")
            if kt:
                row.update({f"model_{k}": v for k, v in kt.items()})
            for k in ("peak_abs_cross_section_cm2", "sigma_abs_peak_cm2", "peak_sigma_abs_cm2"):
                if e.get(k) is not None:
                    row["exp_sigma_abs"] = _val(e[k])
            for k in ("peak_em_cross_section_cm2", "sigma_em_peak_cm2", "peak_sigma_em_cm2"):
                if e.get(k) is not None:
                    row["exp_sigma_em"] = _val(e[k])
            for k in ("fwhm_nm", "emission_fwhm_nm", "em_fwhm_nm"):
                if e.get(k) is not None:
                    row["exp_fwhm_nm"] = _val(e[k])
            taus = [d["tau_ms"].get(pair[0]) for d in r["pumps"].values() if "tau_ms" in d]
            taus = [t for t in taus if t is not None]
            if taus:
                row["model_tau_ms"] = float(taus[0])
            for k in ("lifetime_ms", "tau_ms", "upper_lifetime_ms"):
                if e.get(k) is not None:
                    row["exp_tau_ms"] = _val(e[k])
            if e.get("tau_upper_s") is not None and _val(e["tau_upper_s"]) is not None:
                row["exp_tau_ms"] = float(_val(e["tau_upper_s"])) * 1000.0
            for k in ("lambda_abs_peak_nm", "lambda_em_peak_nm"):
                if e.get(k) is not None:
                    row["exp_" + k] = _val(e[k])
            for k in ("lifetime_us", "tau_us", "upper_lifetime_us"):
                if e.get(k) is not None and _val(e[k]) is not None:
                    row["exp_tau_ms"] = float(_val(e[k])) / 1000.0
        rows.append(row)
        if verbose:
            print(json.dumps({k: v for k, v in row.items() if k != "omega_src"}, default=str))
    return rows
