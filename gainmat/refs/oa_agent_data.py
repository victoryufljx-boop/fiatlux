"""Loader for literature records collected (round 5) by extraction agents from open-access papers:
gainmat/refs/data/oa_records.jsonl (bands, cross sections, lifetimes, Omegas, n), et_records.jsonl
(energy-transfer microparameters) and stark_records.jsonl (Stark levels).  Every record carries its
citation and URL; records failing the checks below are dropped (and counted) rather than repaired.
"""
from __future__ import annotations

import json
import os
import re

D = os.path.join(os.path.dirname(__file__), "data")
FAMILIES = {"silica", "aluminosilicate", "germanosilicate", "phosphosilicate", "phosphate", "germanate",
            "tellurite", "ZBLAN", "fluoroindate", "chalcogenide"}
RANGES = {"em_peak_nm": (300, 5000), "abs_peak_nm": (300, 5000), "em_fwhm_nm": (1, 800),
          "sigma_em_cm2": (1e-23, 3e-19), "sigma_abs_cm2": (1e-23, 3e-19), "tau_rad_ms": (1e-4, 200),
          "tau_ms": (1e-4, 200), "n_refr": (1.3, 3.2)}
_LABELS = {}


def _labels(ion):
    if ion not in _LABELS:
        try:
            from ..forward.ions import get_ion
            _LABELS[ion] = {m["label"] for m in get_ion(ion).multiplets}
        except Exception:
            _LABELS[ion] = set()
    return _LABELS[ion]


def _clean_label(s):
    return str(s).replace(" ", "").replace("(", "").replace(")", "")


def _doi(src):
    m = re.search(r"doi[:\s]*(10\.\S+?)(?:[,;)\s]|$)", str(src), re.I)
    return m.group(1).lower().rstrip(".") if m else None


def _read(name):
    p = os.path.join(D, name)
    out = []
    if os.path.exists(p):
        for line in open(p):
            line = line.strip()
            if line:
                try:
                    out.append(json.loads(line))
                except json.JSONDecodeError:
                    pass
    return out


def load_oa(exclude_dois=()):
    """Benchmark records in the format of oa_benchmark_refs.RECORDS; returns (records, rejected_counts)."""
    ex = {d.lower() for d in exclude_dois}
    keep, rej = [], {}
    for r in _read("oa_records.jsonl"):
        why = None
        ion, fam, q = r.get("ion"), r.get("family"), r.get("q")
        pair = tuple(_clean_label(x) for x in (r.get("pair") or [])[:2])
        if fam not in FAMILIES:
            why = "family"
        elif q not in RANGES and q != "omegas":
            why = "quantity"
        elif len(pair) != 2 or not set(pair) <= _labels(ion):
            why = "labels"
        elif _doi(r.get("src")) in ex:
            why = "duplicate paper"
        elif "ABSTRACT ONLY" in str(r.get("note", "")).upper():
            why = "abstract only"
        else:
            v = r.get("value")
            if q == "omegas":
                if not (isinstance(v, list) and len(v) == 3 and all(isinstance(x, (int, float)) and 0 <= x < 40 for x in v)):
                    why = "value"
            elif not isinstance(v, (int, float)) or not (RANGES[q][0] <= v <= RANGES[q][1]):
                why = "value"
        if why:
            rej[why] = rej.get(why, 0) + 1
            continue
        keep.append(dict(ion=ion, family=fam, proxy=bool(r.get("proxy")), pair=pair, q=q,
                         value=tuple(r["value"]) if q == "omegas" else float(r["value"]),
                         src=str(r.get("src")) + (f" [{r.get('url')}]" if r.get("url") else ""),
                         note=str(r.get("composition", "")) + "; " + str(r.get("note", "")), agent=True))
    return keep, rej


def load_et():
    out = []
    for r in _read("et_records.jsonl"):
        try:
            c = float(r["C_cm6_s"])
        except (KeyError, TypeError, ValueError):
            continue
        if not (1e-46 < c < 1e-34) or r.get("family") not in FAMILIES:
            continue
        dt = tuple(_clean_label(x) for x in r.get("donor_transition") or [])
        at = tuple(_clean_label(x) for x in r.get("acceptor_transition") or [])
        if len(dt) != 2 or len(at) != 2 or not set(dt) <= _labels(r.get("donor_ion")) or not set(at) <= _labels(r.get("acceptor_ion")):
            continue
        out.append(dict(r, donor_transition=dt, acceptor_transition=at, C_cm6_s=c))
    return out


def load_stark():
    return _read("stark_records.jsonl")
