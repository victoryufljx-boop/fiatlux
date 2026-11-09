"""Validation against commercial-fibre data.

1. nLIGHT LIEKKI Er and Yb cross sections (silica, Al-codoped).  The digitised curves live in
   /home/claude/private_data/vendor_liekki.py (source carries a confidentiality notice), so only
   summary metrics are written to the public results.  Skipped if the private file is absent.
2. Le Verre Fluore datasheets (public): lasers demonstrated on their fluoride fibres; is each
   wavelength inside a model gain window for at least one plausible pump?

    python scripts/validate_vendor.py  ->  forward_results/vendor_validation.json
"""
import json
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np

from gainmat.forward.glass import load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ, run, EMAX
from gainmat.refs.vendor_lvf_refs import LVF_APPLICATIONS

load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ, verbose=False)
out = dict(liekki=[], lvf=[])

PRIV = "/home/claude/private_data"
if os.path.exists(os.path.join(PRIV, "vendor_liekki.py")):
    sys.path.insert(0, PRIV)
    from vendor_liekki import LIEKKI_XS as L
    sets = [("Er 1.5 um band", "Er3+", ("4I13/2", "4I15/2"), L["Er"]["wl_1500band_nm"],
             L["Er"]["sigma_abs_1500band_cm2"], L["Er"]["sigma_em_1500band_cm2"]),
            ("Er 980 nm band", "Er3+", ("4I11/2", "4I15/2"), L["Er"]["wl_980band_nm"],
             L["Er"]["sigma_abs_980band_cm2"], None),
            ("Yb", "Yb3+", ("2F5/2", "2F7/2"), L["Yb"]["wl_nm"], L["Yb"]["sigma_abs_cm2"], L["Yb"]["sigma_em_cm2"])]
    cache = {}
    for name, ion, pair, wl, sa, se in sets:
        if ion not in cache:
            cache[ion] = run(ion, "aluminosilicate", n_sites=30, pumps=[])["_spec"]
        sp = cache[ion]
        lam = np.array(wl, float)
        for kind, meas in (("abs", sa), ("em", se)):
            if meas is None:
                continue
            m = np.array([np.nan if v is None else v for v in meas], float)
            ok = np.isfinite(m) & (m > 0)
            if ok.sum() < 3:
                continue
            l, a = lam[ok], m[ok]
            b = np.interp(1e7 / l, sp["nu"], sp["sigma_abs" if kind == "abs" else "sigma_em"][pair])
            row = dict(curve=f"{name} {kind}", ion=ion, host="aluminosilicate",
                       integrated_ratio=float(np.trapezoid(b, l) / np.trapezoid(a, l)),
                       peak_ratio=float(b.max() / a.max()),
                       peak_shift_nm=float(l[b.argmax()] - l[a.argmax()]),
                       shape_rms=float(np.sqrt(np.mean((b / b.max() - a / a.max()) ** 2))),
                       src="nLIGHT LIEKKI cross sections (Application Designer v4.0), via Optiwave")
            out["liekki"].append(row)
            print(f"LIEKKI {row['curve']:22s} integrated x{row['integrated_ratio']:.2f}  peak x{row['peak_ratio']:.2f}  "
                  f"shift {row['peak_shift_nm']:+.0f} nm  shape rms {row['shape_rms']:.2f}")
        if se is not None:
            m_a = np.array([np.nan if v is None else v for v in sa], float)
            m_e = np.array([np.nan if v is None else v for v in se], float)
            ok = np.isfinite(m_a) & np.isfinite(m_e) & (m_a > 0.05 * np.nanmax(m_a)) & (m_e > 0.05 * np.nanmax(m_e))
            l = lam[ok]
            r_m = np.interp(1e7 / l, sp["nu"], sp["sigma_em"][pair]) / np.interp(1e7 / l, sp["nu"], sp["sigma_abs"][pair])
            r_e = m_e[ok] / m_a[ok]
            row = dict(curve=f"{name} sigma_e/sigma_a", ion=ion, host="aluminosilicate",
                       rms_log10=float(np.sqrt(np.mean(np.log10(r_m / r_e) ** 2))),
                       mean_log10=float(np.mean(np.log10(r_m / r_e))),
                       src="nLIGHT LIEKKI cross sections")
            out["liekki"].append(row)
            print(f"LIEKKI {row['curve']:22s} ratio rms {row['rms_log10']:.2f} dec, bias {row['mean_log10']:+.2f} dec")
else:
    print("private LIEKKI data not found: skipped")

# ---------------------------------------------------------------- Le Verre Fluore lasers
PUMPS = {"Er3+": [976.0, (985.0, 1973.0)], "Tm3+": [793.0, 690.0, 1064.0], "Ho3+": [1150.0, 890.0],
         "Dy3+": [1100.0, 1700.0, 2800.0], "Pr3+": [444.0, 1017.0]}
cache = {}
for row in LVF_APPLICATIONS:
    ion = row["dopants"][0] + "3+"
    if ion not in EMAX:
        continue
    for res in row["results"]:
        wl = res.get("wavelength_um"); rng = res.get("wavelength_range_um")
        if isinstance(wl, (list, tuple)):
            targets = [w * 1e3 for w in wl]
        else:
            targets = [wl * 1e3] if wl else ([rng[0] * 1e3, rng[1] * 1e3] if rng else [])
        if not targets:
            continue
        if ion not in cache:
            r = run(ion, "ZBLAN", n_sites=20, pumps=PUMPS.get(ion, [None]))
            cache[ion] = {str(k): v.get("windows", []) for k, v in r["pumps"].items()}
        for t in targets:
            hits = [p for p, ws in cache[ion].items() if any(w["lo_nm"] - 10 <= t <= w["hi_nm"] + 10 for w in ws)]
            near = min([min(abs(t - w["lo_nm"]), abs(t - w["hi_nm"])) for ws in cache[ion].values() for w in ws] or [None])
            e = dict(ion=ion, dopants=row["dopants"], ppm=row["ppm"], target_nm=t, description=res["description"],
                     reproduced=bool(hits), pumps_that_work=hits, nearest_window_edge_nm=None if hits else near,
                     out_of_model_range=t < 1e7 / (EMAX[ion] + 1000), src=row["src"])
            out["lvf"].append(e)
            print(f"LVF {ion:5s} {t:7.0f} nm  {'OK  ' if hits else 'MISS'} {hits}  {res['description']}"
                  + ("" if hits else f"  (nearest edge {near} nm)"))
json.dump(out, open("forward_results/vendor_validation.json", "w"), indent=1, default=float)
print("wrote forward_results/vendor_validation.json")
