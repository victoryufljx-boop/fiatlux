"""Level-by-level validation of the glass forward model against primary-paper numbers:
lifetimes (measured at low concentration, or intrinsic), radiative lifetimes, peak cross
sections, Stark levels, and the wavelengths of demonstrated fibre lasers.

    python scripts/validate_levels.py [--sites 30]  ->  forward_results/level_validation.json, .png
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np

from gainmat.forward.glass import (load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ,
                                   run, key_transitions, FAMILIES, EMAX, omegas_for)
from gainmat.forward.fiber import LocalModel

ap = argparse.ArgumentParser()
ap.add_argument("--sites", type=int, default=30)
ap.add_argument("--out", default="forward_results")
a = ap.parse_args()
load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ, verbose=False)

W92 = "Wetenkamp et al., JNCS 140, 35 (1992) Table 2 (0.1 mol%)"
FA07 = "Faure et al., JNCS 353, 2767 (2007) (<=200 ppm Tm)"
# (ion, family, level, kind, value, unit, src, note)
#   kind: tau = measured lifetime (low concentration / intrinsic); tau_rad = radiative lifetime
LIFETIMES = [
    ("Ho3+", "ZBLAN", "5I7", "tau", 12.01, "ms", W92, ""),
    ("Ho3+", "ZBLAN", "5I6", "tau", 3.5, "ms", W92, ""),
    ("Ho3+", "ZBLAN", "5F5", "tau", 0.29, "ms", W92, ""),
    ("Ho3+", "ZBLAN", "5S2", "tau", 0.30, "ms", W92, "5S2,5F4 thermalised in the paper"),
    ("Er3+", "ZBLAN", "4I13/2", "tau", 8.69, "ms", W92, ""),
    ("Er3+", "ZBLAN", "4I11/2", "tau", 6.71, "ms", W92, ""),
    ("Er3+", "ZBLAN", "4S3/2", "tau", 0.57, "ms", W92, ""),
    ("Nd3+", "ZBLAN", "4F3/2", "tau", 0.46, "ms", W92, ""),
    ("Tm3+", "ZBLAN", "3F4", "tau", 11.2, "ms", "Walsh & Barnes, APB 78, 325 (2004), 0.32 wt%", ""),
    ("Tm3+", "ZBLAN", "3H4", "tau", 1.51, "ms", "Walsh & Barnes 2004, 0.32 wt%", ""),
    ("Pr3+", "ZBLAN", "1G4", "tau", 0.110, "ms", "Ohishi et al., OL 16, 1747 (1991)", ""),
    ("Tm3+", "silica", "3H4", "tau", 0.014, "ms", FA07 + " (pure silica, via Jackson & King 1999)", ""),
    ("Tm3+", "aluminosilicate", "3H4", "tau", 0.037, "ms", FA07, "interpolated to 11.3 mol% AlO3/2"),
    ("Tm3+", "germanosilicate", "3H4", "tau", 0.021, "ms", FA07, "interpolated to 10 mol% GeO2"),
    ("Ho3+", "aluminosilicate", "5I7", "tau", 1.9, "ms", "Le Gouet et al., arXiv:2003.14188, Forster tau0",
     "intrinsic decay; 1/e lifetime 0.7 ms includes Ho-Ho transfer"),
    ("Tm3+", "aluminosilicate", "3F4", "tau", 0.65, "ms", "Agger & Povlsen, OE 14, 50 (2006), 8.4e19 cm^-3",
     "observed; quantum efficiency ~10%: not single-ion"),
    # radiative lifetimes
    ("Tm3+", "aluminosilicate", "3F4", "tau_rad", 6.0, "ms", "Agger & Povlsen 2006 (FL from measured sigma_e)", ""),
    ("Tm3+", "silica", "3F4", "tau_rad", 4.56, "ms", "Walsh & Barnes 2004 (JO)", ""),
    ("Tm3+", "silica", "3H4", "tau_rad", 0.67, "ms", "Walsh & Barnes 2004 (JO)", ""),
    ("Yb3+", "germanosilicate", "2F5/2", "tau_rad", 0.84, "ms", "Pask et al., JSTQE 1, 2 (1995)", ""),
    ("Pr3+", "ZBLAN", "1G4", "tau_rad", 3.24, "ms", "Ohishi 1991 (JO)", ""),
    ("Ho3+", "ZBLAN", "5I7", "tau_rad", 12.65, "ms", W92 + " (JO)", ""),
    ("Ho3+", "ZBLAN", "5I6", "tau_rad", 5.87, "ms", W92 + " (JO)", ""),
    ("Er3+", "ZBLAN", "4I13/2", "tau_rad", 9.52, "ms", W92 + " (JO)", ""),
    ("Er3+", "ZBLAN", "4I11/2", "tau_rad", 8.69, "ms", W92 + " (JO)", ""),
    ("Tm3+", "ZBLAN", "3F4", "tau_rad", 6.56, "ms", W92 + " (JO)", ""),
    ("Tm3+", "ZBLAN", "3H4", "tau_rad", 1.32, "ms", W92 + " (JO)", ""),
    ("Nd3+", "ZBLAN", "4F3/2", "tau_rad", 0.44, "ms", W92 + " (JO)", ""),
]
J12 = "Jackson, Nat. Photon. 6, 423 (2012) Fig. 4 (digitised)"
# (ion, family, upper, lower, kind 'em'|'abs', peak sigma cm^2, peak nm, src)
XS = [
    ("Er3+", "aluminosilicate", "4I13/2", "4I15/2", "em", 5.0e-21, 1533, J12),
    ("Tm3+", "aluminosilicate", "3F4", "3H6", "em", 5.5e-21, 1810, J12),
    ("Ho3+", "aluminosilicate", "5I7", "5I8", "em", 7.0e-21, 2030, J12),
    ("Er3+", "ZBLAN", "4I11/2", "4I13/2", "em", 5.6e-21, 2730, J12),
    ("Ho3+", "ZBLAN", "5I6", "5I7", "em", 5.2e-21, 2850, J12),
    ("Dy3+", "ZBLAN", "6H13/2", "6H15/2", "em", 2.95e-21, 2870, J12),
    ("Tm3+", "aluminosilicate", "3F4", "3H6", "abs", 4.2e-21, 1650, "Agger & Povlsen 2006 (Tm1)"),
    ("Tm3+", "aluminosilicate", "3F4", "3H6", "em", 3.76e-21, 1817, "Agger & Povlsen 2006 (Tm1, digitised)"),
    ("Tm3+", "silica", "3F4", "3H6", "em", 6.0e-21, 1837, "Jackson & King, JLT 17, 948 (1999) Fig. 2"),
    ("Tm3+", "silica", "3F4", "3H6", "abs", 4.47e-21, 1633, "Jackson & King 1999 Fig. 3"),
    ("Ho3+", "aluminosilicate", "5I7", "5I8", "abs", 4.3e-21, 1946, "Le Gouet et al. 2020"),
    ("Ho3+", "aluminosilicate", "5I7", "5I8", "em", 6.0e-21, 1950, "Le Gouet et al. 2020"),
    ("Yb3+", "germanosilicate", "2F5/2", "2F7/2", "em", 2.5e-20, 975, "Pask et al. 1995"),
    ("Pr3+", "ZBLAN", "1G4", "3H5", "em", 3.48e-21, 1322, "Ohishi et al. 1991"),
    ("Ho3+", "ZBLAN", "5I7", "5I8", "em", 0.56e-20, 2030, W92),
    ("Er3+", "ZBLAN", "4I13/2", "4I15/2", "em", 0.42e-20, 1540, W92),
    ("Tm3+", "ZBLAN", "3F4", "3H6", "em", 0.57e-20, 1660, W92 + " (wavelength = mean)"),
    ("Nd3+", "ZBLAN", "4F3/2", "4I11/2", "em", 2.64e-20, 1050, W92),
]
# Stark levels: (ion, family, label, 'lo'|'hi'|'split', value cm^-1, src)
STARK = [
    ("Er3+", "aluminosilicate", "4I13/2", "lo", 6544, "Desurvire & Simpson, OL 15, 547 (1990)"),
    ("Er3+", "aluminosilicate", "4I15/2", "split", 268, "Desurvire & Simpson 1990 (lower bound: 5 of 8 levels)"),
    ("Er3+", "silica", "4I15/2", "split", 445, "Desurvire 1990 Table 1, silicate glasses 395-500"),
    ("Er3+", "silica", "4I13/2", "split", 330, "Desurvire 1990 Table 1, silicate glasses 305-357"),
    ("Yb3+", "germanosilicate", "2F7/2", "split", 1490, "Pask et al. 1995 Fig. 1(a)"),
    ("Yb3+", "germanosilicate", "2F5/2", "split", 1370, "Pask et al. 1995 Fig. 1(a)"),
    ("Yb3+", "germanosilicate", "2F5/2", "lo", 10260, "Pask et al. 1995 Fig. 1(a)"),
    ("Tm3+", "aluminosilicate", "3H4", "lo", 12710, "Faure et al. 2007 Table 2 (GSA peak)"),
    ("Tm3+", "aluminosilicate", "3H5", "hi", 9090, "Faure et al. 2007 Table 2 (10% edge)"),
]

_cache = {}


def model(ion, fam):
    if (ion, fam) not in _cache:
        r = run(ion, fam, n_sites=a.sites, pumps=[])
        f = FAMILIES[fam]
        lm = LocalModel(ion, r["_spec"], omegas_for(ion, fam)[0], f["n"], f["hw_max"],
                        f.get("gap_family"), emax=EMAX[ion])
        _cache[(ion, fam)] = (r, lm)
    return _cache[(ion, fam)]


rows = []
for ion, fam, lev, kind, val, unit, src, note in LIFETIMES:
    r, lm = model(ion, fam)
    if lev not in lm.labels:
        continue
    i = lm.labels.index(lev)
    A = lm.A[i].sum()
    mod = 1e3 / A if kind == "tau_rad" else 1e3 / (A + lm.W[i])
    rows.append(dict(group="lifetime", ion=ion, family=fam, item=lev, kind=kind, measured=val,
                     model=float(mod), unit=unit, ratio=float(mod / val), src=src, note=note))
for ion, fam, u, l, kind, val, lam, src in XS:
    r, lm = model(ion, fam)
    kt = key_transitions(r, [(u, l)]).get(f"{u}->{l}")
    if not kt:
        continue
    s = kt["sigma_em"] if kind == "em" else kt["sigma_abs"]
    pk = kt["em_peak_nm"] if kind == "em" else kt["abs_peak_nm"]
    rows.append(dict(group="cross_section", ion=ion, family=fam, item=f"{u}->{l} {kind}", kind=kind,
                     measured=val, model=float(s), unit="cm2", ratio=float(s / val),
                     peak_meas_nm=lam, peak_model_nm=float(pk), src=src))
for ion, fam, lab, kind, val, src in STARK:
    r, lm = model(ion, fam)
    e = r["_spec"]["edges"].get(lab)
    if e is None:
        continue
    mod = {"lo": e[0], "hi": e[1], "split": r["split"][lab]}[kind]
    rows.append(dict(group="stark", ion=ion, family=fam, item=f"{lab} {kind}", kind=kind, measured=val,
                     model=float(mod), unit="cm-1", diff=float(mod - val), src=src))

# demonstrated fibre lasers (Jackson 2012 Table B1): is the laser wavelength inside a model window?
from gainmat.refs.midir_ho_refs import JACKSON2012
HOST = {"silicate": "aluminosilicate", "ZBLAN": "ZBLAN"}
for d in JACKSON2012["DEMONSTRATED_LASERS"]:
    ion = d["transition"].split()[0] + "3+"
    fam = HOST.get(d["host"])
    if fam is None or ion not in EMAX:
        continue
    lam = d["laser_um"] * 1e3; pump = d["pump_um"] * 1e3
    rr = run(ion, fam, n_sites=max(15, a.sites // 2), pumps=[pump])
    pr = list(rr["pumps"].values())[0]
    wins = pr.get("windows", [])
    hit = any(w["lo_nm"] - 10 <= lam <= w["hi_nm"] + 10 for w in wins)
    near = min([0.0 if hit else min(abs(lam - w["lo_nm"]), abs(lam - w["hi_nm"])) for w in wins] or [None]) \
        if wins else None
    rows.append(dict(group="laser_wavelength", ion=ion, family=fam, item=d["transition"], measured=lam,
                     pump_nm=pump, reproduced=bool(hit), nearest_edge_nm=near, error=pr.get("error"),
                     windows=[(w["lo_nm"], w["hi_nm"]) for w in wins], src="Jackson 2012 Table B1: " + d["primary"]))

json.dump(rows, open(f"{a.out}/level_validation.json", "w"), indent=1, default=float)

# ---------------- summary + figure
for g in ("lifetime", "cross_section", "stark", "laser_wavelength"):
    print(f"--- {g}")
    for x in rows:
        if x["group"] != g:
            continue
        if g == "laser_wavelength":
            print(f'{x["ion"]:5s} {x["family"]:15s} {x["item"]:28s} {x["measured"]:6.0f} nm  pump {x["pump_nm"]:5.0f}  '
                  f'{"OK" if x["reproduced"] else "MISS"} {x.get("error") or ""}')
        elif g == "stark":
            print(f'{x["ion"]:5s} {x["family"]:15s} {x["item"]:14s} meas {x["measured"]:7.0f}  model {x["model"]:7.0f}  '
                  f'diff {x["diff"]:+.0f}')
        else:
            extra = f' peak {x["peak_model_nm"]:.0f}/{x["peak_meas_nm"]:.0f} nm' if g == "cross_section" else ""
            print(f'{x["ion"]:5s} {x["family"]:15s} {x["item"]:20s} {x["kind"]:7s} meas {x["measured"]:9.3g}  '
                  f'model {x["model"]:9.3g}  x{x["ratio"]:.2f}{extra}')

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, 3, figsize=(15, 5))
for ax, g, kinds, lab in ((axs[0], "lifetime", ("tau",), "lifetime (ms), low conc. / intrinsic"),
                          (axs[1], "lifetime", ("tau_rad",), "radiative lifetime (ms)"),
                          (axs[2], "cross_section", ("em", "abs"), "peak cross section (cm$^2$)")):
    pts = [x for x in rows if x["group"] == g and x["kind"] in kinds]
    xs = np.array([p["measured"] for p in pts]); ys = np.array([p["model"] for p in pts])
    lo, hi = min(xs.min(), ys.min()) / 2, max(xs.max(), ys.max()) * 2
    ax.fill_between([lo, hi], [lo / 2, hi / 2], [lo * 2, hi * 2], color="0.9")
    ax.plot([lo, hi], [lo, hi], color="0.5", lw=1)
    for p in pts:
        c = {"ZBLAN": "tab:purple"}.get(p["family"], "tab:orange")
        ax.scatter(p["measured"], p["model"], color=c, edgecolor="k", lw=0.5, zorder=3)
        ax.annotate(f'{p["ion"][:-2]} {p["item"].split()[0]}' + ("" if p["family"] == "ZBLAN" else f' ({p["family"][:5]})'),
                    (p["measured"], p["model"]), fontsize=6.5, xytext=(3, 2), textcoords="offset points")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("measured"); ax.set_ylabel("model"); ax.set_title(lab, fontsize=10)
fig.suptitle("Level-by-level validation (purple = ZBLAN, orange = silica family; grey = within 2x)", fontsize=10)
fig.tight_layout(); fig.savefig(f"{a.out}/level_validation.png", dpi=140)
print("wrote", f"{a.out}/level_validation.json")
