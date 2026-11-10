"""Energy-transfer microparameters predicted from the model's own cross sections (gainmat.forward.transfer)
vs values published for real glasses.  Two kinds of published value:
  overlap    : the authors' Forster-Dexter overlap of THEIR measured spectra (same formula, measured
               spectra) -> tests the model's line shapes and strengths where they matter for transfer
  decay fit  : a microparameter fitted to measured decay curves (Yokota-Tanimoto) -> a direct measurement
Two model variants: given the paper's Omegas and refractive index, and fully predictive (model defaults).

    python scripts/validate_transfer.py -> forward_results/transfer_validation.json
"""
import json, os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np
from gainmat.forward import glass as g
from gainmat.forward.transfer import c_da

g.load_refs(); g.apply_site_refs(); g.calibrate_spm(370.0, g.CARNALL_ER_LAF3_BKQ, verbose=False)
HUANG = "Huang, Liu, Hu & Chen, Sci. Rep. 4, 5053 (2014), Table 4"
HAN = "Han et al., Sci. Rep. 6, 31207 (2016), Table 4 (BGN, bismuth germanate)"
GEB = "Gebavi et al., J. Lumin. 132, 270 (2012) (75TeO2-20ZnO-5Na2O), Tables 6 and 8 (m^6/s x 1e12)"
CAO = "Cao et al., Sci. Rep. 6, 37873 (2016) (Er/Ho Ba-fluorosilicate, n = 1.60, phonon 952 cm^-1), text"
BALAJI = "Balaji et al., Sci. Rep. 6, 29203 (2016) (80TeO2-15(BaF2+BaO)-3La2O3), text"
SONT = "Sontakke, Biswas, Mandal & Annapurna (CGCRI), Nd Ba-Al-metaphosphate, Table 3 (CORE 228425814)"
CASES = [
    # (label, family, donor ion, donor (u,l) emission, acceptor ion, acceptor (u,l) absorption, published cm6/s, kind, src,
    #  paper overrides donor, paper overrides acceptor)
    ("Er 4I13/2 migration", "ZBLAN", "Er3+", ("4I13/2", "4I15/2"), "Er3+", ("4I13/2", "4I15/2"), 4.53e-39, "overlap", HUANG,
     dict(omegas=(3.27, 1.3, 1.75), n=1.499), None),
    ("Er 4I13/2 migration (ZBYA)", "ZBLAN", "Er3+", ("4I13/2", "4I15/2"), "Er3+", ("4I13/2", "4I15/2"), 3.99e-39, "overlap", HUANG,
     dict(omegas=(3.01, 1.03, 1.64), n=1.502), None),
    ("Tm 3H4 migration", "germanate", "Tm3+", ("3H4", "3H6"), "Tm3+", ("3H4", "3H6"), 3.54e-39, "overlap", HAN,
     dict(omegas=(6.26, 0.24, 1.6), n=2.1), None),
    ("Tm cross-relaxation 3H4,3H6->3F4,3F4", "germanate", "Tm3+", ("3H4", "3F4"), "Tm3+", ("3F4", "3H6"), 1.60e-39, "overlap", HAN,
     dict(omegas=(6.26, 0.24, 1.6), n=2.1), None),
    ("Tm 3F4 migration", "tellurite", "Tm3+", ("3F4", "3H6"), "Tm3+", ("3F4", "3H6"), 1.27e-38, "overlap", GEB,
     dict(omegas=(4.94, 2.15, 1.62), n=1.98), None),
    ("Tm 3F4 -> Ho 5I7", "tellurite", "Tm3+", ("3F4", "3H6"), "Ho3+", ("5I7", "5I8"), 5.7e-39, "overlap", GEB,
     dict(omegas=(4.94, 2.15, 1.62), n=1.98), dict(n=1.98)),
    ("Ho 5I7 -> Tm 3F4 (back)", "tellurite", "Ho3+", ("5I7", "5I8"), "Tm3+", ("3F4", "3H6"), 1.5e-40, "overlap", GEB,
     dict(n=1.98), dict(omegas=(4.94, 2.15, 1.62), n=1.98)),
    ("Tm cross-relaxation 3H4,3H6->3F4,3F4", "tellurite", "Tm3+", ("3H4", "3F4"), "Tm3+", ("3F4", "3H6"), 4.93e-39, "decay fit", GEB,
     dict(omegas=(4.94, 2.15, 1.62), n=1.98), None),
    ("Er 4I13/2 -> Ho 5I7", "aluminosilicate", "Er3+", ("4I13/2", "4I15/2"), "Ho3+", ("5I7", "5I8"), 5.52e-40, "overlap", CAO,
     dict(n=1.60, hw_max=952), dict(n=1.60, hw_max=952)),
    ("Er 4I11/2 -> Ho 5I6", "aluminosilicate", "Er3+", ("4I11/2", "4I15/2"), "Ho3+", ("5I6", "5I8"), 0.66e-40, "overlap", CAO,
     dict(n=1.60, hw_max=952), dict(n=1.60, hw_max=952)),
    ("Yb 2F5/2 -> Ho 5I6", "tellurite", "Yb3+", ("2F5/2", "2F7/2"), "Ho3+", ("5I6", "5I8"), 1.38e-40, "overlap", BALAJI,
     None, None),
    ("Nd 4F3/2 migration", "phosphate", "Nd3+", ("4F3/2", "4I9/2"), "Nd3+", ("4F3/2", "4I9/2"), 3.5e-39, "overlap", SONT,
     dict(n=1.552), None),
    ("Nd self-quenching 4F3/2+4I9/2 -> 4I15/2+4I15/2", "phosphate", "Nd3+", ("4F3/2", "4I15/2"), "Nd3+", ("4I15/2", "4I9/2"),
     5.5e-42, "decay fit", SONT, dict(n=1.552), None),
]
_sp = {}


def spec(ion, fam, ov):
    k = (ion, fam, json.dumps(ov, sort_keys=True, default=str))
    if k not in _sp:
        _sp[k] = g.run(ion, fam, n_sites=60, pumps=[], overrides=ov)["_spec"]
    return _sp[k]


rows = []
for lab, fam, iD, pD, iA, pA, meas, kind, src, ovD, ovA in CASES:
    out = dict(process=lab, family=fam, published_cm6_s=meas, kind=kind, src=src)
    for variant, oD, oA in (("paper Omega", ovD, ovA if iA != iD else ovD), ("model defaults", None, None)):
        r = c_da(iD, fam, pD, pA, spec=spec(iD, fam, oD), spec_acceptor=spec(iA, fam, oA), ion_acceptor=iA, overrides=oD)
        out[variant] = dict(model_cm6_s=r["total"], ratio=r["total"] / meas, log10_ratio=float(np.log10(r["total"] / meas)),
                            phonon_fractions=[round(f, 3) for f in r["fractions"]])
    rows.append(out)
    print(f"{lab:42s} {fam:10s} {kind:9s} pub {meas:.2e} | paper-Omega x{out['paper Omega']['ratio']:.2f} | defaults x{out['model defaults']['ratio']:.2f}", flush=True)
summ = {}
for kind in ("overlap", "decay fit"):
    for variant in ("paper Omega", "model defaults"):
        v = [abs(r[variant]["log10_ratio"]) for r in rows if r["kind"] == kind]
        summ[f"{kind} / {variant}"] = dict(n=len(v), median_abs_dex=float(np.median(v)), rms_dex=float(np.sqrt(np.mean(np.square(v)))),
                                           within_2x=float(np.mean(np.array(v) <= np.log10(2))))
print(json.dumps(summ, indent=1))
json.dump(dict(summary=summ, rows=rows), open("forward_results/transfer_validation.json", "w"), indent=1, default=float)
