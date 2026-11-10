"""Macroscopic energy-transfer coefficients (cm^3/s) of Er3+ in ZBLAN predicted from the model's own spectra
(Forster-Dexter microparameters, gainmat.forward.transfer) through the Burshtein hopping (migration-assisted)
model, k = pi (2 pi/3)^(5/2) sqrt(C_DA C_DD) N_Er, vs the rate-equation coefficients measured by Golding,
Jackson, King & Pollnau, PRB 62, 856 (2000), Table I (Er ZBLAN, 2-14e20 cm^-3).  Nothing fitted.

Direct (static Forster) transfer adds k_dir = (4/3) pi^(3/2) sqrt(C_DA / tau_D) (the Inokuti-Hirayama
exp(-gamma sqrt t) term expressed as a rate at t ~ tau_D, independent of N), with tau_D the model's own lifetime of
the donor level; k = k_dir + k_hop.

    python scripts/validate_transfer_rates.py -> forward_results/transfer_rates_validation.json
"""
import sys, os, json, math, numpy as np
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from gainmat.forward import glass as g
from gainmat.forward.transfer import c_da
from gainmat.refs.er_pr_zblan_refs import ENERGY_TRANSFER as ET
g.load_refs(); g.apply_site_refs(); g.calibrate_spm(370.0, g.CARNALL_ER_LAF3_BKQ, verbose=False)
r0 = g.run("Er3+", "ZBLAN", n_sites=60, pumps=[]); sp = r0["_spec"]
from gainmat.forward.fiber import LocalModel
_F = g.FAMILIES["ZBLAN"]; _lm = LocalModel("Er3+", sp, r0["omegas"], _F["n"], _F["hw_max"], _F.get("gap_family"), emax=g.EMAX["Er3+"])
def tau_model(level):
    i = _lm.labels.index(level); return 1.0 / (float(_lm.A[i].sum()) + float(_lm.W[i]))
CASES = [  # name in refs, donor emission, acceptor absorption, migration (donor) transition
    ("ETU1 (Er-Er upconversion from 4I13/2)", ("4I13/2", "4I15/2"), ("4I9/2", "4I13/2"), ("4I13/2", "4I15/2")),
    ("ETU2 (Er-Er upconversion from 4I11/2)", ("4I11/2", "4I15/2"), ("4F7/2", "4I11/2"), ("4I11/2", "4I15/2")),
    ("CR (4S3/2/2H11/2 self-quenching), channel a", None, None, ("4S3/2", "4I15/2")),
]
# CR from the thermalised 4S3/2 + 2H11/2 pair (Golding's lumped W_CR): channel a (2H11/2,4I15/2)->(4I9/2,4I13/2)
# weighted by the 2H11/2 Boltzmann share, plus channel b (4S3/2,4I15/2)->(4I13/2,4I9/2) (phonon-absorbing).
def cr_lumped():
    b = sp["bary"]; dE = b["2H11/2"] - b["4S3/2"]
    fH = 12 * math.exp(-dE / (0.695 * 295)) / (4 + 12 * math.exp(-dE / (0.695 * 295)))
    ca = c_da("Er3+", "ZBLAN", ("2H11/2", "4I9/2"), ("4I13/2", "4I15/2"), spec=sp)["total"]
    cb = c_da("Er3+", "ZBLAN", ("4S3/2", "4I13/2"), ("4I9/2", "4I15/2"), spec=sp)["total"]
    return fH * ca + (1 - fH) * cb, dict(f_2H11=fH, C_a=ca, C_b=cb)
K = math.pi * (2 * math.pi / 3) ** 2.5
rows = []
for name, d, a, mig in CASES:
    ref = next(e for e in ET if e["name"] == name)
    if d is None:
        cda, info = cr_lumped(); print("   CR channels", {k: f"{v:.2e}" for k, v in info.items()})
    else:
        cda = c_da("Er3+", "ZBLAN", d, a, spec=sp)["total"]
    cdd = c_da("Er3+", "ZBLAN", mig, mig, spec=sp)["total"]
    for N, kmeas in sorted(ref["k_cm3_s"].items()):
        tD = tau_model(mig[0])
        k_hop = K * math.sqrt(cda * cdd) * N
        k_dir = (4.0 / 3.0) * math.pi ** 1.5 * math.sqrt(cda / tD)
        k = k_hop + k_dir
        rows.append(dict(process=name, N_cm3=N, C_DA=cda, C_DD=cdd, tau_D_s=tD, k_hop=k_hop, k_direct=k_dir, k_model=k,
                         k_measured=kmeas, ratio=k / kmeas, ratio_hop_only=k_hop / kmeas))
        print(f"{name[:40]:40s} N {N:.1e}  k_hop {k_hop:.2e} k_dir {k_dir:.2e} -> {k:.2e}  measured {kmeas:.2e}  x{k/kmeas:.2f} (hop only x{k_hop/kmeas:.2f})", flush=True)
lr = [abs(math.log10(r["ratio"])) for r in rows]
summ = dict(n=len(rows), median_abs_dex=float(np.median(lr)), within_2x=float(np.mean(np.array(lr) <= math.log10(2))), within_3x=float(np.mean(np.array(lr) <= math.log10(3))))
for proc in {r["process"] for r in rows}:
    v = [abs(math.log10(r["ratio"])) for r in rows if r["process"] == proc]
    summ[proc] = dict(median_abs_dex=float(np.median(v)), ratios=[round(r["ratio"], 2) for r in rows if r["process"] == proc])
print(summ)
json.dump(dict(src="Golding et al., PRB 62, 856 (2000), Table I", model="Burshtein hopping, k = pi(2pi/3)^2.5 sqrt(C_DA C_DD) N", summary=summ, rows=rows),
          open("forward_results/transfer_rates_validation.json", "w"), indent=1, default=float)
