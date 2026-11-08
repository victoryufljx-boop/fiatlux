"""Ce3+:LiYF4 -- embedded [CeF8]5- cluster, SA-CASSCF(1,12)/NEVPT2 over Ce 4f+5d.

* fast : Ce SBKJC 46-e-core ECP (4f,5s,5p,5d,6s in valence; Stevens-Basch-Krauss-
         Jasien-Cundari) + F def2-SVP, no X2C (ECP is scalar-relativistic).
* full : all-electron sfX2C-1e, Ce ANO-RCC-VTZP (shipped with PySCF), F ANO-RCC-VTZP
         (heavy; on Colab) -- or set CE_FULL_BASIS="x2c-tzvpall" for a cheaper AE run.

One SA-CASSCF over all 12 roots (equal weights), CASCI on the SA orbitals,
SC-NEVPT2 per root.  The 12x12 ab initio ligand-field matrix (4f+5d block,
including 4f-5d mixing allowed by the S4 site) is then diagonalised with
spin-orbit zeta_4f = 615, zeta_5d = 1082 cm-1 (van Pieterson et al. 2001/2002).
Absorption = 5d-character Kramers doublets minus the lowest (2F5/2-derived)
doublet.  Emission: vertical lowest-5d -> 2F5/2 and 2F7/2 at the *ground-state*
geometry (no excited-state relaxation -> no Stokes shift predicted; report the
experimental Stokes shift as a separate, documented correction if needed).
"""
from __future__ import annotations

import os
import time
import numpy as np

from . import embedding as emb
from .casci_lf import run_scf, sa_casscf_nevpt2, lf_matrix, add_spin_orbit, HARTREE_CM1
from ..refs.cluster_refs import CE_ZETA

LEVELS = {
    "fast": dict(basis={"Ce": "sbkjc", "F": "def2-svp"}, ecp={"Ce": "sbkjc"},
                 r_pc=12.0, x2c=False, shells=["Ce 4f", "Ce 5d"], max_macro=12),
    "full": dict(basis={"Ce": os.environ.get("CE_FULL_BASIS", "ano-rcc-vtzp"),
                        "F": os.environ.get("F_FULL_BASIS", "ano-rcc-vdzp")}, ecp={},
                 r_pc=20.0, x2c=True, shells=["Ce 4f", "Ce 5d"], max_macro=100),
}


def run(level="fast", verbose=0, log=print, do_nevpt2=True, partial_cap=False, casscf=True):
    p = LEVELS[level]
    t0 = time.time()
    cl = emb.make_cluster("LiYF4", "Y", "Ce", 3, r_qm=2.6, r_cap=4.0, r_pc=p["r_pc"],
                          partial_cap=partial_cap)
    log(f"[Ce] cluster: {cl.summary()}")
    mol, mf = run_scf(cl, p["basis"], p["ecp"], x2c=p["x2c"], verbose=verbose)
    log(f"[Ce] nao={mol.nao} ROHF E={mf.e_tot:.6f} conv={mf.converged} ({time.time()-t0:.1f}s)")
    res = sa_casscf_nevpt2(mol, mf, p["shells"], nroots=12, do_nevpt2=do_nevpt2,
                           verbose=verbose, max_macro=p["max_macro"], casscf=casscf)
    log(f"[Ce] SA-CASSCF {res['t_casscf']:.1f}s  NEVPT2 {res['t_nevpt2']:.1f}s  "
        f"active-orbital metal weights {np.round(res['proj_weights'], 3).tolist()}")
    out = dict(nao=mol.nao, cluster=cl.summary(), casscf_converged=res["converged"],
               basis=p["basis"], x2c=p["x2c"], zeta=dict(f=CE_ZETA["zeta_4f_cm1"],
                                                         d=CE_ZETA["zeta_5d_cm1"]))
    for key in ("e_cas", "e_pt2"):
        if key == "e_pt2" and not do_nevpt2:
            continue
        e = res[key]
        h, blocks, chars, smin = lf_matrix(mol, res, p["shells"], e)
        # scalar (no-SO) roots with 4f/5d character
        sc = (e - e.min()) * HARTREE_CM1
        dchar = chars[:, 1]
        w, v = add_spin_orbit(h, blocks, [CE_ZETA["zeta_4f_cm1"], CE_ZETA["zeta_5d_cm1"]])
        w = w - w[0]
        n = h.shape[0]
        nf = 2 * blocks[0][0] + 1
        dso = np.array([float((np.abs(v[2 * nf:, k]) ** 2).sum()) for k in range(2 * n)])
        kr = w[::2]; kd = dso[::2]                       # one member per Kramers pair
        f_lev = kr[kd < 0.5]; d_lev = kr[kd >= 0.5]
        out[key] = dict(scalar_cm1=sc.tolist(), scalar_5d_char=dchar.tolist(),
                        so_4f_cm1=f_lev.tolist(), so_5d_cm1=d_lev.tolist(),
                        f_span_scalar_cm1=float(np.ptp(np.sort(sc)[:7])),
                        lowdin_min_sv=smin)
        # 2F7/2 reference = mean of the upper 4 4f doublets; emission lowest 5d -> 2F5/2 / 2F7/2
        f5 = f_lev[:3]; f7 = f_lev[3:]
        out[key]["em_to_2F52_cm1"] = float(d_lev[0] - f5[0])
        out[key]["em_to_2F72_cm1"] = float(d_lev[0] - f7.min()) if len(f7) else None
    out["timings"] = dict(total_s=time.time() - t0, casscf_s=res["t_casscf"],
                          nevpt2_s=res["t_nevpt2"])
    return out
