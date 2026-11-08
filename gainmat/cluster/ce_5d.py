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

# Flags (all levels accept overrides as keyword arguments of run()):
#   cap_mode : 'pc' | 'partial' | 'aimp'   Y3+ next-shell treatment (see embedding.py)
#   r_cap    : 4.0 A -> 4 Li + 4 Y ; 5.3 A -> also the second set of 4 Y (8 Y in total)
#   relax    : None | 'breathing' | 'breathing2' (short/long Ce-F sets separately)
#   semicore : True -> extra NEVPT2 with Ce 5s5p frozen, reports the 4f->5d change
LEVELS = {
    # baseline (unchanged from the first version, used by the slow tests)
    "fast": dict(basis={"Ce": "sbkjc", "F": "def2-svp"}, ecp={"Ce": "sbkjc"},
                 r_pc=12.0, x2c=False, shells=["Ce 4f", "Ce 5d"], max_macro=12,
                 cap_mode="pc", r_cap=4.0, relax=None, semicore=False,
                 relax_grid=(0.0, 0.04, 0.08, 0.12)),
    # improved embedding + relaxation, same small basis (fits in ~15 min on 2 cores)
    "fast+": dict(basis={"Ce": "sbkjc", "F": "def2-svp"}, ecp={"Ce": "sbkjc"},
                  r_pc=12.0, x2c=False, shells=["Ce 4f", "Ce 5d"], max_macro=12,
                  cap_mode="aimp", r_cap=5.3, relax="breathing", semicore=True,
                  relax_grid=(0.0, 0.04, 0.08, 0.12)),
    # Colab: proper Ce basis (def2-TZVP, ECP28MWB, 4f/5s5p/5d/6s + g) and diffuse +
    # polarised F (def2-TZVPD), AIMP Y caps, 2-coordinate S4 breathing relaxation.
    "full": dict(basis={"Ce": os.environ.get("CE_FULL_BASIS", "bundled:ce_def2-tzvp.nw"),
                        "F": os.environ.get("F_FULL_BASIS", "def2-tzvpd")}, ecp={},
                 r_pc=20.0, x2c=False, shells=["Ce 4f", "Ce 5d"], max_macro=50,
                 cap_mode="aimp", r_cap=5.3, relax="breathing2", semicore=True,
                 relax_grid=(0.0, 0.05, 0.10)),
    # legacy all-electron route (sfX2C-1e ANO-RCC); no Y AIMP difference in the code path
    "full_ae": dict(basis={"Ce": "ano-rcc-vtzp", "F": "ano-rcc-vdzp"}, ecp={},
                    r_pc=20.0, x2c=True, shells=["Ce 4f", "Ce 5d"], max_macro=100,
                    cap_mode="aimp", r_cap=5.3, relax=None, semicore=False,
                    relax_grid=(0.0, 0.05, 0.10)),
}


def run(level="fast", verbose=0, log=print, do_nevpt2=True, partial_cap=False, casscf=True,
        **flags):
    p = dict(LEVELS[level]); p.update({k: v for k, v in flags.items() if v is not None})
    if partial_cap:
        p["cap_mode"] = "partial"
    t0 = time.time()
    cl = emb.make_cluster("LiYF4", "Y", "Ce", 3, r_qm=2.6, r_cap=p["r_cap"], r_pc=p["r_pc"],
                          cap_mode=p["cap_mode"])
    log(f"[Ce] cluster: {cl.summary()}")
    rel = None
    if p["relax"]:
        from . import relax as _rx
        modes = _rx.breathing_modes(cl, split_sets=(p["relax"] == "breathing2"))
        rel = _rx.relax(cl, p["basis"], p["ecp"], modes, grid=p["relax_grid"], log=log)
        cl = rel.pop("cl_gs"); rel.pop("cl_ex")
    mol, mf = run_scf(cl, p["basis"], p["ecp"], x2c=p["x2c"], verbose=verbose)
    log(f"[Ce] nao={mol.nao} ROHF E={mf.e_tot:.6f} conv={mf.converged} ({time.time()-t0:.1f}s)")
    res = sa_casscf_nevpt2(mol, mf, p["shells"], nroots=12, do_nevpt2=do_nevpt2,
                           verbose=verbose, max_macro=p["max_macro"], casscf=casscf)
    log(f"[Ce] SA-CASSCF {res['t_casscf']:.1f}s  NEVPT2 {res['t_nevpt2']:.1f}s  "
        f"active-orbital metal weights {np.round(res['proj_weights'], 3).tolist()}")
    out = dict(nao=mol.nao, cluster=cl.summary(), casscf_converged=res["converged"],
               basis={k: (v if isinstance(v, str) else "explicit") for k, v in p["basis"].items()},
               x2c=p["x2c"], flags={k: p[k] for k in ("cap_mode", "r_cap", "relax", "semicore")},
               relax=rel, zeta=dict(f=CE_ZETA["zeta_4f_cm1"],
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
    if p["semicore"] and do_nevpt2:
        from .casci_lf import nevpt2_frozen, semicore_indices
        t2 = time.time()
        fz = semicore_indices(mol, res["mo_sa"], res["ncore"])
        _, e_fz = nevpt2_frozen(mf, res["mo_sa"], res["ncore"], res["ncas"], 12, fz)
        h, blocks, chars, _ = lf_matrix(mol, res, p["shells"], e_fz)
        w, v = add_spin_orbit(h, blocks, [CE_ZETA["zeta_4f_cm1"], CE_ZETA["zeta_5d_cm1"]])
        w = w - w[0]; nf = 2 * blocks[0][0] + 1
        dso = np.array([float((np.abs(v[2 * nf:, k]) ** 2).sum()) for k in range(len(w))])
        d_fz = w[::2][dso[::2] >= 0.5]
        out["semicore"] = dict(n_frozen=len(fz), so_5d_frozen_5s5p_cm1=d_fz.tolist(),
                               delta_cm1=(np.array(out["e_pt2"]["so_5d_cm1"]) - d_fz).tolist(),
                               time_s=time.time() - t2)
        log(f"[Ce] 5s5p frozen ({len(fz)} orb): 5d {np.round(d_fz).tolist()}  "
            f"correlated-minus-frozen {np.round(out['semicore']['delta_cm1']).tolist()}")
    if rel is not None:
        key = "e_pt2" if do_nevpt2 else "e_cas"
        st = rel["stokes_scf_cm1"]
        out[key]["em_relaxed_to_2F52_cm1"] = out[key]["em_to_2F52_cm1"] - st
        out[key]["em_relaxed_to_2F72_cm1"] = (out[key]["em_to_2F72_cm1"] - st
                                              if out[key]["em_to_2F72_cm1"] else None)
        out[key]["stokes_cm1"] = st
    out["timings"] = dict(total_s=time.time() - t0, casscf_s=res["t_casscf"],
                          nevpt2_s=res["t_nevpt2"])
    return out
