"""Ti3+:Al2O3 -- embedded [TiO6]9- cluster, SA-CASSCF(1,5)/NEVPT2 on Ti 3d.

Vertical absorption: ground 2T2g-derived root (lowest, C3 splits t2g into a1+e)
to the two 2Eg-derived roots (the e_g pair; in C3 they are degenerate at the
crystal geometry -- the observed 490/560 nm doublet comes from the dynamic
Jahn-Teller effect in the excited state, which a vertical calculation at the
undistorted site cannot produce).  We therefore report the mean 2Eg energy as
the 490-nm band and do NOT predict the shoulder.  Spin-orbit (zeta_3d ~154
cm-1) is added through the ligand-field matrix for completeness; it barely
matters.  Emission/Stokes shift need the JT-relaxed 2Eg minimum (not done).
"""
from __future__ import annotations

import time
import numpy as np

from . import embedding as emb
from .casci_lf import run_scf, sa_casscf_nevpt2, lf_matrix, add_spin_orbit, HARTREE_CM1
from ..refs.cluster_refs import TI_ZETA

# relax: None | 'jt' -> ground-state (2T2g, UHF) and 2Eg (MOM Delta-SCF) minima on the
# (a1g breathing, e_g theta) grid; NEVPT2 is then rerun at both minima.
LEVELS = {
    "fast": dict(basis={"Ti": "def2-svp", "O": "def2-svp"}, r_pc=12.0, x2c=False,
                 relax=None, relax_grid=(-0.04, 0.0, 0.04)),
    "fast+": dict(basis={"Ti": "def2-svp", "O": "def2-svpd"}, r_pc=12.0, x2c=False,
                  relax="jt", relax_grid=(-0.04, 0.0, 0.04)),
    "full": dict(basis={"Ti": "def2-tzvp", "O": "def2-tzvpd"}, r_pc=20.0, x2c=True,
                 relax="jt", relax_grid=(-0.04, 0.0, 0.04)),
}


def _vertical(cl, p, verbose, log):
    mol, mf = run_scf(cl, p["basis"], {}, x2c=p["x2c"], verbose=verbose)
    res = sa_casscf_nevpt2(mol, mf, ["Ti 3d"], nroots=5, do_nevpt2=True, verbose=verbose)
    return mol, mf, res


def run(level="fast", verbose=0, log=print, do_nevpt2=True, **flags):
    p = dict(LEVELS[level]); p.update({k: v for k, v in flags.items() if v is not None})
    t0 = time.time()
    cl = emb.make_cluster("Al2O3", "Al", "Ti", 3, r_qm=2.2, r_cap=3.6, r_pc=p["r_pc"])
    log(f"[Ti] cluster: {cl.summary()}")
    rel = None
    if p["relax"]:
        from . import relax as _rx
        modes = [_rx.breathing_modes(cl)[0] / np.sqrt(6.0), _rx.eg_theta_mode(cl)]
        rel = _rx.relax(cl, p["basis"], {}, modes, grid=p["relax_grid"], log=log, skip_d=2)
        cl_ex = rel.pop("cl_ex"); cl = rel.pop("cl_gs")
        if do_nevpt2:
            _, _, rex = _vertical(cl_ex, p, verbose, log)
            e = (rex["e_pt2"] - rex["e_pt2"].min()) * HARTREE_CM1
            ec = (rex["e_cas"] - rex["e_cas"].min()) * HARTREE_CM1
            rel["at_Q_ex_nevpt2_cm1"] = e.tolist()
            rel["at_Q_ex_casscf_cm1"] = ec.tolist()
            rel["em_nevpt2_cm1"] = float(e[3])          # lower 2Eg component -> ground
            log(f"[Ti] at 2Eg minimum: NEVPT2 roots {np.round(e).tolist()}")
    mol, mf = run_scf(cl, p["basis"], {}, x2c=p["x2c"], verbose=verbose)
    log(f"[Ti] nao={mol.nao} ROHF E={mf.e_tot:.6f} conv={mf.converged} ({time.time()-t0:.1f}s)")
    res = sa_casscf_nevpt2(mol, mf, ["Ti 3d"], nroots=5, do_nevpt2=do_nevpt2, verbose=verbose)
    out = {}
    for key in ("e_cas", "e_pt2"):
        e = res[key]
        out[key + "_cm1"] = ((e - e.min()) * HARTREE_CM1).tolist()
        h, blocks, chars, smin = lf_matrix(mol, res, ["Ti 3d"], e)
        w, _ = add_spin_orbit(h, blocks, [TI_ZETA["zeta_3d_cm1"]])
        w = w[::2] - w[0]                                    # Kramers doublets
        out[key + "_so_cm1"] = w.tolist()
        out[key + "_d_char"] = chars[:, 0].tolist()
    e = np.array(out["e_pt2_cm1"] if do_nevpt2 else out["e_cas_cm1"])
    # t2g-like = 3 lowest, eg-like = 2 highest
    eg = e[3:5]
    out.update(abs_eg_mean_cm1=float(eg.mean()), abs_eg_cm1=eg.tolist(),
               t2g_split_cm1=float(e[2] - e[0]), timings=dict(total_s=time.time() - t0,
               casscf_s=res["t_casscf"], nevpt2_s=res["t_nevpt2"]),
               nao=mol.nao, cluster=cl.summary(), casscf_converged=res["converged"],
               basis=p["basis"], x2c=p["x2c"], relax=rel)
    if rel is not None and "em_nevpt2_cm1" in rel:
        out["em_cm1"] = rel["em_nevpt2_cm1"]
        out["stokes_cm1"] = out["abs_eg_mean_cm1"] - rel["em_nevpt2_cm1"]
    return out
