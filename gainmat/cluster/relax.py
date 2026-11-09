"""Relaxation of the first anion shell inside the frozen embedding (ground and excited
state) -> emission energies and Stokes shift.

Approach (configuration-coordinate, as in Seijo & Barandiaran's a1g / e_g scans):

* Ligand displacements are restricted to a few symmetry-adapted coordinates
  (radial breathing of the whole shell, optionally split into the two S4 Ce-F sets for
  LiYF4, and the e_g(theta) Jahn-Teller coordinate for the Ti3+ octahedron).
  Cations of the embedding (ECP caps, AIMP ions, point charges) stay fixed.
* Ground state: UHF energy of the cluster (4f1 / t2g1) on a small grid, quadratic fit.
* Excited state (lowest 5d / 2Eg): Delta-SCF UHF with the maximum-overlap method
  (MOM): the open-shell electron is moved into the lowest virtual with >50 % dopant-d
  character and the orbitals are re-optimised at every geometry.  Because MOM breaks
  the e_g degeneracy, the Ti 2Eg minimum found this way *is* the static Jahn-Teller
  minimum along e_theta.
* The SCF curves give the geometries Q_gs and Q_ex and the two relaxation energies
  (E_ex(Q_gs)-E_ex(Q_ex) and E_gs(Q_ex)-E_gs(Q_gs)); the vertical energies at Q_gs
  (absorption) and Q_ex (emission) are then recomputed with SA-CASSCF/NEVPT2 when
  ``casscf_at_minima`` is used, or emission = absorption - Stokes(SCF) otherwise.
"""
from __future__ import annotations

import copy
import time
import numpy as np

from . import embedding as emb

HARTREE_CM1 = 219474.6313705


def displaced(cl, coords_disp):
    """Copy of ``cl`` with QM ligand coordinates replaced (dopant stays at origin)."""
    c2 = copy.copy(cl)
    c2.qm_coords = np.array(coords_disp, float)
    if hasattr(cl, "_pc_formal_sum"):
        c2._pc_formal_sum = cl._pc_formal_sum
    return c2


# ------------------------------------------------------------- coordinates
def breathing_modes(cl, split_sets=False, tol=0.02):
    """Radial unit displacements of the ligands.  split_sets: one coordinate per group of
    symmetry-equivalent (same distance) ligands, e.g. LiYF4 short/long Ce-F."""
    X = cl.qm_coords[1:]
    r = np.linalg.norm(X, axis=1)
    u = X / r[:, None]
    if not split_sets:
        return [np.vstack([np.zeros(3), u])]
    groups = []
    for i in np.argsort(r):
        for g in groups:
            if abs(r[g[0]] - r[i]) < tol:
                g.append(i); break
        else:
            groups.append([i])
    modes = []
    for g in groups:
        d = np.zeros_like(X); d[g] = u[g]
        modes.append(np.vstack([np.zeros(3), d]))
    return modes


def eg_theta_mode(cl):
    """Octahedral e_g(theta) radial mode: the most trans-like ligand pair along the local
    z axis moves out by +2 units, the other two pairs in by -1 (normalised later)."""
    X = cl.qm_coords[1:]
    u = X / np.linalg.norm(X, axis=1)[:, None]
    n = len(u)
    used, pairs = set(), []
    for i in range(n):
        if i in used:
            continue
        j = min((k for k in range(n) if k != i and k not in used), key=lambda k: u[i] @ u[k])
        pairs.append((i, j)); used |= {i, j}
    assert len(pairs) == 3, "e_g mode needs an octahedron"
    d = np.zeros_like(X)
    for p, w in zip(pairs, (2.0, -1.0, -1.0)):
        for k in p:
            d[k] = w * u[k]
    return np.vstack([np.zeros(3), d / np.sqrt(12.0)])


# ------------------------------------------------------------- energies
def _dopant_d_weight(mol, C):
    S = mol.intor("int1e_ovlp")
    labs = mol.ao_labels(fmt=False)
    idx = [i for i, (ia, s, sh, c) in enumerate(labs) if ia == 0 and sh.endswith("d")]
    return np.einsum("ik,ij,jk->k", C[idx], S[idx], C)


def scf_energies(cl, basis, ecp, spin=1, excited=True, verbose=0, conv_tol=1e-7, skip_d=0):
    """(E_gs, E_ex, info) for one geometry: UHF ground state and MOM Delta-SCF state with the
    open-shell electron in the lowest dopant-d virtual (after skipping ``skip_d`` d-like
    virtuals: 2 for a d1 ion, whose empty t2g partners come first)."""
    from pyscf import scf, qmmm
    basis, ecp = emb.resolve_basis(basis, ecp)
    mol, mmx, mmq = emb.build_mol(cl, basis, ecp, spin=spin, verbose=verbose)
    mf = qmmm.mm_charge(scf.UHF(mol), mmx, mmq, unit="Angstrom")
    mf = emb.attach_aimp(mf, mol, cl)
    mf.conv_tol = conv_tol
    mf.max_cycle = 150
    mf.level_shift = 0.2
    dm0 = emb.qm_only_guess(mol, len(cl.qm_symbols), spin)
    mf.kernel(dm0=dm0)
    if not mf.converged:
        mf = mf.newton(); mf.kernel(mf.make_rdm1())
    e_gs = mf.e_tot
    if not excited:
        return e_gs, None, {}
    Ca, Cb = mf.mo_coeff
    occa, occb = mf.mo_occ.copy()
    wd = _dopant_d_weight(mol, Ca)
    homo = np.where(occa > 0)[0].max()
    vv = np.where(occa == 0)[0]
    # the (skip_d + 2) most d-like virtuals (Ti d1: 2 empty t2g partners + 2 eg; Ce: 5d),
    # ordered by orbital energy; target = the first after the skipped ones
    virt = sorted(vv[np.argsort(-wd[vv])][:skip_d + 2], key=lambda k: mf.mo_energy[0][k])
    tgt = virt[skip_d]
    occa2 = occa.copy(); occa2[homo] = 0; occa2[tgt] = 1
    mf2 = qmmm.mm_charge(scf.UHF(mol), mmx, mmq, unit="Angstrom")
    mf2 = emb.attach_aimp(mf2, mol, cl)
    mf2 = scf.addons.mom_occ(mf2, mf.mo_coeff, (occa2, occb))
    mf2.conv_tol = conv_tol
    mf2.max_cycle = 80
    mf2.kernel(dm0=mf2.make_rdm1(mf.mo_coeff, (occa2, occb)))
    dE = (mf2.e_tot - e_gs) * HARTREE_CM1
    if not (mf2.converged and 2000 < dE < 80000):
        # MOM collapsed / ran away: retry with damping and no DIIS from the same guess
        mf2.diis = None; mf2.damp = 0.5; mf2.max_cycle = 60
        mf2.kernel(dm0=mf2.make_rdm1(mf.mo_coeff, (occa2, occb)))
        dE = (mf2.e_tot - e_gs) * HARTREE_CM1
    ok = bool(mf2.converged and 2000 < dE < 80000)
    Ca2 = mf2.mo_coeff[0]
    wocc = _dopant_d_weight(mol, Ca2[:, mf2.mo_occ[0] > 0])
    return e_gs, (mf2.e_tot if ok else np.nan), dict(d_weight_ex=float(wocc.max()), conv_ex=ok,
                                 dE_vert_cm1=float((mf2.e_tot - e_gs) * HARTREE_CM1))


def _quad_fit(Q, E):
    """Least-squares quadratic in len(Q[0]) variables -> (Qmin, Emin, coeffs)."""
    Q = np.atleast_2d(Q); E = np.asarray(E, float)
    good = np.isfinite(E); Q, E = Q[good], E[good]
    n = Q.shape[1]
    cols = [np.ones(len(Q))] + [Q[:, i] for i in range(n)] + \
           [Q[:, i] * Q[:, j] for i in range(n) for j in range(i, n)]
    A = np.array(cols).T
    c, *_ = np.linalg.lstsq(A, E, rcond=None)
    g = c[1:1 + n]
    H = np.zeros((n, n)); k = 1 + n
    for i in range(n):
        for j in range(i, n):
            if i == j:
                H[i, i] = 2 * c[k]
            else:
                H[i, j] = H[j, i] = c[k]
            k += 1
    qmin = -np.linalg.solve(H, g)
    emin = c[0] + g @ qmin + 0.5 * qmin @ H @ qmin
    return qmin, float(emin), (c, H)


def relax(cl, basis, ecp, modes, grid=(-0.04, 0.0, 0.04, 0.08), center=None, log=print,
          spin=1, skip_d=0):
    """Scan ``modes`` (list of displacement fields, unit amplitude in Angstrom of radial
    motion) on a product grid around ``center``; fit quadratics to the ground (UHF) and
    MOM excited energies.  Returns dict with Q_gs, Q_ex, relaxation energies (cm-1),
    relaxed clusters and the raw scan."""
    import itertools
    t0 = time.time()
    n = len(modes)
    center = np.zeros(n) if center is None else np.asarray(center, float)
    base = cl.qm_coords.copy()
    pts = [center + np.array(p) for p in itertools.product(grid, repeat=n)]
    Eg, Ee, info = [], [], []
    for q in pts:
        X = base + sum(qi * m for qi, m in zip(q, modes))
        eg, ee, inf = scf_energies(displaced(cl, X), basis, ecp, spin=spin, skip_d=skip_d)
        Eg.append(eg); Ee.append(ee); info.append(inf)
        log(f"[relax] Q={np.round(q, 3).tolist()} E_gs={eg:.6f} E_ex={ee:.6f} "
            f"dE={inf['dE_vert_cm1']:.0f} cm-1 d_w={inf['d_weight_ex']:.2f}")
    pts = np.array(pts)
    qg, eg_min, (cg, _) = _quad_fit(pts, Eg)
    qe, ee_min, (ce, _) = _quad_fit(pts, Ee)

    def ev(c, q):
        q = np.atleast_1d(q)
        v = [1.0] + list(q) + [q[i] * q[j] for i in range(n) for j in range(i, n)]
        return float(np.dot(c, v))
    lam_ex = (ev(ce, qg) - ev(ce, qe)) * HARTREE_CM1      # excited-state relaxation
    lam_gs = (ev(cg, qe) - ev(cg, qg)) * HARTREE_CM1      # ground-state destabilisation
    span = (min(grid), max(grid))
    warn = [f"minimum outside scanned grid (Q={np.round(q, 3).tolist()})"
            for q in (qg, qe) if np.any(q < center + span[0] - 0.02) or np.any(q > center + span[1] + 0.02)]
    mk = lambda q: displaced(cl, base + sum(qi * m for qi, m in zip(q, modes)))
    out = dict(Q_gs=qg.tolist(), Q_ex=qe.tolist(), lambda_ex_cm1=lam_ex, lambda_gs_cm1=lam_gs,
               stokes_scf_cm1=lam_ex + lam_gs,
               dE_scf_abs_cm1=(ev(ce, qg) - ev(cg, qg)) * HARTREE_CM1,
               dE_scf_em_cm1=(ev(ce, qe) - ev(cg, qe)) * HARTREE_CM1,
               scan=dict(Q=pts.tolist(), E_gs=Eg, E_ex=Ee, info=info), warnings=warn,
               time_s=time.time() - t0)
    out["cl_gs"], out["cl_ex"] = mk(qg), mk(qe)
    d0 = np.linalg.norm(base[1:], axis=1)
    out["dist_gs_A"] = np.round(np.linalg.norm(out["cl_gs"].qm_coords[1:], axis=1), 4).tolist()
    out["dist_ex_A"] = np.round(np.linalg.norm(out["cl_ex"].qm_coords[1:], axis=1), 4).tolist()
    out["dist_host_A"] = np.round(d0, 4).tolist()
    log(f"[relax] Q_gs={np.round(qg, 4)} Q_ex={np.round(qe, 4)} lambda_ex={lam_ex:.0f} "
        f"lambda_gs={lam_gs:.0f} Stokes(SCF)={lam_ex + lam_gs:.0f} cm-1 ({out['time_s']:.0f}s)")
    return out
