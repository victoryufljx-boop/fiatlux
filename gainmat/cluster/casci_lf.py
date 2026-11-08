"""Shared machinery: embedded ROHF -> SA-CASSCF(1,n) -> per-root NEVPT2 ->
ab initio ligand-field (AILFT-style) one-electron Hamiltonian -> + spin-orbit.

For a single active electron the CASCI roots are one-electron states
|Psi_i> = sum_a c_ia |phi_a>, so the effective one-electron ligand-field
operator in the active space is h = sum_i E_i |c_i><c_i| (E_i = CASSCF or
NEVPT2 root energies).  It is rotated to pure metal real spherical harmonics
(Lowdin-orthonormalised projection of the active orbitals onto the metal nl
AO shells), and spin-orbit  sum_l zeta_l l.s  is added inside each l block
with literature zeta values (the SCF/CASSCF is scalar-relativistic only).
"""
from __future__ import annotations

import time
import numpy as np

from . import embedding as emb

HARTREE_CM1 = 219474.6313705
HARTREE_EV = 27.211386245988


# ------------------------------------------------------------- angular momentum
def l_matrices_real(l):
    """(Lx, Ly, Lz) in PySCF's real-spherical-harmonic order (m=-l..l)."""
    from pyscf.symm.sph import sph_pure2real
    m = np.arange(-l, l + 1)
    Lz = np.diag(m).astype(complex)
    Lp = np.zeros((2 * l + 1, 2 * l + 1), complex)
    for k in range(2 * l):            # <m+1|L+|m>
        mm = m[k]
        Lp[k + 1, k] = np.sqrt(l * (l + 1) - mm * (mm + 1))
    Lm = Lp.conj().T
    Lx, Ly = (Lp + Lm) / 2, (Lp - Lm) / (2j)
    u = sph_pure2real(l, reorder_p=False) if l == 1 else sph_pure2real(l)
    return [u.conj().T @ L @ u for L in (Lx, Ly, Lz)]


def add_spin_orbit(h_real, blocks, zetas):
    """h_real: (n,n) real one-electron LF matrix (cm-1) in real harmonics, ordered by
    ``blocks`` = [(l, start), ...].  Returns sorted eigenvalues (2n) and eigenvectors."""
    n = h_real.shape[0]
    H = np.kron(h_real.astype(complex), np.eye(2))
    S = [np.array([[0, 1], [1, 0]]) / 2, np.array([[0, -1j], [1j, 0]]) / 2,
         np.array([[1, 0], [0, -1]]) / 2]
    for (l, st), z in zip(blocks, zetas):
        L = l_matrices_real(l)
        P = np.zeros((n, 2 * l + 1)); P[st:st + 2 * l + 1, :] = np.eye(2 * l + 1)
        for La, Sa in zip(L, S):
            H += z * np.kron(P @ La @ P.T, Sa)
    w, v = np.linalg.eigh(H)
    return w, v


def _lowdin(S):
    """Plain symmetric Lowdin orthogonalisation S^-1/2 (keeps AO labels)."""
    w, v = np.linalg.eigh(S)
    return (v / np.sqrt(w)) @ v.T


# ------------------------------------------------------------- SCF / CAS
def run_scf(cl, basis, ecp, x2c=False, spin=1, verbose=0, max_memory=4000, conv_tol=1e-8):
    from pyscf import scf, qmmm
    basis, ecp = emb.resolve_basis(basis, ecp)
    mol, mmx, mmq = emb.build_mol(cl, basis, ecp, spin=spin, verbose=verbose,
                                  max_memory=max_memory)
    mf = scf.ROHF(mol)
    if x2c:
        mf = mf.sfx2c1e()
    mf = qmmm.mm_charge(mf, mmx, mmq, unit="Angstrom")
    mf = emb.attach_aimp(mf, mol, cl)
    mf.conv_tol = conv_tol
    mf.max_cycle = 200
    mf.level_shift = 0.2
    dm0 = emb.qm_only_guess(mol, len(cl.qm_symbols), spin)
    mf.kernel(dm0=dm0)
    if not mf.converged:
        mf = mf.newton()
        mf.kernel(mf.make_rdm1())
    return mol, mf


def metal_ao_indices(mol, shell_label):
    """AO indices of the first contracted shell matching e.g. 'Ce 4f' / 'Ti 3d' on atom 0."""
    labs = mol.ao_labels(fmt=False)
    el, nl = shell_label.split()
    idx = [i for i, (ia, s, sh, comp) in enumerate(labs) if ia == 0 and sh == nl]
    if not idx:      # basis labels shells differently (e.g. ECP bases start at 5d): fall back
        lch = nl[-1]
        cand = sorted({sh for (ia, s, sh, c) in labs if ia == 0 and sh.endswith(lch)},
                      key=lambda x: int(x[:-1]))
        idx = [i for i, (ia, s, sh, comp) in enumerate(labs) if ia == 0 and sh == cand[0]]
    return idx


def atomic_reference(mol, shells, ion_charge=3, spin=1):
    """Pure-l, real-harmonic reference orbitals (AO coefficients on atom 0, padded to
    mol.nao) for ``shells`` (e.g. ['Ce 4f','Ce 5d']).  The radial part is the
    dominant contraction of the corresponding free-ion ROHF orbitals (same basis/ECP),
    so the reference is the actual nl radial function, not a single primitive shell.
    Returns (R[nao, nref], blocks[(l, start)])."""
    from pyscf import gto, scf
    a0 = mol._atom[0][0]
    m = gto.Mole()
    m.atom = [[a0, (0, 0, 0)]]
    m.basis = {a0: mol._basis[a0]}
    m.ecp = {a0: mol._ecp[a0]} if a0 in mol._ecp else {}
    m.charge = ion_charge
    m.spin = spin
    m.verbose = 0
    m.build(parse_arg=False)
    mf = scf.ROHF(m)
    mf.max_cycle = 100
    mf.kernel()
    labs = m.ao_labels(fmt=False)
    R, blocks, st = [], [], 0
    for sh in shells:
        lch = sh.split()[1][-1]
        l = "spdfg".index(lch)
        # AO indices of all shells with this l on atom 0, grouped by contraction
        shl = sorted({x[2] for x in labs if x[2].endswith(lch)}, key=lambda q: int(q[:-1]))
        groups = [[i for i, x in enumerate(labs) if x[2] == q] for q in shl]
        G = np.array(groups)                               # (ncontr, 2l+1)
        # orbitals with this l character: pick the 2l+1 lowest non-doubly-occupied ones
        # whose weight in this l block is largest (4f: SOMO+virtual f; 5d: lowest virtual d)
        S = m.intor("int1e_ovlp")
        w = np.array([(mf.mo_coeff[G.ravel(), k] * (S[G.ravel()] @ mf.mo_coeff[:, k])).sum()
                      for k in range(m.nao)])
        cand = [k for k in np.argsort(mf.mo_energy) if mf.mo_occ[k] < 1.5 and w[k] > 0.9]
        sel = cand[:2 * l + 1]
        A = mf.mo_coeff[G, :][:, :, sel]                   # (ncontr, 2l+1, 2l+1)
        U, sv, Vt = np.linalg.svd(A.reshape(len(shl), -1), full_matrices=False)
        r = U[:, 0]
        for mm in range(2 * l + 1):
            v = np.zeros(mol.nao)
            v[G[:, mm]] = r
            R.append(v)
        blocks.append((l, st)); st += 2 * l + 1
    R = np.array(R).T
    # normalise each in the cluster metric
    Sm = mol.intor("int1e_ovlp")
    R /= np.sqrt(np.einsum("ik,ij,jk->k", R, Sm, R))
    return R, blocks


def select_active(mol, mf, shells):
    """AVAS-like active space: within the span of the non-doubly-occupied ROHF
    orbitals, take the ncas linear combinations with maximal overlap with the
    (symmetric-Lowdin-orthogonalised) metal AOs of ``shells`` (SVD of the
    projector).  Returns (mo_coeff, ncore, ncas, singular values)."""
    S = mol.intor("int1e_ovlp")
    Rf, _ = atomic_reference(mol, shells)
    Rf = Rf @ _lowdin(Rf.T @ S @ Rf)                 # orthonormal references
    occ = mf.mo_occ
    core = np.where(occ > 1.5)[0]
    cand = np.where(occ < 1.5)[0]
    Cc = mf.mo_coeff[:, cand]
    O = Rf.T @ S @ Cc                               # (nref, ncand)
    idx = list(range(Rf.shape[1]))
    U, sv, Vt = np.linalg.svd(O, full_matrices=True)
    ncas = len(idx)
    Crot = Cc @ Vt.T                                # first ncas = active
    mo = np.hstack([mf.mo_coeff[:, core], Crot])
    return mo, len(core), ncas, sv[:ncas]


def sa_casscf_nevpt2(mol, mf, shells, nroots, do_nevpt2=True, verbose=0, frozen=None,
                     log=print, max_macro=100, casscf=True):
    """SA-CASSCF(1, ncas) equal weights over ``nroots`` doublet roots, then CASCI on the
    SA orbitals and strongly-contracted NEVPT2 for each root."""
    from pyscf import mcscf, mrpt
    mo, ncore, ncas, w = select_active(mol, mf, shells)
    t0 = time.time()
    mc = mcscf.CASSCF(mf, ncas, (1, 0))
    mc.verbose = verbose
    if frozen is not None:
        mc.frozen = frozen
    mc.fcisolver.nroots = nroots
    mc = mc.state_average_([1.0 / nroots] * nroots)
    mc.conv_tol = 1e-7
    mc.max_cycle_macro = max_macro
    if casscf:
        mc.kernel(mo)
        mo_sa = mc.mo_coeff
    else:                                  # CASCI on ROHF orbitals (debug/very fast)
        mo_sa = mo
        mc.converged = False
    t_cas = time.time() - t0
    # CASCI on SA orbitals to get per-root CI vectors for NEVPT2
    ci = mcscf.CASCI(mf, ncas, (1, 0))
    ci.verbose = verbose
    ci.fcisolver.nroots = nroots
    ci.kernel(mo_sa)
    e_cas = np.array(ci.e_tot)
    civecs = [np.asarray(c).reshape(ncas, -1)[:, 0] for c in ci.ci]   # (ncas,) one-electron
    e_pt2 = np.zeros(nroots)
    t1 = time.time()
    if do_nevpt2:
        for i in range(nroots):
            e_pt2[i] = mrpt.NEVPT(ci, root=i).kernel()
    t_pt = time.time() - t1
    _, ncore2, _, _ = ncore, ci.ncore, ncas, None
    C_act = mo_sa[:, ci.ncore:ci.ncore + ncas]
    out_extra = {}
    return dict(e_cas=e_cas, e_pt2=e_cas + e_pt2, mo_sa=mo_sa, civecs=np.array(civecs), C_act=C_act,
                ncore=ci.ncore, ncas=ncas, proj_weights=w, t_casscf=t_cas, t_nevpt2=t_pt,
                converged=bool(mc.converged))


def nevpt2_frozen(mf, mo, ncore, ncas, nroots, frozen_idx, verbose=0):
    """SC-NEVPT2 with the orbitals ``frozen_idx`` (indices into the core block of
    ``mo``) left uncorrelated.  PySCF's NEVPT2 has no frozen-core option, so the frozen
    orbitals are folded into an effective one-electron Hamiltonian (h + 2J - K of the
    frozen density, plus their energy as a constant) and removed from the MO set."""
    from pyscf import mcscf, mrpt, scf
    mol = mf.mol
    fz = np.asarray(frozen_idx, int)
    keep = np.setdiff1d(np.arange(mo.shape[1]), fz)
    Cf = mo[:, fz]
    Df = 2 * Cf @ Cf.T
    h0 = mf.get_hcore()
    vj, vk = scf.hf.get_jk(mol, Df)
    veff = vj - 0.5 * vk
    e_fz = float(np.einsum("ij,ji", Df, h0 + 0.5 * veff))
    heff = h0 + veff
    import copy
    mol2 = copy.copy(mol)
    mol2.nelectron = mol.nelectron - 2 * len(fz)
    mf2 = scf.ROHF(mol2)
    mf2.get_hcore = lambda *a, **k: heff
    e_nuc = mf.energy_nuc()
    mf2.energy_nuc = lambda *a, **k: e_nuc + e_fz
    mf2.mo_coeff = mo[:, keep]; mf2.mo_occ = mf.mo_occ[keep] if len(mf.mo_occ) == mo.shape[1] else None
    mf2.mo_energy = np.zeros(len(keep))
    mf2.converged = True
    ci = mcscf.CASCI(mf2, ncas, (1, 0))
    ci.verbose = verbose
    ci.fcisolver.nroots = nroots
    ci.kernel(mo[:, keep])
    e = np.array(ci.e_tot)
    pt = np.array([mrpt.NEVPT(ci, root=i).kernel() for i in range(nroots)])
    return e, e + pt


def semicore_indices(mol, mo, ncore, n_s=1, n_p=3, thresh=0.5):
    """Indices of the outermost dopant (atom 0) s and p core orbitals (Ce 5s, 5p): core
    orbitals with > thresh Mulliken weight on atom-0 s (p) AOs, the last n_s (n_p) of each
    in orbital order (core orbitals are energy ordered)."""
    S = mol.intor("int1e_ovlp")
    labs = mol.ao_labels(fmt=False)
    C = mo[:, :ncore]
    out = []
    for lch, n in (("s", n_s), ("p", n_p)):
        sel = [i for i, (ia, sy, sh, c) in enumerate(labs) if ia == 0 and sh.endswith(lch)]
        w = np.einsum("ik,ij,jk->k", C[sel], S[sel], C)
        out += [int(i) for i in np.where(w > thresh)[0]][-n:]
    return sorted(out)


def lf_matrix(mol, res, shells, energies):
    """Effective one-electron ligand-field matrix (cm-1, relative to mean) in the basis of
    orthonormalised metal real harmonics for ``shells`` (concatenated).  Also returns the
    character of each root in each shell."""
    from pyscf import lo
    S = mol.intor("int1e_ovlp")
    Rf, blocks = atomic_reference(mol, shells)
    Rf = Rf @ _lowdin(Rf.T @ S @ Rf)
    M = res["C_act"].T @ S @ Rf                            # active x metal
    U, s_, Vt = np.linalg.svd(M)
    T = U @ Vt                                              # Lowdin: closest unitary
    c = res["civecs"]                                       # roots x active
    E = (energies - energies.mean()) * HARTREE_CM1
    h_act = (c.T * E) @ c
    h = T.T @ h_act @ T
    h = 0.5 * (h + h.T)
    chars = []
    for i in range(len(E)):
        v = c[i] @ T
        chars.append([float((v[b[1]:b[1] + 2 * b[0] + 1] ** 2).sum()) for b in blocks])
    return h, blocks, np.array(chars), float(np.min(s_))
