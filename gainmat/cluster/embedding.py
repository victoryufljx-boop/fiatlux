"""Embedded-cluster construction for single-open-shell dopant centres.

Cluster layout (all coordinates in Angstrom, dopant at the origin):

  QM      : dopant + first anion shell (full basis, all valence electrons)
  caps    : next cation shell within ``r_cap`` -- each cation is an ECP-only
            pseudo-atom (ECP + one very tight s function that stays empty, i.e.
            no valence electrons) whose ECP core charge equals the formal ionic
            charge (Li+ with a 2-e ECP, Al3+ with a 10-e ECP).  If no such ECP
            exists for an element (Y3+ here), the cation is a bare point charge
            and this is reported in ``Cluster.notes`` (electron-leakage risk).
            Optionally ``partial_cap`` places a smaller-core ECP plus a
            co-located correcting point charge (Zeff - formal) instead.
  PC      : all other ions inside a sphere of radius ``r_pc`` as formal point
            charges; charges in the outer shell ``r_pc - shell <= r < r_pc``
            are adjusted (least squares, with an exact neutrality constraint)
            so that the potential of the finite array reproduces the Ewald
            potential of the infinite crystal at sample points inside the
            QM region.

Point charges enter the SCF via ``pyscf.qmmm.mm_charge``.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np
from scipy.special import erfc

from ..refs.cluster_refs import STRUCTURES

EV_ANG = 14.399645  # e^2/(4 pi eps0) in eV*Angstrom  -> potential in V = EV_ANG * q / r


# ---------------------------------------------------------------- crystals
def build_crystal(key: str):
    """Return (cell[3x3] Angstrom, symbols, frac_coords[N,3], charges[N]) for one unit cell.

    Built from the space group + Wyckoff positions in refs.cluster_refs.STRUCTURES
    with ase.spacegroup (only for symmetry expansion)."""
    from ase.spacegroup import crystal
    s = STRUCTURES[key]
    syms = [x[0] for x in s["sites"]]
    basis = [x[1] for x in s["sites"]]
    qmap = {x[0]: x[2] for x in s["sites"]}
    if s["spacegroup"] == 167:
        cellpar = [s["a"], s["a"], s["c"], 90, 90, 120]
    else:
        cellpar = [s["a"], s["a"], s["c"], 90, 90, 90]
    at = crystal(syms, basis, spacegroup=s["spacegroup"], setting=s["setting"],
                 cellpar=cellpar)
    sy = at.get_chemical_symbols()
    q = np.array([qmap[x] for x in sy], float)
    assert abs(q.sum()) < 1e-8, "unit cell not neutral"
    return np.array(at.cell), sy, at.get_scaled_positions(wrap=True), q


# ---------------------------------------------------------------- Ewald
def ewald_potential(points, cell, frac, q, alpha=None, rcut=None, gtol=1e-10):
    """Electrostatic potential (V) of the infinite periodic point-charge crystal
    at Cartesian ``points`` (Angstrom).  Tin-foil convention (mean potential 0).
    Points must not coincide with an ion (use ``site_potential`` for that)."""
    points = np.atleast_2d(points)
    cart = frac @ cell
    L = np.linalg.norm(cell, axis=1)
    if rcut is None:
        rcut = max(12.0, 1.2 * L.max())
    if alpha is None:
        alpha = 3.5 / rcut          # erfc(3.5) ~ 7e-7
    # real space
    nmax = [int(np.ceil(rcut / np.abs(np.linalg.det(cell)) *
                        np.linalg.norm(np.cross(cell[(i + 1) % 3], cell[(i + 2) % 3])))) + 1
            for i in range(3)]
    T = np.array(list(itertools.product(*[range(-n, n + 1) for n in nmax]))) @ cell
    V = np.zeros(len(points))
    for t in T:
        d = points[:, None, :] - (cart[None, :, :] + t)
        r = np.linalg.norm(d, axis=2)
        m = r < rcut
        V += np.where(m, q[None, :] * erfc(alpha * r) / np.where(m, r, 1.0), 0).sum(1)
    # reciprocal space
    rec = 2 * np.pi * np.linalg.inv(cell).T
    gmax = 2 * alpha * np.sqrt(-np.log(gtol))
    nG = [int(np.ceil(gmax / np.linalg.norm(rec[i]))) + 1 for i in range(3)]
    G = np.array(list(itertools.product(*[range(-n, n + 1) for n in nG]))) @ rec
    g2 = (G ** 2).sum(1)
    G, g2 = G[(g2 > 1e-12) & (g2 < gmax ** 2)], g2[(g2 > 1e-12) & (g2 < gmax ** 2)]
    vol = abs(np.linalg.det(cell))
    sf = np.exp(-1j * (G @ cart.T)) @ q                       # S(G)
    pref = 4 * np.pi / vol * np.exp(-g2 / (4 * alpha ** 2)) / g2
    V += np.real(np.exp(1j * points @ G.T) @ (pref * sf))
    return EV_ANG * V


def site_potential(cell, frac, q, i, **kw):
    """Madelung potential (V) at ion i due to all other ions (Ewald)."""
    cart = frac @ cell
    eps = 1e-4
    # potential at a tiny offset minus the ion's own bare Coulomb term
    offs = np.array([[eps, 0, 0], [-eps, 0, 0], [0, eps, 0], [0, -eps, 0], [0, 0, eps], [0, 0, -eps]])
    v = ewald_potential(cart[i] + offs, cell, frac, q, **kw).mean()
    return v - EV_ANG * q[i] / eps


# ---------------------------------------------------------------- cluster
# ECPs whose core charge equals the formal ionic charge (no valence electrons left)
FULL_CAP_ECP = {"Li": "stuttgart",   # 2-e core  -> Zeff = 1
                "Al": "lanl2dz",     # 10-e core -> Zeff = 3
                "Ca": "stuttgart"}   # 18-e core -> Zeff = 2
PARTIAL_CAP_ECP = {"Y": "sbkjc"}     # 28-e core -> Zeff = 11, needs -8 co-located PC
TIGHT_S = [[0, [1.0e4, 1.0]]]        # dummy basis for ECP-only atoms (never occupied)


@dataclass
class Cluster:
    name: str
    dopant: str
    dopant_charge: int
    qm_symbols: list
    qm_coords: np.ndarray            # Angstrom, dopant first
    qm_formal: np.ndarray
    cap_symbols: list
    cap_coords: np.ndarray
    cap_formal: np.ndarray
    cap_ecp: list                    # ECP name or None (=> treated as point charge)
    pc_coords: np.ndarray
    pc_charges: np.ndarray           # point charges actually passed to PySCF (incl. cap corrections)
    madelung_ewald_V: float = 0.0    # host Madelung potential at the substituted site
    madelung_cluster_V: float = 0.0  # same from the finite cluster (all ions as formal/fitted charges)
    fit_rms_V: float = 0.0
    notes: list = field(default_factory=list)
    aimp_symbols: list = field(default_factory=list)     # frozen-ion (AIMP-like) caps, no basis
    aimp_coords: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))
    aimp_formal: np.ndarray = field(default_factory=lambda: np.zeros(0))

    @property
    def qm_charge(self):
        return int(round(self.qm_formal.sum()))

    def total_charge(self):
        return (self.qm_formal.sum() + self.cap_formal.sum() + self._pc_formal_sum
                + float(np.sum(self.aimp_formal)))

    def summary(self):
        d = np.linalg.norm(self.qm_coords[1:], axis=1)
        return dict(name=self.name, n_qm=len(self.qm_symbols), n_cap=len(self.cap_symbols),
                    n_pc=len(self.pc_charges), n_aimp=len(self.aimp_symbols), qm_charge=self.qm_charge,
                    dopant_ligand_A=sorted(np.round(d, 4).tolist()),
                    total_charge=float(self.total_charge()),
                    madelung_ewald_V=self.madelung_ewald_V,
                    madelung_cluster_V=self.madelung_cluster_V, fit_rms_V=self.fit_rms_V,
                    notes=self.notes)


def make_cluster(key, site, dopant, dopant_charge, r_qm=2.6, r_cap=4.0, r_pc=12.0,
                 shell=2.5, n_fit=400, partial_cap=False, seed=0, cap_mode=None):
    """cap_mode for cations without a full-core ECP: 'pc' (bare point charge, default),
    'partial' (small-core ECP + co-located correcting charge, = partial_cap=True) or
    'aimp' (frozen-ion model potential, see aimp_hcore).

    Build an embedded cluster for ``dopant`` substituting a ``site`` cation of host ``key``."""
    cell, sy, frac, q = build_crystal(key)
    i0 = [i for i, s in enumerate(sy) if s == site][0]
    v_site = site_potential(cell, frac, q, i0)
    # supercell around the site
    cart0 = frac @ cell
    origin = cart0[i0]
    L = np.linalg.norm(cell, axis=1)
    n = [int(np.ceil(r_pc / l)) + 1 for l in L]
    pos, chg, sym = [], [], []
    for t in itertools.product(*[range(-k, k + 1) for k in n]):
        tt = np.array(t) @ cell
        pos.append(cart0 + tt - origin)
        chg.append(q)
        sym += sy
    pos = np.vstack(pos); chg = np.concatenate(chg); sym = np.array(sym)
    r = np.linalg.norm(pos, axis=1)
    keep = r < r_pc
    pos, chg, sym, r = pos[keep], chg[keep], sym[keep], r[keep]
    order = np.argsort(r, kind="stable")
    pos, chg, sym, r = pos[order], chg[order], sym[order], r[order]
    assert r[0] < 1e-6
    anion = chg < 0
    is_qm = (r < 1e-6) | (anion & (r < r_qm))
    is_cap = (~is_qm) & (chg > 0) & (r < r_cap)
    is_bnd = (~is_qm) & (~is_cap) & (r >= r_pc - shell)

    # ---- fit boundary charges to the Ewald potential inside the QM region
    rng = np.random.default_rng(seed)
    rfit = r[is_qm].max() + 0.5
    pts = rng.normal(size=(n_fit, 3))
    pts *= (rfit * rng.random(n_fit) ** (1 / 3) / np.linalg.norm(pts, axis=1))[:, None]
    pts = np.vstack([pts, [[1e-3, 2e-3, -1.5e-3]]])        # include ~dopant site
    target = ewald_potential(pts + origin, cell, frac, q)
    def coul(P, R):
        return EV_ANG / np.linalg.norm(P[:, None, :] - R[None, :, :], axis=2)
    fixed = ~is_bnd
    b = target - coul(pts, pos[fixed]) @ chg[fixed]
    A = coul(pts, pos[is_bnd])
    q0 = chg[is_bnd]
    nb = len(q0)
    mu = 1e-3 * np.trace(A.T @ A) / nb
    # minimise |A(q0+dq) - b|^2 + mu |dq|^2  s.t.  sum(dq) = -sum(all formal)
    K = np.zeros((nb + 1, nb + 1))
    K[:nb, :nb] = A.T @ A + mu * np.eye(nb)
    K[:nb, nb] = K[nb, :nb] = 1.0
    rhs = np.concatenate([A.T @ (b - A @ q0), [-chg.sum()]])
    dq = np.linalg.solve(K, rhs)[:nb]
    chg = chg.copy()
    chg[is_bnd] = q0 + dq
    resid = A @ chg[is_bnd] - b
    fit_rms = float(np.sqrt(np.mean(resid ** 2)))
    # host Madelung potential at the site from the finite array
    v_clu = float(EV_ANG * (chg[1:] / r[1:]).sum())

    notes = []
    if cap_mode is None:
        cap_mode = "partial" if partial_cap else "pc"
    partial_cap = cap_mode == "partial"
    cap_ecp, pc_extra_pos, pc_extra_q = [], [], []
    for s in sym[is_cap]:
        if s in FULL_CAP_ECP:
            cap_ecp.append(FULL_CAP_ECP[s])
        elif cap_mode == "aimp" and s in AIMP_IONS:
            cap_ecp.append("aimp")
        elif partial_cap and s in PARTIAL_CAP_ECP:
            cap_ecp.append(PARTIAL_CAP_ECP[s])
        else:
            cap_ecp.append(None)
    missing = sorted({s for s, e in zip(sym[is_cap], cap_ecp) if e is None})
    if missing:
        notes.append(f"cations {missing} within r_cap={r_cap} A are bare point charges "
                     "(no full-core ECP available): electron-leakage risk")
    qm_sym = list(sym[is_qm]); qm_sym[0] = dopant
    qm_formal = chg[is_qm].copy(); qm_formal[0] = dopant_charge
    capmask_aimp = np.array([e == "aimp" for e in cap_ecp], bool)
    if capmask_aimp.any():
        notes.append(f"{int(capmask_aimp.sum())} cations treated as frozen-ion AIMP caps "
                     f"({sorted(set(sym[is_cap][capmask_aimp]))})")
    capmask_pc = np.array([e is None for e in cap_ecp], bool)
    pc_pos = np.vstack([pos[(~is_qm) & (~is_cap)], pos[is_cap][capmask_pc]])
    pc_q = np.concatenate([chg[(~is_qm) & (~is_cap)], chg[is_cap][capmask_pc]])
    cl = Cluster(name=f"{dopant}{dopant_charge}+:{key}", dopant=dopant,
                 dopant_charge=dopant_charge, qm_symbols=qm_sym, qm_coords=pos[is_qm],
                 qm_formal=qm_formal,
                 cap_symbols=[s for s, e in zip(sym[is_cap], cap_ecp) if e not in (None, "aimp")],
                 cap_coords=pos[is_cap][~capmask_pc & ~capmask_aimp],
                 cap_formal=chg[is_cap][~capmask_pc & ~capmask_aimp],
                 cap_ecp=[e for e in cap_ecp if e not in (None, "aimp")],
                 aimp_symbols=list(sym[is_cap][capmask_aimp]),
                 aimp_coords=pos[is_cap][capmask_aimp], aimp_formal=chg[is_cap][capmask_aimp],
                 pc_coords=pc_pos, pc_charges=pc_q,
                 madelung_ewald_V=float(v_site), madelung_cluster_V=v_clu,
                 fit_rms_V=fit_rms, notes=notes)
    cl._pc_formal_sum = float(pc_q.sum())
    return cl


# ---------------------------------------------------------------- PySCF objects
def _zeff(symbol, ecp):
    from pyscf import gto
    from pyscf.data.elements import charge
    e = gto.basis.load_ecp(ecp, symbol)
    return charge(symbol) - (e[0] if e else 0)


def build_mol(cl: Cluster, basis: dict, ecp: dict | None = None, spin=1, verbose=0,
              max_memory=4000):
    """Return (mol, mm_coords_A, mm_charges).  QM atoms first, then ECP caps
    (labelled "<El>@9" so they get their own basis/ECP)."""
    from pyscf import gto
    from pyscf.data.elements import charge as Z
    atoms, bas, ecps = [], dict(basis), dict(ecp or {})
    for s, x in zip(cl.qm_symbols, cl.qm_coords):
        atoms.append([s, tuple(x)])
    extra_pos, extra_q = [], []
    for s, x, e, qf in zip(cl.cap_symbols, cl.cap_coords, cl.cap_ecp, cl.cap_formal):
        lab = f"{s}@9"
        atoms.append([lab, tuple(x)])
        bas[lab] = TIGHT_S
        ecps[lab] = gto.basis.load_ecp(e, s)
        zf = _zeff(s, e)
        if abs(zf - qf) > 1e-8:          # partial cap: co-located correction charge
            extra_pos.append(x); extra_q.append(qf - zf)
    # electrons: only QM valence electrons
    n_ecp_core_qm = 0
    for s in set(cl.qm_symbols):
        if s in ecps and ecps[s]:
            e = ecps[s] if not isinstance(ecps[s], str) else gto.basis.load_ecp(ecps[s], s)
            if e:
                n_ecp_core_qm += e[0] * cl.qm_symbols.count(s)
    n_el = sum(Z(s) for s in cl.qm_symbols) - n_ecp_core_qm - int(round(cl.qm_formal.sum()))
    mol = gto.Mole()
    mol.atom = atoms
    mol.basis = bas
    mol.ecp = ecps
    mol.unit = "Angstrom"
    mol.verbose = verbose
    mol.max_memory = max_memory
    # total charge so that nelectron == n_el
    mol.charge = 0
    mol.spin = None                     # parity-adaptive for the provisional build
    mol.build(parse_arg=False)
    mol.charge = int(round(mol.atom_charges().sum())) - n_el
    mol.spin = spin
    mol.build(parse_arg=False)
    assert mol.nelectron == n_el
    mm_pos = cl.pc_coords if not extra_pos else np.vstack([cl.pc_coords, extra_pos])
    mm_q = cl.pc_charges if not extra_q else np.concatenate([cl.pc_charges, extra_q])
    return mol, np.asarray(mm_pos), np.asarray(mm_q)


def qm_only_guess(mol, n_qm, spin):
    """Initial ROHF density from the QM atoms alone (MINAO), padded with zeros for caps."""
    from pyscf import gto, scf
    m = gto.Mole()
    m.atom = mol._atom[:n_qm]
    m.unit = "Bohr"
    m.basis = {a[0]: mol._basis[a[0]] for a in mol._atom[:n_qm]}
    m.ecp = {k: v for k, v in mol._ecp.items() if k in m.basis}
    m.charge = 0
    m.spin = None
    m.verbose = 0
    m.build(parse_arg=False)
    m.spin = spin
    m.charge = int(round(m.atom_charges().sum())) - mol.nelectron
    m.build(parse_arg=False)
    dm = scf.ROHF(m).get_init_guess(key="minao")
    nq = m.nao
    out = np.zeros((2, mol.nao, mol.nao))
    out[:, :nq, :nq] = dm
    return out


# ---------------------------------------------------------------- AIMP-like frozen-ion caps
# Seijo & Barandiaran-style ab initio model potential for cations of the next shell
# that have no full-core ECP (Y3+ in LiYF4).  Each such ion is a frozen Y3+ (free-ion
# RHF density, 28-e Stuttgart/def2 ECP + def2-SVP) that contributes to the cluster
# one-electron Hamiltonian, *without* any basis functions on the ion:
#   V = -Zeff/|r-R| + V_ECP + 2 J[rho_ion] - K[rho_ion] + sum_c (-2 eps_c) |c><c|
# Coulomb and exchange are evaluated exactly (4-centre integrals between cluster AO
# pairs and ion AO pairs); the last term is the Huzinaga/AIMP core-projection
# (level-shift) operator that keeps cluster orbitals out of the ion's occupied 4s4p
# shells, i.e. it supplies the Pauli repulsion that a bare point charge lacks.
AIMP_IONS = {"Y": dict(basis="def2-svp", ecp="def2-svp", charge=3)}
_AIMP_CACHE = {}


def _free_ion(symbol):
    if symbol in _AIMP_CACHE:
        return _AIMP_CACHE[symbol]
    from pyscf import gto, scf
    p = AIMP_IONS[symbol]
    m = gto.M(atom=[[symbol, (0, 0, 0)]], basis={symbol: p["basis"]},
              ecp={symbol: p["ecp"]}, charge=p["charge"], spin=0, verbose=0)
    mf = scf.RHF(m)
    mf.conv_tol = 1e-10
    mf.kernel()
    occ = mf.mo_occ > 0
    out = dict(mol=m, dm=mf.make_rdm1(), C=mf.mo_coeff[:, occ], eps=mf.mo_energy[occ],
               zeff=float(m.atom_charges()[0]), e=mf.e_tot)
    _AIMP_CACHE[symbol] = out
    return out


def aimp_hcore(mol, symbols, coords_A, shift_scale=1.0):
    """One-electron embedding matrix (Hartree, mol AO basis) for frozen ions at
    ``coords_A`` (Angstrom).  Returns (h_emb, info)."""
    from pyscf import gto
    from pyscf.scf import jk as _jk
    BOHR = 0.52917721092
    h = np.zeros((mol.nao, mol.nao))
    has_ecp = bool(mol._ecpbas.size) if hasattr(mol, "_ecpbas") else False
    ecp0 = mol.intor("ECPscalar") if has_ecp else 0.0
    for s, x in zip(symbols, coords_A):
        ion = _free_ion(s)
        m0 = ion["mol"]
        my = gto.M(atom=[[s, tuple(x)]], basis={s: AIMP_IONS[s]["basis"]},
                   ecp={s: AIMP_IONS[s]["ecp"]}, charge=AIMP_IONS[s]["charge"], spin=0,
                   unit="Angstrom", verbose=0)
        # nuclear (ECP core-screened) attraction
        with mol.with_rinv_origin(np.asarray(x) / BOHR):
            h -= ion["zeff"] * mol.intor("int1e_rinv")
        # semilocal ECP of the ion acting on cluster AOs
        cm = gto.conc_mol(mol, my)
        h += cm.intor("ECPscalar")[:mol.nao, :mol.nao] - ecp0
        # Coulomb and exchange of the frozen ion density
        D = ion["dm"]
        vj = _jk.get_jk((mol, mol, my, my), D, scripts="ijkl,lk->ij", intor="int2e")
        vk = _jk.get_jk((mol, my, my, mol), D, scripts="ijkl,jk->il", intor="int2e")
        h += vj - 0.5 * vk
        # core projection (level shift -2 eps_c)
        Sx = gto.intor_cross("int1e_ovlp", mol, my)
        SC = Sx @ ion["C"]
        h += shift_scale * (SC * (-2.0 * ion["eps"])) @ SC.T
    return 0.5 * (h + h.T), dict(n_ions=len(symbols))


def attach_aimp(mf, mol, cl):
    """Add the AIMP embedding of ``cl.aimp_symbols`` to ``mf`` (hcore + nuclear energy)."""
    if not len(getattr(cl, "aimp_symbols", [])):
        return mf
    BOHR = 0.52917721092
    h_emb, _ = aimp_hcore(mol, cl.aimp_symbols, cl.aimp_coords)
    h0 = mf.get_hcore
    e0 = mf.energy_nuc
    # nuclear repulsion of QM nuclei/caps with the ions (net charge q_ion at long range)
    Z = mol.atom_charges()
    R = mol.atom_coords()
    enuc = 0.0
    for x, qf in zip(cl.aimp_coords, cl.aimp_formal):
        d = np.linalg.norm(R - np.asarray(x) / BOHR, axis=1)
        enuc += float((Z * qf / d).sum())
    mf.get_hcore = lambda *a, **k: h0(*a, **k) + h_emb
    mf.energy_nuc = lambda *a, **k: e0(*a, **k) + enuc
    mf._aimp_h = h_emb
    return mf


# ---------------------------------------------------------------- bundled basis files
def load_bundled(fn, symbol):
    """(basis, ecp) for ``symbol`` from an NWChem-format file in cluster/basis_data
    (downloaded once from basis_set_exchange so Colab needs no network)."""
    import os
    from pyscf import gto
    txt = open(os.path.join(os.path.dirname(__file__), "basis_data", fn)).read()
    b_part = txt.split("BASIS", 1)[1].split("END", 1)[0]
    b_part = "\n".join(b_part.splitlines()[1:])
    bas = gto.basis.parse(b_part, symb=symbol)
    ecp = None
    if "\nECP" in txt:
        e_part = txt.split("\nECP", 1)[1].split("END", 1)[0]
        ecp = gto.basis.parse_ecp(e_part, symb=symbol)
    return bas, ecp


def resolve_basis(basis: dict, ecp: dict):
    """Expand 'bundled:<file>' entries into explicit basis/ECP data."""
    basis, ecp = dict(basis), dict(ecp or {})
    for el, b in list(basis.items()):
        if isinstance(b, str) and b.startswith("bundled:"):
            bb, ee = load_bundled(b.split(":", 1)[1], el)
            basis[el] = bb
            if ee is not None:
                ecp[el] = ee
    return basis, ecp
