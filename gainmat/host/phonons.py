"""Phonons from a MACE foundation model + phonopy.

DEFAULT (since the Nov 2026 bake-off, BAKEOFF_PHONONS.md): MACE-OMAT-0 medium,
experimental cell (positions relaxed), LO-TO NAC where Born charges are known.
Legacy: MACE-MP-0 "medium" (Batatia et al., arXiv:2401.00096), loaded through
`mace.calculators.mace_mp(model="medium")`.  It is downloaded automatically
from the ACEsuite GitHub release on first use (~20 MB) and is the most
widely-validated MACE foundation model for phonons (e.g. Loew et al.,
"Universal MLIPs are ready for phonons", arXiv:2412.16551, ~10-20 cm^-1 MAE on
average frequencies, larger errors for the top optical modes in oxides).
Newer models ("medium-mpa-0", "medium-omat-0") are typically better for phonons;
pass model="medium-mpa-0" to use one (it downloads the same way).

Outputs
  phonon_max_cm1 : highest frequency on the band path + Gamma-centred DOS mesh,
                   INCLUDING the highest LO mode when NAC is on (out["nac"]).
                   phonon_max_TO_cm1 = max over mesh/path/Gamma(no direction); with NAC
                   the mesh already contains near-Gamma LO branches.
  phonon_eff_cm1 : DOS-based effective cutoff = frequency below which
                   `eff_fraction` (default 0.99) of the integrated phonon DOS lies.
                   This is a robust "cutoff phonon" for the energy-gap law that
                   ignores a few stray modes.
"""
from __future__ import annotations

import time
import numpy as np

THZ_TO_CM1 = 33.35641


def get_device(device: str | None = None) -> str:
    if device:
        return device
    try:
        import torch
        return "cuda" if torch.cuda.is_available() else "cpu"
    except Exception:
        return "cpu"


_CALC_CACHE: dict = {}


def get_mace_calc(model: str = "medium", device: str | None = None,
                  default_dtype: str = "float64"):
    """MACE foundation-model calculator (float64 recommended for phonons)."""
    from mace.calculators import mace_mp
    device = get_device(device)
    key = (model, device, default_dtype)
    if key not in _CALC_CACHE:
        _CALC_CACHE[key] = mace_mp(model=model, device=device,
                                   default_dtype=default_dtype)
    return _CALC_CACHE[key]


# --------------------------------------------------------------------------
# Universal-MLIP registry (see BAKEOFF_PHONONS.md).  Names are bake-off keys;
# legacy MACE names ("medium", "medium-mpa-0", ...) are still accepted.
# --------------------------------------------------------------------------
MLIP_MODELS = {
    "mace-mp-0-medium": ("mace", "medium"),
    "mace-mpa-0-medium": ("mace", "medium-mpa-0"),
    "mace-omat-0-medium": ("mace", "medium-omat-0"),
    "7net-mf-ompa": ("sevenn", ("7net-mf-ompa", "mpa")),
    "7net-0": ("sevenn", ("7net-0", None)),
    "orb-v3-con-omat": ("orb", "orb_v3_conservative_inf_omat"),
    "mattersim-v1-1m": ("mattersim", "MatterSim-v1.0.0-1M.pth"),
    "mattersim-v1-5m": ("mattersim", "MatterSim-v1.0.0-5M.pth"),
}
LEGACY_MACE = {"medium": "mace-mp-0-medium", "medium-mpa-0": "mace-mpa-0-medium",
               "medium-omat-0": "mace-omat-0-medium"}
DEFAULT_MODEL = "mace-omat-0-medium"   # bake-off winner (BAKEOFF_PHONONS.md)
LEGACY_MODEL = "mace-mp-0-medium"   # pre-bake-off default (MACE-MP-0 medium)


def get_calculator(name: str = DEFAULT_MODEL, device: str | None = None):
    """ASE calculator for a bake-off model key (cached)."""
    name = LEGACY_MACE.get(name, name)
    device = get_device(device)
    key = (name, device)
    if key in _CALC_CACHE:
        return _CALC_CACHE[key]
    fam, arg = MLIP_MODELS[name]
    if fam == "mace":
        calc = get_mace_calc(arg, device=device)
    elif fam == "sevenn":
        from sevenn.calculator import SevenNetCalculator
        mdl, modal = arg
        kw = dict(modal=modal) if modal else {}
        calc = SevenNetCalculator(mdl, device=device, **kw)
    elif fam == "orb":
        from orb_models.forcefield import pretrained
        from orb_models.forcefield.calculator import ORBCalculator
        orbff = getattr(pretrained, arg)(device=device, precision="float32-highest")
        calc = ORBCalculator(orbff, device=device)
    elif fam == "mattersim":
        from mattersim.forcefield import MatterSimCalculator
        calc = MatterSimCalculator(load_path=arg, device=device)
    else:
        raise KeyError(name)
    _CALC_CACHE[key] = calc
    return calc


def gamma_irreps(ph, tol: float = 1e-3):
    """Gamma-point modes labelled by irrep: list of dict(freq_cm1, irrep, deg).
    Acoustic modes (|f| < 5 cm-1, lowest 3) are dropped."""
    ph.set_irreps([0, 0, 0], degeneracy_tolerance=tol)
    ir = ph.irreps
    out = []
    for i, ds in enumerate(ir._degenerate_sets):
        lab = ir._ir_labels[i] if ir._ir_labels is not None else None
        f = float(ir._freqs[ds[0]]) * THZ_TO_CM1
        out.append(dict(freq_cm1=round(f, 1), irrep=lab, deg=len(ds)))
    # drop the acoustic triplet
    nac, keep = 0, []
    for m in sorted(out, key=lambda m: abs(m["freq_cm1"])):
        if nac < 3 and abs(m["freq_cm1"]) < 5:
            nac += m["deg"]
            continue
        keep.append(m)
    return sorted(keep, key=lambda m: m["freq_cm1"])


def nac_params_for(atoms, born_ref: dict, eps_inf):
    """phonopy nac_params from per-element Born charges (scalar or 3x3) and
    eps_inf (scalar, (perp, par) pair, or 3x3).  Charges are made neutral
    (acoustic sum rule) by phonopy's symmetrisation.  Units: eV/A forces ->
    phonopy default NAC factor (Hartree*Bohr = 14.40 eV*A)."""
    born = []
    for s in atoms.get_chemical_symbols():
        z = np.asarray(born_ref[s], dtype=float)
        born.append(np.eye(3) * z if z.ndim == 0 else z)
    born = np.array(born)
    born -= born.mean(axis=0)[None]           # enforce charge neutrality
    e = np.asarray(eps_inf, dtype=float)
    if e.ndim == 0:
        eps = np.eye(3) * e
    elif e.shape == (2,):                      # uniaxial (perp, par); c || z
        eps = np.diag([e[0], e[0], e[1]])
    else:
        eps = e
    return dict(born=born, dielectric=eps, factor=14.399645)


def lo_frequencies(ph, directions=((0, 0, 1), (1, 0, 0))) -> dict:
    """Gamma frequencies (cm-1) with NAC for q -> 0 along Cartesian directions.
    Directions are converted to reduced reciprocal coordinates."""
    cell = np.asarray(ph.primitive.cell)
    res = {}
    for d in directions:
        qred = cell @ np.asarray(d, float)    # q_cart . a_i (reduced coords)
        ph.run_qpoints([[0, 0, 0]], nac_q_direction=qred)
        f = np.sort(ph.get_qpoints_dict()["frequencies"][0]) * THZ_TO_CM1
        res["".join(map(str, d))] = np.round(f, 1).tolist()
    return res


def relax(atoms, calc, fmax: float = 1e-3, steps: int = 300, cell: bool = True):
    """Relax positions (and cell, keeping symmetry) with FIRE/BFGS."""
    from ase.optimize import BFGS
    from ase.constraints import FixSymmetry
    try:
        from ase.filters import FrechetCellFilter as CellFilter
    except Exception:  # older ase
        from ase.constraints import ExpCellFilter as CellFilter
    a = atoms.copy()
    a.calc = calc
    try:
        a.set_constraint(FixSymmetry(a))
    except Exception:
        pass
    obj = CellFilter(a) if cell else a
    opt = BFGS(obj, logfile=None)
    opt.run(fmax=fmax, steps=steps)
    a.set_constraint()
    return a, opt.get_number_of_steps()


def _supercell_matrix(atoms, min_len: float):
    """Diagonal supercell with every lattice vector length >= min_len (A)."""
    lens = atoms.cell.lengths()
    return np.diag([max(1, int(np.ceil(min_len / L))) for L in lens])


def compute_phonons(atoms, calc=None, model: str = DEFAULT_MODEL, device=None,
                    supercell_min_len: float = 10.0, displacement: float = 0.01,
                    mesh=(20, 20, 20), relax_cell: bool = False,
                    eff_fraction: float = 0.99, fmax: float = 1e-3,
                    band_npoints: int = 51, nac: bool | dict = True,
                    return_objects: bool = False) -> dict:
    """Finite-displacement phonons (phonopy) with a universal-MLIP calculator.

    Defaults = bake-off winner: MACE-OMAT-0 medium, cell fixed at the input
    (experimental) lattice with positions relaxed (relax_cell=False), LO-TO
    non-analytic correction from gainmat.refs.phonon_refs when Born charges +
    eps_inf are known for atoms.info['structure_key'] (or pass nac=dict(born=
    {el: Z}, eps_inf=...)).  Old behaviour: model="mace-mp-0-medium",
    relax_cell=True, nac=False.  Use relax_cell=True for hypothetical hosts
    without an experimental lattice (expect +2..7 % volume and softer modes).
    """
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    t0 = time.time()
    user_calc = calc is not None
    if calc is None:
        calc = get_calculator(model, device=device)
    a0 = atoms.copy()
    a, nsteps = relax(a0, calc, fmax=fmax, cell=relax_cell)
    label = LEGACY_MACE.get(model, model) if not user_calc else type(calc).__name__
    out = dict(model=label, device=get_device(device), relax_steps=nsteps,
               cell_protocol="relaxed" if relax_cell else "expcell",
               vol_change_pct=100 * (a.get_volume() / a0.get_volume() - 1))

    unit = PhonopyAtoms(symbols=a.get_chemical_symbols(), cell=a.cell[:],
                        scaled_positions=a.get_scaled_positions())
    sc = _supercell_matrix(a, supercell_min_len)
    ph = Phonopy(unit, supercell_matrix=sc, primitive_matrix=None)
    ph.generate_displacements(distance=displacement)
    forces = []
    from ase import Atoms
    for s in ph.supercells_with_displacements:
        sa = Atoms(s.symbols, cell=s.cell, scaled_positions=s.scaled_positions, pbc=True)
        sa.calc = calc
        f = sa.get_forces().copy()
        forces.append(f - f.mean(axis=0))
    ph.forces = np.array(forces)
    ph.produce_force_constants()
    ph.symmetrize_force_constants()
    # Gamma TO modes with irreps (before NAC)
    try:
        out["gamma_modes"] = gamma_irreps(ph)
    except Exception as e:
        out["irreps_error"] = repr(e)
    # LO-TO non-analytic correction
    nacp = None
    if isinstance(nac, dict):
        nacp = nac_params_for(a, nac["born"], nac["eps_inf"])
        out["nac_source"] = "user"
    elif nac:
        from gainmat.refs.phonon_refs import PHONON_REFS
        ref = PHONON_REFS.get(atoms.info.get("structure_key", ""), {})
        if ref.get("born") and ref.get("eps_inf"):
            nacp = nac_params_for(a, ref["born"], ref["eps_inf"])
            out["nac_source"] = "gainmat.refs.phonon_refs"
    out["nac"] = nacp is not None
    if nacp is not None:
        ph.nac_params = nacp
        out["gamma_nac_cm1"] = lo_frequencies(ph)
    out["supercell"] = np.diag(sc).tolist()
    out["n_displacements"] = len(forces)

    # DOS on a mesh
    ph.run_mesh(list(mesh), with_eigenvectors=False, is_gamma_center=True)
    freqs_mesh = ph.get_mesh_dict()["frequencies"] * THZ_TO_CM1
    weights = ph.get_mesh_dict()["weights"]
    ph.run_total_dos(sigma=None, freq_pitch=0.05)
    dos = ph.get_total_dos_dict()
    f_dos = np.asarray(dos["frequency_points"]) * THZ_TO_CM1
    g_dos = np.asarray(dos["total_dos"])

    # band path (Gamma always included)
    fmax_path = None
    try:
        from phonopy.phonon.band_structure import get_band_qpoints_by_seekpath  # noqa
        ph.auto_band_structure(npoints=band_npoints, with_eigenvectors=False)
        bs = ph.get_band_structure_dict()
        fmax_path = float(max(np.max(f) for f in bs["frequencies"]) * THZ_TO_CM1)
        fmin_path = float(min(np.min(f) for f in bs["frequencies"]) * THZ_TO_CM1)
        out["band_min_cm1"] = fmin_path
    except Exception as e:  # seekpath not installed etc.
        out["band_error"] = repr(e)
    # Gamma
    ph.run_qpoints([[0, 0, 0]])
    fg = ph.get_qpoints_dict()["frequencies"][0] * THZ_TO_CM1
    out["gamma_cm1"] = np.round(np.sort(fg), 1).tolist()

    fm = float(freqs_mesh.max())
    lo = max((max(v) for v in out.get("gamma_nac_cm1", {}).values()), default=0.0)
    out["phonon_max_TO_cm1"] = max(fm, fmax_path or 0.0, float(fg.max()))
    # includes LO when NAC is on (LO is the true top of the spectrum)
    out["phonon_max_cm1"] = max(out["phonon_max_TO_cm1"], lo)
    out["phonon_eff_cm1"] = dos_cutoff(freqs_mesh, weights, eff_fraction)
    out["min_freq_cm1"] = float(freqs_mesh.min())   # < -10 => instability
    out["dos"] = dict(freq_cm1=f_dos.tolist(), dos=g_dos.tolist())
    out["t_s"] = round(time.time() - t0, 1)
    if return_objects:
        out["_phonopy"] = ph
        out["_atoms"] = a
    return out


def dos_cutoff(freqs, weights=None, fraction=0.99) -> float:
    """Frequency below which `fraction` of modes lie (mesh-weighted)."""
    f = np.asarray(freqs)
    w = np.ones_like(f) if weights is None else \
        np.repeat(np.asarray(weights)[:, None], f.shape[1], axis=1)
    f = f.ravel(); w = w.ravel()
    o = np.argsort(f)
    c = np.cumsum(w[o].astype(float)); c /= c[-1]
    return float(f[o][np.searchsorted(c, fraction)])


def compute_phonons_gpaw(atoms, ecut: float = 500, kpts=(4, 4, 4),
                         supercell=(2, 2, 2), displacement: float = 0.01,
                         mesh=(16, 16, 16)) -> dict:
    """Optional DFT cross-check (PBE, GPAW PW, phonopy finite displacements)
    on the GIVEN geometry (no relaxation).  Intended for CaF2 (3-atom cell,
    2x2x2 supercell = 24 atoms, 2 displacements).  Expensive: run on Colab only.
    k-points refer to the supercell."""
    from gpaw import GPAW, PW, FermiDirac
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    from ase import Atoms
    unit = PhonopyAtoms(symbols=atoms.get_chemical_symbols(), cell=atoms.cell[:],
                        scaled_positions=atoms.get_scaled_positions())
    ph = Phonopy(unit, supercell_matrix=np.diag(supercell))
    ph.generate_displacements(distance=displacement)
    forces = []
    for s in ph.supercells_with_displacements:
        sa = Atoms(s.symbols, cell=s.cell, scaled_positions=s.scaled_positions, pbc=True)
        sa.calc = GPAW(mode=PW(ecut), xc="PBE", kpts={"size": kpts, "gamma": True},
                       occupations=FermiDirac(0.01), symmetry="off", txt=None)
        f = sa.get_forces()
        forces.append(f - f.mean(axis=0))
    ph.forces = np.array(forces)
    ph.produce_force_constants()
    ph.symmetrize_force_constants()
    ph.run_mesh(list(mesh))
    fr = ph.get_mesh_dict()["frequencies"] * THZ_TO_CM1
    ph.run_qpoints([[0, 0, 0]])
    fg = ph.get_qpoints_dict()["frequencies"][0] * THZ_TO_CM1
    return dict(method="GPAW-PBE finite displacement", ecut=ecut,
                phonon_max_cm1=float(max(fr.max(), fg.max())),
                phonon_eff_cm1=dos_cutoff(fr, ph.get_mesh_dict()["weights"]),
                gamma_cm1=np.round(np.sort(fg), 1).tolist())
