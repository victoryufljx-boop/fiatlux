"""Phonons from a MACE foundation model + phonopy.

Model: MACE-MP-0 "medium" (Batatia et al., arXiv:2401.00096), loaded through
`mace.calculators.mace_mp(model="medium")`.  It is downloaded automatically
from the ACEsuite GitHub release on first use (~20 MB) and is the most
widely-validated MACE foundation model for phonons (e.g. Loew et al.,
"Universal MLIPs are ready for phonons", arXiv:2412.16551, ~10-20 cm^-1 MAE on
average frequencies, larger errors for the top optical modes in oxides).
Newer models ("medium-mpa-0", "medium-omat-0") are typically better for phonons;
pass model="medium-mpa-0" to use one (it downloads the same way).

Outputs
  phonon_max_cm1 : highest frequency on the band path + DOS mesh (the top
                   optical mode at the TO frequency: phonopy without NAC gives TO frequencies).
                   NOTE: without Born charges (non-analytic correction) the
                   LO-TO splitting at Gamma is missing, so this is a TO-like max.
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


def compute_phonons(atoms, calc=None, model: str = "medium", device=None,
                    supercell_min_len: float = 10.0, displacement: float = 0.01,
                    mesh=(20, 20, 20), relax_cell: bool = True,
                    eff_fraction: float = 0.99, fmax: float = 1e-3,
                    band_npoints: int = 51, return_objects: bool = False) -> dict:
    """Relax + finite-displacement phonons (phonopy) with a MACE calculator."""
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    t0 = time.time()
    user_calc = calc is not None
    if calc is None:
        calc = get_mace_calc(model=model, device=device)
    a0 = atoms.copy()
    if relax_cell:
        a, nsteps = relax(a0, calc, fmax=fmax)
    else:
        a, nsteps = a0, 0
    label = f"mace_mp:{model}" if not user_calc else type(calc).__name__
    out = dict(model=label, device=get_device(device),
               relax_steps=nsteps,
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
        forces.append(sa.get_forces().copy())
    ph.forces = np.array(forces)
    ph.produce_force_constants()
    ph.symmetrize_force_constants()
    out["supercell"] = np.diag(sc).tolist()
    out["n_displacements"] = len(forces)

    # DOS on a mesh
    ph.run_mesh(list(mesh), with_eigenvectors=False)
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
    out["phonon_max_cm1"] = max(fm, fmax_path or 0.0, float(fg.max()))
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
