"""Phonon method bake-off driver (see BAKEOFF_PHONONS.md).

    python -m gainmat.host.bakeoff MODEL [HOST ...] --out results/bakeoff

For each host and each cell protocol ("relaxed" = full symmetry-constrained
cell + position relaxation; "expcell" = experimental lattice, positions relaxed)
computes Gamma-point phonons with phonopy, labels them by irrep, and writes one
JSON per (model, host, protocol).  Gamma frequencies from finite displacements
are exact for any supercell (the periodic-image sum *is* the q=0 dynamical
matrix), so the primitive cell is used unless --sc is given; that keeps the
bake-off cheap on CPU.
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

from gainmat.host.phonons import (relax, gamma_irreps, get_calculator,
                                  nac_params_for, lo_frequencies)
from gainmat.refs.phonon_refs import PHONON_REFS
from gainmat.structures import get_structure

HOSTS = ["CaF2", "Al2O3", "SiO2-alpha-quartz", "LiYF4"]


def gamma_run(atoms, calc, relax_cell, fmax=1e-3, sc_min_len=0.0, host=None):
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    from ase import Atoms
    t0 = time.time()
    a, nsteps = relax(atoms, calc, fmax=fmax, cell=relax_cell, steps=500)
    vol_pct = 100 * (a.get_volume() / atoms.get_volume() - 1)
    unit = PhonopyAtoms(symbols=a.get_chemical_symbols(), cell=a.cell[:],
                        scaled_positions=a.get_scaled_positions())
    lens = a.cell.lengths()
    sc = np.diag([max(1, int(np.ceil(sc_min_len / L))) for L in lens])
    ph = Phonopy(unit, supercell_matrix=sc, primitive_matrix=None)
    ph.generate_displacements(distance=0.01)
    forces = []
    for s in ph.supercells_with_displacements:
        sa = Atoms(s.symbols, cell=s.cell, scaled_positions=s.scaled_positions,
                   pbc=True)
        sa.calc = calc
        f = sa.get_forces().copy()
        forces.append(f - f.mean(axis=0))
    ph.forces = np.array(forces)
    ph.produce_force_constants()
    ph.symmetrize_force_constants()
    modes = gamma_irreps(ph)
    out = dict(vol_pct=vol_pct, relax_steps=nsteps,
               cellpar=[round(x, 4) for x in a.cell.cellpar()],
               n_disp=len(forces), modes=modes)
    ref = PHONON_REFS.get(host, {})
    if ref.get("born") and ref.get("eps_inf"):
        ph.nac_params = nac_params_for(a, ref["born"], ref["eps_inf"])
        out["nac_gamma"] = lo_frequencies(ph)
    out["t_s"] = round(time.time() - t0, 1)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("model")
    p.add_argument("hosts", nargs="*", default=HOSTS)
    p.add_argument("--out", default="results/bakeoff")
    p.add_argument("--protocols", default="relaxed,expcell")
    p.add_argument("--device", default="cpu")
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    calc = get_calculator(a.model, device=a.device)
    t_load = time.time() - t0
    for host in a.hosts:
        for prot in a.protocols.split(","):
            fn = os.path.join(a.out, f"{a.model}__{host}__{prot}.json")
            if os.path.exists(fn):
                continue
            try:
                r = gamma_run(get_structure(host), calc,
                              relax_cell=(prot == "relaxed"), host=host)
            except Exception as e:  # keep going
                r = dict(error=repr(e))
            r.update(model=a.model, host=host, protocol=prot, t_load=t_load)
            with open(fn, "w") as f:
                json.dump(r, f, indent=1)
            print(a.model, host, prot, r.get("vol_pct"), r.get("t_s"),
                  r.get("error", ""), flush=True)


if __name__ == "__main__":
    main()
