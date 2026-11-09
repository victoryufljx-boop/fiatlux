"""Er3+ site geometry in fiber glasses from fixed-density melt-quench MD with a
MACE foundation model.  Wall-clock-budgeted: step counts are set from a
measured step time so a run always finishes inside `--budget-min`.

usage: python md_sites/run_md_sites.py <glass> --budget-min 20 [--threads 1]
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

# composition as {formula unit: count}, experimental density (g/cm3), T_melt
GLASSES = {
    # 64 SiO2 + 2 Er2O3 (3 mol% Er2O3 -> ~6 at% Er per cation; heavy but needed
    # for statistics) -- charge balanced by stoichiometry, no Al co-doping.
    "Er-silica": dict(units={"SiO2": 64, "Er2O3": 2}, density=2.30,
                      T_melt=4000.0, ligands=("O",)),
    # 53ZrF4-20BaF2-4LaF3-3AlF3-20NaF scaled to ~200 atoms; both La -> Er
    "Er-ZBLAN": dict(units={"ZrF4": 27, "BaF2": 10, "ErF3": 2, "AlF3": 2,
                            "NaF": 10}, density=4.33, T_melt=1800.0,
                     ligands=("F",)),
    # 75TeO2-20ZnO-5Na2O (TZN), density ~5.3; 1 ErTe swap via Er2O3
    "Er-tellurite": dict(units={"TeO2": 45, "ZnO": 12, "Na2O": 3, "Er2O3": 1},
                         density=5.30, T_melt=1600.0, ligands=("O",)),
}


def build_cell(name, seed=0, dmin=1.7):
    from ase import Atoms
    from ase.formula import Formula
    from ase.data import atomic_masses, atomic_numbers
    g = GLASSES[name]
    syms = []
    for u, c in g["units"].items():
        syms += list(Formula(u * c))
    mass = sum(atomic_masses[atomic_numbers[s]] for s in syms)
    L = (mass * 1.66053906660 / g["density"]) ** (1 / 3)
    rng = np.random.default_rng(seed)
    syms = [syms[k] for k in rng.permutation(len(syms))]
    pos = []
    for s in syms:
        for _ in range(50000):
            p = rng.random(3) * L
            if not pos:
                break
            d = np.array(pos) - p
            d -= L * np.round(d / L)
            if np.min(np.linalg.norm(d, axis=1)) > dmin:
                break
        pos.append(p)
    return Atoms(syms, positions=pos, cell=[L] * 3, pbc=True)


def pick_calc(threads):
    import torch
    torch.set_num_threads(threads)
    from mace.calculators import mace_mp
    from gainmat.host.phonons import DEFAULT_MODEL, MLIP_MODELS
    fam, arg = MLIP_MODELS[DEFAULT_MODEL]
    local = os.path.expanduser(f"~/.cache/mace/{DEFAULT_MODEL}.model")
    try:
        calc = mace_mp(model=local if os.path.exists(local) else arg,
                       device="cpu", default_dtype="float32")
        zs = set(int(z) for z in calc.z_table.zs) if hasattr(calc, "z_table") else \
            set(int(z) for z in calc.models[0].atomic_numbers)
        assert {68, 70, 69, 59} <= zs, "model lacks Er/Yb/Tm/Pr"
        return calc, DEFAULT_MODEL
    except Exception as e:
        print("OMAT model unavailable/unsuitable:", repr(e), "-> MACE-MP-0 medium")
        return mace_mp(model="medium", device="cpu", default_dtype="float32"), \
            "mace-mp-0-medium"


def er_sites(frames, ligands, rcut=None):
    """First-shell analysis around every Er over `frames` (list of Atoms)."""
    # Er-ligand RDF to place cutoff at its first minimum
    ds = []
    for a in frames:
        er = [i for i, s in enumerate(a.get_chemical_symbols()) if s == "Er"]
        lig = [i for i, s in enumerate(a.get_chemical_symbols()) if s in ligands]
        for i in er:
            ds += list(a.get_distances(i, lig, mic=True))
    ds = np.array(ds)
    h, e = np.histogram(ds, bins=np.arange(1.5, 5.0, 0.05))
    c = 0.5 * (e[1:] + e[:-1])
    g = h / (c ** 2)
    gs = np.convolve(g, np.ones(5) / 5, mode="same")
    ipk = np.argmax(gs * (c < 2.8))
    if rcut is None:
        win = (c > c[ipk] + 0.2) & (c < 3.6)
        rcut = float(c[win][np.argmin(gs[win])]) if win.any() else 3.0
    sites, cns, rs = [], [], []
    for fi, a in enumerate(frames):
        syms = a.get_chemical_symbols()
        for i, s in enumerate(syms):
            if s != "Er":
                continue
            D = a.get_distances(i, list(range(len(a))), mic=True, vector=True)
            r = np.linalg.norm(D, axis=1)
            nb = [j for j in np.argsort(r) if j != i and r[j] < rcut
                  and syms[j] in ligands]
            sites.append(dict(frame=fi, er_index=i, ligands=[
                [syms[j], *[round(float(x), 4) for x in D[j]]] for j in nb]))
            cns.append(len(nb)); rs += [float(r[j]) for j in nb]
    hist = {int(k): int(v) for k, v in zip(*np.unique(cns, return_counts=True))}
    return sites, dict(CN_mean=float(np.mean(cns)), CN_hist=hist,
                       R_mean=float(np.mean(rs)), R_std=float(np.std(rs)),
                       rcut=rcut, rdf_peak=float(c[ipk]),
                       rdf_r=c.round(3).tolist(), rdf_g=gs.round(4).tolist())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("glass")
    ap.add_argument("--budget-min", type=float, default=20)
    ap.add_argument("--threads", type=int, default=os.cpu_count())
    ap.add_argument("--dt", type=float, default=1.0)
    ap.add_argument("--nframes", type=int, default=20)
    a = ap.parse_args()
    from ase import units
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase.optimize import FIRE
    from gainmat.host.glass import vdos_from_velocities
    t0 = time.time()
    g = GLASSES[a.glass]
    calc, mname = pick_calc(a.threads)
    atoms = build_cell(a.glass)
    atoms.calc = calc
    FIRE(atoms, logfile=None).run(fmax=1.0, steps=40)
    tsetup = time.time() - t0
    # measure step time
    dt = a.dt * units.fs
    MaxwellBoltzmannDistribution(atoms, temperature_K=g["T_melt"],
                                 rng=np.random.default_rng(0))
    md = Langevin(atoms, dt, temperature_K=g["T_melt"], friction=0.02 / units.fs)
    ts = time.time(); md.run(10); tstep = (time.time() - ts) / 10
    avail = a.budget_min * 60 - (time.time() - t0) - 30
    nsteps = max(200, int(avail / tstep))
    # split: melt 20 %, quench 55 %, production 25 %
    nmelt, nq = int(0.2 * nsteps), int(0.55 * nsteps)
    nprod = nsteps - nmelt - nq
    print(f"[{a.glass}] {len(atoms)} atoms, {mname}, step {tstep:.3f}s, "
          f"melt {nmelt} quench {nq} prod {nprod}", flush=True)
    md.run(nmelt)
    nstage = 40
    for T in np.linspace(g["T_melt"], 300.0, nstage):
        md.set_temperature(temperature_K=T)
        md.run(max(1, nq // nstage))
    print(f"[{a.glass}] quenched T={atoms.get_temperature():.0f} K "
          f"t={time.time()-t0:.0f}s", flush=True)
    md = Langevin(atoms, dt, temperature_K=300.0, friction=0.002 / units.fs)
    every_v = 2
    every_f = max(every_v, (nprod // a.nframes) // every_v * every_v)
    frames, vel = [], []
    for k in range(nprod // every_v):
        md.run(every_v)
        vel.append(atoms.get_velocities().copy())
        if (k * every_v) % every_f == 0 and len(frames) < a.nframes:
            f = atoms.copy(); f.calc = None; frames.append(f)
    sites, summ = er_sites(frames, g["ligands"])
    vd = vdos_from_velocities(np.array(vel), atoms.get_masses(), dt_fs=a.dt * every_v)
    rho = atoms.get_masses().sum() * 1.66053906660 / atoms.get_volume()
    out = dict(
        glass=a.glass, model=mname, n_atoms=len(atoms),
        composition_units=g["units"], density=round(float(rho), 4),
        protocol=dict(ensemble="NVT Langevin, fixed volume (experimental density)",
                      T_melt_K=g["T_melt"], T_final_K=300.0, timestep_fs=a.dt,
                      melt_ps=nmelt * a.dt / 1000, quench_ps=nq * a.dt / 1000,
                      production_ps=nprod * a.dt / 1000,
                      quench_rate_K_per_ps=(g["T_melt"] - 300) / (nq * a.dt / 1000),
                      n_frames=len(frames), threads=a.threads,
                      step_time_s=round(tstep, 4), setup_s=round(tsetup, 1),
                      wall_s=round(time.time() - t0, 1)),
        er_sites=sites, summary=summ,
        vdos=dict(cutoff99_cm1=vd["cutoff_cm1"], top_peak_cm1=vd["top_peak_cm1"],
                  note="99% cumulative mass-weighted VDOS; short production -> coarse"))
    os.makedirs("md_sites", exist_ok=True)
    with open(f"md_sites/{a.glass}.json", "w") as fh:
        json.dump(out, fh, indent=1)
    from ase.io import write
    write(f"md_sites/{a.glass}_final.extxyz", atoms)
    print(f"[{a.glass}] done: CN {summ['CN_mean']:.2f} R {summ['R_mean']:.3f}"
          f"+-{summ['R_std']:.3f} rcut {summ['rcut']} wall {time.time()-t0:.0f}s",
          flush=True)


if __name__ == "__main__":
    main()
