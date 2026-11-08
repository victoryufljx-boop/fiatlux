"""Melt-quench MD for glasses with the MACE foundation model (ASE), and the
vibrational DOS from the velocity autocorrelation function (VACF).

Protocol (defaults sized for an A100; use `preset="smoke"` on CPU):
  1. random packing of the composition at the experimental glass density
     (minimum-distance rejection), short FIRE pre-relaxation
  2. NVT Langevin melt at T_melt
  3. linear quench T_melt -> 300 K (NVT; optional NPT tail via ase NPTBerendsen
     to let the density relax)
  4. NVE/low-friction production at 300 K, record velocities
  5. VDOS = |FFT(VACF)|  -> cutoff = frequency below which `eff_fraction` of
     the mass-weighted VDOS lies, plus the highest-frequency peak position.

Quench rates reachable in MD (1e12-1e14 K/s) are 10+ orders of magnitude
faster than lab rates; the high-frequency stretch band (which is what sets the
multiphonon cutoff) is local and fairly insensitive to this, but ring
statistics / density are not.  Treat glass results as semi-quantitative.
"""
from __future__ import annotations

import time
import numpy as np

THZ_TO_CM1 = 33.35641
CM1_PER_INV_FS = 1e15 / 2.99792458e10   # 1/fs -> cm^-1

# composition: {formula unit: mole fraction}, density g/cm3, T_melt K
GLASSES = {
    "SiO2-glass": dict(units={"SiO2": 1.0}, density=2.20, T_melt=4000.0),
    # ZBLAN 53ZrF4-20BaF2-4LaF3-3AlF3-20NaF, density 4.33 g/cm3
    # (Zhu & Peyghambarian, Adv. OptoElectron. 2010, 501956)
    "ZBLAN": dict(units={"ZrF4": 0.53, "BaF2": 0.20, "LaF3": 0.04,
                         "AlF3": 0.03, "NaF": 0.20},
                  density=4.33, T_melt=1500.0),
}

PRESETS = {
    # A100: ~600 atoms, 10 ps melt, 4000 K -> 300 K over ~20 ps (~1.85e14 K/s),
    # 10 ps production.
    "full": dict(n_target=600, dt_fs=1.0, melt_ps=10.0, quench_ps=20.0,
                 npt_ps=5.0, prod_ps=10.0, sample_every=2),
    "medium": dict(n_target=300, dt_fs=1.0, melt_ps=5.0, quench_ps=10.0,
                   npt_ps=0.0, prod_ps=5.0, sample_every=2),
    # CPU smoke test: proves the pipeline runs; numbers are not converged.
    "smoke": dict(n_target=24, dt_fs=1.0, melt_ps=0.03, quench_ps=0.03,
                  npt_ps=0.0, prod_ps=0.1, sample_every=1),
}


def _integer_units(units: dict, n_target: int):
    """Choose integer formula-unit counts approximating mole fractions."""
    from ase.formula import Formula
    natoms_per = {u: len(Formula(u)) for u in units}
    mean = sum(units[u] * natoms_per[u] for u in units)
    scale = n_target / mean
    counts = {u: max(1, int(round(units[u] * scale))) for u in units}
    return counts


def random_glass_cell(name: str, n_target: int = 300, seed: int = 0,
                      dmin: float = 1.6):
    from ase import Atoms
    from ase.formula import Formula
    from ase.data import atomic_masses, atomic_numbers
    g = GLASSES[name]
    counts = _integer_units(g["units"], n_target)
    syms = []
    for u, c in counts.items():
        syms += list(Formula(u * c))
    mass = sum(atomic_masses[atomic_numbers[s]] for s in syms)
    vol = mass * 1.66053906660 / g["density"]
    L = vol ** (1 / 3)
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(syms))
    pos = []
    for k in order:
        for _ in range(20000):
            p = rng.random(3) * L
            if not pos:
                break
            d = np.array(pos) - p
            d -= L * np.round(d / L)
            if np.min(np.linalg.norm(d, axis=1)) > dmin:
                break
        pos.append(p)
    atoms = Atoms([syms[k] for k in order], positions=pos, cell=[L] * 3, pbc=True)
    atoms.info["counts"] = counts
    return atoms


def vdos_from_velocities(vel, masses, dt_fs: float, eff_fraction: float = 0.99,
                         smooth_cm1: float = 10.0):
    """Mass-weighted VDOS from velocities (nframes, natoms, 3)."""
    v = np.asarray(vel) * np.sqrt(np.asarray(masses))[None, :, None]
    n = v.shape[0]
    # Wiener-Khinchin: power spectrum of velocities == FT of VACF
    win = np.hanning(n)[:, None, None]
    spec = np.abs(np.fft.rfft(v * win, axis=0)) ** 2
    dos = spec.sum(axis=(1, 2))
    freq = np.fft.rfftfreq(n, d=dt_fs) * CM1_PER_INV_FS
    if smooth_cm1 and len(freq) > 3:
        df = freq[1] - freq[0]
        s = max(1, int(round(smooth_cm1 / df)))
        k = np.exp(-0.5 * (np.arange(-3 * s, 3 * s + 1) / s) ** 2)
        dos = np.convolve(dos, k / k.sum(), mode="same")
    dos[0] = 0.0
    c = np.cumsum(dos); c /= c[-1]
    f_eff = float(freq[np.searchsorted(c, eff_fraction)])
    # highest-frequency peak above half of the median band level
    hi = freq > 0.5 * f_eff
    peak = float(freq[hi][np.argmax(dos[hi])]) if hi.any() else None
    return dict(freq_cm1=freq, dos=dos, cutoff_cm1=f_eff, top_peak_cm1=peak)


def melt_quench(name: str, preset: str = "smoke", calc=None, model: str = "medium",
                device=None, seed: int = 0, T_melt: float | None = None,
                T_final: float = 300.0, eff_fraction: float = 0.99,
                traj: str | None = None, log=print, **overrides) -> dict:
    from ase import units
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase.optimize import FIRE
    from .phonons import get_mace_calc, get_device
    p = dict(PRESETS[preset]); p.update(overrides)
    g = GLASSES[name]
    T_melt = T_melt or g["T_melt"]
    t0 = time.time()
    if calc is None:
        calc = get_mace_calc(model=model, device=device, default_dtype="float32")
    atoms = random_glass_cell(name, p["n_target"], seed=seed)
    atoms.calc = calc
    FIRE(atoms, logfile=None).run(fmax=0.5, steps=20 if preset == "smoke" else 200)
    dt = p["dt_fs"] * units.fs
    MaxwellBoltzmannDistribution(atoms, temperature_K=T_melt,
                                 rng=np.random.default_rng(seed))
    md = Langevin(atoms, dt, temperature_K=T_melt, friction=0.02 / units.fs)
    nmelt = int(p["melt_ps"] * 1000 / p["dt_fs"])
    md.run(nmelt)
    log(f"[{name}] melt done ({nmelt} steps, T={atoms.get_temperature():.0f} K)")
    # linear quench in 50 stages
    nq = int(p["quench_ps"] * 1000 / p["dt_fs"])
    nstage = max(1, min(50, nq))
    for k, T in enumerate(np.linspace(T_melt, T_final, nstage)):
        md.set_temperature(temperature_K=T)
        md.run(max(1, nq // nstage))
    log(f"[{name}] quench done (T={atoms.get_temperature():.0f} K)")
    if p.get("npt_ps", 0) > 0:
        try:
            from ase.md.nptberendsen import NPTBerendsen
            npt = NPTBerendsen(atoms, dt, temperature_K=T_final, pressure_au=0.0,
                               taut=100 * units.fs, taup=1000 * units.fs,
                               compressibility_au=1e-6 / units.bar)
            npt.run(int(p["npt_ps"] * 1000 / p["dt_fs"]))
        except Exception as e:
            log(f"[{name}] NPT skipped: {e!r}")
    # production: weak friction (keeps T, barely perturbs spectrum)
    md = Langevin(atoms, dt, temperature_K=T_final, friction=1e-4 / units.fs)
    nprod = int(p["prod_ps"] * 1000 / p["dt_fs"])
    every = p["sample_every"]
    vel = []
    for _ in range(nprod // every):
        md.run(every)
        vel.append(atoms.get_velocities().copy())
    vd = vdos_from_velocities(np.array(vel), atoms.get_masses(),
                              dt_fs=p["dt_fs"] * every, eff_fraction=eff_fraction)
    if traj:
        atoms.write(traj)
    rho = atoms.get_masses().sum() * 1.66053906660 / atoms.get_volume()
    return dict(glass=name, preset=preset, n_atoms=len(atoms),
                counts=atoms.info.get("counts"), T_melt=T_melt,
                quench_rate_K_per_ps=(T_melt - T_final) / max(p["quench_ps"], 1e-9),
                density=float(rho), phonon_eff_cm1=vd["cutoff_cm1"],
                phonon_max_cm1=vd["top_peak_cm1"],
                vdos=dict(freq_cm1=vd["freq_cm1"].tolist(), dos=vd["dos"].tolist()),
                model=f"mace_mp:{model}", device=get_device(device),
                t_s=round(time.time() - t0, 1), _atoms=atoms)
