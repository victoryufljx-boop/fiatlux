"""Run the GPAW host steps (gap, eps_inf) under MPI from a notebook or script.

Why: inside a Jupyter/Colab kernel GPAW runs on one core (MPI is compiled in
but the kernel is a single rank, and OpenMP is off).  CaF2 at ecut 500 then
takes ~270 s and the 80-atom YAG cell takes hours.  This module writes a small
per-host driver script, runs it with

    mpiexec -n <ncpu> --allow-run-as-root --oversubscribe python <script>

with GPAW_MPI_BACKEND=cgpaw (required since GPAW 25: without it GPAW refuses to
run under an MPI launcher unless started as `gpaw python`) and OMP_NUM_THREADS=1,
and collects the JSON the script writes on rank 0.

    from gainmat.host.mpi_run import run_gpaw_mpi
    r = run_gpaw_mpi("CaF2", level="full")             # gap only
    r = run_gpaw_mpi("CaF2", level="full", do_eps=True)  # + RPA eps_inf
    r["gap"]["gap_eV"], r["eps"]["n_inf"]

The result dicts are the same as gap.compute_gap / gap.compute_eps_inf.
`run_host(..., gpaw_nproc=N)` uses this automatically.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

PKG_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

_DRIVER = r'''
import json, os, sys, time, traceback, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, {pkg_root!r})
from ase.io import read
from gpaw.mpi import world
from gainmat.host.gap import compute_gap, compute_eps_inf
atoms = read({traj!r})
gap_kw = {gap_kw!r}
eps_kw = {eps_kw!r}
out = dict(nproc=world.size, errors={{}})
t0 = time.time()
if gap_kw is not None:
    try:
        out["gap"] = compute_gap(atoms, txt={gap_txt!r}, **gap_kw)
    except Exception as e:
        out["errors"]["gap"] = "%s: %s" % (type(e).__name__, e)
        out["traceback_gap"] = traceback.format_exc(limit=4)
if eps_kw is not None:
    try:
        g = out.get("gap", {{}})
        if eps_kw.pop("scissor_from_gap", False) and g.get("gap_eV") and g.get("gap_pbe_eV"):
            eps_kw["scissor_eV"] = g["gap_eV"] - g["gap_pbe_eV"]
        out["eps"] = compute_eps_inf(atoms, workdir={workdir!r}, txt={eps_txt!r}, **eps_kw)
    except Exception as e:
        out["errors"]["eps"] = "%s: %s" % (type(e).__name__, e)
        out["traceback_eps"] = traceback.format_exc(limit=4)
out["t_wall_s"] = round(time.time() - t0, 1)
if world.rank == 0:
    with open({out_json!r}, "w") as f:
        json.dump(out, f, default=float)
'''


def mpi_env(nthreads: int = 1) -> dict:
    env = dict(os.environ)
    env.update(GPAW_MPI_BACKEND="cgpaw", OMP_NUM_THREADS=str(nthreads),
               OPENBLAS_NUM_THREADS=str(nthreads), MKL_NUM_THREADS=str(nthreads),
               OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1",
               PYTHONPATH=PKG_ROOT + os.pathsep + env.get("PYTHONPATH", ""))
    return env


def mpi_command(script: str, nproc: int | None = None, python: str | None = None) -> list:
    nproc = nproc or os.cpu_count() or 1
    mpiexec = shutil.which("mpiexec") or shutil.which("mpirun")
    py = python or sys.executable
    if not mpiexec or nproc == 1:
        return [py, script]
    return [mpiexec, "-n", str(nproc), "--allow-run-as-root", "--oversubscribe", py, script]


def write_driver(atoms, workdir: str, gap_kw: dict | None, eps_kw: dict | None,
                 tag: str = "host") -> tuple[str, str]:
    """Write <workdir>/<tag>.traj and <workdir>/run_<tag>.py; return (script, out_json)."""
    from ase.io import write
    os.makedirs(workdir, exist_ok=True)
    traj = os.path.join(workdir, f"{tag}.traj")
    write(traj, atoms)
    out_json = os.path.join(workdir, f"{tag}_gpaw.json")
    script = os.path.join(workdir, f"run_{tag}.py")
    with open(script, "w") as f:
        f.write(_DRIVER.format(pkg_root=PKG_ROOT, traj=traj, gap_kw=gap_kw, eps_kw=eps_kw,
                               gap_txt=os.path.join(workdir, f"{tag}_gap.txt"),
                               eps_txt=os.path.join(workdir, f"{tag}_eps.txt"),
                               workdir=workdir, out_json=out_json))
    return script, out_json


def run_gpaw_mpi(key_or_atoms, level: str = "full", nproc: int | None = None,
                 do_gap: bool = True, do_eps: bool = False, gap_kwargs: dict | None = None,
                 eps_kwargs: dict | None = None, workdir: str | None = None,
                 timeout: float | None = None, verbose: bool = True) -> dict:
    """Run gap (and optionally eps_inf) for one host under mpiexec; return the JSON dict."""
    from .run_host import LEVELS, KPTS_OVERRIDE, BAND_EDGES_OVERRIDE
    if isinstance(key_or_atoms, str):
        from ..structures import get_structure
        key, atoms = key_or_atoms, get_structure(key_or_atoms)
    else:
        key, atoms = key_or_atoms.get_chemical_formula(), key_or_atoms
    cfg = LEVELS[level]
    gk = None
    if do_gap:
        gk = dict(cfg["gap"])
        if key in KPTS_OVERRIDE:
            gk["kpts"] = KPTS_OVERRIDE[key][level]
        if key in BAND_EDGES_OVERRIDE:
            gk["band_edges"] = BAND_EDGES_OVERRIDE[key]
        gk.update(gap_kwargs or {})
        if "kpts" in gk and gk["kpts"] is not None:
            gk["kpts"] = tuple(gk["kpts"])
    ek = None
    if do_eps:
        ek = dict(cfg["eps"] or LEVELS["full"]["eps"])
        ek["scissor_from_gap"] = do_gap
        ek.update(eps_kwargs or {})
    nproc = nproc or os.cpu_count() or 1
    if gk is not None and nproc > 1 and "parallel" not in gk:
        # GPAW builds without ScaLAPACK (pip/Colab) can only parallelise over k-points:
        # any domain or band split trips `assert c == 1` in the BLACS matrix code.
        # So use one rank per irreducible k-point (at most nproc of them).
        nibz = n_ibz_kpoints(atoms, gk)
        nproc = max(1, min(nproc, nibz))
        gk["parallel"] = {"kpt": nproc, "domain": 1, "band": 1}
    workdir = workdir or tempfile.mkdtemp(prefix=f"gpaw_{key}_")
    script, out_json = write_driver(atoms, workdir, gk, ek, tag=key)
    cmd = mpi_command(script, nproc)
    if verbose:
        print(f"[{key}] {' '.join(cmd)}", flush=True)
    t0 = time.time()
    p = subprocess.run(cmd, env=mpi_env(), cwd=workdir, capture_output=True, text=True,
                       timeout=timeout)
    if not os.path.exists(out_json):
        raise RuntimeError(f"GPAW MPI run failed (exit {p.returncode}):\n"
                           f"{p.stdout[-2000:]}\n{p.stderr[-3000:]}")
    with open(out_json) as f:
        out = json.load(f)
    out["cmd"] = " ".join(cmd)
    out["workdir"] = workdir
    out["t_launch_s"] = round(time.time() - t0, 1)
    if verbose:
        g = out.get("gap", {})
        print(f"[{key}] {out['nproc']} ranks, {out['t_launch_s']} s: "
              f"gap({g.get('gap_method')})={g.get('gap_eV')} PBE={g.get('gap_pbe_eV')} "
              f"errors={out['errors'] or None}", flush=True)
    gerr = (out.get("errors") or {}).get("gap", "")
    if nproc > 1 and gerr and any(w in gerr for w in ("AssertionError", "BadParallelization",
                                                      "parallelize", "ranks")):
        if verbose:
            print(f"[{key}] parallel layout refused by GPAW ({gerr[:80]}); retrying serially",
                  flush=True)
        gk2 = {k: v for k, v in (gap_kwargs or {}).items() if k != "parallel"}
        gk2["parallel"] = {"kpt": 1, "domain": 1, "band": 1}
        return run_gpaw_mpi(key_or_atoms, level=level, nproc=1, do_gap=do_gap, do_eps=do_eps,
                            gap_kwargs=gk2, eps_kwargs=eps_kwargs, timeout=timeout,
                            verbose=verbose)
    return out


def n_ibz_kpoints(atoms, gap_kw: dict) -> int:
    """Number of irreducible k-points of the Gamma-centred SCF mesh compute_gap will use
    (spglib, with time reversal; matches GPAW's reduction for these crystals)."""
    import numpy as np
    import spglib
    from .gap import _kmesh
    mesh = gap_kw.get("kpts") or _kmesh(atoms, gap_kw.get("kdens", 3.0))
    cell = (atoms.cell[:], atoms.get_scaled_positions(), atoms.numbers)
    mapping, _ = spglib.get_ir_reciprocal_mesh(tuple(int(m) for m in mesh), cell,
                                               is_shift=[0, 0, 0])
    return int(len(np.unique(mapping)))


if __name__ == "__main__":   # python -m gainmat.host.mpi_run CaF2 smoke [nproc]
    k = sys.argv[1]
    lv = sys.argv[2] if len(sys.argv) > 2 else "smoke"
    n = int(sys.argv[3]) if len(sys.argv) > 3 else None
    print(json.dumps(run_gpaw_mpi(k, level=lv, nproc=n), indent=1, default=float))
