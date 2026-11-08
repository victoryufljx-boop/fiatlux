"""Host-tier driver: run_host(key, level) -> schema.HostProps.

Every sub-step is wrapped so that a failure (missing package, OOM, SCF not
converging) is recorded in HostProps.method["errors"] and the rest still runs.

Levels
  "smoke" : CPU seconds-minutes; tiny k-mesh / cutoff; proves the code runs.
  "fast"  : CPU-friendly-ish; PBE+GLLB-SC on a modest mesh, MACE phonons.
  "full"  : Colab A100 settings; converged-ish meshes, path band edges,
            RPA eps_inf (optional), glass MD with ~600 atoms.
"""
from __future__ import annotations

import traceback

from ..schema import HostProps, ev_to_nm

CRYSTAL_KEYS = ["CaF2", "LiYF4", "Al2O3", "Y3Al5O12", "SiO2-alpha-quartz"]
GLASS_KEYS = ["SiO2-glass", "ZBLAN"]

LEVELS = {
    "smoke": dict(gap=dict(ecut=300, kdens=1.5, band_edges="mesh", method="gllbsc"),
                  phon=dict(supercell_min_len=7.0, mesh=(8, 8, 8)),
                  eps=None, glass_preset="smoke"),
    "fast": dict(gap=dict(ecut=500, kdens=2.5, band_edges="mesh", method="gllbsc"),
                 phon=dict(supercell_min_len=9.0, mesh=(16, 16, 16)),
                 eps=None, glass_preset="medium"),
    "full": dict(gap=dict(ecut=500, kdens=3.5, band_edges="path", method="gllbsc"),
                 phon=dict(supercell_min_len=10.0, mesh=(24, 24, 24)),
                 eps=dict(ecut=500, kdens=4.0, ecut_resp=50.0, nbands_factor=6),
                 glass_preset="full"),
}

# YAG primitive cell (80 atoms, 12 A) is expensive in DFT: cap its k-mesh.
KPTS_OVERRIDE = {"Y3Al5O12": {"smoke": (1, 1, 1), "fast": (2, 2, 2), "full": (2, 2, 2)}}


def _err(method, step, e):
    method.setdefault("errors", {})[step] = f"{type(e).__name__}: {e}"
    method.setdefault("tracebacks", {})[step] = traceback.format_exc(limit=3)


def _store_gap(hp, g):
    hp.gap_pbe_eV = g.get("gap_pbe_eV")
    hp.gap_gllbsc_eV = g.get("gap_gllbsc_eV")
    hp.gap_eV = g.get("gap_eV")
    hp.gap_method = g.get("gap_method")
    if hp.gap_eV:
        hp.uv_edge_nm = ev_to_nm(hp.gap_eV)


def run_host(key: str, level: str = "fast", do_gap: bool = True,
             do_phonons: bool = True, do_eps: bool | None = None,
             mace_model: str = "medium", device=None, gap_kwargs=None,
             phon_kwargs=None, verbose: bool = True,
             gap_method: str | None = None,
             gpaw_nproc: int | None = None) -> HostProps:
    """gap_method: "gllbsc" (default, see BAKEOFF_GAPS.md), "pbe", "tb09", "hse06".
    gpaw_nproc: if > 1, run the GPAW steps (gap, eps_inf) in a subprocess under
    `mpiexec -n gpaw_nproc` (gainmat.host.mpi_run) instead of in this process.
    Use this in notebooks: the kernel itself is a single serial rank.
    """
    cfg = LEVELS[level]
    hp = HostProps()
    m = hp.method
    m["level"] = level
    m["key"] = key
    log = print if verbose else (lambda *a, **k: None)

    if key in GLASS_KEYS:
        if do_phonons:
            try:
                from .glass import melt_quench
                r = melt_quench(key, preset=cfg["glass_preset"], model=mace_model,
                                device=device, log=log)
                hp.phonon_eff_cm1 = r["phonon_eff_cm1"]
                hp.phonon_max_cm1 = r["phonon_max_cm1"]
                m["phonons"] = {k: v for k, v in r.items()
                                if k not in ("_atoms", "vdos")}
                m["vdos"] = r["vdos"]
            except Exception as e:
                _err(m, "glass_md", e)
        m["note"] = "glass: band gap / n not computed in this tier (use crystal proxy or DFT on MD snapshots)"
        return hp

    from ..structures import get_structure
    atoms = get_structure(key, primitive=True)
    m["n_atoms"] = len(atoms)

    if do_eps is None:
        do_eps = cfg["eps"] is not None
    gk = dict(cfg["gap"])
    if key in KPTS_OVERRIDE:
        gk["kpts"] = KPTS_OVERRIDE[key][level]
    if gap_method:
        gk["method"] = gap_method
    gk.update(gap_kwargs or {})

    mpi_res = None
    if gpaw_nproc and gpaw_nproc > 1 and (do_gap or do_eps):
        try:
            from .mpi_run import run_gpaw_mpi
            log(f"[{key}] GPAW under MPI x{gpaw_nproc}: gap {gk} eps={do_eps}")
            mpi_res = run_gpaw_mpi(key, level=level, nproc=gpaw_nproc, do_gap=do_gap,
                                   do_eps=do_eps, gap_kwargs=gk, verbose=verbose)
            m["mpi"] = dict(cmd=mpi_res["cmd"], nproc=mpi_res["nproc"],
                            t_launch_s=mpi_res["t_launch_s"])
            for step, msg in mpi_res.get("errors", {}).items():
                m.setdefault("errors", {})[step] = msg
        except Exception as e:
            _err(m, "gpaw_mpi", e)
            mpi_res = {}

    if do_gap:
        try:
            if mpi_res is not None:
                g = mpi_res.get("gap")
                if g is None:
                    raise RuntimeError("no gap from MPI run")
            else:
                from .gap import compute_gap
                log(f"[{key}] GPAW gap {gk}")
                g = compute_gap(atoms, **gk)
            _store_gap(hp, g)
            m["gap"] = dict(code="GPAW", **g)
        except Exception as e:
            _err(m, "gap", e)

    if do_eps:
        try:
            if mpi_res is not None:
                r = mpi_res.get("eps")
                if r is None:
                    raise RuntimeError("no eps from MPI run")
            else:
                from .gap import compute_eps_inf
                ek = dict(cfg["eps"] or LEVELS["full"]["eps"])
                scis = None
                if hp.gap_pbe_eV and hp.gap_eV:
                    scis = hp.gap_eV - hp.gap_pbe_eV
                r = compute_eps_inf(atoms, scissor_eV=scis, **ek)
            hp.eps_inf = r["eps_inf"]
            hp.n_inf = r["n_inf"]
            m["eps"] = r
        except Exception as e:
            _err(m, "eps_inf", e)

    if do_phonons:
        try:
            from .phonons import compute_phonons
            pk = dict(cfg["phon"]); pk.update(phon_kwargs or {})
            log(f"[{key}] MACE phonons {pk}")
            r = compute_phonons(atoms, model=mace_model, device=device, **pk)
            hp.phonon_max_cm1 = r["phonon_max_cm1"]
            hp.phonon_eff_cm1 = r["phonon_eff_cm1"]
            m["phonons"] = {k: v for k, v in r.items() if k != "dos"}
            m["phonon_dos"] = r["dos"]
        except Exception as e:
            _err(m, "phonons", e)
    return hp


def compare_to_experiment(key: str, hp: HostProps) -> list:
    """Rows (quantity, computed, experiment, %err, exp kind, source)."""
    from ..refs.host_refs import HOST_REFS
    ref = HOST_REFS.get(key, {})
    rows = []

    def add(q, comp, rk):
        r = ref.get(rk)
        if r is None:
            return
        ex = r["value"]
        err = None if (comp is None or ex in (None, 0)) else 100 * (comp - ex) / ex
        rows.append(dict(host=key, quantity=q, computed=comp, experiment=ex,
                         err_pct=err, exp_kind=r["kind"], source=r["source"],
                         ref_verified=r["verified"]))
    add("gap_PBE_eV", hp.gap_pbe_eV, "gap_fund_eV")
    add("gap_GLLBSC_eV", hp.gap_gllbsc_eV, "gap_fund_eV")
    add("gap_GLLBSC_vs_exciton_eV", hp.gap_gllbsc_eV, "gap_exciton_eV")
    if hp.gap_method and hp.gap_method not in ("gllbsc", "pbe"):
        add(f"gap_{hp.gap_method.upper()}_eV", hp.gap_eV, "gap_fund_eV")
    add("n_inf", hp.n_inf, "n")
    add("phonon_max_cm1 (vs highest LO/optical)", hp.phonon_max_cm1, "phonon_LO_cm1")
    add("phonon_eff_cm1 (vs energy-gap-law value)", hp.phonon_eff_cm1, "phonon_max_cm1")
    return rows
