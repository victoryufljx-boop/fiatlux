"""Electronic band gap (GLLB-SC default; PBE, TB09, HSE06 selectable) and optical dielectric constant (RPA).

API checked against the installed GPAW 26.7.0 source and the GPAW tutorial
"Calculating band gap using the GLLB-sc functional"
(https://gpaw.readthedocs.io/tutorialsexercises/electronic/band_gap/band_gap.html):

    calc = GPAW(mode=PW(ecut), xc='GLLBSC', kpts=..., occupations=FermiDirac(0.01))
    atoms.get_potential_energy()
    homo, lumo = calc.get_homo_lumo()
    response = calc.hamiltonian.xc.response
    dxc_pot = response.calculate_discontinuity_potential(homo, lumo)
    KS_gap, dxc = response.calculate_discontinuity(dxc_pot)
    QP_gap = KS_gap + dxc

The VBM and CBM must lie on the k-mesh: we use Gamma-centred meshes with even
divisions, which contain Gamma, X, L, ... for the cubic/tetragonal/hexagonal hosts here.
For a production run, use `band_edges="path"` (non-self-consistent band
structure along the high-symmetry path) as in the tutorial's second recipe.

Expected accuracy (literature, e.g. Castelli et al. Energy Environ. Sci. 5,
5814 (2012); Tran et al. PRB 2018): PBE underestimates gaps of wide-gap
ionics by 30-45 %; GLLB-SC is usually within ~10-15 % of the experimental
fundamental gap (it tends to overshoot slightly for some fluorides).
"""
from __future__ import annotations

import time
import numpy as np


def _gpaw_available():
    try:
        import gpaw  # noqa: F401
        return True
    except Exception:
        return False


def _kmesh(atoms, kdens):
    """Gamma-centred even mesh with ~kdens points per reciprocal A (x 2pi)."""
    rec = np.linalg.norm(atoms.cell.reciprocal(), axis=1)  # 1/A (no 2pi)
    k = [max(2, int(np.ceil(kdens * r))) for r in rec]
    k = [x + (x % 2) for x in k]   # even -> contains zone-boundary points
    return tuple(k)


GAP_METHODS = ("gllbsc", "pbe", "tb09", "hse06")
DEFAULT_GAP_METHOD = "gllbsc"


def compute_gap(atoms, ecut: float = 450.0, kpts=None, kdens: float = 3.0,
                method: str | tuple | list = DEFAULT_GAP_METHOD,
                do_pbe: bool | None = None, do_gllbsc: bool | None = None,
                band_edges: str = "mesh", txt: str | None = None,
                parallel: dict | None = None, nbands_extra: int = 8) -> dict:
    """Band gap with GPAW plane waves.

    method : "gllbsc" (default; best accuracy per cost for wide-gap fluorides
             and oxides, see BAKEOFF_GAPS.md), "pbe", "tb09" (TB-mBJ, ~5-10x
             the cost, underestimates by 1-2 eV with GPAW PAW setups),
             "hse06" (non-self-consistent HSE06@PBE at the band-edge k-points;
             underestimates wide gaps by 1-3 eV).  A list runs several.
             PBE is always run alongside (it is cheap and gives the scissor for
             eps_inf) unless do_pbe=False.
    do_pbe / do_gllbsc : legacy switches (override `method` if given).
    band_edges : "mesh" (gap on the SCF mesh) or "path" (non-SCF band
                 structure along the high-symmetry path; PBE and GLLB-SC only).
    Returns dict with gap_eV (from the primary method), gap_method,
    gap_pbe_eV, gap_gllbsc_eV, gap_gllbsc_KS_eV, dxc_eV, gap_tb09_eV,
    gap_hse06_eV, timings ...
    """
    from gpaw import GPAW, PW, FermiDirac
    methods = [method] if isinstance(method, str) else list(method)
    methods = [m.lower().replace("-", "").replace("_", "") for m in methods]
    if do_gllbsc is not None:
        methods = [m for m in methods if m != "gllbsc"] + (["gllbsc"] if do_gllbsc else [])
    for m_ in methods:
        if m_ not in GAP_METHODS:
            raise ValueError(f"unknown gap method {m_!r}; choose from {GAP_METHODS}")
    if do_pbe is None:
        do_pbe = True
    primary = methods[0] if methods else "pbe"
    out = dict(ecut=ecut, band_edges=band_edges, gap_method=primary)
    if kpts is None:
        kpts = _kmesh(atoms, kdens)
    out["kpts"] = list(kpts)
    kp = {"size": tuple(kpts), "gamma": True}
    common = dict(mode=PW(ecut), kpts=kp, occupations=FermiDirac(0.01),
                  txt=txt, convergence={"bands": "occupied"})
    if parallel:
        common["parallel"] = parallel

    pbe_calc = None
    if do_pbe or "pbe" in methods or "hse06" in methods:
        t0 = time.time()
        a = atoms.copy()
        calc = GPAW(xc="PBE", **common)
        a.calc = calc
        a.get_potential_energy()
        pbe_calc = calc
        if band_edges == "path":
            calc = _path_calc(calc, a)
        homo, lumo = calc.get_homo_lumo()
        out["gap_pbe_eV"] = float(lumo - homo)
        out["t_pbe_s"] = round(time.time() - t0, 1)

    if "gllbsc" in methods:
        t0 = time.time()
        a = atoms.copy()
        calc = GPAW(xc="GLLBSC", **common)
        a.calc = calc
        a.get_potential_energy()
        response = calc.hamiltonian.xc.response
        if band_edges == "path":
            bs = _path_calc(calc, a)
            homo, lumo = bs.get_homo_lumo()
            dxc_pot = response.calculate_discontinuity_potential(homo, lumo)
            ks, dxc = bs.hamiltonian.xc.response.calculate_discontinuity(dxc_pot)
        else:
            homo, lumo = calc.get_homo_lumo()
            dxc_pot = response.calculate_discontinuity_potential(homo, lumo)
            ks, dxc = response.calculate_discontinuity(dxc_pot)
        out["gap_gllbsc_KS_eV"] = float(ks)
        out["dxc_eV"] = float(dxc)
        out["gap_gllbsc_eV"] = float(ks + dxc)
        out["t_gllbsc_s"] = round(time.time() - t0, 1)

    if "tb09" in methods:
        from gpaw import Mixer
        t0 = time.time()
        a = atoms.copy()
        c2 = dict(common)
        calc = GPAW(xc={"name": "TB09", "stencil": 1}, mixer=Mixer(0.5, 5, 50.0),
                    nbands=-nbands_extra, **c2)
        a.calc = calc
        a.get_potential_energy()
        homo, lumo = calc.get_homo_lumo()
        out["gap_tb09_eV"] = float(lumo - homo)
        try:
            out["tb09_c"] = float(calc.hamiltonian.xc.c)
        except Exception:
            pass
        out["t_tb09_s"] = round(time.time() - t0, 1)

    if "hse06" in methods:
        t0 = time.time()
        out["gap_hse06_eV"] = _hse06_nsc_gap(pbe_calc)
        out["t_hse06_s"] = round(time.time() - t0, 1)

    out["gap_eV"] = out.get(f"gap_{primary}_eV")
    return out


def _hse06_nsc_gap(calc) -> float:
    """Non-self-consistent HSE06@PBE gap at the PBE VBM/CBM k-points (+Gamma)."""
    from gpaw.hybrids.eigenvalues import non_self_consistent_eigenvalues as nsc
    nk = len(calc.get_ibz_k_points())
    nv = int(round(calc.get_number_of_electrons())) // 2
    e = np.array([calc.get_eigenvalues(kpt=k) for k in range(nk)])
    ks = sorted({int(np.argmax(e[:, nv - 1])), int(np.argmin(e[:, nv])), 0})
    edft, vdft, vhyb = nsc(calc, "HSE06", n1=nv - 1, n2=nv + 1, kpt_indices=ks)
    eh = edft - vdft + vhyb
    return float(eh[0, :, 1].min() - eh[0, :, 0].max())


def _path_calc(calc, atoms, npoints: int = 60):
    path = atoms.cell.bandpath(npoints=npoints)
    nocc = int(round(calc.get_number_of_electrons())) // 2
    return calc.fixed_density(nbands=nocc + 8, kpts=path, symmetry="off",
                              convergence={"bands": nocc + 4}, txt=None)


# ---------------------------------------------------------------- eps_inf
def compute_eps_inf(atoms, ecut: float = 450.0, kpts=None, kdens: float = 4.0,
                    ecut_resp: float = 50.0, nbands_factor: int = 6,
                    scissor_eV: float | None = None, directions=None,
                    direction: str | None = None,
                    workdir: str = ".", txt: str | None = None) -> dict:
    """RPA macroscopic dielectric constant (with local fields) via gpaw.response.

    PBE ground state -> fixed-density diagonalisation with nbands_factor x nocc
    bands -> DielectricFunction(...).get_macroscopic_dielectric_constant().

    Bake-off result (BAKEOFF_GAPS.md): the *unshifted* PBE RPA+LFE value is the
    best estimate of eps_inf for wide-gap fluorides/oxides (within +0-7 % of
    experiment for LiF, CaF2, Al2O3, quartz; +6 % for MgO).  Adding a scissor
    (`eshift`, conduction bands shifted up by scissor_eV, e.g. gap_gllbsc -
    gap_pbe) without renormalising the momentum matrix elements lowers eps_inf
    by 10-15 % and makes it worse, so `eps_inf` is always the unshifted value;
    the scissored one is returned as `eps_inf_rpa_scissor` for reference.
    GPAW has no DFPT or periodic finite-field (Berry-phase field) eps_inf.

    directions : e.g. ["x"] (cubic) or ["x", "z"] (uniaxial: x = ordinary,
                 z = extraordinary for hexagonal cells with c || z).  Default:
                 ["x"] for cubic cells, ["x", "y", "z"] otherwise.
                 eps_inf is the trace average; per-direction values are kept.
    Convergence (LiF): nbands 6->8 x nocc changes eps_inf by < 0.2 %.
    """
    import os
    from gpaw import GPAW, PW, FermiDirac
    from gpaw.mpi import world
    from gpaw.response.df import DielectricFunction
    if kpts is None:
        kpts = _kmesh(atoms, kdens)
    if direction is not None:
        directions = [direction]
    if directions is None:
        cp = atoms.cell.cellpar()
        cubic = np.allclose(cp[:3], cp[0], rtol=1e-3) and np.allclose(cp[3:], cp[3], atol=0.1)
        directions = ["x"] if cubic else ["x", "y", "z"]
    a = atoms.copy()
    calc = GPAW(mode=PW(ecut), xc="PBE", kpts={"size": tuple(kpts), "gamma": True},
                occupations=FermiDirac(0.01), txt=txt)
    a.calc = calc
    a.get_potential_energy()
    nocc = int(round(calc.get_number_of_electrons())) // 2
    nb = int(nbands_factor * nocc)
    calc = calc.fixed_density(nbands=nb + 8, convergence={"bands": nb}, txt=txt)
    gpw = os.path.join(workdir, "eps_gs.gpw")
    calc.write(gpw, mode="all")
    world.barrier()
    out = dict(kpts=list(kpts), ecut_resp=ecut_resp, nbands=nb, directions=list(directions))

    def _eps(shift):
        res = {}
        for d in directions:
            df = DielectricFunction(gpw, frequencies=[0.0], ecut=ecut_resp, nbands=nb,
                                    intraband=False, hilbert=False, eta=1e-4,
                                    eshift=shift, txt=txt)
            e_nlf, e_lf = df.get_macroscopic_dielectric_constant(direction=d)
            res[d] = (float(np.real(e_nlf)), float(np.real(e_lf)))
        return res

    r0 = _eps(None)
    out["eps_inf_rpa_dir"] = {d: v[1] for d, v in r0.items()}
    out["eps_inf_rpa_nlfe"] = float(np.mean([v[0] for v in r0.values()]))
    out["eps_inf_rpa"] = float(np.mean([v[1] for v in r0.values()]))
    if scissor_eV:
        r1 = _eps(scissor_eV)
        out["eps_inf_rpa_scissor_dir"] = {d: v[1] for d, v in r1.items()}
        out["eps_inf_rpa_scissor"] = float(np.mean([v[1] for v in r1.values()]))
        out["scissor_eV"] = scissor_eV
    out["eps_inf"] = out["eps_inf_rpa"]
    out["n_inf"] = float(np.sqrt(out["eps_inf"]))
    out["n_inf_dir"] = {d: float(np.sqrt(v)) for d, v in out["eps_inf_rpa_dir"].items()}
    return out
