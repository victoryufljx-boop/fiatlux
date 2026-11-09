"""Fast smoke tests for the host tier (2-core CPU box, a few minutes total).

    pytest -q tests/test_host.py -s
Tests skip (not fail) when an optional heavy dependency is missing.
"""
import importlib.util
import os
import shutil
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gainmat.structures import STRUCTURES, verify_structure, get_structure  # noqa: E402
from gainmat.refs.host_refs import HOST_REFS  # noqa: E402

HAVE_GPAW = importlib.util.find_spec("gpaw") is not None
HAVE_MACE = importlib.util.find_spec("mace") is not None


def _mace_ok():
    if not HAVE_MACE:
        return False
    try:
        import torch  # noqa: F401
        from mace.calculators import mace_mp  # noqa: F401
        return True
    except Exception:
        return False


@pytest.mark.parametrize("key", list(STRUCTURES))
def test_structures(key):
    r = verify_structure(key)
    print(r)
    expected_prim = {"CaF2": 3, "LiYF4": 12, "Al2O3": 10, "Y3Al5O12": 80,
                     "SiO2-alpha-quartz": 9}
    assert r["n_prim"] == expected_prim[key]


def test_refs_complete():
    for host in ["CaF2", "LiYF4", "Al2O3", "Y3Al5O12", "SiO2-alpha-quartz",
                 "SiO2-glass", "ZBLAN"]:
        assert host in HOST_REFS
        for q in ("gap_fund_eV", "n", "phonon_max_cm1"):
            r = HOST_REFS[host][q]
            assert r["source"] and "verified" in r


@pytest.mark.skipif(not HAVE_GPAW, reason="gpaw not installed")
def test_caf2_gap_smoke():
    from gainmat.host.gap import compute_gap
    g = compute_gap(get_structure("CaF2"), ecut=300, kpts=(2, 2, 2))
    print("CaF2 smoke gap:", g)
    assert g["gap_method"] == "gllbsc" and g["gap_eV"] == g["gap_gllbsc_eV"]
    assert 5.0 < g["gap_pbe_eV"] < 9.0          # PBE ~7 eV
    assert g["gap_gllbsc_eV"] > g["gap_pbe_eV"] + 2.0   # GLLB-SC opens it
    assert 9.0 < g["gap_gllbsc_eV"] < 15.0       # exp 12.1 eV (fundamental)


def test_gap_defaults():
    """Bake-off winner (BAKEOFF_GAPS.md): GLLB-SC is the default everywhere."""
    from gainmat.host.run_host import LEVELS
    for lv in LEVELS.values():
        assert lv["gap"]["method"] == "gllbsc"
        assert lv["gap"]["ecut"] >= 300
    assert LEVELS["full"]["gap"]["ecut"] >= 500
    if HAVE_GPAW:
        from gainmat.host import gap
        assert gap.DEFAULT_GAP_METHOD == "gllbsc"
        assert set(gap.GAP_METHODS) == {"gllbsc", "pbe", "tb09", "hse06"}
        with pytest.raises(ValueError):
            gap.compute_gap(get_structure("CaF2"), method="b3lyp")


def test_mpi_command():
    from gainmat.host.mpi_run import mpi_command, mpi_env
    env = mpi_env()
    assert env["GPAW_MPI_BACKEND"] == "cgpaw" and env["OMP_NUM_THREADS"] == "1"
    cmd = mpi_command("x.py", nproc=12)
    if shutil.which("mpiexec") or shutil.which("mpirun"):
        assert cmd[1:5] == ["-n", "12", "--allow-run-as-root", "--oversubscribe"]
    assert cmd[-1] == "x.py"
    assert mpi_command("x.py", nproc=1)[-1] == "x.py"


@pytest.mark.skipif(not HAVE_GPAW or not shutil.which("mpiexec"),
                    reason="gpaw or mpiexec missing")
def test_caf2_gap_mpi_smoke(tmp_path):
    from gainmat.host.mpi_run import run_gpaw_mpi
    n = min(2, os.cpu_count() or 1)
    r = run_gpaw_mpi("CaF2", level="smoke", nproc=n, workdir=str(tmp_path),
                     gap_kwargs=dict(ecut=300, kpts=(2, 2, 2)))
    assert not r["errors"], r["errors"]
    assert r["nproc"] == n
    assert 9.0 < r["gap"]["gap_eV"] < 15.0


@pytest.mark.skipif(not _mace_ok(), reason="mace/torch not importable")
def test_caf2_mace_phonons():
    from gainmat.host.phonons import compute_phonons
    r = compute_phonons(get_structure("CaF2"), device="cpu",
                        supercell_min_len=7.0, mesh=(8, 8, 8))
    print("CaF2 MACE phonons: max %.1f eff %.1f cm-1, Gamma %s" %
          (r["phonon_max_cm1"], r["phonon_eff_cm1"], r["gamma_cm1"]))
    assert r["min_freq_cm1"] > -10.0           # dynamically stable
    # default: MACE-OMAT-0, exp cell, NAC -> max includes LO (exp 463 cm^-1)
    assert r["nac"] and r["cell_protocol"] == "expcell"
    assert 420 < r["phonon_max_cm1"] < 520
    modes = {m["irrep"]: m["freq_cm1"] for m in r["gamma_modes"]}
    assert abs(modes["T2g"] - 322) < 25 and abs(modes["T1u"] - 257) < 30


@pytest.mark.skipif(not _mace_ok(), reason="mace/torch not importable")
def test_caf2_legacy_protocol_selectable():
    from gainmat.host.phonons import compute_phonons
    r = compute_phonons(get_structure("CaF2"), device="cpu", model="medium",
                        relax_cell=True, nac=False, supercell_min_len=7.0,
                        mesh=(6, 6, 6))
    assert not r["nac"] and r["model"] == "mace-mp-0-medium"
    assert r["phonon_max_cm1"] == r["phonon_max_TO_cm1"] < 420


def test_phonon_refs_and_matching():
    from gainmat.refs.phonon_refs import PHONON_REFS, match_modes
    for h, ref in PHONON_REFS.items():
        assert ref["source"] and ref["url"]
    comp = [dict(freq_cm1=300.0, irrep="T2g", deg=3), dict(freq_cm1=250.0, irrep="T1u", deg=3)]
    pairs = match_modes(comp, {"T2g": [322], "T1u": [257]})
    assert pairs == [("T2g", 322, 300.0), ("T1u", 257, 250.0)]


@pytest.mark.skipif(not _mace_ok() or not os.environ.get("GAINMAT_SLOW"),
                    reason="set GAINMAT_SLOW=1 (GPAW+MACE end-to-end)")
def test_run_host_caf2_smoke():
    from gainmat.host.run_host import run_host, compare_to_experiment
    hp = run_host("CaF2", level="smoke", do_gap=HAVE_GPAW, device="cpu")
    assert not hp.method.get("errors"), hp.method.get("errors")
    rows = compare_to_experiment("CaF2", hp)
    for r in rows:
        print(r["quantity"], r["computed"], r["experiment"], r["err_pct"])
    assert hp.phonon_max_cm1 is not None


@pytest.mark.skipif(not _mace_ok() or not os.environ.get("GAINMAT_SLOW"),
                    reason="glass MD smoke test is slow on CPU; set GAINMAT_SLOW=1")
def test_glass_smoke():
    from gainmat.host.glass import melt_quench
    r = melt_quench("SiO2-glass", preset="smoke", device="cpu", n_target=24,
                    log=lambda *a: None)
    print("SiO2 glass smoke:", r["n_atoms"], r["phonon_eff_cm1"], r["phonon_max_cm1"])
    assert r["phonon_eff_cm1"] > 0


@pytest.mark.skipif(not HAVE_GPAW, reason="gpaw not installed")
def test_lif_eps_inf_smoke(tmp_path):
    """Unshifted PBE RPA+LFE is the reported eps_inf; the scissor lowers it."""
    from ase.build import bulk
    from gainmat.host.gap import compute_eps_inf
    r = compute_eps_inf(bulk("LiF", "rocksalt", a=4.0263), ecut=300, kpts=(4, 4, 4),
                        ecut_resp=30, nbands_factor=4, scissor_eV=5.0,
                        workdir=str(tmp_path))
    print("LiF eps smoke:", r)
    assert r["directions"] == ["x"]
    assert r["eps_inf"] == r["eps_inf_rpa"]
    assert 1.6 < r["eps_inf"] < 2.3            # exp 1.92-1.93
    assert r["eps_inf_rpa_scissor"] < r["eps_inf"]
