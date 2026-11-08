"""Fast smoke tests for the host tier (2-core CPU box, a few minutes total).

    pytest -q tests/test_host.py -s
Tests skip (not fail) when an optional heavy dependency is missing.
"""
import importlib.util
import os
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
    assert 5.0 < g["gap_pbe_eV"] < 9.0          # PBE ~7 eV
    assert g["gap_gllbsc_eV"] > g["gap_pbe_eV"] + 2.0   # GLLB-SC opens it
    assert 9.0 < g["gap_gllbsc_eV"] < 15.0       # exp 12.1 eV (fundamental)


@pytest.mark.skipif(not _mace_ok(), reason="mace/torch not importable")
def test_caf2_mace_phonons():
    from gainmat.host.phonons import compute_phonons
    r = compute_phonons(get_structure("CaF2"), device="cpu",
                        supercell_min_len=7.0, mesh=(8, 8, 8))
    print("CaF2 MACE phonons: max %.1f eff %.1f cm-1, Gamma %s" %
          (r["phonon_max_cm1"], r["phonon_eff_cm1"], r["gamma_cm1"]))
    assert r["min_freq_cm1"] > -10.0           # dynamically stable
    # exp: Raman T2g 322 cm^-1, TO 257 cm^-1; no NAC -> max is TO-like (~330-400)
    assert 250 < r["phonon_max_cm1"] < 550


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
