"""Tests for gainmat.cluster (embedding geometry/electrostatics + fast CASSCF).

The CASSCF tests take a few minutes on 2 cores; deselect with  -m "not slow".
"""
import numpy as np
import pytest

from gainmat.cluster import embedding as emb
from gainmat.cluster.casci_lf import add_spin_orbit, HARTREE_CM1


@pytest.fixture(scope="module")
def ce_cluster():
    return emb.make_cluster("LiYF4", "Y", "Ce", 3, r_qm=2.6, r_cap=4.0, r_pc=12.0)


@pytest.fixture(scope="module")
def ti_cluster():
    return emb.make_cluster("Al2O3", "Al", "Ti", 3, r_qm=2.2, r_cap=3.6, r_pc=12.0)


def test_ce_geometry(ce_cluster):
    d = np.sort(np.linalg.norm(ce_cluster.qm_coords[1:], axis=1))
    assert ce_cluster.qm_symbols.count("F") == 8
    np.testing.assert_allclose(d[:4], 2.246, atol=0.01)     # Y-F short
    np.testing.assert_allclose(d[4:], 2.293, atol=0.01)     # Y-F long
    assert ce_cluster.qm_charge == -5


def test_ti_geometry(ti_cluster):
    d = np.sort(np.linalg.norm(ti_cluster.qm_coords[1:], axis=1))
    assert ti_cluster.qm_symbols.count("O") == 6
    np.testing.assert_allclose(d[:3], 1.855, atol=0.01)
    np.testing.assert_allclose(d[3:], 1.971, atol=0.01)
    assert ti_cluster.qm_charge == -9
    assert len(ti_cluster.cap_symbols) > 0 and set(ti_cluster.cap_symbols) == {"Al"}


@pytest.mark.parametrize("fx", ["ce_cluster", "ti_cluster"])
def test_neutral_and_madelung(fx, request):
    cl = request.getfixturevalue(fx)
    assert abs(cl.total_charge()) < 1e-8
    # potential at the substituted site from the finite array reproduces Ewald within 0.05 V
    assert abs(cl.madelung_cluster_V - cl.madelung_ewald_V) < 0.05
    assert cl.fit_rms_V < 0.1


def test_ewald_alpha_independent():
    cell, sy, frac, q = emb.build_crystal("Al2O3")
    i = sy.index("Al")
    v1 = emb.site_potential(cell, frac, q, i, alpha=0.2, rcut=18.0)
    v2 = emb.site_potential(cell, frac, q, i, alpha=0.35, rcut=12.0)
    assert abs(v1 - v2) < 1e-3
    assert -40 < v1 < -33          # Al site in corundum ~ -36.6 V


def test_ewald_vs_nacl_madelung():
    """Rock-salt check: Madelung constant 1.747565."""
    a = 5.64
    cell = np.eye(3) * a
    frac = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5],
                     [.5, 0, 0], [0, .5, 0], [0, 0, .5], [.5, .5, .5]])
    q = np.array([1, 1, 1, 1, -1, -1, -1, -1.])
    v = emb.site_potential(cell, frac, q, 0)
    assert abs(-v * (a / 2) / emb.EV_ANG - 1.747565) < 1e-4


def test_spin_orbit_free_ion():
    """Free d1 / f1 ion: SO splits into j=l-1/2 and j=l+1/2 with gaps zeta*(2l+1)/2."""
    for l, z in [(2, 1000.0), (3, 600.0)]:
        n = 2 * l + 1
        w, _ = add_spin_orbit(np.zeros((n, n)), [(l, 0)], [z])
        w = np.round(w, 6)
        assert np.isclose(w.max() - w.min(), z * (2 * l + 1) / 2)
        assert (np.isclose(w, w.min()).sum(), np.isclose(w, w.max()).sum()) == (2 * l, 2 * l + 2)


@pytest.mark.slow
def test_ce_fast_casscf():
    from gainmat.cluster import ce_5d
    o = ce_5d.run("fast", do_nevpt2=False, log=lambda *a: None)
    r = o["e_cas"]
    sc = np.array(r["scalar_cm1"]); dch = np.array(r["scalar_5d_char"])
    f = sc[dch < 0.5]; d = sc[dch >= 0.5]
    assert len(f) == 7 and len(d) == 5
    assert np.ptp(f) < 0.5 * 8065.5                  # 4f manifold < 0.5 eV
    assert d.min() - f.max() > 2.0 * 8065.5          # 5d several eV above
    assert len(r["so_4f_cm1"]) == 7 and len(r["so_5d_cm1"]) == 5


@pytest.mark.slow
def test_ti_fast_casscf():
    from gainmat.cluster import ti_3d
    o = ti_3d.run("fast", do_nevpt2=False, log=lambda *a: None)
    e = np.array(o["e_cas_cm1"])
    assert len(e) == 5
    assert e[2] - e[0] < 2000                        # t2g-like triplet
    assert 12000 < e[3] < 35000                      # 10Dq in the visible


# ---------------------------------------------------------------- improvements (AIMP, basis, relax)
def test_aimp_cluster_layout():
    cl = emb.make_cluster("LiYF4", "Y", "Ce", 3, r_qm=2.6, r_cap=5.3, r_pc=12.0, cap_mode="aimp")
    assert cl.aimp_symbols == ["Y"] * 8                  # all 8 nearest Y3+ are AIMP ions
    assert abs(cl.total_charge()) < 1e-8
    assert not any(s == "Y" for s in cl.cap_symbols)


def test_aimp_is_repulsive_near_ion():
    """Frozen Y3+ potential: net +3 attraction far away, Pauli repulsion close in."""
    from pyscf import gto
    x = np.array([0.0, 0.0, 0.0])
    m = gto.M(atom=[["H", (1.0, 0, 0)]], basis={"H": [[0, [1.0, 1.0]]]}, spin=1, verbose=0)
    h, _ = emb.aimp_hcore(m, ["Y"], [x])
    with m.with_rinv_origin(x / 0.52917721092):
        v_pc = -3 * m.intor("int1e_rinv")
    assert h[0, 0] > v_pc[0, 0] + 0.05
    m2 = gto.M(atom=[["H", (8.0, 0, 0)]], basis={"H": [[0, [1.0, 1.0]]]}, spin=1, verbose=0)
    h2, _ = emb.aimp_hcore(m2, ["Y"], [x])
    with m2.with_rinv_origin(x):
        assert abs(h2[0, 0] + 3 * m2.intor("int1e_rinv")[0, 0]) < 1e-3


def test_bundled_ce_basis():
    bas, ecp = emb.load_bundled("ce_def2-svp.nw", "Ce")
    assert ecp[0] == 28 and len(bas) > 10


def test_relax_modes():
    from gainmat.cluster import relax
    ti = emb.make_cluster("Al2O3", "Al", "Ti", 3, r_qm=2.2, r_cap=3.6, r_pc=12.0)
    m = relax.eg_theta_mode(ti)
    assert np.isclose(np.linalg.norm(m), 1.0) and abs(m[0]).sum() == 0
    ce = emb.make_cluster("LiYF4", "Y", "Ce", 3, r_qm=2.6, r_cap=4.0, r_pc=12.0)
    assert len(relax.breathing_modes(ce, split_sets=True)) == 2
    q, e, _ = relax._quad_fit(np.array([[-1.], [0.], [1.], [2.]]), [1.0, 0.0, 1.0, 4.0])
    assert abs(q[0]) < 1e-8 and abs(e) < 1e-8


@pytest.mark.slow
def test_ti_jt_relax():
    from gainmat.cluster import ti_3d
    o = ti_3d.run("fast+", log=lambda *a: None)
    assert 0 < o["stokes_cm1"] < 15000
    assert o["em_cm1"] < o["abs_eg_mean_cm1"]
