"""Fast checks of the glass forward model (site-ensemble lineshapes + rate equations)."""
import numpy as np
import pytest

from gainmat.forward.ions import get_ion
from gainmat.forward.spectra import random_sites, ensemble_spectra, nu_grid, spm_bkq
from gainmat.forward.glass import run, ir_loss_dB_m

SPM = {"O": dict(R0=2.3, B2=400., B4=350., B6=150., t2=5., t4=6., t6=6.)}


def test_kramers_and_stark_counts():
    m = get_ion("Er3+")
    E, W, par = m.stark(spm_bkq(random_sites(1, "O", 6, 2.27, 0.05, seed=3)[0], SPM))
    g = m.mult_index("4I15/2")
    e = np.sort(E[par == g])
    assert len(e) == 16 and np.allclose(e[0::2], e[1::2], atol=1e-6)   # 8 Kramers doublets


def test_integrated_cross_section_matches_judd_ofelt_independent_of_sites():
    # the integrated absorption must not depend on the crystal field (JO sum rule)
    om, n = (4.23, 1.04, 0.61), 1.45
    nu = nu_grid(1300, 1800, 1.0)
    ints = []
    for ang in (2.0, 20.0):
        r = ensemble_spectra("Er3+", random_sites(20, "O", 6.5, 2.27, 0.07, ang, seed=1), SPM, om, n,
                             pairs=[("4I13/2", "4I15/2")], nu=nu, T=5000.0)
        ints.append(np.trapezoid(r["sigma_abs"][("4I13/2", "4I15/2")], nu))
    assert ints[0] == pytest.approx(ints[1], rel=0.05)


def test_more_disorder_broadens():
    om, n, nu = (4.23, 1.04, 0.61), 1.45, nu_grid(1400, 1700, 1.0)
    from gainmat.forward.spectra import fwhm_nm
    w = []
    for ang, rs in ((3.0, 0.01), (20.0, 0.12)):
        r = ensemble_spectra("Er3+", random_sites(30, "O", 6.5, 2.27, rs, ang, seed=2), SPM, om, n,
                             pairs=[("4I13/2", "4I15/2")], nu=nu)
        w.append(fwhm_nm(nu, r["sigma_em"][("4I13/2", "4I15/2")]))
    assert w[1] > w[0]


def test_ir_edge_ordering():
    assert ir_loss_dB_m(2500, 1100) > 1 > ir_loss_dB_m(2500, 600)


def test_er_silica_c_band_and_no_2p7():
    r = run("Er3+", "silica", n_sites=15)
    assert any(w["lo_nm"] < 1535 < w["hi_nm"] for w in r["windows"])
    assert not any(w["lo_nm"] < 2750 < w["hi_nm"] for w in r["windows"])


def test_pr_oband_dead_in_silica_alive_in_zblan():
    hit = lambda r: any(w["lo_nm"] < 1310 < w["hi_nm"] for w in r["pumps"][1017.0]["windows"])
    assert not hit(run("Pr3+", "silica", n_sites=12))
    assert hit(run("Pr3+", "ZBLAN", n_sites=12))


def test_family_matching_prefers_exact():
    from gainmat.forward.validate import _family
    assert _family("aluminosilicate") == "aluminosilicate"
    assert _family("phosphosilicate") == "phosphosilicate"
    assert _family("silica") == "silica"
