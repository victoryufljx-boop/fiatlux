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
    # with Stark-edge multiphonon gaps the 1G4 lifetime depends on the crystal-field strength,
    # so use the production setup (measured site geometry + calibrated SPM scale)
    from gainmat.forward.glass import load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ
    load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ)
    hit = lambda r: any(w["lo_nm"] < 1310 < w["hi_nm"] for w in r["pumps"][1017.0]["windows"])
    assert not hit(run("Pr3+", "silica", n_sites=12))
    assert hit(run("Pr3+", "ZBLAN", n_sites=12))


def test_family_matching_prefers_exact():
    from gainmat.forward.validate import _family
    assert _family("aluminosilicate") == "aluminosilicate"
    assert _family("phosphosilicate") == "phosphosilicate"
    assert _family("silica") == "silica"


def test_er_dual_pump_opens_3p5um_band_in_zblan():
    # Henderson-Sapir et al., Opt. Lett. 39, 493 (2014): 985 + 1973 nm pumping, lasing 3.47-3.60 um
    from gainmat.forward.glass import load_refs, apply_site_refs, calibrate_spm, CARNALL_ER_LAF3_BKQ
    load_refs(); apply_site_refs(); calibrate_spm(370.0, CARNALL_ER_LAF3_BKQ)
    r = run("Er3+", "ZBLAN", n_sites=12, pumps=[980.0, (985.0, 1973.0)])
    in_band = lambda d: any(w["lo_nm"] < 3500 < w["hi_nm"] for w in d["windows"])
    assert not in_band(r["pumps"][980.0])
    assert in_band(r["pumps"]["985+1973"])


def test_stark_resolved_intensities_keep_judd_ofelt_sums():
    # odd-rank crystal-field intensities redistribute strength among Stark components but the
    # multiplet sum must stay the Judd-Ofelt value of that site
    from gainmat.forward.spectra import intensity_B, site_omegas, aniso_S
    from gainmat.forward.glass import spm_for
    m = get_ion("Yb3+")
    lig = random_sites(1, "O", 6.0, 2.28, 0.1, 10.0, seed=3)[0]
    E, W, par = m.stark(spm_bkq(lig, {"O": spm_for("Yb3+", "O")}))
    Uw = {lam: [W.conj().T @ Uq @ W for Uq in m.Uq[lam]] for lam in (2, 4, 6)}
    B = intensity_B(lig, {2: 0.7, 4: 1.3, 6: -0.5}); om = site_omegas(B)
    Sa = aniso_S(Uw, B, {2: 1, 4: 1, 6: 1})
    Si = sum(om[lam] * sum(np.abs(X) ** 2 for X in Uw[lam]) for lam in (2, 4, 6))
    a = np.where(par == m.mult_index("2F5/2"))[0]; b = np.where(par == m.mult_index("2F7/2"))[0]
    assert abs(Sa[np.ix_(a, b)].sum() / Si[np.ix_(a, b)].sum() - 1) < 0.02
    assert np.allclose(Sa, Sa.T, atol=1e-12 * Sa.max())


def test_measured_absorption_mode_mccumber_and_rate_scale():
    nu = np.linspace(9500, 11200, 400); lam = 1e7 / nu
    sig = 2e-20 * np.exp(-0.5 * ((nu - 10250) / 60) ** 2) + 5e-21 * np.exp(-0.5 * ((nu - 10600) / 200) ** 2)
    r = run("Yb3+", "germanosilicate", n_sites=8, pumps=[], overrides=dict(measured_absorption={("2F5/2", "2F7/2"): (lam, sig)}))
    s = r["_spec"]; p = ("2F5/2", "2F7/2")
    assert "2F5/2->2F7/2" in s["measured_bands"] and s["A_scale"][p] > 0
    i = np.argmin(abs(s["nu"] - 10250))
    eps = s["measured_bands"]["2F5/2->2F7/2"]["eps_cm1"]
    assert np.isclose(s["sigma_em"][p][i] / s["sigma_abs"][p][i], np.exp((eps - s["nu"][i]) / (0.6950348 * 295)), rtol=1e-6)


def test_judd_ofelt_fit_round_trip():
    from gainmat.forward.jo_fit import fit_omegas
    pairs = [("3F4", "3H6"), ("3H5", "3H6"), ("3H4", "3H6"), ("3F3", "3H6")]
    om_true = (6.23, 1.91, 1.36)
    r = run("Tm3+", "silica", n_sites=8, pumps=[], seed=5, overrides=dict(omegas=om_true))
    s = r["_spec"]; lam = 1e7 / s["nu"]
    om, info = fit_omegas("Tm3+", "silica", {p: (lam, s["sigma_abs"][p]) for p in pairs}, n_sites=8)
    assert info["rank"] == 3
    assert np.allclose(om, om_true, rtol=0.08)


def test_level_anchor_matches_measured_levels():
    """Measured-level mode: the site-mean Stark levels equal the anchors, scatter is kept."""
    from gainmat.forward.glass import spm_for
    from gainmat.forward.spectra import ensemble_spectra, random_sites
    sites = random_sites(20, "O", 6.0, 2.26, 0.08, 13, seed=3)
    anc = {"2F7/2": dict(levels=[0, 425, 530, 865]), "2F5/2": dict(levels=[0, 630], E0=10235.0)}
    sp = ensemble_spectra("Yb3+", sites, {"O": spm_for("Yb3+", "O")}, (1.0, 1.0, 1.0), 1.45,
                          pairs=[("2F5/2", "2F7/2")], level_anchor=anc)
    lo7, hi7 = sp["edges"]["2F7/2"]; lo5, hi5 = sp["edges"]["2F5/2"]
    assert abs((hi7 - lo7) - 865) < 5.0      # per-site map: mean of mapped != map of mean (few cm^-1)
    assert abs((lo5 - lo7) - 10235) < 5.0
    assert sp["split"]["2F7/2"] > 0


def test_predict_returns_bands_uncertainty_and_provenance():
    from gainmat.forward.predict import predict
    p = predict("Er3+", "aluminosilicate", n_sites=10)
    b = p["bands"][0]
    assert b["transition"] == "4I13/2->4I15/2" and 1500 < b["abs_peak_nm"] < 1560
    assert b["uncertainty"]["sigma_peak_log10"] is not None
    assert p["provenance"]["stark"].startswith("measured levels")
    assert any(l["level"] == "4I13/2" and 5 < l["tau_rad_ms"] < 20 for l in p["levels"])


def test_composition_front_end_families():
    from gainmat.forward.composition import family_for
    assert family_for({"SiO2": 92, "Al2O3": 6, "Yb2O3": 2})[0] == "aluminosilicate"
    assert family_for("70TeO2-15K2O-15Nb2O5")[0] == "tellurite"
    assert family_for("53ZrF4-20BaF2-4LaF3-3AlF3-20NaF")[0] == "ZBLAN"
    assert family_for("60GeO2-15Ga2O3-10BaO-15Na2O")[0] == "germanate"
    assert family_for("40P2O5-20Al2O3-35Na2O-5K2O")[0] == "phosphate"
    assert family_for({"SiO2": 99.0, "GeO2": 1.0})[0] == "germanosilicate"
