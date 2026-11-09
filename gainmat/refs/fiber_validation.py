"""Measured properties of known rare-earth-doped fiber gain media, used to
validate the forward model (Stark-resolved spectra + rate equations).

Each VALIDATION entry is a dict with fixed keys; every numeric field is a
dict(value, src, verified[, unit, note]). value=None means no sourced number
was obtained in this build (always verified=False). Wavelengths in nm, cross
sections in cm^2, lifetimes in s. `verified=True` means the number was read in
the cited source or in a page/paper that quotes it (named in src) during the
2026-11-08 build. Do not add numbers without a source.

Blocked during build (permission timeouts): rp-photonics.com pages.
"""
from __future__ import annotations


def _v(value, src, verified=False, unit=None, note=None):
    d = dict(value=value, src=src, verified=verified)
    if unit:
        d["unit"] = unit
    if note:
        d["note"] = note
    return d


_NA = _v(None, "not found in this build")


def _entry(name, ion, host, transition, band_nm=_NA, sigma_abs_peak=_NA,
           lambda_abs_peak_nm=_NA, sigma_em_peak=_NA, lambda_em_peak_nm=_NA,
           em_fwhm_nm=_NA, tau_upper_s=_NA, extra=None, notes=""):
    d = dict(name=name, ion=ion, host=host, transition=transition,
             band_nm=band_nm, sigma_abs_peak_cm2=sigma_abs_peak,
             lambda_abs_peak_nm=lambda_abs_peak_nm,
             sigma_em_peak_cm2=sigma_em_peak,
             lambda_em_peak_nm=lambda_em_peak_nm, em_fwhm_nm=em_fwhm_nm,
             tau_upper_s=tau_upper_s, notes=notes)
    d["extra"] = extra or {}
    return d


_MELK = ("M.A. Melkumov, I.A. Bufetov, K.S. Kravtsov, A.V. Shubin, E.M. Dianov, "
         "Preprint No. 5, GPI RAS (2004), arXiv:1502.02885; short version "
         "Quantum Electron. 34, 843 (2004)")
_PASCH = ("R. Paschotta, J. Nilsson, P.R. Barber, J.E. Caplen, A.C. Tropper, "
          "D.C. Hanna, 'Lifetime quenching in Yb doped fibres' (Southampton "
          "eprint 78060)")
_KNALL = ("J. Knall, M. Engholm, T. Boilard, M. Bernier, M. Digonnet, "
          "'A radiation-balanced silica fiber amplifier', arXiv:2103.02698")
_GILES = ("C.R. Giles & E. Desurvire, J. Lightwave Technol. 9, 271 (1991), "
          "Table I (via Optiwave modelling.pdf reproduction)")
_TM_REV = ("Tm-silica parameter survey in arXiv:1511.02266 (mode-instability "
           "thresholds in Tm fiber amplifiers pumped at 790 nm)")
_LEGOUET = "J. Le Gouet et al., arXiv:2003.14188 (Ho-doped silica fiber)"
_PAYNE = ("D.N. Payne, R.I. Laming, D.W. Hewak, 'Towards pump-efficient 1.3 um "
          "fibre amplifiers', ORC Southampton (eprint 76975)")
_NTT335 = "US patent 5,351,335 (NTT, 1994)"
_MAL = ("Numerical modeling of 3.5 um dual-wavelength pumped Er mid-IR fiber "
        "lasers, arXiv:1606.07842, Tables II-III (lifetimes from Bogdanov "
        "[24])")

VALIDATION = [
    # ---------------------------------------------------------------- Yb
    _entry(
        "Yb:aluminosilicate", "Yb3+", "aluminosilicate", "2F5/2->2F7/2",
        band_nm=_v((976, 1200), _PASCH, True, "nm",
                   "'Gain is available between ~976nm and ~1200nm' "
                   "(germanosilicate/aluminosilicate fibres)"),
        sigma_abs_peak=_v(2.7e-20, _MELK, True, "cm^2", "2.7 pm^2"),
        lambda_abs_peak_nm=_v(976, _MELK, True, "nm"),
        sigma_em_peak=_v(3.0e-20, _MELK, True, "cm^2",
                         "3.0 pm^2 at the zero-line; ratio em/abs 1.1"),
        lambda_em_peak_nm=_v(976, _MELK, False, "nm",
                             "peak emission is at the zero-phonon line; "
                             "exact nm not separately stated"),
        em_fwhm_nm=_v(7.7, _MELK, True, "nm", "FWHM of main (976 nm) peak"),
        tau_upper_s=_v(0.83e-3, _MELK, True, "s"),
        extra=dict(
            tau_range_s=_v((0.72e-3, 1.35e-3), _PASCH, True, "s",
                           "six fibres, Table 1"),
            tau_rad_s_radbal=_v(765e-6, _KNALL, True, "s",
                                "2.52 wt% Yb, 2.0 wt% Al"),
            mean_fluor_wavelength_nm=_v(1003.9, _KNALL, True, "nm"),
            sigma_em_1030nm=_v(None, "not found numerically", False, "cm^2"),
        ),
        notes="Stark levels in site_refs.CF_CHECKS. Phosphosilicate variant: "
              "sigma_abs 1.4 pm^2 @974.5 nm, sigma_em 1.5 pm^2, tau 1.45 ms, "
              "FWHM 4.8 nm (Melkumov 2004, verified).",
    ),
    _entry(
        "Yb:phosphosilicate", "Yb3+", "phosphosilicate", "2F5/2->2F7/2",
        sigma_abs_peak=_v(1.4e-20, _MELK, True, "cm^2"),
        lambda_abs_peak_nm=_v(974.5, _MELK, True, "nm"),
        sigma_em_peak=_v(1.5e-20, _MELK, True, "cm^2"),
        em_fwhm_nm=_v(4.8, _MELK, True, "nm"),
        tau_upper_s=_v(1.45e-3, _MELK, True, "s"),
    ),
    # ---------------------------------------------------------------- Er 1.5
    _entry(
        "Er:aluminosilicate C/L band", "Er3+", "aluminosilicate",
        "4I13/2->4I15/2",
        band_nm=_v((1530, 1610), "commonly quoted C+L band edges", False, "nm"),
        sigma_abs_peak=_v(None, "Barnes et al., IEEE JQE 27, 1004 (1991) "
                          "(numbers not on abstract page)", False, "cm^2"),
        lambda_abs_peak_nm=_v(1530, "commonly quoted", False, "nm"),
        sigma_em_peak=_v(None, "Barnes et al. 1991 (not accessed)", False),
        tau_upper_s=_v(10e-3, _GILES, True, "s",
                       "metastable lifetime used for Al-silicate fibre"),
        extra=dict(
            Leff_silicate_nm=_v(34.8, "quoted in F. Huang et al., Sci. Rep. "
                                "4, 5053 (2014)", True, "nm",
                                "effective linewidth, silicate glass"),
            tau_mg_silica_s=_v((11.7e-3, 6.7e-3), "Dussardier, Blanc, "
                               "Peterka, arXiv:1201.1113", True, "s",
                               "fibres A/B (0.1 / 1 mol/l Mg); B FWHM 44 nm"),
        ),
        notes="Barnes 1991: Fuchtbauer-Ladenburg vs saturation methods give "
              "significantly different cross sections, esp. their ratio.",
    ),
    _entry(
        "Er:tellurite L band", "Er3+", "tellurite", "4I13/2->4I15/2",
        band_nm=_v((1561, 1611), "A. Mori et al., J. Lightwave Technol. 20, "
                   "822 (2002)", True, "nm",
                   "25.3 dB gain, 0.6 dB excursion, NF<6 dB"),
        notes="Broadband Er tellurite cross sections/lifetime not obtained.",
    ),
    _entry(
        "Er:phosphate", "Er3+", "phosphate", "4I13/2->4I15/2",
        extra=dict(Leff_nm=_v(46.0, "quoted in F. Huang et al., Sci. Rep. 4, "
                              "5053 (2014)", True, "nm",
                              "effective linewidth, phosphate glass")),
    ),
    _entry(
        "Er:ZBLAN 1.5um", "Er3+", "ZBLAN", "4I13/2->4I15/2",
        tau_upper_s=_v(6.44e-3, "F. Huang et al., Sci. Rep. 4, 5053 (2014), "
                       "Table 2", True, "s", "RADIATIVE (Judd-Ofelt) lifetime"),
        extra=dict(JO_Omega_2_4_6=_v((3.27e-20, 1.3e-20, 1.75e-20),
                   "Huang et al. 2014", True, "cm^2")),
    ),
    # ---------------------------------------------------------------- Tm
    _entry(
        "Tm:silica 2um", "Tm3+", "aluminosilicate", "3F4->3H6",
        band_nm=_v((1800, 2050), "task spec; laser 1886-1929 nm in Cajzl et "
                   "al. Fibers 6, 66 (2018)", False, "nm"),
        sigma_em_peak=_v((4.0e-21, 6.0e-21), _TM_REV, True, "cm^2",
                         "literature peaks: Jackson&King 6e-25 m^2, "
                         "Walsh&Barnes 4.6e-25, Turri 4e-25 (300 K)"),
        sigma_abs_peak=_v(None, "Agger & Povlsen, Opt. Express 14, 50 (2006) "
                          "-- numbers not in abstract", False, "cm^2",
                          "peak em/abs ratio 'significantly below 1'"),
        tau_upper_s=_v(650e-6, "Agger & Povlsen 2006 via " + _TM_REV, True,
                       "s", "others: 335, 430, 635 us; 756 us in 11.7 mol% "
                            "Al2O3 (Cajzl 2018)"),
        extra=dict(
            tau_rad_s=_v((4.56e-3, 6.0e-3), _TM_REV, True, "s"),
            pump_sigma_abs_790=_v((4.5e-21, 1.0e-20), _TM_REV, True, "cm^2"),
        ),
    ),
    _entry(
        "Tm:ZBLAN S band", "Tm3+", "ZBLAN", "3H4->3F4",
        band_nm=_v((1480, 1510), "S. Aozasa et al., J. Lightwave Technol. 25, "
                   "2108 (2007)", True, "nm",
                   ">26 dB gain; host not stated in abstract (NTT TDFA, "
                   "fluoride assumed)"),
        tau_upper_s=_v(1.5e-3, "US patent 6,667,257 text ('observed in low "
                       "phonon fluorozirconate')", True, "s", "approximate"),
        extra=dict(
            tau_3H4_silica_s=_v(14e-6, "Dussardier et al. arXiv:1201.1113 "
                                "(citing Faure 2007)", True, "s",
                                "rises to 50 us at 17.4 mol% AlO3/2"),
            tau_3H4_TmAl_silica_s=_v(20e-6, "US 6,667,257", True, "s"),
            gap_3H4_3H5_cm1=_v(4400, "US 6,667,257", True, "cm^-1"),
        ),
        notes="Quenching in silica: 3H4 tau ~14-58 us vs ~1.5 ms in ZBLAN.",
    ),
    # ---------------------------------------------------------------- Pr
    _entry(
        "Pr:ZBLAN O band", "Pr3+", "ZBLAN", "1G4->3H5",
        band_nm=_v((1250, 1340), "US patent 5,309,452", True, "nm",
                   "stated operating range; ESA 1G4->1D2 limits >1310 nm"),
        sigma_em_peak=_v(3.48e-21, "Y. Ohishi et al., Opt. Lett. 16, 1747 (1991), p. 1747 and Table 1 "
                         "(the earlier 3.5e-22 from US patent 5,351,335 was a factor-10 slip)", True, "cm^2"),
        lambda_em_peak_nm=_v(1322, _NTT335, True, "nm",
                             "US 5,309,452 shows peak ~1340 nm"),
        tau_upper_s=_v(110e-6, _PAYNE, True, "s"),
        extra=dict(
            eta_ZBLAN=_v(0.04, _PAYNE, True, None,
                         "~3% per Ohishi et al. Electron. Lett. 27 (1991) "
                         "via US 5,351,335"),
            tau_GLS_s=_v(290e-6, _PAYNE, True, "s"),
            eta_GLS=_v(0.53, _PAYNE, True),
            tau_CdFCl_s=_v(325e-6, _PAYNE, True, "s"),
            eta_CdFCl=_v(0.12, _PAYNE, True),
            tau_silica=_v(None, "severely quenched (US 5,309,452), no number",
                          False),
            tau_fluoroindate=_v(None, "not found", False),
            branching_1G4_3H5=_v(0.64, "US 5,309,452", True),
            net_gain_dB=_v(30.1, "Y. Ohishi et al., IEEE PTL 3, 715 (1991)",
                           True, "dB", "at 1309 nm, pump 1017 nm, Pr<1000 ppm"),
        ),
    ),
    # ---------------------------------------------------------------- mid-IR
    _entry(
        "Er:ZBLAN 2.8um", "Er3+", "ZBLAN", "4I11/2->4I13/2",
        band_nm=_v(None, "lasing 2710-2783 nm reported; full band not found",
                   False),
        sigma_em_peak=_v(7e-21, "F. Huang et al., Sci. Rep. 4, 5053 (2014)",
                         True, "cm^2", "'above 7e-21' for all fluoride glasses"),
        lambda_em_peak_nm=_v(2708, "Huang et al. 2014", True, "nm"),
        tau_upper_s=_v(6.9e-3, "Li, Zhang, Wang, Luo (Xiamen U.) 3-um fiber "
                       "laser review, core.ac.uk 343508871", True, "s",
                       "7.9 ms in " + _MAL),
        extra=dict(tau_lower_4I13_2_s=_v(9e-3, "Li et al. review", True, "s",
                                         "9.9 ms in arXiv:1606.07842")),
        notes="Self-terminating: tau(4I13/2) > tau(4I11/2).",
    ),
    _entry(
        "Ho:ZBLAN 2.9um", "Ho3+", "ZBLAN", "5I6->5I7",
        tau_upper_s=_v(3.5e-3, "Li et al. Xiamen review (core.ac.uk "
                       "343508871); also 3.5 ms in arXiv:1902.06945", True,
                       "s"),
        extra=dict(tau_lower_5I7_s=_v(12e-3, "Li et al. review", True, "s")),
    ),
    _entry(
        "Ho:silica 2.1um", "Ho3+", "silica", "5I7->5I8",
        band_nm=_v((2019, 2163), "S.D. Jackson & Y. Li, Electron. Lett. 40 "
                   "(2004) doi:10.1049/el:20046979", True, "nm",
                   "laser tuning range"),
        sigma_abs_peak=_v(4.3e-21, _LEGOUET, True, "cm^2",
                          "(4.3+-0.2)e-25 m^2; lit. 2.9-8.5e-25 m^2"),
        lambda_abs_peak_nm=_v(1946, _LEGOUET, True, "nm"),
        sigma_em_peak=_v(6.0e-21, _LEGOUET, True, "cm^2",
                         "direct gain measurement; em/abs peak ratio 1.40"),
        tau_upper_s=_v(0.7e-3, _LEGOUET, True, "s"),
    ),
    _entry(
        "Dy:ZBLAN 3um", "Dy3+", "ZBLAN", "6H13/2->6H15/2",
        band_nm=_v((2807, 3380), "M.R. Majewski et al., arXiv:1802.07413",
                   True, "nm", "CW tuning range, 1.7 um pumping"),
        lambda_em_peak_nm=_v(2950, "arXiv:1802.07413 (free-running gain peak)",
                             True, "nm"),
        tau_upper_s=_v(650e-6, "Gomes et al. via arXiv:1902.06945", True, "s"),
        extra=dict(
            sigma_abs_2825=_v(3.9e-21, "arXiv:1902.06945 Table 1 (model)",
                              True, "cm^2", "units printed as m^-3 (typo)"),
            sigma_em_2825=_v(2.6e-21, "arXiv:1902.06945 Table 1", True, "cm^2"),
            sigma_abs_3100=_v(2.8e-22, "arXiv:1902.06945 Table 1", True,
                              "cm^2"),
            sigma_em_3100=_v(6.8e-22, "arXiv:1902.06945 Table 1", True,
                             "cm^2"),
        ),
        notes="Dy:fluoroindate: no numbers found.",
    ),
    _entry(
        "Er:ZBLAN 3.5um", "Er3+", "ZBLAN", "4F9/2->4I9/2",
        band_nm=_v((3330, 3780), "Henderson-Sapir et al., Opt. Lett. 41, 1676 "
                   "(2016)", True, "nm", "450 nm tuning, longest 3.78 um"),
        sigma_em_peak=_v(1.2e-21, _MAL, True, "cm^2", "12e-26 m^2 @3470 nm"),
        lambda_em_peak_nm=_v(3470, _MAL, False, "nm",
                             "laser wavelength in model, not peak per se"),
        tau_upper_s=_v(177e-6, _MAL, True, "s"),
        extra=dict(
            tau_lower_4I9_2_s=_v(8e-6, _MAL, True, "s",
                                 "~5 us in Henderson-Sapir 2016"),
            tau_4I11_2_s=_v(7.9e-3, _MAL, True, "s"),
        ),
    ),
    _entry(
        "Er:fluoroindate 3.3um", "Er3+", "fluoroindate", "4F9/2->4I9/2",
        band_nm=_v((3100, 3850), "H. He et al., Opt. Lett. 46, 1057 (2021)",
                   True, "nm", "dual-pumped 976/1973 nm emission"),
        lambda_em_peak_nm=_v(3300, "He et al. 2021", True, "nm"),
        extra=dict(
            tau_4I11_2_fluoroindate_s=_v(7.09e-3, "J. She et al., J. Rare "
                                         "Earths (2021) doi:10.1016/j.jre."
                                         "2021.05.011", True, "s",
                                         "2 mol% Er; 6.49 ms ZBLAN"),
            sigma_em_2700_fluoroindate=_v(6.95e-21, "She et al. 2021", True,
                                          "cm^2"),
        ),
    ),
]
