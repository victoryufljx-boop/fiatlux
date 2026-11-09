"""Fibre spectroscopy reference values, read from the papers themselves (page / table cited).

[DS90]  E. Desurvire, J. R. Simpson, "Evaluation of 4I15/2 and 4I13/2 Stark-level energies in
        erbium-doped aluminosilicate glass fibers", Opt. Lett. 15, 547-549 (1990).
[OH91]  Y. Ohishi, T. Kanamori, T. Kitagawa, S. Takahashi, E. Snitzer, G. H. Sigel Jr.,
        "Pr3+-doped fluoride fiber amplifier operating at 1.31 um", Opt. Lett. 16, 1747-1749 (1991).
[HS14]  O. Henderson-Sapir, J. Munch, D. J. Ottaway, "Mid-infrared fiber lasers at and beyond
        3.5 um using dual-wavelength pumping", Opt. Lett. 39, 493-496 (2014).
[TAO15] M. Tao, F. Wang, Z. Wang, P. Yang, G. Feng, X. Ye, "Theoretical analysis and comparison of
        Tm-doped fiber lasers with different pump bands", Proc. SPIE 9255, 92551O (2015).
[FA07]  B. Faure, W. Blanc, B. Dussardier, G. Monnom, "Improvement of the Tm3+:3H4 level lifetime
        in silica optical fibres by lowering the local phonon energy", J. Non-Cryst. Solids 353,
        2767 (2007); read from the preprint arXiv:1003.0166.

Only numbers are stored here; the PDFs are not redistributed.
"""

# ------------------------------------------------------------------------------------------ DS90
DESURVIRE1990 = dict(
    src="[DS90]",
    host="Al2O3:GeO2:SiO2 MCVD fibre, ~7.5 mol% GeO2 + ~7.5 mol% Al2O3 (p. 547)",
    N_Er_cm3=2.8e18, tau_4I13_2_s=10e-3,                      # p. 547
    # Table 1 (p. 549), "This paper" column, cm^-1 relative to the lowest level of each multiplet;
    # ranges are as printed.  Only 5 of the 8 4I15/2 and 6 of the 7 4I13/2 Kramers doublets were
    # identified, so the total splittings are LOWER bounds.
    stark_4I15_2_cm1=[0, (26, 34), (51, 59), (125, 133), (201, 268)],
    stark_4I13_2_cm1=[(0, 8), (96, 104), (163, 171), (222, 230)],
    E_4I13_2_lowest_cm1=(6540, 6548),                          # Fig. 4/5 and p. 547
    total_splitting_cm1={"4I15/2": 268, "4I13/2": 230},        # p. 549 text
    # Table 1, other Er:silica glasses (refs 10, 12, 13 of DS90): top of 4I15/2 and 4I13/2
    other_silicates_total_splitting_cm1={"4I15/2": (395, 500), "4I13/2": (305, 357),
                                         "E_4I13_2_lowest": (6510, 6540)},
    homogeneous_fwhm_295K_cm1=49, inhomogeneous_fwhm_cm1=49,   # p. 549 (ref 11 of DS90)
    homogeneous_fwhm_4K_nm=0.1, zpl_um=1.531,                  # p. 547
)

# ------------------------------------------------------------------------------------------ OH91
OHISHI1991 = dict(
    src="[OH91]",
    host="ZBLAN (ZrF4-BaF2-LaF3-AlF3-NaF) for spectroscopy; fibre ZrF4-BaF2-LaF3-YF3-AlF3-LiF-NaF",
    Pr_wppm=500,                                               # p. 1748
    Omega_1e20_cm2=(1.431, 4.221, 4.872),                      # p. 1747 (Omega_6 printed '4.872 x 10-2', a slip)
    tau_meas_1G4_s=110e-6, tau_rad_1G4_s=3.24e-3, eta_1G4=0.034,   # pp. 1747-1748
    sigma_em_1G4_3H5_cm2=3.48e-21,                             # p. 1747 and Table 1 (FL from mean lambda, FWHM)
    sigma_em_1G4_3H6_cm2=1.96e-21,                             # Table 1
    fluorescence_peak_um=1.322, gain_peak_um=1.31,             # pp. 1747-1748
    hw_max_cm1=580,                                            # p. 1747
    # Table 1: 1G4 -> (lower): (frequency cm^-1, A s^-1, branching ratio)
    transitions_1G4={"3F4": (2913, 11.8, 0.038), "3F3": (3344, 2.7, 0.009), "3F2": (4696, 1.8, 0.006),
                     "3H6": (5466, 76.3, 0.248), "3H5": (7564, 197.0, 0.639), "3H4": (9833, 18.5, 0.060)},
    esa="1G4 -> 1D2, oscillator strength 1.27e-6 vs 2.32e-6 for 1G4 -> 3H5; ESA peaks on the long-"
        "wavelength side of the 1.3 um band, shifting the gain peak to 1.31 um (p. 1748)",
    amplifier=dict(pump_nm=1017, pump_mW=180, length_m=8, core_um=4, dn_pct=0.6,
                   internal_gain_dB=7.1, net_gain_dB=5.2, net_gain_band_um=(1.28, 1.335),
                   gain_coeff_dB_per_mW=0.03),                 # p. 1748
)

# ------------------------------------------------------------------------------------------ HS14
HENDERSON_SAPIR2014 = dict(
    src="[HS14]",
    host="Er:ZBLAN (IR-Photonics), 1 mol% Er, 18 cm, core 10 um, NA 0.15 (p. 494)",
    scheme="dual-wavelength pumping: P1 985 nm 4I15/2->4I11/2 fills the long-lived 4I11/2; P2 1973 nm "
           "4I11/2->4F9/2 (ESA); laser 4F9/2->4I9/2; 4I9/2 decays to 4I11/2 by multiphonon (p. 494)",
    pump1_nm=985.0, pump2_nm=1973.0,
    P1_incident_mW=194, P1_absorbed_mW=140,                    # p. 495
    P2_incident_max_mW=2400,                                   # Fig. 3
    tau_4I11_2_s=6.8e-3,                                       # p. 494 (citing Pollnau 1997)
    laser_lines_um=(3.472, 3.538, 3.604),                      # p. 495
    fluorescence_band_um=(3.2, 3.9),                           # Fig. 5 and p. 495 (tuning range claim)
    output_cw_mW=260, slope_eff_incident_P2=0.254, slope_eff_absorbed_P2=0.523,
    note="No cross sections are given; ESA cross sections must come from elsewhere "
         "(e.g. the Malouf et al. model, arXiv:1606.07842, already cited in fiber_validation.py).",
)

# ----------------------------------------------------------------------------------------- TAO15
TAO2015 = dict(
    src="[TAO15] Table 1 (modelling parameters, taken from its refs 8, 10-12, incl. Jackson & King 1999)",
    host="Tm:silica LMA fibre, 25 um core",
    N_Tm_cm3=2.91e20,
    tau_3F4_s=334.7e-6, tau_3H4_s=14.2e-6, tau_3H5_s=7e-9,
    sigma_em_signal_cm2=2.5e-21, sigma_abs_signal_cm2=1.0e-22,
    sigma_abs_pump_cm2={793: 5e-21, 1600: 1.52e-21, 1900: 2.33e-22},
    k_CR_3101_cm3_s=1.41e-19,
    note="Model inputs, not new measurements; useful as a second Tm:silica 3F4 lifetime at high "
         "concentration (0.335 ms at 2.9e20 cm^-3 vs 0.42 ms at 1.1e20 in Walsh & Barnes 2004).",
)

# ------------------------------------------------------------------------------------------ FA07
FAURE2007 = dict(
    src="[FA07]",
    Tm_conc_mol_ppm="30-200 (Table 1) -- low enough to exclude Tm-Tm transfer when 3H4 is pumped (p. 8)",
    tau_rad_3H4_s=650e-6,                                      # p. 9 (Digonnet JO value)
    # 1/e lifetimes of 3H4 (us): pure silica from Jackson & King 1999; others measured (Fig. 5, text)
    tau_3H4_us={"SiO2": 14, "Al 17.4 mol% AlO3/2": 50, "Ge 20 mol% GeO2": 28, "P 8 mol% PO5/2": 9},
    tau_vs_Al="linear in AlO3/2 content from 14 us (0) to 50 us (17.4 mol%) (p. 8)",
    Al_series_mol_pct_AlO3_2=(5.6, 6.4, 7.8, 10.4, 13.0, 17.4),   # Table 1 (lifetimes only in Fig. 5)
    # Table 2: 3H6->3H4 GSA peak (~ lowest 3H4 Stark level), 3H5 10%-edge (~ highest 3H5 Stark level)
    gap_3H4_3H5_cm1={"Al": (12710, 9090, 3620), "P": (12670, 8940, 3730), "Ge": (12760, 8990, 3770),
                     "SiO2 (fitted to 14 us)": (None, None, 3693)},
    gap_error_cm1=40,
    # Table 3: van Dijk-type law W = W0 exp[-alpha (dE - 2 Ep)] for pure oxides (from refs 28, 29)
    vdS_params={"SiO2": dict(Ep=1100, W0=7.8e7, alpha=4.7e-3),
                "AlO3/2": dict(Ep=870, W0=6.2e7, alpha=4.7e-3, note="averages, not measured"),
                "GeO2": dict(Ep=900, W0=6.1e7, alpha=4.6e-3),
                "PO5/2": dict(Ep=1320, W0=7.6e7, alpha=4.7e-3)},
    decays="non-exponential: Al-poor sites ~10 us, Al-rich sites ~100 us (p. 11)",
)


def _w(tau_us, tau_rad_us=650.0):
    return 1e6 / tau_us - 1e6 / tau_rad_us


# Multiphonon rates of Tm 3H4 -> 3H5 implied by [FA07] (s^-1), with the measured edge gaps.
FAURE2007_W = {
    "SiO2": dict(W=_w(14), dE=3693),
    "Al 17.4": dict(W=_w(50), dE=3620),
    "Ge 20": dict(W=_w(28), dE=3770),
    "P 8": dict(W=_w(9), dE=3730),
}
