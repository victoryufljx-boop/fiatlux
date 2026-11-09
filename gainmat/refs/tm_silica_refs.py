"""Numerical data extracted from three papers on Tm3+-doped silica fibre.

Full citations
--------------
JACKSON1999 : S. D. Jackson and T. A. King, "Theoretical modeling of Tm-doped
              silica fiber lasers," J. Lightwave Technol. 17(5), 948-956 (1999).
JACKSON2003 : S. D. Jackson and S. Mossman, "Efficiency dependence on the Tm3+
              and Al3+ concentrations for Tm3+-doped silica double-clad fiber
              lasers," Appl. Opt. 42(15), 2702-2707 (2003).
AGGER2006   : S. D. Agger and J. H. Povlsen, "Emission and absorption cross
              section of thulium doped silica fibers," Opt. Express 14(1),
              50-57 (2006).

Page numbers ("p.") are the journal page numbers printed on the paper.

READ THIS FIRST: level-label convention
----------------------------------------
Jackson & King 1999 (p. 949, "we have arbitrarily chosen the zero-spin orbit
convention") and Jackson & Mossman 2003 (Fig. 1, text) use the OLD labels in
which 3H4 and 3F4 are swapped relative to modern usage (Agger 2006 uses the
modern labels):

    energy (approx.)       OLD label (Jackson)   MODERN label (Agger, this file)
    ~5600 cm^-1  (1.8 um)  "3H4"                 3F4   <- 2 um upper laser level
    ~8300 cm^-1  (1.2 um)  "3H5"                 3H5
    ~12600 cm^-1 (790 nm)  "3F4"                 3H4   <- 790 nm pump level

Every key/comment in THIS FILE uses MODERN labels.  Where a quote from a
Jackson paper is given, the original old label is shown as e.g. 'old "3H4"'.
Consequently the "3F4 lifetime ~20-25 us" quoted in Jackson 2003 is the
MODERN 3H4 lifetime, and the cross-relaxation "3F4,3H6 -> 3H4,3H4" in Jackson
1999/2003 is the MODERN 3H4+3H6 -> 3F4+3F4 ("two-for-one") process.

What is NOT in these papers (so you do not go looking)
------------------------------------------------------
* Jackson & Mossman 2003 contains NO lifetime measurements and NO
  lifetime-vs-Tm or lifetime-vs-Al series.  Its concentration series is a
  LASER-PERFORMANCE series (slope efficiency / threshold vs Tm wt% and Al:Tm
  ratio).  The only lifetimes it mentions are literature values for the modern
  3H4 level (~20 us germanosilicate, ~25 us aluminosilicate, from Lincoln et
  al., J. Lumin. 50, 297 (1991)), p. 2705.  No energy-transfer coefficients.
* None of the three papers gives Judd-Ofelt Omega parameters, branching ratios,
  Stark-level energies, or a numerical refractive index.  (Jackson 1999 assumed
  Tm:silica energy levels coincide with Tm:YAG, ref. [15] Gruber 1989, p. 951;
  Agger's Eq. (9) contains n but no value is stated.)
* Agger 2006 did NOT use McCumber or Fuchtbauer-Ladenburg scaling to obtain
  sigma_e; it explicitly rejected McCumber (p. 51, p. 54) and used the gain-tilt
  method anchored at one measured eta = sigma_e/sigma_a.  F-L (their Eq. 9) is
  used only afterwards, to derive tau_rad from the scaled sigma_e.
* Agger gives no core diameter, NA, or cut-off wavelength for Tm1/Tm2.
* The only ion-ion coefficients (cm^3/s) in the three papers are in Jackson
  1999 Table I, and those were BORROWED from Tm:YAG / Tm,Ho:YAG literature and
  from the authors' earlier modelling, not measured in silica.

Unit conversions used
---------------------
* wt% (or ppm by weight) of Tm -> ions/cm^3:
      N = w * rho_glass * N_A / M_Tm,   rho = 2.20 g/cm^3, M_Tm = 168.934 g/mol
      => 1 wt% = 7.843e19 cm^-3 ;  1000 ppm wt = 7.843e18 cm^-3.
  This is the convention the authors themselves use: Jackson 1999 states
  3000 ppm = 2.35e25 m^-3 (Fig. 4 caption) and 800 ppm = 6.3e24 m^-3 (p. 952),
  both reproduced by this formula; Jackson 2003's wt% values are consistent
  with alpha_core/sigma_abs (Table 1, Eq. 3) under the same formula (checked:
  26 m^-1 / 4.8e-25 m^2 = 5.4e19 cm^-3 vs 0.7 wt% -> 5.5e19 cm^-3).
* 1 m^3/s = 1e6 cm^3/s ; 1e-25 m^2 = 1e-21 cm^2 ; 1e25 m^-3 = 1e19 cm^-3.
* dB/m -> 1/m : multiply by ln(10)/10 = 0.2303.

Typos / inconsistencies found and how they are handled
------------------------------------------------------
1. Jackson 2003 p. 2704: "sigma_abs = 4.8e-25 m^-2"  -> should be m^2.
2. Agger 2006 p. 53: "sigma_a = 4.3 +/- 0.1e-25 m2 for Tm1" (second value) ->
   must be Tm2 (first value 3.5e-25 is Tm1; confirmed by digitizing Fig. 4/5:
   sigma_a(1600 nm) = 3.50 (Tm1) and 4.30 (Tm2)).
3. Agger 2006 p. 55: peak sigma_e "3.5 and 3.9e-25 m2" is listed in the same
   order as sigma_a "4.2 and 4.5" (Tm1, Tm2), but digitizing Figs. 4/5 gives
   sigma_e,peak = 3.76 (Tm1) and 3.45 (Tm2).  So Tm1 ~3.8-3.9, Tm2 ~3.5; the
   text order of the emission values is reversed (or ambiguous).
4. Agger Table 1 lists sigma_a,peak "4.2-4.4" whereas text p. 55 and Fig. 3
   give 4.2 (Tm1) and 4.5 (Tm2).
5. Agger Eq. (5) misprints (1 - x sigma_as); should read (1 - x) sigma_as.
6. Agger Fig. 3/Fig. 1 axis label "lambda (mn)" -> nm.  "4H6" on p. 57 -> 3H4.
7. Agger: eta(1600 nm) stated 0.15 +/- 0.02 for Tm2 (p. 53), but the plotted
   Tm2 curves give sigma_e/sigma_a(1600) ~ 0.31 (Tm1: plotted 0.12 vs stated
   0.10, consistent).  The plotted Tm2 sigma_e is noisy below 1650 nm.
8. Agger: tau_rad from F-L with my digitized sigma_e (n = 1.45) is 6.2 ms
   (Tm1; paper 6.0 ms, consistent) but 5.6 ms (Tm2; paper 6.6 ms).  The Tm2
   paper value is not reproduced from the plotted curve; use with caution.
9. Jackson 1999 Fig. 3 caption says the RIGHT inset is 3H4->3F2,3 ESA, but the
   insets themselves are labelled the other way round: the LEFT inset is
   labelled old "3H4 to 3F2,3" (E1, units 1e-25 m^2) and the RIGHT inset old
   "3F4 to 1G4" (E2, units 1e-26 m^2).  The inset labels are used here.
10. Jackson 2003 ref. 3 (Hayward et al.) printed as Electron. Lett. 36 (2002);
    vol. 36 is 2000.  Irrelevant to numbers.

Digitization
------------
Values marked 'digitized' were read from the PDF figures in this work.
Agger Figs. 3-5 are vector colour plots; they were rasterized at 300 dpi and
the red/blue curves extracted automatically by colour, with axes calibrated on
the tick marks (estimated uncertainty: +/-3 nm, +/-0.05e-25 m^2, larger for
noisy Tm2 sigma_e).  Jackson 1999 Figs. 2-3 and Jackson 2003 Figs. 3-5 were
read by eye from 300-400 dpi rasters with pixel calibration on the axis ticks
(uncertainties stated per block).
"""

# ---------------------------------------------------------------------------
# helpers / constants used for conversions (documented above)
# ---------------------------------------------------------------------------
RHO_SILICA_g_cm3 = 2.20
M_TM_g_mol = 168.934
N_A = 6.02214e23
WT_PCT_TO_CM3 = RHO_SILICA_g_cm3 * N_A / M_TM_g_mol / 100.0  # 7.843e19 cm^-3 per wt%


def wt_pct_to_cm3(w):
    """Tm wt% (of element, in silica, rho=2.2) -> ions/cm^3."""
    return w * 1e-2 * RHO_SILICA_g_cm3 * N_A / M_TM_g_mol


# ===========================================================================
# JACKSON & KING 1999 (J. Lightwave Technol. 17, 948)
# ===========================================================================
JACKSON1999 = {
    "citation": "S. D. Jackson and T. A. King, J. Lightwave Technol. 17(5), 948-956 (1999)",
    "type": "rate-equation model of core-pumped CW Tm:silica fibre laser; spectra "
            "measured on a Tm-doped silica PREFORM (composition not stated)",
    "label_note": "paper uses OLD labels (old '3H4' = modern 3F4, old '3F4' = modern 3H4); "
                  "keys below are MODERN",

    # ---------------- lifetimes ----------------
    "lifetimes": {
        "src": "Fig. 1 (p. 949): 'calculated overall lifetimes (in us)' from radiative rates of "
               "Peng & Izumitani [9] (Tm in silica glass) + multiphonon rates by Reisfeld-Eyal "
               "method [10], with Tm:YAG level energies [15]. NOT measured; no concentration "
               "or composition attached.",
        "tau_3F4_us_calc": 334.7,      # Fig.1 N1, old "3H4"
        "tau_3H5_us_calc": 0.007,      # Fig.1 N2
        "tau_3H4_us_calc": 14.2,       # Fig.1 N3, old "3F4"
        "tau_3F23_us_calc": 0.0004,    # Fig.1 N4 (3F2 and 3F3 lumped)
        "tau_1G4_us_calc": 783.8,      # Fig.1 N5
        "tau_1D2_us_calc": 128.1,      # Fig.1 N6
        "tau_3P0_us_calc": 284.9,      # Fig.1 N7
        "tau_3F4_literature_range_us": (200.0, 500.0),
        "tau_3F4_literature_src": "p. 954: 'the 3H4 [old; = modern 3F4] energy level lifetime has "
                                  "been reported to be at values between 200 and 500 us'; Fig. 6 "
                                  "models threshold vs this lifetime over 200-500 us.",
        "tau_3F4_rad_ms_inferred_by_Agger": 4.2,
        "tau_3F4_rad_src": "NOT stated in Jackson 1999; Agger 2006 Table 1 lists 4.2 ms for ref [11] "
                           "(= Jackson 1999) in italics, i.e. 'estimated from relevant data or graphs'.",
        "tau_3H4_ZBLAN_ms": 1.5,
        "tau_3H4_ZBLAN_src": "p. 954: '~1.5 ms [old] 3F4 energy level lifetime for ZBLAN' (modern 3H4).",
        "comparison_note": "p. 951: for each Table II comparison the [modern 3F4] lifetime was set to "
                           "the value reported in that experiment (values not listed).",
    },

    # ---------------- concentration conversions used in the paper ----------------
    "concentrations": {
        "src": "Fig. 4 caption p. 953; p. 951-952",
        "conversion_check": "3000 ppm wt = 2.35e25 m^-3 stated (Fig. 4 caption) -> implies rho=2.2 g/cm3, "
                            "M=168.93: 3000e-6*2.2*6.022e23/168.93 = 2.35e19 cm^-3",
        "default_model_ppm_wt": 3000, "default_model_N_cm3": 2.35e19,
        "list": [
            # (ppm wt, N cm^-3, context)
            (300, wt_pct_to_cm3(0.03), "low end of concentration scan (p. 954)"),
            (320, wt_pct_to_cm3(0.032), "reported conc. of Yamamoto [8] in-band laser (p. 952)"),
            (800, 6.3e18, "adjusted conc. for [8]; paper states 6.3e24 m^-3 (p. 952)"),
            (3000, 2.35e19, "default model value (Fig. 4 caption)"),
            (3700, 2.9e19, "adjusted conc. for Hanna [7] from 1.6 dB/m at 1064 nm; 2.9e25 m^-3 stated (p. 951)"),
            (4700, wt_pct_to_cm3(0.47), "adjusted conc. for Hanna [4] from 2 dB/m at 1064 nm (p. 952)"),
            (6000, wt_pct_to_cm3(0.6), "model (Fig. 5)"),
            (30000, wt_pct_to_cm3(3.0), "high end of concentration scan (p. 954)"),
        ],
        "al_ge_content": "not stated (preform composition not given)",
    },

    # ---------------- ion-ion energy transfer (Table I, p. 951; Eqs. 6, 20) ----------------
    "energy_transfer": {
        "src": "Table I (p. 951) and Eqs. (6a,b), (20a-d) (p. 949, 951). Values in m^3/s in paper; "
               "converted to cm^3/s (x1e6). Indices: 0=3H6, 1=modern 3F4, 2=3H5, 3=modern 3H4, "
               "4=3F2,3, 5=1G4, 6=1D2, 7=3P0. k_ijkl: ions in i,j -> k,l.",
        "rate_form": "CR1 = k3101*N3*N0 - k1310*N1^2 ; CR2 = k2101*N2*N0 - k1012*N1^2 ; "
                     "CR_I = k5631 N5 N3 ; CR_II = k5751 N5^2 ; CR_III = k5313 N5 N1 ; CR_IV = k5203 N5 N0",
        "k3101_cross_relaxation_3H4_3H6_to_3F4_3F4": {
            "value_cm3_s": 1.8e-16, "value_m3_s": 1.8e-22,
            "normalised_at": "Tm concentration 5.75 wt.% (from Falconieri et al. [12], Tm,Ho:YAG crystal)",
            "normalised_at_N_cm3_if_silica_rho2.2": wt_pct_to_cm3(5.75),  # 4.51e20, only if read as silica wt%
            "concentration_dependence": "'Given functional dependence on Tm3+ concentration according to "
                                        "[11]' (Petrin et al., Opt. Mater. 1, 111 (1992), Tm,Ho:YAG). The "
                                        "functional form itself is NOT given in the paper.",
            "src": "Table I, p. 951",
        },
        "k1310_inverse_3F4_3F4_to_3H4_3H6_upconversion": {
            "value": "0.084 * k3101", "value_cm3_s_at_norm": 0.084 * 1.8e-16,  # 1.51e-17
            "basis": "ratio reverse/forward CR determined for Tm(6%),Ho(0.5%):YAG, Armagan et al. [13] (p. 951)",
            "src": "Table I, p. 951",
        },
        "k2101_3H5_3H6_to_3F4_3F4": {"value_cm3_s": 3e-18, "value_m3_s": 3e-24,
                                     "basis": "ref [13] Armagan (Tm,Ho:YAG)", "src": "Table I"},
        "k1012_3F4_3F4_to_3H5_3H6": {"value": "0.5 * k2101", "value_cm3_s": 1.5e-18,
                                     "basis": "factor 0.5 from Rustad & Stenersen [14]", "src": "Table I, p. 951"},
        "k5631_1G4_3H4_to_1D2_3F4": {"value_cm3_s": 6.0e-17, "value_m3_s": 6.0e-23, "basis": "ref [1] Jackson&King JQE 1998", "src": "Table I"},
        "k5751_1G4_1G4_to_3P0_3F4": {"value_cm3_s": 1.2e-17, "value_m3_s": 1.2e-23, "basis": "ref [1]", "src": "Table I"},
        "k5313_1G4_3F4_to_3H4_3F23": {"value_cm3_s": 1.6e-17, "value_m3_s": 1.6e-23, "basis": "ref [1]", "src": "Table I"},
        "k5203_1G4_3H6_to_3H5_3H4": {"value_cm3_s": 1.0e-17, "value_m3_s": 1.0e-23, "basis": "ref [1]", "src": "Table I"},
        "conclusion": "p. 954/956: varying Tm from 300 to 30000 ppm wt (50% OC, 790 nm pump) changed modelled "
                      "slope only 39 -> 41.5% and threshold 75 -> 89 mW; authors concluded ion-ion processes "
                      "involving modern 3H4 have 'practically a negligible effect' because of its short "
                      "lifetime. (Contradicted experimentally by Jackson 2003.)",
    },

    # ---------------- cross sections ----------------
    "cross_sections": {
        "src": "Fig. 2 (emission, p. 949) and Fig. 3 (GSA + ESA, p. 950), preform of Tm-doped silica; "
               "peak values digitized here (+/-0.1e-25 m^2, +/-5 nm) unless noted",
        "sigma_e_3F4_3H6_peak_cm2": 6.0e-21, "sigma_e_3F4_3H6_peak_nm": 1837,      # digitized, Fig. 2
        "sigma_e_3F4_3H6_at_1945nm_cm2": 4.6e-21,                                   # digitized +/-0.3e-21
        "sigma_e_3F4_3H6_at_1630nm_cm2": 0.6e-21,                                   # digitized +/-0.2e-21
        "sigma_e_3F4_3H6_band_nm": (1570, 2200),
        "sigma_e_3F4_peak_quoted_by_Agger_cm2": 6.0e-21,  # Agger Table 1, ref [11]
        "sigma_a_3H6_3F4_peak_cm2": 4.47e-21, "sigma_a_3H6_3F4_peak_nm": 1633,     # digitized, Fig.3 old "3H4"
        "sigma_a_3H6_3F4_FWHM_nm": 150,   # stated p. 953
        "sigma_a_3H6_3H5_peak_cm2": 4.6e-21, "sigma_a_3H6_3H5_peak_nm": 1209,      # digitized; shoulder 1.45e-21 @ ~1175 nm
        "sigma_a_3H6_3H4_peak_cm2": 8.6e-21, "sigma_a_3H6_3H4_peak_nm": 790,       # digitized, old "3F4"; text: peak at 790 nm
        "sigma_a_3H6_3F23_peak_cm2": 3.2e-21, "sigma_a_3H6_3F23_peak_nm": 681,     # digitized
        "sigma_a_1064nm_cm2": 1.26e-22,   # stated p. 951: 'accurately measured' 1.26e-26 m^2
        "ESA_3F4_to_3F23_peak_cm2": 3.85e-21, "ESA_3F4_to_3F23_peak_nm": 1068,     # digitized, inset 'E1', band ~980-1155 nm
        "ESA_3H4_to_1G4_peak_cm2": 5.7e-22, "ESA_3H4_to_1G4_peak_nm": 1133,       # digitized, inset 'E2' (1e-26 m^2 units), band ~1050-1187 nm
        "ESA_note": "inset labels old '3H4 to 3F2,3' (= modern 3F4->3F2,3) and old '3F4 to 1G4' "
                    "(= modern 3H4->1G4); caption swaps left/right (see typo 9).",
        "background_loss_pump_per_m": 1.2e-2,   # p. 949 (model input, delta_p)
        "background_loss_laser_per_m": 2.3e-3,  # p. 949 (model input, delta_s)
        "judd_ofelt": None, "branching_ratios": None, "stark_levels": "assumed = Tm:YAG [15]; none listed",
        "refractive_index": None,
    },

    # ---------------- validation targets (Table II, p. 951; text p. 951-955) ----------------
    "laser_results": {
        "src": "Table II (p. 951): measured (literature) vs calculated slope & threshold",
        "table_II": [
            # pump band (modern), literature ref, measured slope %, measured threshold mW, model slope %, model thr mW
            {"pump_level": "3H4 (~790 nm)", "ref": "[16] Barnes & Townsend, Electron. Lett. 26, 746 (1990)",
             "meas_slope_pct": 17, "meas_thr_mW": 4, "calc_slope_pct": 18.5, "calc_thr_mW": 5.5},
            {"pump_level": "3H4 (~790 nm)", "ref": "[2] Myslinski et al., Opt. Eng. 32, 2025 (1993)",
             "meas_slope_pct": 30, "meas_thr_mW": 36, "calc_slope_pct": 35, "calc_thr_mW": 24},
            {"pump_level": "3H5 (1064 nm)", "ref": "[7] Hanna et al., Electron. Lett. 25, 1365 (1989)",
             "meas_slope_pct": 30, "meas_thr_mW": 60, "calc_slope_pct": 19, "calc_thr_mW": 30,
             "note": "Tm adjusted to 3700 ppm (2.9e19 cm^-3) to match 1.6 dB/m at 1064 nm; with 5% OC model "
                     "gives 30% / 40 mW (p. 952)"},
            {"pump_level": "3H5 (1064 nm)", "ref": "[4] Hanna et al., Opt. Commun. 80, 52 (1990)",
             "meas_slope_pct": 37, "meas_thr_mW": 600, "calc_slope_pct": 33, "calc_thr_mW": 240,
             "note": "Tm adjusted to 4700 ppm to match 2 dB/m at 1064 nm"},
            {"pump_level": "3F4 in-band (1570 nm)", "ref": "[8] Yamamoto et al., Electron. Lett. 30, 220 (1994)",
             "meas_slope_pct": 71, "meas_thr_mW": 95, "calc_slope_pct": 67, "calc_thr_mW": 100,
             "note": "reported 320 ppm gives only 33% pump absorption in 4.0 m; model used 800 ppm "
                     "(6.3e18 cm^-3, 90% absorbed in 4.0 m); co-propagating output in experiment"},
        ],
        "model_config_Fig4": "core 10 um, Tm 3000 ppm (2.35e19 cm^-3), 10% OC, alpha_p*L = 10, "
                             "lambda_laser fixed 1945 nm; L = 0.5 m (790 band), 0.924 m (1.2 um band), "
                             "0.953 m (1.6 um band) (Fig. 4 caption, p. 953)",
        "model_max_slope_pct": {"790nm_band": 40, "1210nm_band": 57, "1630nm_band": 84},
        "model_details": "790 band: max 40% for 785-795 nm, threshold 17-18 mW (p. 952); 1.2 um band: 57%, "
                         "22 mW (p. 952); 1.6 um band: 83% @1700 nm, 84% @1720 nm (p. 953); Stokes limit "
                         "87% for 1700 nm pump.",
        "own_double_clad_laser": "p. 955 (ref [3] Jackson&King OL 23, 1462 (1998)): 790 nm diode cladding "
                                 "pumped, 1.58 m fibre, 100%-50% mirrors: 5.4 W out, 31% slope, 3 W incident "
                                 "threshold, 22% optical-optical.",
    },
    "quenching_OH": "No OH data. Model treats modern 3H4 decay as multiphonon (14.2 us). Ion-ion "
                    "processes deemed negligible for silica (p. 954, 956).",
}


# ===========================================================================
# JACKSON & MOSSMAN 2003 (Appl. Opt. 42, 2702)
# ===========================================================================
_J03_FIBRES = [
    # Table 1, p. 2703 ; Al:Tm ratio is that of the SOLUTION (Al not measured, p. 2703)
    {"fibre": 1, "Al_Tm_ratio": 5.0, "area_ratio_core_clad": 0.0063, "alpha_eff_per_m": 0.16,
     "alpha_core_per_m": 26.0, "Tm_wt_pct": 0.7, "maker": "Nufern"},
    {"fibre": 2, "Al_Tm_ratio": 2.0, "area_ratio_core_clad": 0.0074, "alpha_eff_per_m": 0.34,
     "alpha_core_per_m": 48.0, "Tm_wt_pct": 1.2, "maker": "Nufern"},
    {"fibre": 3, "Al_Tm_ratio": 9.0, "area_ratio_core_clad": 0.0032, "alpha_eff_per_m": 0.16,
     "alpha_core_per_m": 50.0, "Tm_wt_pct": 1.3, "maker": "OFTC (Univ. Sydney)"},
    {"fibre": 4, "Al_Tm_ratio": 10.0, "area_ratio_core_clad": 0.0039, "alpha_eff_per_m": 0.22,
     "alpha_core_per_m": 58.0, "Tm_wt_pct": 1.6, "maker": "OFTC (Univ. Sydney)"},
]
for _f in _J03_FIBRES:
    _f["N_Tm_cm3_from_wt"] = round(wt_pct_to_cm3(_f["Tm_wt_pct"]), -17)
    _f["N_Tm_cm3_from_alpha_core"] = round(_f["alpha_core_per_m"] / 4.8e-25 * 1e-6, -17)  # Eq. (3)
    _f["N_Al_cm3_est"] = round(_f["Al_Tm_ratio"] * _f["N_Tm_cm3_from_wt"], -18)
    # core diameter estimate from area ratio with ~300 um (full-circle) cladding; truncation ignored
    _f["core_diam_um_est_from_area_ratio"] = round(300.0 * _f["area_ratio_core_clad"] ** 0.5, 1)
    _f["src"] = ("Table 1, p. 2703; N, N_Al, core-diameter estimates derived in this work "
                 "(core estimates for fibres 3,4 fall below the stated ~20-25 um, so the area ratio or "
                 "cladding size is uncertain; prefer the stated range)")

JACKSON2003 = {
    "citation": "S. D. Jackson and S. Mossman, Appl. Opt. 42(15), 2702-2707 (2003)",
    "label_note": "paper uses OLD labels (old '3H4' = modern 3F4 upper laser level; old '3F4' = modern "
                  "3H4 pump level); keys below are MODERN",
    "lifetimes": {
        "measured_in_this_paper": None,
        "note": "NO lifetimes measured; no lifetime-vs-concentration data in this paper.",
        "tau_3H4_us_literature_germanosilicate": 20.0,
        "tau_3H4_us_literature_aluminosilicate": 25.0,
        "src": "p. 2705: 'measured lifetimes of ~20 us in germanosilicate fibers and of ~25 us in "
               "aluminosilicate fibers [ref 21: Lincoln et al., J. Lumin. 50, 297 (1991)], even for Tm3+ "
               "concentrations in the former fibers of as much as ~1.2 wt.%' (old label '3F4' = modern 3H4). "
               "1.2 wt% -> 9.4e19 cm^-3.",
    },
    "fibres": _J03_FIBRES,
    "fibre_common": {
        "src": "p. 2703",
        "NA": 0.19, "core_diam_um": (20, 25), "clad_diam_um": 300, "clad_shape": "truncated circle (D-shape)",
        "coating_index": 1.41, "codopants": "Al (Al:Tm in soak solution, Table 1) + Ge (index raiser, varied "
                                            "to keep NA ~0.19); MCVD + solution doping",
        "background_loss_dB_per_km_at_1990nm": 400,
        "background_loss_note": "small-signal background ~400 dB/km at 1.99 um, comprising Tm 3H6->3F4 "
                                "reabsorption, OH absorption and host (GeO2, Al2O3) loss; proportions not "
                                "determined (p. 2703). ~5 dB background loss at optimal lengths with 805 nm "
                                "pumping (p. 2706).",
        "sigma_abs_pump_cm2": 4.8e-21,
        "sigma_abs_src": "p. 2704, Eq. (3): pump-spectrum-weighted absorption cross section for the ~805 nm "
                         "diode, computed from Jackson&King 1999 spectrum; '~1/2 of the value at 790 nm'; "
                         "printed unit m^-2 is a typo for m^2.",
        "conc_method": "cutback alpha_eff -> alpha_core = alpha_eff/area_ratio (Eq. 1,2) -> N0 = alpha_core/"
                       "sigma_abs (Eq. 3)",
    },
    "laser_setup": {
        "src": "p. 2703-2704, Fig. 2",
        "pump": "FLS-Iron GaAlAs diode, 50 W max, centre ~805 nm (15 nm longer than 790 nm peak), "
                "cladding pumped, f=30 mm launch lens",
        "cavity": "dichroic mirror (HR 1.9-2.1 um, HT 800 nm) butted at pump end; Fresnel reflection at output "
                  "end (value not given, ~3.5% for silica)",
        "lasing_wavelength": "~2 um (not specified precisely; monochromator in setup)",
        "stokes_limit_pct_805nm": (39, 41),  # p. 2706
        "slope_basis": "slope efficiency w.r.t. LAUNCHED pump power",
    },
    "laser_results_summary": {
        "src": "text p. 2704-2705; Fig. 3 inset, Fig. 4 inset",
        # fibre: (max slope %, min launched threshold W)
        "fibre1_0.7wt_AlTm5": {"max_slope_pct": 21, "min_thr_W": 6.7},
        "fibre2_1.2wt_AlTm2": {"max_slope_pct": 8, "min_thr_W": 10.7,
                               "note": "'~8%' and '1.6-fold increase in P_launch at threshold' vs fibre 1 "
                                       "(1.6*6.7 = 10.7 W, derived); not plotted"},
        "fibre3_1.3wt_AlTm9": {"max_slope_pct": 41, "min_thr_W": 4.4},
        "fibre4_1.6wt_AlTm10": {"max_slope_pct": 46, "min_thr_W": 4.6},
        "max_slope_vs_conc": {"Tm_wt_pct": [0.7, 1.3, 1.6],
                              "N_Tm_cm3": [wt_pct_to_cm3(w) for w in (0.7, 1.3, 1.6)],
                              "max_slope_pct": [21, 41, 46], "src": "Fig. 3 inset + text"},
        "min_thr_vs_conc": {"Tm_wt_pct": [0.7, 1.3, 1.6], "min_thr_W": [6.7, 4.4, 4.6],
                            "src": "Fig. 4 inset + text"},
    },
    # digitized curves; L*alpha_eff dimensionless. Fibre length L = (L*alpha_eff)/alpha_eff (Table 1).
    "fig3_slope_vs_L_alpha_eff_digitized": {
        "src": "Fig. 3, p. 2704; digitized by eye at 400 dpi; uncertainty +/-0.04 in L*alpha, +/-0.7 % slope",
        "fibre4_1.6wt": [(0.285, 8.3), (0.46, 15.6), (0.64, 20.9), (0.81, 25.3), (0.985, 28.8), (1.16, 32.9),
                         (1.33, 33.6), (1.50, 36.4), (1.67, 36.8), (1.84, 38.3), (2.02, 44.6), (2.19, 46.1),
                         (2.37, 42.8), (2.46, 43.7), (2.54, 43.2)],
        "fibre3_1.3wt": [(0.50, 16.9), (0.66, 19.3), (0.83, 22.9), (0.985, 26.5), (1.15, 26.9), (1.31, 32.7),
                         (1.47, 32.3), (1.63, 34.4), (1.79, 36.4), (1.95, 37.1), (2.11, 35.6), (2.27, 36.2),
                         (2.43, 40.1), (2.60, 41.3), (2.75, 36.2), (2.92, 36.2), (3.08, 34.3), (3.24, 33.0),
                         (3.40, 31.3), (3.56, 30.5), (3.73, 32.9), (3.89, 31.9), (4.05, 30.2), (4.21, 28.0),
                         (4.37, 27.9)],
        "fibre1_0.7wt": [(0.19, 2.8), (0.28, 3.7), (0.36, 5.7), (0.44, 7.3), (0.53, 8.8), (0.62, 9.6),
                         (0.71, 10.7), (0.79, 11.8), (0.88, 12.3), (0.97, 13.2), (1.06, 14.5), (1.14, 15.5),
                         (1.23, 16.6), (1.32, 16.4), (1.40, 18.2), (1.56, 18.4), (1.64, 18.0), (1.73, 19.1),
                         (1.82, 20.6), (1.91, 20.9), (1.99, 20.9), (2.08, 20.1), (2.17, 19.0), (2.25, 18.7),
                         (2.34, 18.5), (2.43, 18.3)],
        "units": "(L*alpha_eff, slope_pct)",
    },
    "fig4_threshold_vs_L_alpha_eff_digitized": {
        "src": "Fig. 4, p. 2705; digitized by eye at 400 dpi; uncertainty +/-0.04 in L*alpha, +/-0.07 W",
        "fibre4_1.6wt": [(0.28, 5.71), (0.46, 5.30), (0.63, 4.80), (0.80, 4.65), (0.98, 4.61), (1.15, 4.66),
                         (1.32, 4.83), (1.49, 4.99), (1.66, 4.96), (1.83, 5.14), (2.01, 5.36), (2.18, 5.23),
                         (2.36, 5.23), (2.44, 5.31), (2.53, 5.48)],
        "fibre3_1.3wt": [(0.50, 4.94), (0.66, 4.50), (0.82, 4.39), (0.99, 4.48), (1.15, 4.57), (1.31, 4.75),
                         (1.47, 4.85), (1.63, 4.85), (1.79, 5.22), (1.95, 5.29), (2.11, 5.23), (2.27, 5.31),
                         (2.43, 5.38), (2.59, 5.95), (2.75, 5.94), (2.92, 6.25), (3.08, 5.98), (3.25, 6.25),
                         (3.41, 6.31), (3.57, 6.55), (3.74, 6.91), (3.90, 6.90), (4.06, 7.12), (4.22, 7.08),
                         (4.38, 7.42), (4.54, 7.65), (4.71, 7.68), (4.79, 7.91)],
        "fibre1_0.7wt": [(0.18, 10.29), (0.27, 9.66), (0.36, 8.97), (0.44, 9.23), (0.53, 8.12), (0.62, 7.59),
                         (0.70, 7.13), (0.79, 7.21), (0.88, 7.02), (0.97, 6.90), (1.06, 6.98), (1.14, 7.37),
                         (1.23, 7.04), (1.31, 7.05), (1.40, 6.95), (1.48, 6.74), (1.57, 6.97), (1.65, 6.97),
                         (1.74, 7.08), (1.82, 7.45), (1.91, 7.35), (1.99, 7.63), (2.08, 7.60), (2.17, 8.21),
                         (2.25, 8.15), (2.34, 8.28), (2.43, 8.85)],
        "units": "(L*alpha_eff, P_launch_threshold_W)",
    },
    "fig5_output_vs_launched_pump_digitized": {
        "src": "Fig. 5, p. 2705 (fibre lengths giving largest slope); digitized by eye at 400 dpi, "
               "+/-0.07 W in pump, +/-0.02 W in output. Line fits from plotted lines.",
        "fibre4_1.6wt": {"points_W": [(5.68, 0.26), (6.46, 0.57), (7.26, 0.88), (8.05, 1.26), (8.86, 1.68),
                                      (9.64, 2.04), (10.42, 2.40)],
                         "fit_slope": 0.46, "fit_threshold_W": 5.24,
                         "L_alpha_eff_at_max_slope": 2.2, "L_m_est": 10.0},
        "fibre3_1.3wt": {"points_W": [(6.46, 0.28), (7.26, 0.49), (8.05, 0.84), (8.85, 1.17), (9.64, 1.57),
                                      (10.41, 1.89), (11.21, 2.15)],
                         "fit_slope": 0.415, "fit_threshold_W": 5.99,
                         "L_alpha_eff_at_max_slope": 2.6, "L_m_est": 16.0},
        "fibre1_0.7wt": {"points_W": [(8.05, 0.15), (8.94, 0.34), (9.82, 0.53), (10.69, 0.72), (11.58, 0.90),
                                      (12.46, 1.07), (14.22, 1.48), (16.86, 2.02)],
                         "fit_slope": 0.213, "fit_threshold_W": 7.41,
                         "L_alpha_eff_at_max_slope": 2.0, "L_m_est": 12.5},
        "L_note": "L_m_est = (L*alpha_eff at max slope, Fig. 3)/alpha_eff (Table 1) -- derived, not stated. "
                  "Min-threshold lengths: fibre4 ~1.0/0.22 = 4.5 m, fibre3 ~0.82/0.16 = 5.1 m, fibre1 "
                  "~1.48/0.16 = 9.3 m (derived from Fig. 4).",
    },
    "energy_transfer": {
        "coefficients": None,
        "processes_named": {
            "cross_relaxation": "modern 3H4+3H6 -> 3F4+3F4 (old '3F4,3H6 -> 3H4,3H4'); invoked to explain slope "
                                "> Stokes limit (39-41% at 805 nm) at >=1.3 wt% Tm (p. 2706)",
            "energy_migration": "modern 3H4+3H6 -> 3H6+3H4 (old '3F4,3H6 -> 3H6,3F4'), proposed to enhance CR; "
                                "argued by analogy with dipole-quadrupole / quadrupole-quadrupole transfer "
                                "in Tm:fluoroindogallate (refs 22,23) (p. 2705-2706)",
            "ETU_quenching": "modern 3F4+3F4 -> 3H4+3H6 (old '3H4,3H4 -> 3F4,3H6'), endothermic; attributed to "
                             "Tm clustering at low Al:Tm; evidenced (ref 26, Jackson&King JOSA B 16, 2178 "
                             "(1999)) by 1.5-1.55 um (modern 3H4->3F4) fluorescence under 1.98 um pumping of "
                             "heavily doped fibre (p. 2706)",
        },
        "src": "Abstract; pp. 2702, 2705-2706",
    },
    "quenching_OH": {
        "OH": "OH absorption listed as one (unquantified) contributor to ~400 dB/km background at 1.99 um "
              "(p. 2703). No OH concentration given; no OH quenching discussed.",
        "clustering": "Recommendation: Al:Tm > 10 for Tm > 1 wt% to reduce clustering-induced ETU quenching of "
                      "modern 3F4 and/or scattering (abstract, p. 2706). Al and Tm both 6-coordinated; "
                      "Al-O-Tm bonds assumed (p. 2706).",
    },
}


# ===========================================================================
# AGGER & POVLSEN 2006 (Opt. Express 14, 50)
# ===========================================================================
_AG_GRID_nm = [1450, 1475, 1500, 1525, 1550, 1575, 1600, 1625, 1650, 1675, 1700, 1725, 1750, 1775,
               1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000, 2025, 2050]

AGGER2006 = {
    "citation": "S. D. Agger and J. H. Povlsen, Opt. Express 14(1), 50-57 (2006)",
    "label_note": "modern labels (3F4 = 2 um upper level, 3H4 = 786 nm pump level)",
    "fibres": {
        "src": "p. 51 (sec. 2), p. 56 (sec. 2.4), p. 57 (conclusion)",
        "Tm1": {"codopants": "Al + La ('aluminum and lanthanide', Al/La-silica), OFS Fitel Denmark custom; "
                             "chosen to optimise radiative efficiency and minimise clustering",
                "N_Tm_peak_m3": 8.4e25, "N_Tm_peak_cm3": 8.4e19,
                "Tm_wt_pct_equiv": round(8.4e19 / WT_PCT_TO_CM3, 2),   # 1.07 wt% (this work, rho=2.2)
                "N_src": "p. 56, Eq. (10) rho0 = GSA/(sigma_a*Gamma); peak conc., profile assumed ~ index profile",
                "Al_La_content": "not given"},
        "Tm2": {"codopants": "Al/Ge", "N_Tm_peak_m3": 6.1e25, "N_Tm_peak_cm3": 6.1e19,
                "Tm_wt_pct_equiv": round(6.1e19 / WT_PCT_TO_CM3, 2),   # 0.78 wt%
                "N_src": "p. 56", "Al_Ge_content": "not given"},
        "geometry": "core diameter, NA, index profile values NOT given (index profile measured, used in FEM "
                    "mode solver); ratio of confinement factors Tm2/Tm1 ~0.7 (p. 56, ambiguous direction)",
        "Er_check_fibre": "same design as Tm1 with Er ~11x lower (fab.); measured sigma_a(1470)=1.24e-21 cm^2 "
                          "using eta=0.22, GSA=0.63 /m -> N_Er = 8.7e18 cm^-3 = Tm1/9.5 (p. 56)",
    },
    "lifetimes": {
        "src": "p. 52 (Fig. 1 and text), p. 55-57",
        "conditions": "fluorescence decay after step pump at 786 nm (into 3H4), 2 mm long fibre, detected "
                      "3F4 fluorescence, double-exponential least-squares fit; room temperature (implied)",
        "tau_3F4_fl_us": {"Tm1": 650.0, "Tm2": 560.0},
        "tau_3F4_rad_ms": {"Tm1": 6.0, "Tm2": 6.6},
        "tau_rad_method": "Eq. (9) Fuchtbauer-Ladenburg: 1/tau_r = 8 pi n^2 c Int(sigma_e/lambda^4 dlambda) "
                          "using the gain-tilt-scaled sigma_e (n value not stated). Digitized-curve check "
                          "(n=1.45): Tm1 6.2 ms (ok), Tm2 5.6 ms (paper 6.6, not reproduced).",
        "quantum_efficiency_3F4": 0.10,  # 'only around 10%' p. 56 (650/6000 = 0.108; 560/6600 = 0.085)
        "tau_3H4_us_Tm1_fast_component": 12.0,
        "tau_3H4_src": "Fig. 1 legend: 2xEXP fit tau1 = 12 us, tau2 = 650 us (Tm1). p. 57: 3H4 lifetime not "
                       "measured directly; fast component 'not very accurate'; literature 14-20 us.",
        "tau_3F4_literature_us": (300, 600),   # p. 52, ref [7] Digonnet
        "max_inversion_786nm": 0.97,  # p. 57, (1+tau32/tau21)^-1
        "conc_dependence": "none measured (two fibres only)",
    },
    "cross_sections": {
        "src": "pp. 53-56, Figs. 3-5, Table 1 (p. 56)",
        "method": "sigma_a(1600) from saturated-fluorescence fit (Eq. 2) using tau_fl and eta_p(1600) from "
                  "ratio of saturated fluorescence for 786 vs 1600 nm pumping (assumes eta=0 at 786 nm); "
                  "sigma_a(lambda) by GSA ratio with FEM confinement factor (Eq. 6); sigma_e by gain-tilt "
                  "(Eq. 7-8) anchored at eta(1740) from small-signal gain vs GSA. NO McCumber scaling used "
                  "(rejected: broad, largely homogeneous line; sigma_a tiny beyond 1.9 um).",
        "eta_1600nm": {"Tm1": (0.10, 0.02), "Tm2": (0.15, 0.02)},          # p. 53
        "sigma_a_1600nm_cm2": {"Tm1": (3.5e-21, 0.1e-21), "Tm2": (4.3e-21, 0.1e-21)},  # p. 53 (typo fixed)
        "eta_1740nm": {"Tm1": (1.20, 0.02), "Tm2": (1.25, 0.02)},          # p. 54
        "sigma_a_peak_cm2_text": {"Tm1": 4.2e-21, "Tm2": 4.5e-21},          # p. 55
        "sigma_a_peak_cm2_table1": "4.2-4.4e-21",
        "sigma_a_peak_nm_digitized": {"Tm1": 1649, "Tm2": 1617},           # +/-5 nm; Tm2 flat top 1610-1630
        "sigma_a_peak_cm2_digitized": {"Tm1": 4.18e-21, "Tm2": 4.48e-21},
        "sigma_e_peak_cm2_text": "3.5 and 3.9e-21 (p. 55); Table 1 '3.5-3.9'",
        "sigma_e_peak_cm2_digitized": {"Tm1": 3.76e-21, "Tm2": 3.45e-21},  # see typo 3
        "sigma_e_peak_nm_digitized": {"Tm1": 1817, "Tm2": "broad plateau ~1740-1830 (noisy; max ~1760)"},
        "sigma_e_over_sigma_a_peak_ratio": "<1 (abstract): ~0.90 (Tm1), ~0.77 (Tm2) from digitized peaks",
        "ESA_note": "signal ESA present below 1550 nm and above 2065 nm: gain-tilt sigma_e invalid there (p. 55)",
        "ETU_note": "faint blue (~470 nm) emission under 786 nm pumping attributed to 3H4+3H4 -> 3F4+1G4 ETU; "
                    "judged weak (1850 nm fluorescence never decreased with pump up to >500 mW) (p. 57). "
                    "No coefficient given.",
        "spectra_digitized_1e-25_m2": {
            "src": "Figs. 4 (Tm1) and 5 (Tm2), p. 55; automated colour digitization at 300 dpi, +/-3 nm, "
                   "+/-0.05e-25 m^2 (Tm2 sigma_e noise ~+/-0.1); None = curve absent; values ~0 beyond "
                   "~1950 nm for sigma_a are at noise floor",
            "lambda_nm": _AG_GRID_nm,
            "Tm1_sigma_a": [0.02, 0.05, 0.19, 0.51, 1.20, 2.37, 3.50, 4.09, 4.18, 4.10, 3.79, 3.25, 2.50, 1.79,
                            1.24, 0.86, 0.56, 0.33, 0.18, 0.09, 0.03, 0.02, 0.01, 0.02, 0.04],
            "Tm1_sigma_e": [None, None, None, None, None, None, 0.43, 0.89, 1.27, 1.91, 2.61, 3.15, 3.49, 3.63,
                            3.73, 3.75, 3.60, 3.27, 2.85, 2.49, 2.14, 1.75, 1.29, 0.73, 0.26],
            "Tm2_sigma_a": [0.05, 0.12, 0.32, 0.86, 1.92, 3.29, 4.29, 4.47, 4.39, 4.17, 3.72, 3.09, 2.37, 1.73,
                            1.29, 0.93, 0.64, 0.45, 0.31, 0.23, 0.19, 0.13, 0.11, 0.08, 0.11],
            "Tm2_sigma_e": [None, None, None, None, 0.09, 0.68, 1.33, 1.86, 2.36, 2.89, 3.26, 3.34, 3.39, 3.38,
                            3.36, 3.36, 3.14, 2.80, 2.39, 2.06, 1.73, 1.36, 0.97, 0.64, 0.04],
            "checks": "digitized Tm1 sigma_a(1600)=3.50, eta(1740)=1.20; Tm2 sigma_a(1600)=4.30, eta(1740)=1.27 "
                      "-- match text values",
        },
        "judd_ofelt": None, "branching_ratios": None, "stark_levels": None, "refractive_index": None,
    },
    "table1_literature_comparison": {
        "src": "Table 1, p. 56 (italic = estimated by Agger from graphs/data)",
        "rows": [
            {"medium": "SM fiber (this work)", "sigma_a_peak_1e-25m2": "4.2-4.4", "sigma_e_peak_1e-25m2": "3.5-3.9", "tau_r_ms": "6.0-6.6"},
            {"medium": "Fiber", "ref": "Jackson&King 1999 [11]", "sigma_a_peak_1e-25m2": 4.5, "sigma_e_peak_1e-25m2": 6.0, "tau_r_ms": 4.2, "italic": "tau_r"},
            {"medium": "Bulk", "ref": "Zou&Toratani JNCS 195, 113 (1996) [12]", "sigma_e_peak_1e-25m2": 6.1, "tau_r_ms": 6.3, "note": "McCumber + JO"},
            {"medium": "Preform", "ref": "Walsh&Barnes APB 78, 325 (2004) [13]", "sigma_a_peak_1e-25m2": 4.3, "sigma_e_peak_1e-25m2": 4.6, "tau_r_ms": 4.56, "italic": "sigma_a", "note": "JO-scaled"},
            {"medium": "Bulk", "ref": "Gandy et al. JAP 36, 3030 [14]", "sigma_e_peak_1e-25m2": 4.0},
        ],
    },
    "laser_results": None,
    "quenching_OH": "No OH data. Low radiative QE (~10%) attributed to multiphonon decay (high silica phonon "
                    "energy, small 3F4-3H6 gap) (p. 56). Clustering believed small (long lifetimes, modest "
                    "conc.) (p. 57).",
}


if __name__ == "__main__":
    print("1 wt%% Tm in silica = %.3e cm^-3" % WT_PCT_TO_CM3)
    for f in JACKSON2003["fibres"]:
        print(f["fibre"], f["Tm_wt_pct"], "%.2e" % f["N_Tm_cm3_from_wt"], "%.2e" % f["N_Tm_cm3_from_alpha_core"],
              "%.2e" % f["N_Al_cm3_est"], f["core_diam_um_est_from_area_ratio"])
