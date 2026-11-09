"""
Reference data extracted from two papers, for a forward model of rare-earth-doped
fibre gain (gain windows, lifetimes, cross sections) and for validation against
demonstrated fibre lasers.

Sources
-------
[J12]  S. D. Jackson, "Towards high-power mid-infrared emission from a fibre laser",
       Nature Photonics 6, 423-431 (2012), doi:10.1038/nphoton.2012.149  (review).
       Page numbers below are journal pages (p.423 = first page of the PDF).
[LG20] J. Le Gouet, F. Gustave, P. Bourdon, T. Robin, A. Laurent, B. Cadier,
       "Realization and simulation of high power holmium doped fiber laser for
       long-range transmission", arXiv:2003.14188v1 [physics.optics] (31 Mar 2020)
       (later published by OSA, (c) 2021 notice on p.1).  ONERA / iXblue.
       Section / figure / table numbers refer to the arXiv v1 PDF.

Caveats (read before using)
---------------------------
* Jackson 2012 is a REVIEW. It contains NO tables of level lifetimes, branching
  ratios, absorption cross sections or energy-transfer rate parameters for any ion.
  The only table is Table B1 (demonstrated lasers, p.424).  The only quantitative
  spectroscopic data are the emission cross-section spectra in Fig. 4 (p.426),
  which have been DIGITIZED here (coarse, +/-0.15e-25 m^2, +/-0.01 um; values for
  noisy traces -- Ho 5I7->5I8, Ho 5I6->5I7 -- are the upper envelope of the noise
  in each 0.02 um bin and therefore biased high by up to ~0.3e-25 m^2).  Fig. 4 does
  not state measurement temperature, concentrations or the primary source of the
  spectra.  Fig. 4 groups: "oxide glass" = aluminosilicate, "fluoride" = ZBLAN.
* Values in Table B1 are reproduced as printed.  Known internal inconsistencies:
  - Fig. 2 plots Ho3+ [25] at ~5 W (2.9 um), Table B1 gives 2.5 W at 2.86 um.
  - Fig. 2 plots Er3+ [24] at 2.7 um, Table B1 gives 2.8 um (text: "24 W at 2.7 um", p.428,
    and "24 W at wavelengths of around 2.8 um", p.426).
  - Row "Tm3+, Ho3+ / ZBLAN / 1.94 um / 3F4->3H6" (ref 33, Eichhorn & Jackson 2008,
    a Tm silica-vs-fluoride comparison) lists Ho as co-dopant although the transition
    is Tm 3F4->3H6; treat the dopant entry as suspect and check the primary paper.
  - Pump/slope/power should be checked against primaries before being used as hard
    validation targets (slope-efficiency definition -- absorbed vs launched -- is
    not stated in Table B1).
* Fig. 2 (p.424) is stated in text (p.425) to show the "maximum CW output power";
  all Table B1 entries are therefore taken as CW.  The 3.95 um Ho:ZBLAN laser
  (ref 29) required cryogenic cooling (Fig. 2 caption, p.426).
* Le Gouet 2020: the model is a TWO-LEVEL (5I8-5I7 only) homogeneous model with an
  effective lifetime; ETU is NOT modelled, PIQ is modelled only as a constant
  non-saturable loss.  No 5I6 lifetime is reported anywhere in the paper.  Host
  composition is given only qualitatively (Ho + Al solution-doped MCVD silica);
  no Al/Ho ratio is given.  Cross-section curves below are DIGITIZED from Fig. 1
  (left) of LG20 (+/-0.1e-25 m^2, +/-2 nm); they reproduce the textual anchor values
  (sigma_a,peak 4.3 @1946 nm, sigma_e,peak 6.0, sigma_a/sigma_e(1908 nm) = 1.9/1.6).
* LG20 text says the single-exponential fit covers "the first part of the
  fluorescence signal (until about 4 us)"; Fig. 2 shows the exponential fit curve
  extending to 4 ms, so "4 us" is almost certainly a typo for 4 ms.
* LG20 threshold is NOT stated; the value here is a linear extrapolation of points
  digitized from Fig. 6 (left).
"""

# =============================================================================
# JACKSON 2012
# =============================================================================

JACKSON2012 = {
    "citation": "S. D. Jackson, Nat. Photonics 6, 423-431 (2012), doi:10.1038/nphoton.2012.149",
    "note_no_spectroscopic_tables": (
        "The review contains no tables of lifetimes, branching ratios, absorption cross "
        "sections or energy-transfer parameters for Er, Ho, Dy, Pr, Tm or any other ion. "
        "Only Table B1 (lasers, p.424) and Fig. 4 (emission cross sections, p.426)."
    ),

    # ------------------------------------------------------------------ (b) lasers
    # Table B1, p.424: "Characteristics of infrared fibre lasers with emission
    # wavelengths >= 1.5 um."  Columns reproduced as printed. Reference numbers are
    # Jackson-2012 reference numbers; full primary citations from the reference list.
    "DEMONSTRATED_LASERS": [
        {"dopants": "Er3+,Yb3+", "host": "silicate", "pump_um": 0.975, "laser_um": 1.5,
         "transition": "Er 4I13/2->4I15/2", "power_max_W": 297, "slope_eff_pct": 19,
         "mode": "CW", "ref_no": 21,
         "primary": "Y. Jeong et al., IEEE JSTQE 13, 573-579 (2007) (Er:Yb large-core, 297 W CW)",
         "source": "Table B1 p.424; Fig.2 p.424 (Er3+[21], ~300 W at ~1.55 um)"},
        {"dopants": "Tm3+,Ho3+ (as printed; see caveat)", "host": "ZBLAN", "pump_um": 0.792,
         "laser_um": 1.94, "transition": "Tm 3F4->3H6", "power_max_W": 20, "slope_eff_pct": 49,
         "mode": "CW", "ref_no": 33,
         "primary": "M. Eichhorn & S. D. Jackson, Appl. Phys. B 90, 35-41 (2008)",
         "source": "Table B1 p.424; text p.428 ('20 W at 1.94 um (ref. 33)')"},
        {"dopants": "Tm3+", "host": "silicate", "pump_um": 0.793, "laser_um": 2.05,
         "transition": "Tm 3F4->3H6", "power_max_W": 1050, "slope_eff_pct": 53,
         "mode": "CW", "ref_no": 22,
         "primary": "T. Ehrenreich et al., '1 kW, all-glass Tm:fiber laser', SPIE Photonics West 2010, conf. 7850",
         "source": "Table B1 p.424; Fig.2 (Tm3+[22], ~1 kW at 2.0 um)"},
        {"dopants": "Tm3+,Ho3+", "host": "silicate", "pump_um": 0.793, "laser_um": 2.1,
         "transition": "Ho 5I7->5I8", "power_max_W": 83, "slope_eff_pct": 42,
         "mode": "CW", "ref_no": 34,
         "primary": "S. D. Jackson et al., Opt. Lett. 32, 241-243 (2007) (83 W Ho silica)",
         "source": "Table B1 p.424"},
        {"dopants": "Ho3+", "host": "silicate", "pump_um": 1.950, "laser_um": 2.14,
         "transition": "Ho 5I7->5I8", "power_max_W": 140, "slope_eff_pct": 55,
         "mode": "CW", "ref_no": 23, "pump_source": "Tm3+-doped silica fibre lasers (resonant cladding pumping)",
         "primary": "A. V. Hemming et al., Proc. SPIE 8237, 82371J (2012)",
         "source": "Table B1 p.424; Fig.2 (Ho3+[23], ~140 W at 2.1 um); p.429"},
        {"dopants": "Tm3+", "host": "ZBLAN", "pump_um": 1.064, "laser_um": 2.31,
         "transition": "Tm 3H4->3H5", "power_max_W": 0.15, "slope_eff_pct": 8,
         "mode": "CW", "ref_no": 35, "pump_scheme": "upconversion pumping at 1.064 um",
         "primary": "R. M. El-Agmy & N. M. Al-Hosiny, Electron. Lett. 46, 936-937 (2010)",
         "source": "Table B1 p.424"},
        {"dopants": "Er3+", "host": "ZBLAN", "pump_um": 0.975, "laser_um": 2.8,
         "transition": "Er 4I11/2->4I13/2", "power_max_W": 24, "slope_eff_pct": 13,
         "mode": "CW", "ref_no": 24, "cooling": "liquid-cooled",
         "primary": "S. Tokita et al., Opt. Lett. 34, 3062-3064 (2009) (liquid-cooled 24 W Er:ZBLAN)",
         "source": "Table B1 p.424; Fig.2 plots it at 2.7 um / ~24 W; text p.426, p.428"},
        {"dopants": "Ho3+,Pr3+", "host": "ZBLAN", "pump_um": 1.1, "laser_um": 2.86,
         "transition": "Ho 5I6->5I7", "power_max_W": 2.5, "slope_eff_pct": 29,
         "mode": "CW", "ref_no": 25, "beam": "single transverse mode",
         "primary": "S. D. Jackson, Opt. Lett. 29, 334-336 (2004)",
         "source": "Table B1 p.424 (Fig.2 plots ~5 W at 2.9 um -- inconsistent)"},
        {"dopants": "Dy3+", "host": "ZBLAN", "pump_um": 1.1, "laser_um": 2.9,
         "transition": "Dy 6H13/2->6H15/2", "power_max_W": 0.275, "slope_eff_pct": 4.5,
         "mode": "CW", "ref_no": 36,
         "primary": "S. D. Jackson, Appl. Phys. Lett. 83, 1316-1318 (2003)",
         "source": "Table B1 p.424"},
        {"dopants": "Ho3+", "host": "ZBLAN", "pump_um": 1.15, "laser_um": 3.002,
         "transition": "Ho 5I6->5I7 (cascade with 5I7->5I8 at ~2.1 um)", "power_max_W": 0.77,
         "slope_eff_pct": 12.4, "mode": "CW", "ref_no": 26, "pump_source": "diode, 1.15 um",
         "primary": "J. Li, D. D. Hudson & S. D. Jackson, Opt. Lett. 36, 3642-3644 (2011)",
         "source": "Table B1 p.424; Fig.2 (Ho3+[26] ~0.8 W at 3.0 um); text p.426, p.429"},
        {"dopants": "Ho3+", "host": "ZBLAN", "pump_um": 0.532, "laser_um": 3.22,
         "transition": "Ho 5S2->5F5", "power_max_W": 0.011, "slope_eff_pct": 2.8,
         "mode": "CW (room temperature)", "ref_no": 27,
         "primary": "C. Carbonnier, H. Tobben & U. B. Unrau, Electron. Lett. 34, 893-894 (1998)",
         "source": "Table B1 p.424"},
        {"dopants": "Er3+", "host": "ZBLAN", "pump_um": 0.653, "laser_um": 3.45,
         "transition": "Er 4F9/2->4I9/2", "power_max_W": 0.008, "slope_eff_pct": 3,
         "mode": "CW (room temperature)", "ref_no": 28,
         "primary": "H. Tobben, Electron. Lett. 28, 1361-1363 (1992)",
         "source": "Table B1 p.424"},
        {"dopants": "Ho3+", "host": "ZBLAN", "pump_um": 0.89, "laser_um": 3.95,
         "transition": "Ho 5I5->5I6", "power_max_W": 0.011, "slope_eff_pct": 3.7,
         "mode": "CW, CRYOGENICALLY cooled", "ref_no": 29,
         "primary": "J. Schneider, C. Carbonnier & U. B. Unrau, Appl. Opt. 36, 8595-8600 (1997)",
         "source": "Table B1 p.424; Fig.2 caption (cryogenic); p.426 (longest fibre-laser wavelength)"},
    ],

    # Additional laser results mentioned only in text / Fig. 2 (secondary targets)
    "OTHER_LASER_RESULTS_IN_TEXT": [
        {"desc": "Yb3+ silica, ~1.07 um, 10 kW single-mode", "ref_no": 20,
         "source": "Fig.2 p.424 (digitized: 1.07 um, 1.0e4 W)", "primary": "V. Fomin et al., 10 kW single-mode fiber laser (2010)"},
        {"desc": "Tm3+ silica tuning range 1.86-2.09 um (cladding pumped)", "ref_no": 43,
         "tuning_um": (1.86, 2.09), "source": "p.426", "primary": "W. A. Clarkson et al., Opt. Lett. 27, 1989 (2002)"},
        {"desc": "Longest silicate-glass fibre-laser wavelength 2.188 um (Tm:silica)", "ref_no": 63,
         "laser_um": 2.188, "source": "p.427", "primary": "Z. S. Sacks et al., Proc. SPIE 6453, 645320 (2007)"},
        {"desc": "Er:ZBLAN tunable 2.71-2.88 um, stable 10 W", "ref_no": 48,
         "tuning_um": (2.71, 2.88), "power_W": 10, "source": "p.426", "primary": "S. Tokita et al., Opt. Lett. 35, 3943 (2010)"},
        {"desc": "Er:ZBLAN all-fibre 20 W passively cooled single-mode at 2.8 um (FBGs)", "ref_no": 49,
         "power_W": 20, "laser_um": 2.8, "source": "p.426 (cited for Bragg gratings)", "primary": "D. Faucher et al., Opt. Lett. 36, 1104 (2011)"},
        {"desc": "Heavily Er-doped ZBLAN all-fibre laser at 2.94 um, slope 35.6% (> Stokes limit, ETU recycling)",
         "ref_no": 51, "slope_eff_pct": 35.6, "laser_um": 2.94, "source": "p.426", "primary": "D. Faucher et al., Opt. Lett. 34, 3313 (2009)"},
        {"desc": "Er:ZBLAN Q-switched 2.8 um: 90 ns, 0.9 kW peak (12 W avg in title of ref)", "ref_no": 99,
         "source": "p.428", "primary": "S. Tokita et al., Opt. Lett. 36, 2812 (2011)"},
        {"desc": "Er:ZBLAN gain-switched: 307 ns, 68 W peak", "ref_no": 100, "source": "p.428",
         "primary": "M. Gorjan et al., Opt. Lett. 36, 1923 (2011)"},
        {"desc": "Ho3+/Pr3+ fluoride Q-switched (single transition, 2.87 um) and cascaded Q-switched Ho fluoride: 70 ns pulses",
         "ref_no": (101, 102), "source": "p.428",
         "primary": "T. Hu et al., Opt. Lett. 37, 2145 (2012); J. Li et al., Opt. Lett. 37, 2208 (2012)"},
        {"desc": "Ho3+-doped fluoride fibre emission range 2.8-3.02 um", "ref_no": 57, "source": "p.426",
         "primary": "L. Wetenkamp et al., J. Non-Cryst. Solids 140, 35 (1992)"},
        {"desc": "Tm fluoride fibre amplifier 5 kW peak", "ref_no": 103, "source": "p.428",
         "primary": "M. Eichhorn, Opt. Lett. 30, 3329 (2005)"},
        {"desc": "Shortest 2 um fibre pulse 108 fs (Tm)", "ref_no": 92, "source": "p.428"},
        {"desc": "Directly diode-pumped Ho silica fibre lasers (1.15 um diodes)", "ref_no": 47, "source": "p.426",
         "primary": "S. D. Jackson, F. Bugge, G. Erbert, Opt. Lett. 32, 2496 (2007)"},
        {"desc": "Nd3+ chalcogenide (Ga-La-S) fibre laser at 1.08 um -- only RE chalcogenide fibre laser", "ref_no": 85,
         "source": "p.428"},
    ],

    # Fig. 2, p.424 -- digitized output power vs wavelength (log scale; +/-5% in power)
    "FIG2_DIGITIZED_um_W": {
        "Yb[20]": (1.07, 1.0e4), "Er[21]": (1.55, 299), "Tm[22]": (2.00, 1.0e3),
        "Ho[23]": (2.10, 141), "Er[24]": (2.70, 24), "Ho[25]": (2.90, 5.1),
        "Ho[26]": (3.00, 0.78), "Ho[27]": (3.22, 0.011), "Er[28]": (3.45, 0.0081),
        "Ho[29]": (3.95, 0.011),
        "_source": "Fig. 2 p.424, digitized by marker centroid; DIGITIZED",
    },

    # ------------------------------------------- (a) spectroscopic data (Fig. 4 only)
    # Emission cross sections, units 1e-25 m^2 (= 1e-21 cm^2).  DIGITIZED from Fig. 4 p.426.
    "FIG4_EMISSION_XS_PEAKS": {
        "_units": "peak_um, sigma_e_peak_1e-25m2; DIGITIZED from Fig.4 p.426",
        "Er_4I13/2-4I15/2_aluminosilicate": {"peak_um": 1.533, "sigma_e_peak_1e-25m2": 5.0,
                                              "span_um": (1.48, 1.60), "shoulder": "~3.3 at 1.55 um"},
        "Tm_3F4-3H6_aluminosilicate": {"peak_um": 1.81, "sigma_e_peak_1e-25m2": 5.5,
                                        "span_um": (1.55, 2.19)},
        "Ho_5I7-5I8_aluminosilicate": {"peak_um": 2.03, "sigma_e_peak_1e-25m2": 7.1,
                                        "secondary_peak": "~2.9 at 1.95 um", "span_um": (1.90, 2.20),
                                        "note": "noisy trace; peak is noise envelope (true ~6.8-7.0)"},
        "Er_4I11/2-4I13/2_ZBLAN": {"peak_um": 2.73, "sigma_e_peak_1e-25m2": 5.6,
                                   "shoulder": "~4.7 at 2.76-2.80 um", "span_um": (2.60, 2.94)},
        "Ho_5I6-5I7_ZBLAN": {"peak_um": 2.85, "sigma_e_peak_1e-25m2": 5.2,
                             "span_um": (2.60, 3.20), "note": "noisy trace"},
        "Dy_6H13/2-6H15/2_ZBLAN": {"peak_um": 2.87, "sigma_e_peak_1e-25m2": 2.95,
                                   "span_um": (2.44, 3.17)},
    },
    # Coarse curves (wavelength_um, sigma_e_1e-25m2), 0.02 um bins, upper envelope. DIGITIZED.
    "FIG4_EMISSION_XS_CURVES_um_1e-25m2": {
        "_source": "Fig. 4 p.426; DIGITIZED (colour extraction, axis calibrated on ticks)",
        "Ho_5I7-5I8_aluminosilicate": [
            (1.90, 0.49), (1.92, 0.96), (1.94, 2.81), (1.96, 2.89), (1.98, 2.80), (2.00, 5.95),
            (2.02, 7.15), (2.04, 7.04), (2.06, 6.74), (2.08, 5.82), (2.10, 4.18), (2.12, 2.31),
            (2.14, 1.34), (2.16, 0.68), (2.18, 0.60), (2.20, 0.46)],
        "Tm_3F4-3H6_aluminosilicate": [
            (1.56, 0.26), (1.60, 0.85), (1.64, 1.92), (1.68, 3.24), (1.72, 4.58), (1.76, 5.33),
            (1.80, 5.52), (1.84, 5.46), (1.88, 5.02), (1.92, 4.51), (1.96, 3.86), (2.00, 2.94),
            (2.04, 1.92), (2.08, 1.10), (2.12, 0.63), (2.16, 0.35), (2.20, 0.21)],
        "Er_4I13/2-4I15/2_aluminosilicate": [
            (1.48, 0.97), (1.50, 1.76), (1.52, 5.04), (1.54, 5.04), (1.56, 3.23), (1.58, 1.97),
            (1.60, 0.89)],
        "Er_4I11/2-4I13/2_ZBLAN": [
            (2.62, 0.53), (2.64, 1.08), (2.66, 1.58), (2.68, 2.19), (2.70, 3.34), (2.72, 5.61),
            (2.74, 5.61), (2.76, 4.70), (2.78, 4.66), (2.80, 4.49), (2.82, 3.85), (2.84, 2.44),
            (2.86, 1.56), (2.88, 1.20), (2.90, 0.83), (2.92, 0.59), (2.94, 0.49)],
        "Ho_5I6-5I7_ZBLAN": [
            (2.70, 0.42), (2.74, 0.96), (2.78, 1.69), (2.80, 3.44), (2.82, 4.55), (2.84, 5.12),
            (2.86, 5.19), (2.88, 4.87), (2.90, 3.01), (2.92, 2.10), (2.94, 1.79), (2.96, 1.29),
            (2.98, 0.88), (3.00, 0.87), (3.04, 0.60), (3.08, 0.35), (3.12, 0.35), (3.16, 0.16),
            (3.20, 0.26)],
        "Dy_6H13/2-6H15/2_ZBLAN": [
            (2.50, 0.08), (2.60, 0.16), (2.64, 0.30), (2.70, 0.64), (2.74, 1.12), (2.78, 1.94),
            (2.82, 2.69), (2.86, 2.96), (2.88, 2.98), (2.92, 2.77), (2.96, 2.40), (3.00, 1.85),
            (3.04, 1.31), (3.08, 0.95), (3.12, 0.64), (3.16, 0.49)],
    },

    # Fig. 3 p.425: laser transitions and primary lasing wavelengths (energy-level diagram)
    "FIG3_TRANSITIONS_um": {
        "Tm": {"3F4->3H6": 2.0, "3H4->3H5": 2.3, "pump_ESA": "3H5/3F4 -> 1G4 (pump ESA at 0.79 um, blue fluorescence)"},
        "Er": {"4I13/2->4I15/2": 1.55, "4I11/2->4I13/2": 2.75, "4F9/2->4I9/2": 3.45,
               "pump_ESA": "4I13/2->4F9/2 and 4I11/2->4F7/2 (pump ESA)"},
        "Ho": {"5I7->5I8": 2.1, "5I6->5I7": 2.85, "5I5->5I6": 3.95, "5S2->5F5": 3.22},
        "Dy": {"6H13/2->6H15/2": 2.9},
        "_source": "Fig. 3 p.425",
    },

    # ------------------------------------------------- (c) physical mechanisms
    "MECHANISMS": [
        {"ion": "Tm", "host": "silicate", "process": "cross relaxation 3H4+3H6 -> 3F4+3F4 ('two-for-one')",
         "numbers": "resonant in silicate; Tm concentration > 2.5 wt% with cluster-reducing co-dopants "
                    "can nearly double the slope efficiency (ref 44: Jackson, Opt. Commun. 230, 197 (2004))",
         "source": "p.426, Fig. 5a p.427"},
        {"ion": "Tm", "host": "silicate", "process": "ETU 3F4+3F4 -> 3H4+3H6 (reverse of CR)",
         "numbers": "'can occur with a large rate parameter in a clustered system' (no value given)",
         "source": "Fig. 5 caption p.427; p.426 (cluster-reducing co-dopants mitigate gain-lowering ETU)"},
        {"ion": "Tm", "host": "silicate", "process": "pump ESA at 0.79 um populating 1G4 (blue fluorescence)",
         "numbers": "possible photodarkening; 30-40% of absorbed pump becomes heat in present systems",
         "source": "p.429"},
        {"ion": "Tm", "host": "silicate", "process": "in-band pumping alternative",
         "numbers": "low-concentration (<2.5 mol%) Tm pumped at 1.53-1.56 um improves efficiency (refs 126,127)",
         "source": "p.429"},
        {"ion": "Tm->Ho", "host": "silicate", "process": "Tm sensitisation of Ho",
         "numbers": "limited by ETU loss and reversible energy transfer between Tm 3F4 and Ho 5I7 (no rates given)",
         "source": "p.426"},
        {"ion": "Ho", "host": "silicate", "process": "resonant (in-band) pumping of 5I7 by Tm fibre lasers",
         "numbers": "quantum defect < 7%; 'negligible energy-transfer processes between Ho3+ cations' for "
                    "low-concentration (<1 mol%) Ho silicate; long 5I7 lifetime exploited (no value)",
         "source": "p.426, p.429"},
        {"ion": "Er", "host": "ZBLAN", "process": "ETU 4I13/2+4I13/2 -> 4I9/2+4I15/2 recycling lower laser level",
         "numbers": "with high Er concentration and <100 dB/km background loss gives slope 35.6% (> Stokes limit) "
                    "(ref 51); energy-transfer rate parameters may be >= 10x smaller in fibres than in bulk glass (ref 52, "
                    "Gorjan et al., IEEE JQE 47, 262 (2011)); upper laser level 4I11/2 lifetime shorter than lower 4I13/2 "
                    "(self-terminating) -- no lifetime values given",
         "source": "p.426, Fig. 5b p.427"},
        {"ion": "Er", "host": "ZBLAN", "process": "pump ESA from upper laser level",
         "numbers": "causes roll-off in calculated (ref 56) and measured (ref 51) output power",
         "source": "p.426"},
        {"ion": "Er", "host": "ZBLAN", "process": "cascade lasing 2.8 um + 1.5 um (refs 53-55)",
         "numbers": "manages heat of single-transition Er systems", "source": "p.426"},
        {"ion": "Ho", "host": "ZBLAN", "process": "5I6->5I7 laser diode-pumped at 1.15 um",
         "numbers": "reduced pump ESA and higher Stokes limit vs Er; cascade 2.9 um + 2.1 um raises power (ref 26)",
         "source": "p.426"},
        {"ion": "Ho/Pr", "host": "ZBLAN", "process": "Pr3+ deactivation co-doping",
         "numbers": "only implicit: Table B1 row Ho3+,Pr3+ 2.86 um (ref 25) and Q-switched Ho3+Pr3+ (ref 101); "
                    "no deactivation rates or Pr concentrations given in this review",
         "source": "Table B1 p.424; p.428"},
        {"ion": "Er", "host": "tellurite", "process": "energy transfer Er -> OH- impurities",
         "numbers": "together with low radiative efficiency of upper level suppresses 3 um lasing", "source": "p.428"},
        {"ion": "all RE", "host": "silica", "process": "clustering in inter-network (percolation) regions",
         "numbers": "Al2O3/P2O5 co-doping improves solubility (raises Rayleigh loss)", "source": "p.425, p.427"},
        {"ion": "general", "host": "-", "process": "phonon broadening homogeneous; site-to-site field "
         "variation inhomogeneous & T-independent; nephelauxetic red shift with covalency", "source": "p.425-426"},
    ],

    # ------------------------------------------------- (d) host glasses / phonons
    "HOSTS": {
        "silica": {
            "max_phonon_energy_cm-1": 1100, "src_phonon": "p.427 (ref 60, Layne et al. PRB 16, 10 (1977))",
            "bandgap_eV": 9, "rayleigh_loss_min_dB_per_km": 0.15, "damage_threshold_MW_per_cm2": 500,
            "background_loss_dilute_RE_dB_per_km": "<2 (for <0.1 mol% RE)", "src": "p.423, p.427",
            "longest_laser_um": 2.188, "note": "lasing beyond 2.188 um not expected in most silicates (p.427)",
            "gain_bandwidth_THz_up_to": 50,
        },
        "ZBLAN": {
            "composition_mol_pct": {"ZrF4": 53, "BaF2": 20, "LaF3": 4, "AlF3": 3, "NaF": 20},
            "max_phonon_energy_cm-1": 565, "src_phonon": "p.427 (ref 69, Almeida & Mackenzie, JCP 74, 5954 (1981))",
            "multiphonon_rule": "RT fluorescence for energy gaps > 2825 cm^-1 (lambda < 3.5 um) (~5 phonons)",
            "refractive_index": 1.49,
            "transmission_window_um": (0.2, 4.5), "transmission_window_def": "attenuation < 200 dB/km",
            "min_loss_dB_per_km": {"2.59um": 0.65, "2.9um(OH band)": "<20"},
            "OH_stretch_um": 2.9,
            "damage_threshold_MW_per_cm2": 25, "damage_condition": "10 ms pulses at 2.8 um (as printed)",
            "RE_site": "RE3+ substitutes La3+, higher doping without clustering",
            "Er_ZBLAN_fibre_background_loss_dB_per_km": "<100 (high-quality double-clad)",
            "src": "p.426-428",
        },
        "germanate": {"max_phonon_energy_cm-1": 900, "src": "p.428 (ref 74)"},
        "sulphide": {"max_phonon_energy_cm-1": (350, 425), "src": "p.428 (ref 80)",
                     "note": "fluorescence up to 4.3 um measured in sulphide (Dy:GLS, ref 84)"},
        "selenide": {"max_phonon_energy_cm-1": (250, 300), "src": "p.428 (ref 81)"},
        "telluride": {"max_phonon_energy_cm-1": (150, 200), "src": "p.428 (ref 82)",
                      "transmission_up_to_um": 15},
        "general_statement": "longest fibre-laser wavelength always shorter than max transmissible wavelength "
                             "because the phonon density of states sets radiative efficiencies (p.424)",
    },
}


# =============================================================================
# LE GOUET 2020
# =============================================================================

LEGOUET2020 = {
    "citation": "J. Le Gouet et al., arXiv:2003.14188v1 (2020), 'Realization and simulation of high power "
                "holmium doped fiber laser for long-range transmission'",

    "host": {
        "glass": "Al-codoped silica (MCVD + solution doping with water-based solution of Ho and Al ions; "
                 "chlorination drying). Al/Ho ratio and mol% NOT given.",
        "OH_routine_ppb": "<50 (porous silica layer)",
        "source": "Sec. 2.4 p.5",
    },

    # ---------------------------------------------------------- fibres
    "fibres": {
        "fiber1_SMF": {"type": "iXblue HDF-PM-8-130 (PM, single-clad)", "core_diam_um": 8, "NA": 0.17,
                       "Ho_conc_m-3": 3.7e25, "Ho_conc_err_m-3": 0.6e25, "Ho_conc_cm-3": 3.7e19,
                       "conc_method": "SEM-WDS", "abs_peak_dB_per_m": 57, "abs_peak_err": 1,
                       "abs_peak_nm": 1946, "mode_diam_um_at_1946": 9.4, "overlap_at_1946": 0.76,
                       "source": "Sec. 2.1 p.2-3"},
        "fiber2_SMF": {"type": "HDF-PM-8-130, higher doping", "core_diam_um": 8, "NA": 0.21,
                       "Ho_conc_m-3": 27e25, "Ho_conc_err_m-3": 1e25, "Ho_conc_cm-3": 2.7e20,
                       "abs_dB_per_m_at_1950": 425, "abs_err": 5, "overlap_at_1946": 0.85,
                       "source": "Sec. 2.1 p.3"},
        "laser_fibre_3CF_LMA": {"type": "IXF-3CF-Ho-O-20-250-300-0.08 triple-clad", "core_diam_um": 20,
                                "core_NA": 0.08, "pedestal_Ge_diam_um": 64,
                                "pump_clad_um": 250, "pump_clad_shape": "octagonal, flat-to-flat",
                                "F_doped_clad_diam_um": 300, "pump_NA": 0.22, "coating_OD_um": 450,
                                "Ho_conc_m-3": 5.3e25, "Ho_conc_cm-3": 5.3e19,
                                "conc_method": "estimated from Ho oxide quantity in preform",
                                "core_abs_peak_dB_per_m": 86, "core_abs_err": 3, "abs_peak_nm": 1946,
                                "overlap": 0.84, "abs_peak_calc_dB_per_m": 87,
                                "source": "Sec. 3.1 p.7-8; Table 1 p.10"},
    },

    # ---------------------------------------------------------- cross sections
    "cross_sections": {
        "sigma_abs_peak_m2": 4.3e-25, "sigma_abs_peak_err_m2": 0.2e-25, "abs_peak_nm": 1946,
        "sigma_abs_peak_fiber1_m2": 4.8e-25, "sigma_abs_peak_fiber1_err_m2": 0.8e-25,
        "sigma_em_peak_m2": 6.0e-25, "sigma_em_peak_err_m2": 0.4e-25, "em_peak_nm_digitized": 1950,
        "ratio_em_over_abs_peak": 1.40, "ratio_err": 0.06,
        "method": "absolute: fluorescence shape (pumped at 1150 nm into 5I6, L<=2 cm) scaled by small-signal gain "
                  "58+/-3 dB/m at 1995 nm (1.16 dB over 20 mm); S_peak/S(1995)=1.37. NOT McCumber.",
        "at_1908nm_m2": {"sigma_abs": 1.9e-25, "sigma_em": 1.6e-25},
        "McCumber_note": "McCumber (300 K, crossing at 1946 nm) underestimates near peak (~4.5e-25 at 1950, dip ~3.2e-25 "
                         "at ~1980, ~4.05e-25 at ~2015 nm -- read from Fig.1, DIGITIZED) but overlaps measured "
                         "emission beyond ~2100 nm. Zero-phonon absorption linewidth ~60 nm vs kT-equivalent ~80 nm.",
        "literature_sigma_abs_1950_m2": {"Simakov2013[18]": 2.9e-25, "Wang2015[19]": 4.65e-25,
                                         "Simakov2016[20]": 5.8e-25, "Ryabochkina2015[21]": 8.5e-25},
        "abs_ratio_1950_over_1150": "~2 (from Kurkov 2010 [17]); alpha(1150) ~29 dB/m for fiber #1",
        "5I8-5I6_transition_nm": "~1150 (from Jackson 2007 [23])",
        "source": "Sec. 2.1-2.2 p.2-4; Sec. 2.5 p.6; Fig. 1",
    },
    # DIGITIZED from Fig. 1 (left), p.3: (lambda_nm, sigma_abs_1e-25m2, sigma_em_1e-25m2). Measured curves.
    "XS_CURVE_DIGITIZED_nm_abs_em_1e-25m2": [
        (1760, 0.07, 0.01), (1770, 0.08, 0.01), (1780, 0.12, 0.03), (1790, 0.17, 0.03),
        (1800, 0.23, 0.07), (1810, 0.32, 0.08), (1820, 0.41, 0.14), (1830, 0.54, 0.19),
        (1840, 0.70, 0.25), (1850, 0.87, 0.36), (1860, 1.03, 0.49), (1870, 1.24, 0.59),
        (1880, 1.39, 0.88), (1890, 1.55, 1.11), (1900, 1.69, 1.29), (1910, 1.99, 1.72),
        (1920, 2.53, 2.35), (1930, 3.35, 3.99), (1940, 4.12, 5.43), (1946, 4.30, 5.97),
        (1950, 4.25, 6.01), (1960, 3.55, 5.37), (1970, 2.52, 4.31), (1980, 2.11, 4.07),
        (1990, 1.93, 4.35), (2000, 1.91, 4.69), (2010, 1.81, 4.97), (2020, 1.62, 4.97),
        (2030, 1.40, 4.80), (2040, 1.19, 4.48), (2050, 1.01, 4.10), (2060, 0.84, 3.78),
        (2070, 0.71, 3.45), (2080, 0.59, 3.16), (2090, 0.46, 2.71), (2100, 0.36, 2.19),
        (2110, 0.25, 1.67), (2120, 0.17, 1.29), (2130, 0.11, 1.01), (2140, 0.07, 0.69),
        (2150, 0.05, 0.51), (2160, 0.03, 0.41), (2170, 0.01, 0.29), (2180, 0.01, 0.17),
        (2190, 0.01, 0.12),
    ],
    "XS_CURVE_note": "DIGITIZED from Fig.1 left (1750-2200 nm plotted); +/-0.1e-25 m^2, +/-2 nm. Fig. 1 plots only "
                     "to 2200 nm. Emission sigma at 2120 nm (laser line) ~1.3e-25 m^2, absorption ~0.17e-25 m^2.",

    # ---------------------------------------------------------- lifetimes
    "lifetimes": {
        "5I7_single_exp_measured_ms": 0.7, "5I7_single_exp_err_ms": 0.1,
        "5I7_fit_window": "first part of decay (text says 'until about 4 us' -- Fig. 2 shows 4 ms; likely typo)",
        "5I7_Forster_tau0_ms": 1.9, "5I7_Forster_gamma_s-0.5": 34,
        "Forster_form": "P(t) = P0 + A*exp(-gamma*sqrt(t) - t/tau0)  (Eq. 5)",
        "5I7_used_in_model_ms": 0.7,
        "measurement": "fiber #1, 10 cm, pumped into 5I6 by 1150 nm Yb fibre laser (<=1 W), 100 Hz AOM modulation, "
                       "detection response ~10 us",
        "agreement": "0.7 ms 'agrees with the value generally reported' -- Bowman et al. 2013 [24], Hanna et al. 1989 [33]",
        "5I6_ms": None, "5I6_note": "not reported",
        "source": "Sec. 2.3 p.4-5; Fig. 2; Table 1 p.10",
    },

    # ---------------------------------------------------------- energy transfer / PIQ
    "energy_transfer": {
        "ETU_coefficient": None,
        "ETU_note": "Non-exponential decay attributed to dipolar energy transfer between adjacent Ho3+ "
                    "(ref [31] Gomes et al. 2016); Forster gamma = 34 s^-1/2 is the only ET parameter. "
                    "ETU not included in rate equations (effective lifetime used).",
        "PIQ_pair_fraction_2k_fiber1": 0.15, "PIQ_2k_err": 0.01,
        "PIQ_alpha_small_signal_1908_dB_per_m": 25.2, "PIQ_alpha_small_signal_err": 0.3,
        "PIQ_alpha_nonsat_1908_dB_per_m": 1.3, "PIQ_alpha_nonsat_err": 0.1,
        "PIQ_formula": "alpha_ns = 2k * alpha_Ho * [1 - (sa+se)/(2sa+se)]  (Eq. 6, from Myslinski 1997 [36])",
        "PIQ_method": "transmission vs input power at 1908 nm (Tm fibre laser probe), cut-back (Boivin 2010 [35])",
        "PIQ_2k_used_in_model_laser_fibre": 0.23,
        "PIQ_2k_scaling_assumption": "2k proportional to peak absorption: 15% @57 dB/m -> 23% @86 dB/m",
        "PIQ_literature_2k": "30% inferred by Wang et al. 2018 [39] for Hemming 2013 fibre [6] (70 dB/m peak abs)",
        "Psat_formula": "Psat = h nu A_core / ((sa+se) Gamma tau0)",
        "Psat_pump_W": 0.040, "Psat_signal_W": 0.260,
        "PIQ_model_treatment": "constant non-saturable loss; added to core signal loss; for clad pump weighted by A_core/A_clad",
        "source": "Sec. 2.5 p.6; Sec. 3.3 p.9-10",
    },

    # ---------------------------------------------------------- background loss
    "background_loss": {
        "OH_peak_1380nm_dB_per_km": 2.5, "OH_conc_ppm": 0.04,
        "OH_coeff_dB_per_km_per_ppm_1380": 60, "OH_coeff_source": "Humbach et al. 1996 [14]",
        "OH_contribution_at_laser_dB_per_km": "<1",
        "OH_literature_ppm": {"Simakov2013[18]": 1, "Friebele2014[7]": 0.14},
        "zeta_core_signal_2120_dB_per_km_model": 40,
        "zeta_clad_pump_1940_dB_per_km_model": 4,
        "zeta_note": "multiphonon (SiO4) contributions used as FIT parameters to match output signal / residual pump; "
                     "'same order as' Izawa 1977 [40], Garrett & Todd 1982 [41]",
        "source": "Sec. 2.4 p.5 (Fig. 3, 233 m cut-back); Sec. 3.3 p.10; Table 1",
    },

    "model": {
        "type": "two-level homogeneous 5I8-5I7 rate + propagation equations, effective overlap (Giles & Desurvire [15])",
        "Table1": {"core_diam_um": 20, "clad_diam_um": 250, "NA_core": 0.08, "tau0_ms": 0.7,
                   "alpha_Ho_peak_1946_dB_per_m": 86, "2k": 0.23,
                   "zeta_clad_pump_dB_per_km": 4, "zeta_core_signal_dB_per_km": 40},
        "optimal_length_m": 16,
        "parametric_study": "Fig. 7: max slope efficiency vs [OH] 0.01-10 ppm for 2k=20/25/30%. "
                            "1->0.1 ppm OH gains <3%; same gain from 2k 30->25%. DIGITIZED (2k=20%): 0.51@0.01, "
                            "0.505@0.1, 0.48@1, 0.435@3, 0.32@10 ppm; (2k=30%): 0.45, 0.45, 0.43, 0.39, 0.295.",
        "replication_Hemming2013": "with 2k=30%: Ps=183 W at Pp=400 W, slope 47% (agrees with exp. [6])",
        "source": "Sec. 3.3, Sec. 4 p.9-11; Table 1; Fig. 7",
    },

    # ---------------------------------------------------------- laser results
    "laser_results": [
        {"config": "clad-pumped oscillator, FBG HR (~99%, 2 nm 3-dB BW) + flat-cleave Fresnel (~3%) output coupler",
         "fibre_length_m": 16, "laser_nm": 2120, "linewidth_nm": "~1 (FBG-defined)",
         "pump_nm": 1940, "pump": "2x Tm fibre lasers, up to 120 W each, single-mode, <1 nm linewidth, via 2x1 "
                                  "combiner 105/125/0.22; 93+/-1% transmission to active fibre",
         "P_out_max_W": 90, "P_pump_coupled_max_W": 200, "slope_eff_vs_coupled_pct": 49,
         "threshold_coupled_W_extrapolated": 13,
         "threshold_note": "not stated; linear fit of Fig. 6 points (DIGITIZED), slope 0.492",
         "mode": "CW", "cooling": "small fan",
         "fig6_signal_digitized_Wpump_Wsig": [(25.7, 6.8), (50.8, 18.3), (76.3, 30.7), (101.6, 42.8),
                                              (125.3, 54.7), (144.9, 65.2), (164.4, 75.0), (180.2, 82.1),
                                              (195.3, 89.5)],
         "fig6_residual_pump_digitized_Wpump_W": [(15.6, 2.2), (25.8, 4.3), (50.8, 8.6), (76.0, 12.8),
                                                  (101.3, 17.0), (125.1, 21.4), (164.5, 29.1), (180.4, 31.6),
                                                  (195.0, 33.9)],
         "note": "'Up to a pump power of 400 W, the laser efficiency was slightly lower than 50%' (Sec. 3.2)",
         "source": "Abstract; Sec. 3.1-3.2 p.7-9; Fig. 6"},
        {"config": "same, shorter fibre", "fibre_length_m": 10, "laser_nm": 2120, "pump_nm": 1940,
         "P_out_W": 70, "P_pump_coupled_W": 170, "slope_eff_vs_coupled_pct": 43, "mode": "CW",
         "source": "Sec. 3.2 p.9"},
    ],
    "comparison_laser": {
        "Hemming2013[6]": "monolithic clad-pumped Ho silica, ~400 W at ~2120 nm for ~1 kW pump, ~50% efficiency; "
                          "Nufern fibre clad 112 um (Simakov 2013 [18])",
        "source": "Sec. 1, Sec. 3.2",
    },
}
