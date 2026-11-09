"""Reference data taken from two papers on rare-earth-doped glass and fibre.

[W92] L. Wetenkamp, G.F. West and H. Toebben, "Optical properties of rare earth-doped
      ZBLAN glasses", J. Non-Cryst. Solids 140 (1992) 35-40, North-Holland.
      Journal pages 35-40 = PDF pages 1-6. The PDF is an image scan with an OCR text
      layer. Every number below was checked against the page images (300 dpi).

[P95] H.M. Pask, R.J. Carman, D.C. Hanna, A.C. Tropper, C.J. Mackechnie, P.R. Barber
      and J.M. Dawes, "Ytterbium-doped silica fiber lasers: versatile sources for the
      1-1.2 um region", IEEE J. Sel. Top. Quantum Electron. 1(1) (April 1995) 2-13.
      Journal page N = PDF page N-1 (12 PDF pages, journal pp. 2-13). This is a scan with
      no text layer. Text values were read by eye from 200/400-dpi renderings and
      cross-checked against tesseract OCR. OCR errors were corrected by eye, for example
      "~840 y1's" -> "~840 us", "0.43£0.1 pm?" -> "0.43+-0.1 pm^2",
      "5504100 ppm" -> "550+-100 ppm".

Conventions
-----------
* Units are in the key names. A cross section in "1e-24m2" equals 1e-20 cm^2,
  which equals 1 pm^2.
* Each value or block carries a "cite" string giving the paper, journal page,
  table/figure/equation, and where useful the PDF page.
* Values marked DIGITIZED were read off figures by pixel analysis of 400-dpi renderings,
  after locating the axis frame. Typical uncertainty is about +-0.02e-24 m^2 (about one
  line width) on the Pask cross-section curves away from the 975 nm line. Near the 975 nm
  line (965-985 nm) the curves are nearly vertical and the absorption (dashed) and
  emission (solid) traces overlap, so values there are uncertain by about +-0.4e-24 m^2.
  Fig. 2 critical-power values are good to about +-10 % in log space. Energy levels read
  from the W92 Fig. 1 diagram are good to about +-300 cm^-1, because the bands are drawn as
  shaded regions.

Caveats - Wetenkamp 1992
------------------------
1. The paper does NOT give the refractive index, density, or ion number densities. It gives
   only the mol% composition and the dopant mol%. Absorption was normalised to dB/cm/mol%.
2. The authors flag no ion's JO fit as unreliable. They say every fit has a small rms
   deviation delta (Table 1). They do say that for Pr3+ and Dy3+ the IR radiative
   transitions are "very strongly influenced by multiphonon transitions" (p.38). The Pr3+
   ZBLAN set (Omega2 = 0.84e-20 cm^2) should be used with care: Pr3+ JO fits are well known
   to be anomalous, and the ZBLA literature value in the same table has Omega2 = 0.06e-20
   cm^2 with delta = 221e-8. This warning is ours, not the authors'.
3. Table 2 has some internal inconsistencies (values reproduced as printed):
   a. Ho 5I6: beta*eta is printed as 6 % (->5I7) and 5 % (->5I8). From beta and
      tau_m/tau_c = 3.50/5.87 = 0.596 you would get about 5 % and about 54 %. The ->5I8
      value looks like a misprint.
   b. Ho 5S2,5F4 -> 5I7: beta*eta is printed as 27 %. 0.34*0.30/0.34 gives 30 %.
   c. Er 4I11/2: A = 21 + 86 s^-1 gives beta = 0.20/0.80 and tau_rad = 9.3 ms. The table
      prints beta = 0.13/0.87 and tau_c = 8.69 ms. Note that 8.69 ms is the same number as
      the measured 4I13/2 lifetime.
   d. Er 4S3/2: the listed A values sum to 1157 s^-1 (0.86 ms), but tau_c is printed as
      0.70 ms. The text notes that 4S3/2 is thermally coupled to 2H11/2, which may account
      for the difference (not stated).
   e. Er 4I13/2: A = 100 s^-1 (10.0 ms), but tau_c is printed as 9.52 ms. A may be the ED
      part only, with MD omitted. Not stated.
4. The captions of Figs. 2 and 3 appear SWAPPED in the printed paper. The plot above the
   "Fig. 2" caption (p.38) compares LaF3, ZBLA, ZBLAN, YAG and SiO2. The plot above the
   "Fig. 3" caption (p.40) shows the measured ZBLAN multiphonon points with the fitted
   line. The text on p.38 cites "Fig. 2" as the regression of the ZBLAN data.
5. Text misprints on p.39: "5I7 ground state" should be 5I8, and "4S3/2 -> 5I9/2" should be
   4I9/2.
6. The multiphonon rates were measured on 0.1 mol% Ho- and Er-doped ZBLAN FIBRES. The
   tau_m values in Table 2 are those low-concentration values. The text says
   (p.40) that at higher concentration the Er 4S3/2, 4F9/2 and Ho 5S2, 5F4, 5F5 lifetimes
   can fall by up to 20 % through energy transfer.

Caveats - Pask 1995
-------------------
1. Only a GERMANOSILICATE host is reported: Fig. 1(b) and all experimental fibres (MCVD
   germanosilicate cores). No aluminosilicate data appear in this paper.
2. The Stark energies in Fig. 1(a) are labelled "approximate" by the authors. They are not
   explicitly tied to a host, but the paper is about germanosilicate fibres.
3. The absorption spectrum was SCALED so that sigma_a equals sigma_e at the 975 nm peak.
   The authors note that a McCumber analysis suggests sigma_a(peak) may really be 5-10 %
   larger than sigma_e(peak) (p.6).
4. The cross sections are tabulated only as figure curves. Every curve value in this
   module is DIGITIZED. The peak wavelength reads as 974 +-1 nm on the digitized axis.
   The text says 975 nm.
5. The Fig. 6 caption says "pump wavelength 840 nm", but the text (p.7) says the same
   experiment used a 1047 nm pump (dichroic mirror transmitting 1047 nm). The Fig. 10
   caption gives a fibre length of 100 m; the text (p.9) says 90 m. Both are reproduced
   as printed.
6. The "Yb concentration ~550 ppm" in the experiments is an EFFECTIVE concentration
   (concentration times pump/dopant overlap) from cutback measurements, 550 +-100 ppm.
   The paper does not say whether ppm is by weight or by mole.
"""

# ---------------------------------------------------------------------------
# Wetenkamp, West & Toebben 1992 (ZBLAN)
# ---------------------------------------------------------------------------
WETENKAMP1992 = {
    "citation": ("L. Wetenkamp, G.F. West, H. Toebben, J. Non-Cryst. Solids 140 (1992) "
                 "35-40"),
    "glass": {
        "composition_mol_pct": {"ZrF4": 52, "BaF2": 19, "AlF3": 4, "NaF": 20,
                                "LaF3": "5 - x", "REF3": "x"},
        "composition_note": "RE fluoride substitutes for LaF3 (52ZrF4-19BaF2-4AlF3-20NaF-(5-x)LaF3 + x REF3)",
        "dopant_x_mol_pct": {
            "Pr": [3], "Dy": [3], "Tm": [4], "Nd": [0.1, 0.5, 4],
            "Er": [0.1, 0.5, 1, 3, 4], "Ho": [0.1, 0.5, 1, 3, 4],
        },
        "preparation": ("high-purity dried fluorides, Pt crucibles, dry N2; cast, annealed "
                        "at Tg, polished with water-free lubricant"),
        "refractive_index": None,   # NOT GIVEN in paper
        "density_g_cm3": None,      # NOT GIVEN in paper
        "ion_density_cm-3": None,   # NOT GIVEN; absorption normalised as dB/cm/mol%
        "cite": "W92 p.35, sec. 2.1 (PDF p.1)",
    },
    "measurement": {
        "absorption_spectrometer": "Cary 2300, 0.2-3 um; alpha = 10*D/(d*C) in dB/cm/mol%",
        "fluorescence_range_um": (0.5, 3.0),
        "lifetime_excitation_nm": 488,
        "lifetime_accuracy_pct": 3,     # "better than +-3%"
        "temperature": "room temperature",
        "JO_method": "Shinn et al., Phys. Rev. B 27 (1983) 6635 [ref 5]",
        "cite": "W92 p.36, sec. 2.2 (PDF p.2)",
    },
    # Table 1. Omega in units of 1e-20 cm^2 (table prints x1e-22 cm^2; e.g. 84e-22 = 0.84e-20).
    # delta = rms deviation of the fit, units 1e-8 (oscillator strength).
    "judd_ofelt_ZBLAN_this_work": {
        "units": "Omega in 1e-20 cm^2 (Table 1 prints integers in 1e-22 cm^2); delta in 1e-8",
        "Pr3+": {"Omega2": 0.84, "Omega4": 4.79, "Omega6": 9.13, "delta_1e-8": 9},
        "Nd3+": {"Omega2": 2.10, "Omega4": 3.71, "Omega6": 4.62, "delta_1e-8": 3.5},
        "Dy3+": {"Omega2": 1.86, "Omega4": 1.42, "Omega6": 2.37, "delta_1e-8": 2},
        "Ho3+": {"Omega2": 2.30, "Omega4": 2.30, "Omega6": 1.71, "delta_1e-8": 5},
        "Er3+": {"Omega2": 2.91, "Omega4": 1.78, "Omega6": 1.00, "delta_1e-8": 3},
        "Tm3+": {"Omega2": 2.70, "Omega4": 1.59, "Omega6": 1.09, "delta_1e-8": 5},
        "reliability_flag_by_authors": None,  # none flagged; see module caveat 2
        "cite": "W92 p.37, Table 1 (PDF p.3)",
    },
    "judd_ofelt_literature_in_table1": {
        "units": "Omega in 1e-20 cm^2; delta in 1e-8 (None = not given)",
        "Pr3+_ZBLA_ref9":   {"Omega2": 0.06, "Omega4": 5.05, "Omega6": 6.92, "delta_1e-8": 221},
        "Nd3+_ZBN_ref10":   {"Omega2": 1.95, "Omega4": 3.65, "Omega6": 4.17, "delta_1e-8": 1.5},
        "Nd3+_ZBLA_ref11":  {"Omega2": 1.10, "Omega4": 3.80, "Omega6": 5.53, "delta_1e-8": None},
        "Dy3+_ZBLA_ref12":  {"Omega2": 3.22, "Omega4": 1.35, "Omega6": 2.38, "delta_1e-8": 8},
        "Dy3+_ZBLALi_ref13": {"Omega2": 2.70, "Omega4": 1.80, "Omega6": 2.00, "delta_1e-8": 11},
        "Ho3+_ZBLA_ref14":  {"Omega2": 2.28, "Omega4": 2.08, "Omega6": 1.73, "delta_1e-8": 23},
        "Er3+_ZBLAN_ref15": {"Omega2": 2.97, "Omega4": 1.23, "Omega6": 1.17, "delta_1e-8": None},
        "Er3+_ZBLA_ref5":   {"Omega2": 2.54, "Omega4": 1.39, "Omega6": 0.96, "delta_1e-8": 14},
        "Tm3+_ZBLA_ref16":  {"Omega2": 2.80, "Omega4": 1.90, "Omega6": 1.00, "delta_1e-8": 30},
        "cite": "W92 p.37, Table 1 (PDF p.3)",
    },
    # Table 2: all values as printed. Columns: lambda_nm (emission peak), A_s-1 (JO radiative
    # rate), beta (JO branching ratio), dlambda_nm (normalised fluorescence FWHM),
    # sigma_p (peak induced-emission cross section, 1e-20 cm^2), tau_c_ms (JO radiative
    # lifetime of the initial level), tau_m_ms (measured), beta_eta_pct (beta * tau_m/tau_c).
    # None = "-" in the table.
    "table2_emission_ZBLAN": {
        "columns": ["lambda_nm", "A_s-1", "beta", "dlambda_nm", "sigma_p_1e-20cm2",
                    "tau_c_ms", "tau_m_ms", "beta_eta_pct"],
        "Nd3+": {
            "4F3/2->4I13/2": [1320, 227, 0.10, 56, 0.74, 0.44, 0.46, 10],
            "4F3/2->4I11/2": [1050, 1135, 0.50, 32, 2.64, 0.44, 0.46, 50],
            "4F3/2->4I9/2":  [885, 915, 0.40, 44, 0.78, 0.44, 0.46, 40],
        },
        "Ho3+": {
            "5I7->5I8":       [2030, 79, 1.00, 145, 0.56, 12.65, 12.01, 95],
            "5I6->5I7":       [2850, 21, 0.09, 88, 0.96, 5.87, 3.50, 6],
            "5I6->5I8":       [1149, 150, 0.91, 48, 0.47, 5.87, 3.50, 5],   # sic; see caveat 3a
            "5F5->5I8":       [650, 1536, 0.78, 34, 0.86, 0.50, 0.29, 45],
            "5S2,5F4->5I5":   [1355, 59, 0.02, 42, 0.28, 0.34, 0.30, 2],
            "5S2,5F4->5I6":   [1010, 189, 0.06, 29, 0.40, 0.34, 0.30, 5],
            "5S2,5F4->5I7":   [740, 892, 0.34, 21, 0.75, 0.34, 0.30, 27],  # sic; caveat 3b
            "5S2,5F4->5I8":   [545, 1794, 0.61, 24, 0.39, 0.34, 0.30, 54],
        },
        "Er3+": {
            "4I13/2->4I15/2": [1540, 100, 1.00, 82, 0.42, 9.52, 8.69, 91],
            "4I11/2->4I13/2": [2720, 21, 0.13, 123, 0.57, 8.69, 6.71, 10],  # caveat 3c
            "4I11/2->4I15/2": [980, 86, 0.87, 38, 0.13, 8.69, 6.71, 67],
            "4S3/2->4I9/2":   [1660, 46, 0.04, None, None, 0.70, 0.57, None],
            "4S3/2->4I11/2":  [1131, 25, 0.02, None, None, 0.70, 0.57, None],
            "4S3/2->4I13/2":  [850, 316, 0.27, 31, 0.31, 0.70, 0.57, 22],
            "4S3/2->4I15/2":  [540, 770, 0.67, None, None, 0.70, 0.57, None],
        },
        "Tm3+": {
            "3F4->3H6": [1660, 152, 1.00, 177, 0.57, 6.56, 5.64, 86],
            "3H4->3H5": [2270, 16, 0.02, None, None, 1.32, None, None],
            "3H4->3F4": [1500, 55, 0.08, None, None, 1.32, None, None],
            "3H4->3H6": [805, 665, 0.90, 36, 0.46, 1.32, None, None],
        },
        "cite": "W92 p.37, Table 2 (PDF p.3); tau_m from 0.1 mol% Ho/Er fibres (p.38)",
    },
    # Lifetimes collected per level (from Table 2)
    "level_lifetimes_ms": {
        "Nd3+ 4F3/2":     {"tau_JO": 0.44, "tau_meas": 0.46},
        "Ho3+ 5I7":       {"tau_JO": 12.65, "tau_meas": 12.01},
        "Ho3+ 5I6":       {"tau_JO": 5.87, "tau_meas": 3.50},
        "Ho3+ 5F5":       {"tau_JO": 0.50, "tau_meas": 0.29},
        "Ho3+ 5S2,5F4":   {"tau_JO": 0.34, "tau_meas": 0.30},
        "Er3+ 4I13/2":    {"tau_JO": 9.52, "tau_meas": 8.69},
        "Er3+ 4I11/2":    {"tau_JO": 8.69, "tau_meas": 6.71},
        "Er3+ 4S3/2":     {"tau_JO": 0.70, "tau_meas": 0.57},
        "Tm3+ 3F4":       {"tau_JO": 6.56, "tau_meas": 5.64},
        "Tm3+ 3H4":       {"tau_JO": 1.32, "tau_meas": None},
        "note": ("text p.39 quotes 'considerably longer lifetime (11 ms) of the 5I7'. "
                 "Lifetimes of Ho 5I6, 5I7 and Er 4I11/2, 4I13/2 do not change up to 4 mol%. "
                 "Ho 5F4/5S2 and Er 4S3/2 drop sharply already at 0.5 mol% (cross relaxation)."),
        "cite": "W92 p.37 Table 2; p.39 text (PDF pp.3,5)",
    },
    "multiphonon": {
        "law": "W_MP = C * exp(-alpha * dE)",
        "C_s-1": 1.99e5,
        "alpha_cm": 0.0021,
        "dE_definition": "energy gap to next-lower level, cm^-1",
        "derived_from": "W_MP = 1/tau_m - 1/tau_c, 0.1 mol% Ho- and Er-doped ZBLAN fibres, 300 K",
        "claim": "MP rates for gaps < 3200 cm^-1 are lower in ZBLAN (300 K) than in ZBLA (77 K), SiO2, YAG",
        "cite": "W92 p.38 text (PDF p.4); fit shown in figure captioned Fig. 3 on p.40 (captions swapped, caveat 4)",
        # DIGITIZED from the ZBLAN data plot (printed above the 'Fig. 3' caption, p.40, PDF p.6).
        # Each entry is (level, dE_cm-1, W_MP_s-1). Where Table 2 gives both lifetimes, the
        # value of 1/tau_m - 1/tau_c computed from Table 2 is shown alongside.
        "data_points_DIGITIZED": [
            ("Er 4I9/2",          2140, 1700, None),
            ("Ho 5F5",            2280, 1530, 1448),
            ("Er 4F9/2",          2840, 460, None),
            ("Ho 5F4,5S2",        3030, 400, 392),
            ("Er 4S3/2,2H11/2",   3100, 345, 326),
            ("Ho 5I6",            3540, 115, 115),
            ("Er 4I11/2",         3740, 43, 34),
            ("Ho 5I7",            5120, 4.5, 4.2),
        ],
        "data_points_columns": ["level", "dE_cm-1", "W_MP_s-1_digitized", "W_MP_s-1_from_table2"],
        "host_comparison_lines": ("figure printed above 'Fig. 2' caption (p.38) shows straight "
                                  "lines for 1 LaF3, 2 ZBLA, 3 ZBLAN, 4 YAG, 5 SiO2 over "
                                  "dE ~2000-6000 cm^-1; the line parameters for the other hosts "
                                  "are not given"),
    },
    "esa_gsa_oscillator_strengths_1e-8": {
        "Ho3+ 530nm ESA 5I7->5G5": 280,
        "Ho3+ 530nm GSA 5I8->5F4": 300,
        "Ho3+ 640nm ESA ->5F3": 153,
        "Ho3+ 640nm GSA 5I8->5F5": 275,
        "Ho3+ 750nm GSA ->5F4": 2,     # printed '5I7 ground state' (sic)
        "Ho3+ 750nm ESA ->5F4": 70,
        "Er3+ 640nm ESA 4I13/2->4F5/2": 108,
        "Er3+ 790nm ESA 4I13/2->2H11/2": 30,
        "Er3+ 836nm ESA 4I13/2->4S3/2": 42,
        "Tm3+ 750nm ESA 3H5->1G4": 134,
        "cite": "W92 p.39 text (PDF p.5)",
    },
    "cross_relaxation_noted": {
        "Ho3+": ["5F4/5S2->5I4 + 5I8->5I7", "5F4/5S2->5I7 + 5I8->5I4"],
        "Er3+": ["4S3/2->4I9/2 + 4I15/2->4I13/2", "4S3/2->4I13/2 + 4I15/2->4I9/2"],
        "cite": "W92 p.39 (PDF p.5)",
    },
    "absorption_band_positions": {
        "note": ("No table of absorption band positions is given. Spectra are in ref [6] "
                 "(Wetenkamp PhD thesis). Fig. 1 (p.36) is a level diagram. Shaded bands "
                 "mark absorption > 0.5 dB/cm/mol%. Values below are DIGITIZED band "
                 "centres from Fig. 1, +-300 cm^-1."),
        "units": "cm^-1",
        "Pr3+": {"3H5": 2150, "3H6": 4300, "3F2": 5200, "3F3": 6350, "3F4": 7100,
                 "1G4": 9850, "1D2": 17100, "3P0": 20800, "3P1": 21500, "3P2": 22700},
        "Nd3+": {"4I11/2": 2000, "4I13/2": 4050, "4I15/2": 5900, "4F3/2": 11650,
                 "4F5/2": 12650, "4F7/2": 13450, "4F9/2": 14750, "2H11/2": 16050,
                 "4G5/2": 17400, "4G7/2": 19000, "4G9/2": 19650, "2P1/2": 23500,
                 "2D5/2": 23950, "2P3/2": 26300},
        "Dy3+": {"6H13/2": 3650, "6H11/2": 5950, "6F11/2": 7800, "6F9/2": 9150,
                 "6H5/2": 10350, "6F7/2": 11100, "6F5/2": 12550, "6F3/2": 13300,
                 "6F1/2": 13850, "4F9/2": 21200, "4I15/2": 22200, "4G11/2": 23500},
        "Ho3+": {"5I7": 5200, "5I6": 8750, "5I5": 11300, "5I4": 13300, "5F5": 15600,
                 "5F4,5S2": 18600, "5F3": 20600, "5F2,3K8": 21200, "5G6": 22300,
                 "5G5": 24150, "5G4": 26050, "3H6": 27700},
        "Er3+": {"4I13/2": 6650, "4I11/2": 10400, "4I9/2": 12650, "4F9/2": 15400,
                 "4S3/2": 18550, "2H11/2": 19300, "4F7/2": 20600, "4F5/2": 22300,
                 "4F3/2": 22750, "2H9/2": 24750, "4G11/2": 26400},
        "Tm3+": {"3F4": 5900, "3H5": 8550, "3H4": 12800, "3F3": 14650, "3F2": 15200,
                 "1G4": 21600, "1D2": 28000},
        "cite": "W92 p.36, Fig. 1 (PDF p.2), DIGITIZED",
    },
    "transition_wavelengths_discussed_um": {
        "Pr3+ 1G4": [3.4, 2.9, 2.1, 1.7, 1.3],
        "Dy3+ 6H11, 6H13/2 (three-level)": [1.7, 2.8],
        "Ho3+": {"5I7": [2.02], "5I6": [1.15, 2.85], "5S2,5F4": [0.55, 0.74, 1.0, 1.35]},
        "Er3+": {"4I13/2": [1.5], "4I11/2": [1.0, 2.7], "4S3/2": [0.54, 0.85, 1.01, 1.35]},
        "Tm3+": {"3H4": [0.8, 1.5, 2.3], "3F4": [1.8], "1G4->3F4_nm": 644, "1G4->3H6_nm": 480},
        "cite": "W92 pp.38-40 (PDF pp.4-6)",
    },
}

# ---------------------------------------------------------------------------
# Pask et al. 1995 (Yb3+ in germanosilicate fibre)
# ---------------------------------------------------------------------------

# Fig. 1(b): germanosilicate host. DIGITIZED (400 dpi, PDF p.2 / journal p.3).
# Units 1e-24 m^2 (= 1e-20 cm^2 = pm^2).
# None = below digitization resolution (< ~0.02).
_PASK_XS_WL_NM = [850, 860, 870, 880, 890, 900, 910, 912, 920, 930, 940, 950, 960,
                  970, 974, 980, 985, 990, 1000, 1010, 1020, 1030, 1035, 1040, 1050,
                  1060, 1070, 1080, 1090, 1100, 1110, 1120, 1130, 1140, 1150]
_PASK_SIGMA_ABS = [0.07, 0.14, 0.22, 0.38, 0.60, 0.75, 0.82, 0.83, 0.71, 0.47, 0.29, 0.21,
                   0.23, 1.0, 2.48, 0.7, 0.15, 0.09, 0.045, 0.045, 0.045, 0.04, 0.035,
                   0.025, 0.02, 0.01, None, None, None, None, None, None, None, None, None]
_PASK_SIGMA_EM = [0.02, 0.02, 0.02, 0.02, 0.02, 0.03, 0.035, 0.04, 0.045, 0.055, 0.045,
                  0.06, 0.11, 1.45, 2.48, 0.9, 0.29, 0.20, 0.19, 0.26, 0.415, 0.57, 0.59,
                  0.57, 0.51, 0.38, 0.305, 0.29, 0.26, 0.225, 0.205, 0.177, 0.146, 0.12,
                  0.095]

PASK1995 = {
    "citation": ("H.M. Pask, R.J. Carman, D.C. Hanna, A.C. Tropper, C.J. Mackechnie, "
                 "P.R. Barber, J.M. Dawes, IEEE JSTQE 1(1) (1995) 2-13"),
    "host": "germanosilicate (MCVD) silica fibre; no aluminosilicate data in this paper",
    "stark_levels_cm-1": {
        "2F7/2": {"a": 0, "b": 600, "c": 1060, "d": 1490},
        "2F5/2": {"e": 10260, "f": 11000, "g": 11630},
        "note": "'Approximate energies in wavenumbers above ground energy' (caption)",
        "cite": "P95 p.3, Fig. 1(a) (PDF p.2)",
    },
    "spectral_features": {
        "A_975nm": "a<->e zero-phonon line; narrow peak in both absorption and emission",
        "B_~910nm": "absorption a->f,g (digitized peak 912 nm, 0.83e-24 m^2)",
        "C_shoulder_>1000nm": "absorption from b (~6% of level-a population at RT, Boltzmann)",
        "D_~1035nm": "emission e->b,c,d; tail to 1200 nm (digitized peak 1035 nm, 0.59e-24 m^2)",
        "E": "weak emission from f (thermal population); no lasing observed on it",
        "absorption_band_nm": (800, 1064),   # abstract: '~800 nm to ~1064 nm'
        "emission_band_nm": (970, 1200),     # abstract: '~970 nm to ~1200 nm'
        "cite": "P95 p.2 abstract; p.3 sec. II (PDF pp.1-2)",
    },
    "cross_sections": {
        "units": "1e-24 m^2 (= 1e-20 cm^2 = pm^2)",
        "peak_975_sigma_abs": 2.5,   # scaled equal to emission at peak (method)
        "peak_975_sigma_em": 2.5,    # digitized max 2.48 at 974 +-1 nm
        "peak_975_wavelength_nm": 975,
        "peak_975_FWHM_nm_DIGITIZED": 10,   # emission trace, half-max crossings ~969 and ~979 nm
        "abs_peak_B_nm_DIGITIZED": 912, "abs_peak_B_sigma_DIGITIZED": 0.83,
        "em_peak_D_nm_DIGITIZED": 1035, "em_peak_D_sigma_DIGITIZED": 0.59,
        "sum_at_1017nm_from_spectra": 0.39,
        "sum_at_1017nm_measured_pump_saturation": 0.43, "sum_at_1017nm_measured_err": 0.1,
        "mccumber_note": ("absorption scaled to equal emission at the 975 nm peak; McCumber "
                          "analysis suggests sigma_a(peak) may be 5-10% larger; profiles "
                          "agree with McCumber sigma_a = sigma_e exp(eps - h nu) using "
                          "eps = 975 nm peak energy"),
        "methods": ("emission: sidelight fluorescence (840 nm pump, integrating sphere) "
                    "scaled by lambda^5 and normalised via 1/tau = (8 pi n^2/c^2) "
                    "int nu^2 sigma_em dnu (eq. 10); absorption: white-light cutback"),
        "cite": "P95 p.6 sec. IV (PDF p.5); Fig. 1(b) p.3",
        "curve_DIGITIZED": {
            "wavelength_nm": _PASK_XS_WL_NM,
            "sigma_abs_1e-24m2": _PASK_SIGMA_ABS,
            "sigma_em_1e-24m2": _PASK_SIGMA_EM,
            "notes": ("DIGITIZED from Fig. 1(b) (germanosilicate). +-0.02 away from the "
                      "975 line. Values at 970 and 980 nm are only good to +-0.4, because the "
                      "traces are near-vertical and overlap. Values of 0.02 or less are at "
                      "the baseline/axis-line resolution. Absorption beyond ~1060 nm is "
                      "below resolution. Emission keeps falling to ~0.03 at 1180-1190 nm."),
            "extra_emission_points_nm_sigma": [(1160, 0.067), (1170, 0.059), (1180, 0.031),
                                               (1190, 0.034)],
            "cite": "P95 p.3, Fig. 1(b) (PDF p.2), DIGITIZED",
        },
    },
    "upper_level_lifetime": {
        "tau_us": 840,
        "notes": ("2F5/2 sidelight fluorescence decay, single exponential over the first three "
                  "e-folds, taken as the radiative lifetime (gap >10000 cm^-1 vs max phonon "
                  "~1100 cm^-1 in silica). The pump-saturation measurement at 1017 nm supports "
                  "a predominantly radiative lifetime."),
        "max_phonon_energy_silica_cm-1": 1100,
        "cite": "P95 p.6 sec. IV (PDF p.5)",
    },
    "model_equations": {
        "critical_pump_power": "P_cr = (A h nu_p / tau) / (sigma_el*sigma_ap/sigma_al - sigma_ep)",
        "pump_saturation_power": "P_s = h nu_p A / ((sigma_ep + sigma_ap) tau phi_p)  (eq. 6)",
        "pump_propagation": "ln(Pp(z)/Pp(0)) + (Pp(z)-Pp(0))/Ps = -N sigma_ap z  (eq. 5)",
        "gain_dB": "4.34*[phi_p (sigma_el+sigma_al) tau P_a/(A h nu_p) - N sigma_al l]  (eq. 9)",
        "phi_p": 1.0,
        "cladding_pumping": "replace sigma_ap, sigma_ep by sigma/x (x = inner-cladding/core area ratio)",
        "assumptions": "uniform intensity over core area A; uniform doping; steady state; "
                       "homogeneous broadening; negligible background loss",
        "cite": "P95 pp.3-6, eqs. (1)-(9) (PDF pp.2-5)",
    },
    # Fig. 2: critical pump power (transparency) vs signal wavelength. Core diameter 3.75 um
    # (right axis gives Pcr/A in mW/um^2; A ~ 11 um^2). DIGITIZED (PDF p.3 / journal p.4).
    "critical_power_fig2_DIGITIZED": {
        "core_diameter_um": 3.75,
        "units": "mW (left axis); divide by A ~= 11.0 um^2 for intensity in mW/um^2",
        "pump_840nm": {960: 290, 980: 114, 1000: 46, 1010: 30, 1020: 17, 1030: 8.6,
                       1040: 4.1, 1050: 2.4, 1060: 1.9, 1070: 1.5, 1080: 0.97,
                       1090: 0.50, 1100: 0.26, 1110: 0.15, 1120: 0.10},
        "pump_900nm": {920: 179, 940: 35, 960: 11, 980: 4.5, 990: 2.4, 1000: 1.8,
                       1010: 1.2, 1020: 0.67, 1030: 0.34, 1040: 0.16, 1050: 0.10},
        "pump_975nm": {980: 35, 990: 1.4, 1000: 0.85, 1010: 0.46, 1020: 0.22, 1030: 0.10},
        "pump_1047nm": {1050: 42, 1060: 12, 1070: 6.8, 1080: 3.0, 1090: 1.24, 1100: 0.58,
                        1110: 0.33, 1120: 0.21, 1130: 0.13},
        "pump_1064nm": {1070: 52, 1080: 10, 1090: 3.4, 1100: 1.5, 1110: 0.81, 1120: 0.50,
                        1130: 0.32, 1140: 0.19},
        "notes": ("The 975 nm pump curve exists only for lambda_l > ~977 nm; it is near-vertical "
                  "at 980 nm (very uncertain). Each curve falls to 0.1 mW at roughly 1030 (975), "
                  "1050 (900), 1120 (840), 1135 (1047) and 1150 nm (1064)."),
        "cite": "P95 p.4, Fig. 2 (PDF p.3), DIGITIZED",
    },
    "fibres": {
        "core_pumped_standard": {"host": "germanosilicate", "NA": 0.17, "Yb_ppm": 550,
                                 "Yb_ppm_err": 100, "core_diameter_um": [3.0, 3.75],
                                 "background_loss_dB_per_km_max": 10,
                                 "launch_efficiency_pct": {"Ti:sapphire": 50, "Nd:YLF": 70},
                                 "pump_abs_at_1064nm_dB_per_m": 0.1,
                                 "cite": "P95 pp.5-7, sec. III-V (PDF pp.4-6)"},
        "cladding_pumped_expt": {"inner_core_diameter_um": 4.25,
                                 "inner_cladding_diameter_um": 12.75,
                                 "NA_inner": 0.16, "NA_outer": 0.15, "Yb_ppm": 500,
                                 "doping": "solution doping, MCVD germanosilicate",
                                 "cite": "P95 p.10 sec. VII-A (PDF p.9)"},
        "cladding_pumped_model": {"inner_core_diameter_um": 6.8, "outer_diameter_um": 120,
                                  "NA_inner": 0.13, "NA_outer": 0.35, "V_at_1.14um": 2.405,
                                  "area_ratio_x": 311, "Yb_ppm": 1000, "length_m": 15,
                                  "pump_nm": 975,
                                  "cite": "P95 p.11 sec. VII-B, Fig. 12 (PDF p.10)"},
        "single_frequency": {"Yb_ppm": 2000, "NA": 0.28, "cutoff_nm": 1000,
                             "cite": "P95 pp.11-12 sec. VIII (PDF pp.10-11)"},
    },
    "laser_results": {
        "summary": {"pump_range_nm": (840, 1064), "free_running_range_nm": (1035, 1115),
                    "slope_eff_max_vs_absorbed_pct": 90,
                    "demonstrated_overall_range_nm": (975, 1180),
                    "grating_tuning_pct": 1,
                    "cite": "P95 p.6 sec. V; p.5; p.12 sec. IX (PDF pp.5,4,11)"},
        "fig4_threshold_vs_length_840nm_pump": {
            "fibre": "3.75 um, NA 0.17, 550 ppm",
            "model_min_threshold_mW_DIGITIZED": 30, "model_opt_length_m_DIGITIZED": 6.5,
            "meas_min_threshold_mW_DIGITIZED": 45, "meas_at_length_m_DIGITIZED": 5,
            "meas_threshold_28m_mW_DIGITIZED": 63, "meas_threshold_90m_mW_DIGITIZED": 88,
            "short_fibre_lasing_nm": 974,  # L < ~1.1 m, thresholds 240-350 mW (digitized)
            "lasing_nm_vs_length_DIGITIZED": {"1-5 m": "1037-1049", "28 m": 1058, "90 m": 1098},
            "model_vs_meas": "calculated thresholds within 30% of measured",
            "cite": "P95 p.7, Fig. 4(a),(b) (PDF p.6); text p.6",
        },
        "fig5_free_running_850nm_pump": {
            "lengths_m": [8, 28, 90], "feedback": "Fresnel, both cleaved ends",
            "slope_eff_vs_launched_pct": (70, 80),
            "lambda_nm_long_fibre": 1090,
            "with_HR_output_mirror_lambda_nm": {"8m": (1038, 1050), "28m": (1058, 1060),
                                                "90m": (1090, 1106)},
            "cite": "P95 p.7, Fig. 5 and text (PDF p.6)",
        },
        "fig6_1102nm": {"length_m": 90, "pump_nm": 1047, "threshold_mW_launched": 30,
                        "slope_eff_vs_launched_pct_min": 90,
                        "mirror": "launch end R>98% at 1100 nm, transmits 1047 nm",
                        "note": "Fig. 6 caption says pump 840 nm; text says 1047 nm",
                        "cite": "P95 pp.7-8, Fig. 6 (PDF pp.6-7)"},
        "high_power": {"P_out_W": 2, "lambda_nm": 1090, "pump": "CW lamp-pumped Nd:YAG 1064 nm",
                       "length_m": 90, "cite": "P95 p.7 (PDF p.6)"},
        "1020nm_gratings": {"length_m": 10, "core_um": 3.0, "pump_nm": 840,
                            "R_in_pct": 95, "R_out_pct": 30,
                            "output_coupling_loss_roundtrip_dB": 5.4,
                            "splice_loss_each_pct": 5, "splice_loss_each_dB": 0.22,
                            "single_pass_gain_threshold_dB": 3.1,
                            "free_running_1050nm_single_pass_gain_dB": 15,
                            "discrimination_dB": 12,
                            "threshold_mW_launched": 55,
                            "P_out_mW_at_510mW_launched": 240,
                            "unabsorbed_pump_mW_at_510mW": 30,
                            "slope_eff_vs_absorbed_pct": 60,
                            "quantum_limit_pct": 82,
                            "model_note": "Fig. 7: 3.0 um fibre, length optimised for 1020 nm gain",
                            "discrimination_for_1020_vs_1040_at_50mW_dB": 5,
                            "cite": "P95 pp.5,8, Figs. 7-8 (PDF pp.4,7)"},
        "pr_zblan_pump_application": {"length_m": 5, "core_um": 3.0,
                                      "outputs_mW": {"491nm": 6.5, "520nm": 18, "635nm": 55},
                                      "combined_incident_pump_mW": 380,
                                      "cite": "P95 p.9 (PDF p.8)"},
        "1140nm_gratings": {"length_m": 90, "core_um": 3.75, "pump_nm": 1047,
                            "pump_source": "diode-pumped Nd:YLF",
                            "R_grating_pct": 50, "R_mirror_pct": 99,
                            "discrimination_vs_1106nm_dB_per_pass": 4,
                            "threshold_mW_launched": 6,
                            "P_out_mW_at_500mW_launched": 330,
                            "slope_eff_vs_launched_pct": 66, "quantum_limit_pct": 92,
                            "note": "Fig. 10 caption gives length 100 m",
                            "application": "pumped Tm:ZBLAN upconversion laser, 30 mW at 480 nm",
                            "cite": "P95 p.9, Figs. 9-10 (PDF p.8)"},
        "cladding_pumped_1042nm": {"length_m": 5, "pump": "Ti:sapphire ~975 nm (MOPA simulation), x5 objective, ~12 um spot",
                                   "threshold_mW_incident": 75,
                                   "slope_eff_vs_incident_pct": 80,
                                   "launch_efficiency_pct": 90,
                                   "objective_transmission_pct": 95,
                                   "quantum_limit_pct": 93.8,
                                   "cite": "P95 p.10, Fig. 11 (PDF p.9)"},
        "cladding_pumped_1020nm": {"length_m": 4, "R_launch_pct": 62, "R_output_pct": 37,
                                   "threshold_mW_incident": 32,
                                   "slope_eff_output_end_pct": 42,
                                   "overall_conversion_both_ends_pct": 70,
                                   "cite": "P95 p.10 (PDF p.9)"},
        "cladding_pumped_diode_model": {"threshold_mW": 10, "cavity_loss_dB": 3,
                                        "slope_eff_vs_launched_pct_min": 70,
                                        "P_out_mW_at_1100nm": 350,
                                        "launch_efficiency_assumed_pct": 50,
                                        "diode": "SDL-6360 at 975 nm",
                                        "coverable_range_nm": (1050, 1150),
                                        "cite": "P95 p.11, Fig. 12 (PDF p.10)"},
        "single_frequency_1083nm": {"cavity_length_cm": 5.6, "R_launch_pct": 92,
                                    "R_output_pct": 60, "grating_bw_nm": 0.25,
                                    "pump_nm": 970, "threshold_mW_incident": 12,
                                    "laser_FSR_GHz": 1.8, "slope_eff_vs_launched_pct": 6,
                                    "cite": "P95 p.12, Fig. 13 (PDF p.11)"},
        "core_pumped_overall_conversion_pct": (40, 60),
    },
    "three_vs_four_level": {
        "975nm": "truly three-level (emission into lowest Stark level a)",
        "long_wavelength": ("approaches four-level at longer wavelengths (terminal levels c, d "
                            "nearly empty); level b population ~6% of a at RT"),
        "quantum_defect_vs_NdYAG": "about 3x smaller (Yb:YAG vs Nd:YAG pumped at 800 nm, lasing at 1064 nm)",
        "ASE_limit_single_pass_gain_dB": (30, 40),
        "transparency": "see critical_power_fig2_DIGITIZED; P_cr is independent of Yb concentration and fibre length",
        "example": "for lasing at 975 nm a 900 nm pump is preferable to an 840 nm pump (Fig. 2)",
        "cite": "P95 pp.2-5, sec. II (PDF pp.1-4)",
    },
}
