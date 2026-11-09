"""Reference data extracted from

    C. R. Giles and E. Desurvire, "Modeling erbium-doped fiber amplifiers,"
    J. Lightwave Technol., vol. 9, no. 2, pp. 271-283, Feb. 1991.
    doi:10.1109/50.65886  (Invited paper, AT&T Bell Laboratories)

Source: the 13-page scanned PDF (Acrobat Capture OCR text layer). Numbers in
tables and text were checked against 250-400 dpi page renders. Curves were
digitized from 400-dpi renders: axes calibrated from tick-mark pixel positions,
curve centres taken from dark-pixel runs per column, assigned to curves by hand,
and checked by overlaying the points on the scan.

READ FIRST - WHAT IS MEASURED AND WHAT IS MODELLED
* This is a modelling/review paper. The ONLY measured data in it are the
  absorption alpha(lambda) and full-inversion gain g*(lambda) spectra of two
  fibres in Fig. 2 (white-light + monochromator loss; gain with strong 980 or
  528 nm pump). The cross-section axis in Fig. 2 is computed with the
  Ladenburg-Fuchtbauer equation; the authors themselves say it is "not very
  accurate".
* Every gain-vs-pump, gain-vs-length, saturation, ASE-spectrum, efficiency and
  n_sp curve (Figs. 4-12) is a CALCULATION by the authors, using a spatial
  model, a spectral (overlap-integral) model, or the Saleh analytical model.
  Use them as model-to-model benchmarks, not as experimental validation.
  The paper reports no measured noise figure and no measured output saturation
  power.
* Fibre parameters in Tables I/II are model inputs. The paper says Table I is
  "similar to" a reported high-efficiency fibre (Zyskind et al., ref. [41]).
  The fibres in Fig. 2 are not given core radius, NA, Er concentration,
  lifetime or length.

KNOWN INCONSISTENCIES IN THE PAPER (copied as printed, not corrected)
* Fig. 8 is said in the text (p. 278) to use "the same fiber used in the
  spatial model to obtain Fig. 4(b)" with zeta = 1.5e15 m^-1 s^-1. But the
  Fig. 8 caption gives alpha_p=3, g*_p=1.1, alpha_s=3.3, g*_s=4.8 dB/m and
  zeta=3e15 m^-1 s^-1. Table I gives alpha_p=1.6, g*_p=0.5, alpha_s=2.6,
  g*_s=3.6 dB/m. The caption also swaps the descriptions of panels (b) and (c):
  in the figure, (b) is gain vs length and (c) is gain vs output power.
* Fig. 7: the inset gives the gains at 1531 and 1555 nm (Table II also says
  signal #2 = 1555 nm), but the text says 1531 and 1550 nm, with 21.6 dB and
  21.9 dB filtered small-signal gains. The inset says 21.4 and 21.2 dB.
* Fig. 4 text: "for P_p > 5 mW, there is no substantial difference in the
  small-signal gain ... 0.6 um <= b <= 1.8 um", yet the digitized Fig. 4(c)
  (b = 1.8 um) is about 3.5 dB below (a)/(b) at 10 mW.

DIGITIZATION UNCERTAINTY (1 sigma, rough)
* Fig. 2: about +/-0.3 nm in wavelength. About +/-0.08 dB/m (Fig. 2a) and
  +/-0.04 dB/m (Fig. 2b) where the curves are smooth, and about 2x that on the
  steep 1525-1538 nm flanks. The Fig. 2a absorption peak touches the top frame
  (10 dB/m), so it may be clipped.
* Fig. 4/5/8 gain curves: +/-0.4 dB on smooth parts, +/-1.5 dB on the steep
  low-pump rise. Pump-axis values are good to about +/-0.05 mW.
* Fig. 6/7 ASE levels: +/-0.7 dB. Curves are dense and labels overlap them.
* Fig. 8(c) y-axis: the tick marks and the dashed 0-dB line disagree by about
  0.5 dB near 0 dB; the line is used as zero, so allow +/-0.6 dB.
* Figs. 9 and 10 are 3-D mesh plots; only rough values are given.

Unit conventions: key suffixes give the units. dB/m is the natural-log
coefficient times 4.343. "digitized": True marks values read off figures;
anything else is copied from text or tables. "derived" marks values computed
here from the paper's own numbers.
"""

GILES1991 = {
    "citation": ("C. R. Giles and E. Desurvire, 'Modeling erbium-doped fiber amplifiers', "
                 "J. Lightwave Technol. 9(2), 271-283 (1991), doi:10.1109/50.65886"),

    # ------------------------------------------------------------------
    # 1. FIBRES
    # ------------------------------------------------------------------
    "fibres": {
        # ---- Table I (p. 275): spatial-model fibre (Figs. 4, 5); also used for Fig. 8 per text ----
        "tableI_spatial_model_fibre": {
            "cite": "Table I, p. 275 (model input; 'similar to' high-efficiency fibre of ref. [41])",
            "host": "Al-silicate (Fig. 4 caption: 'Al-silicate, Er3+-doped fiber amplifier')",
            "index_core_radius_um": 1.2,
            "index_step_delta_n": 0.035,
            "er_core_radius_b_um": 1.2,        # uniform Er disk; Fig. 4 also uses b = 0.6 and 1.8 um
            "amplifier_length_m": 18.8,         # Fig. 4
            "metastable_lifetime_ms": 10.0,
            "ase_bandwidth_nm": {"fig4": 10.0, "fig5": 0.0},   # text p.274: 1250 GHz (10 nm) noise BW
            "pump_wavelength_nm": 1480.0,
            "signal_wavelength_nm": 1550.0,
            "alpha_pump_dB_per_m": 1.6,
            "gstar_pump_dB_per_m": 0.5,
            "alpha_signal_dB_per_m": 2.6,
            "gstar_signal_dB_per_m": 3.6,
            "excess_loss_dB_per_m": {"fig4": 0.03, "fig5": 0.0},   # both pump and signal
            "sigma_abs_signal_m2": 3e-25,
            "pump_input_mW": {"fig4": (0.0, 10.0), "fig5": (0.0, 5.0)},
            "signal_input": {"fig4_nW": 100.0, "fig5_uW": 50.0},
            "overlap_pump_Gamma_p": 0.62,       # text p.274, for a = b = 1.2 um
            "overlap_signal_Gamma_s": 0.58,     # text p.274
            "zeta_m_inv_s_inv": 1.5e15,         # text p.278 ('same fiber ... Fig. 4(b)')
            "derived": {
                "NA_approx": 0.319,  # sqrt(2*n*dn) with n~1.45; NOT stated in paper
                "Er_density_m3": 3.44e24,  # alpha_s/(sigma_as*Gamma_s) = 0.5987/(3e-25*0.58)
                "zeta_check_m_inv_s_inv": 1.56e15,  # pi b^2 n_t / tau, consistent with 1.5e15
                "note": "derived here, not printed in the paper",
            },
        },

        # ---- Table II (p. 277): spectral-model fibre (Figs. 6, 7) ----
        "tableII_spectral_model_fibre": {
            "cite": "Table II, p. 277 (model input)",
            "host": ("alumino-silicate; text p.276: 'higher pump threshold than the fiber used in the spatial "
                     "models'. alpha(lambda), g*(lambda) spectra used are not stated (probably Al-type like Fig. 2b; unconfirmed)"),
            "zeta_m_inv_s_inv": 4.2e15,           # common to Fig. 6 and Fig. 7
            "fig6": {"length_m": 10.0, "pump_nm": 1480.0, "pump_mW": 50.0,
                     "signal1_nm": 1545.0, "signal2_nm": 1555.0,
                     "signal_input_each_dBm": [-20.0, -10.0],    # text p.276
                     "ase_resolution_GHz": 125.0, "ase_resolution_nm": 1.0},
            "fig7": {"length_m": 5.0, "pump_nm": 980.0, "pump_mW": 100.0,
                     "signal1_nm": 1531.0, "signal2_nm": 1555.0,   # text says 1550; see caveats
                     "filter": {"type": "absorptive Lorentzian notch, mid-amplifier",
                                "centre_nm": 1531.0, "fwhm_3dB_nm": 15.0,
                                "min_transmission_fraction": 0.054}},
        },

        # ---- Fig. 8 caption (p. 278): analytical (Saleh) model parameters ----
        "fig8_caption_params": {
            "cite": "Fig. 8 caption, p. 278 (conflicts with Table I / text; see docstring)",
            "pump_nm": 1480.0, "signal_nm": 1550.0,
            "alpha_pump_dB_per_m": 3.0, "gstar_pump_dB_per_m": 1.1,
            "alpha_signal_dB_per_m": 3.3, "gstar_signal_dB_per_m": 4.8,
            "zeta_m_inv_s_inv": 3e15,
            "text_says_zeta_m_inv_s_inv": 1.5e15,
            "signal_input_nW_fig8a": 100.0,  # text p.278
            "fig8c_fig8d_length_m": 15.0,    # text p.278 ('15-m long amplifier')
        },

        # ---- Fig. 2 fibres (p. 272): MEASURED spectra; no other parameters given ----
        "fig2a_Ge_silicate": {
            "cite": "Fig. 2(a), p. 272 - MEASURED (white-light loss; gain at full inversion w/ 980 or 528 nm pump)",
            "codopants": "Ge only",
            "sigma_axis_m2_per_dB_per_m": 0.717e-25,  # digitized right-axis scale (LF-computed)
            "derived_Gamma_nt_m3": 3.21e24,             # (ln10/10)/scale; implied by axis scaling
        },
        "fig2b_Al_Ge_silicate": {
            "cite": "Fig. 2(b), p. 272 - MEASURED",
            "codopants": "Al + Ge",
            "sigma_axis_m2_per_dB_per_m": 1.747e-25,
            "derived_Gamma_nt_m3": 1.318e24,
        },
    },

    # ------------------------------------------------------------------
    # Giles parameter spectra (MEASURED), Fig. 2. Solid = absorption alpha, dashed = gain g*.
    # Range covered by the figure: 1470-1570 nm only (nothing at 1450-1470 or 1570-1620).
    # sigma (m^2) = value_dB_per_m * sigma_axis_m2_per_dB_per_m (LF-scaled, low accuracy).
    # ------------------------------------------------------------------
    "spectra_fig2": {
        "digitized": True,
        "cite": "Fig. 2, p. 272",
        "uncertainty": "+/-0.3 nm; +/-0.08 dB/m (a), +/-0.04 dB/m (b); x2 on 1525-1538 nm flanks",
        "wavelength_nm": [1471, 1475, 1480, 1485, 1490, 1495, 1500, 1505, 1510, 1515, 1520, 1522, 1524,
                          1526, 1528, 1530, 1531, 1532, 1533, 1534, 1535, 1536, 1537, 1538, 1540, 1542,
                          1544, 1546, 1548, 1550, 1552, 1554, 1556, 1558, 1560, 1562, 1564, 1566, 1568,
                          1569],
        "a_Ge_silicate": {
            "alpha_dB_per_m": [2.00, 2.35, 2.64, 2.88, 3.10, 3.33, 3.58, 3.85, 4.16, 4.61, 5.21, 5.58, 6.10,
                               6.85, 7.72, 8.60, 9.04, 9.52, 9.93, 10.0, 9.82, 9.23, 7.89, 6.82, 5.5, 4.98,
                               4.82, 4.64, 4.57, 4.44, 4.27, 3.99, 3.68, 3.29, 2.85, 2.45, 2.08, 1.77, 1.48,
                               1.33],
            "gstar_dB_per_m": [0.47, 0.58, 0.73, 0.87, 1.05, 1.21, 1.43, 1.74, 2.10, 2.59, 3.22, 3.52, 4.01,
                               5.01, 6.42, 7.66, 8.3, 8.74, 9.13, 9.19, 8.81, 8.22, 7.37, 6.82, 5.5, 4.98,
                               4.82, 5.05, 5.22, 5.32, 5.30, 5.14, 4.88, 4.55, 4.19, 3.82, 3.47, 3.17, 2.75,
                               2.61],
            "notes": ("alpha peak ~10 dB/m at 1533.5-1534 nm (touches frame: possibly clipped); g* peak ~9.2 dB/m "
                      "at 1534 nm; curves cross ~1538 and ~1544 nm; g* secondary peak 5.3 dB/m at ~1550 nm. "
                      "1540-1544 nm values are both curves overlapping."),
        },
        "b_Al_Ge_silicate": {
            "alpha_dB_per_m": [1.00, 1.06, 1.16, 1.29, 1.38, 1.46, 1.53, 1.60, 1.70, 1.86, 2.17, 2.36, 2.59,
                               2.94, 3.32, 3.51, 3.60, 3.68, 3.75, 3.69, 3.43, 3.20, 2.99, 2.74, 2.40, 2.22,
                               2.09, 2.03, 1.95, 1.88, 1.80, 1.73, 1.66, 1.55, 1.42, 1.24, 1.02, 0.93, 0.89,
                               0.88],
            "gstar_dB_per_m": [0.23, 0.27, 0.33, 0.39, 0.46, 0.55, 0.65, 0.78, 0.93, 1.13, 1.42, 1.60, 1.88,
                               2.30, 2.75, 3.13, 3.29, 3.37, 3.40, 3.30, 3.07, 2.81, 2.60, 2.39, 2.17, 2.05,
                               2.05, 2.03, 2.06, 2.07, 2.06, 2.03, 1.96, 1.86, 1.75, 1.62, 1.47, 1.31, 1.16,
                               1.10],
            "notes": "alpha peak 3.75 dB/m at 1533 nm; g* peak 3.40 dB/m at 1533 nm; crossing ~1545-1546 nm.",
        },
    },

    # ------------------------------------------------------------------
    # 2. AMPLIFIER RESULTS - ALL MODEL CALCULATIONS (no measured amplifier data in paper)
    # ------------------------------------------------------------------
    "fig4_spatial_model_gain_vs_pump": {
        "model": True, "digitized": True,
        "cite": "Fig. 4(a)-(c), p. 275; Table I fibre, L=18.8 m, 1480 nm pump, 1550 nm 100 nW signal, ASE 10 nm",
        "uncertainty": "+/-0.4 dB above 3 mW; +/-1.5 dB on the steep rise; pump +/-0.05 mW",
        "axes_note": ("shared y-axis: GAIN dB (-10..40); gain coeff dB/mW = (gain_axis+10)/5; "
                      "P_ASE dBm = gain_axis-30. ASE (-) and (+) taken as backward/forward (u_k = -1/+1)."),
        "b0p6um": {
            "pump_mW": [0.96, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 9.8],
            "signal_gain_dB": [-9.9, -8.7, 0.3, 7.6, 13.4, 17.3, 22.0, 24.5, 26.05, 27.4, 28.9, 29.9, 30.8,
                               31.45, 32.1, 32.5],
            "ase_backward_dBm": {"pump_mW": [4, 6, 8, 9, 9.8], "P_dBm": [-5.1, -2.8, -1.1, -0.6, -0.3]},
            "ase_forward_dBm": {"pump_mW": [8, 9, 9.8], "P_dBm": [-1.6, -1.1, -0.7]},
            "peak_gain_coeff_dB_per_mW": 8.1,   # text p.274 (dotted curve peaks ~3 mW)
        },
        "b1p2um": {
            "pump_mW": [1.27, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 9.8],
            "signal_gain_dB": [-10.0, -5.0, 1.0, 5.8, 10.4, 13.9, 19.0, 22.2, 24.6, 27.1, 28.9, 30.2, 31.1,
                               32.1, 32.6],
            "ase_backward_dBm": {"pump_mW": [6, 7, 8, 9, 9.8], "P_dBm": [-3.3, -2.3, -1.3, -0.5, -0.2]},
            "ase_forward_dBm": {"pump_mW": [6, 7, 8, 9, 9.8], "P_dBm": [-4.2, -3.1, -2.15, -1.3, -0.8]},
            "peak_gain_coeff_dB_per_mW_approx": 6.2,   # digitized dotted-curve peak at ~3 mW
        },
        "b1p8um": {
            "pump_mW": [1.88, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 9.8],
            "signal_gain_dB": [-9.8, -8.0, -0.7, 6.2, 11.5, 14.8, 17.85, 20.1, 23.15, 25.4, 27.0, 28.2, 29.0],
            "ase_backward_dBm": {"pump_mW": [5, 6, 7, 8, 9, 9.8], "P_dBm": [-10.6, -8.2, -6.5, -5.05, -4.0, -3.5]},
            "ase_forward_dBm": {"pump_mW": [5, 6, 7, 8, 9, 9.8], "P_dBm": [-12.8, -9.6, -7.4, -5.85, -4.8, -4.1]},
            "peak_gain_coeff_dB_per_mW_approx": 4.0,
        },
        "pump_curve_note": ("an unlabelled 'PUMP' solid curve is also plotted (gain-axis value -10 at threshold, "
                            "rising to about -3 dB at 10 mW for b=0.6/1.2 um and about -4 dB for 1.8 um). It is "
                            "probably residual pump transmission in dB, but the paper does not say; not tabulated."),
    },

    "fig5_optimum_length_vs_pump": {
        "model": True, "digitized": True,
        "cite": "Fig. 5, p. 276; Table I fibre, no ASE, no excess loss, 50 uW signal at 1550, 1480 pump, length optimized",
        "uncertainty": "+/-0.3 dB gain, +/-0.4 m length",
        "pump_mW": [0.5, 0.6, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0, 4.25],
        "signal_gain_dB": [0.0, 0.08, 0.42, 1.72, 3.16, 4.6, 6.05, 7.5, 8.7, 9.72, 11.6, 13.1, 14.38, 14.75],
        "optimum_length_m": [0.0, 1.1, 2.8, 4.9, 6.9, 8.6, 10.3, 11.6, 13.0, 14.2, 15.9, 17.6, 19.1, 19.8],
        "gain_coeff_dB_per_mW": {"pump_mW": [2.5, 4.0, 4.5, 4.9], "value": [3.89, 3.53, 3.33, 3.18]},
        "notes": ("gain = 0, L_opt = 0 below transparency threshold ~0.5 mW. Gain and length curves leave the plot "
                  "at the top (15 dB / 20 m) above ~4.35 mW. From gain coeff (= G/Pp): G ~15.0 dB at 4.5 mW and "
                  "~15.6 dB at 4.9 mW (derived). Solid 'PUMP' curve (probably pump transmission in dB, unstated): "
                  "0 at 0.5, -3.1 at 1, -6.2 at 2, -8.1 at 3, -9.4 at 4, -10 at 5 mW."),
    },

    "fig6_two_channel_ase_spectra": {
        "model": True, "digitized": True,
        "cite": "Fig. 6, p. 277; Table II (L=10 m, 1480 nm 50 mW pump, signals 1545/1555 nm), 1-nm ASE resolution",
        "signal_gain_dB_printed": {"Pin_-20dBm": {"1545nm": 24.5, "1555nm": 24.3},
                                   "Pin_-10dBm": {"1545nm": 19.9, "1555nm": 20.5}},  # figure labels and text p.276
        "text_notes": "compression of 4.6 dB at 1545 and 3.8 dB at 1555 nm; backward ASE ~50% > forward when saturated",
        "uncertainty": "+/-0.7 dB",
        "wavelength_nm": [1520, 1525, 1528, 1531.5, 1535, 1540, 1550, 1560, 1564, 1566],
        "Pin_-20dBm_backward_dBm": [-28.0, -20.8, -15.7, -13.6, -17.0, -19.3, -18.9, -20.0, -22.4, -24.0],
        "Pin_-20dBm_forward_dBm": [-28.6, -21.6, -16.4, -14.1, -17.4, -19.7, -19.1, -20.0, -22.7, -24.0],
        "Pin_-10dBm_backward_dBm": [-30.4, -24.2, -20.4, -18.7, -21.1, -22.6, -22.0, -22.5, -24.7, -25.9],
        "Pin_-10dBm_forward_dBm": [-32.6, -27.2, -23.7, -21.4, -23.7, -24.4, -23.3, -23.5, -25.4, -26.6],
        "notes": "plot spans 1510-1570 nm, -35..-10 dBm; ASE peak at ~1531.5 nm; spikes at 1545/1555 are signals.",
    },

    "fig7_gain_equalization_filter": {
        "model": True,
        "cite": "Fig. 7 and text p. 277; Table II (L=5 m, 980 nm 100 mW pump)",
        "inset_gain_dB": {"no_filter": {"1531nm": 31.7, "1555nm": 21.3},
                          "filter": {"1531nm": 21.4, "1555nm": 21.2}},
        "text_gain_dB_with_filter": {"1531nm": 21.6, "1550nm": 21.9},
        "text_notes": "unfiltered gain difference > 10 dB; noise spectral density only 1 dB higher at 1531 than 1550 nm with filter",
        "ase_digitized": {
            "digitized": True, "uncertainty": "+/-0.7 dB",
            "wavelength_nm": [1520, 1525, 1528, 1530, 1538, 1545, 1550, 1560, 1566],
            "no_filter_backward_dBm": [-26.9, -19.6, -14.5, -12.9, -20.0, -21.5, -22.4, -24.2, -28.1],
            "no_filter_forward_dBm": [-27.4, -20.5, -15.1, -13.5, -20.6, -21.8, -22.6, -24.5, -28.1],
            "filter_backward_dBm": [-28.8, -23.9, -21.5, -21.4, -24.4, -22.9, -22.5, -24.5, -28.1],
            "filter_forward_dBm": [-29.3, -24.8, -22.5, -22.6, -25.1, -23.2, -22.8, -24.9, -28.1],
            "notes": "solid/dashed assumed forward/backward as in Fig. 6 (caption does not say)",
        },
    },

    "fig8a_gain_vs_pump_lengths": {
        "model": True, "digitized": True,
        "cite": "Fig. 8(a), p. 278; Saleh analytical model, 1480 pump, 1550 signal 100 nW, no ASE",
        "params_note": "use fig8_caption_params, or Table I with zeta = 1.5e15 per the text (conflict)",
        "uncertainty": "+/-0.4 dB smooth parts, +/-1.5 dB where curves cross / are steep",
        "pump_mW": [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.2, 7.0,
                    8.0, 9.0, 9.9],
        "L5m_gain_dB": [-8.75, -4.15, -0.74, 1.25, 2.56, 4.09, 5.11, 5.3, 6.14, 6.48, 6.99, 7.5, 7.84, 8.07,
                        8.41, 8.92, 9.15, 9.43, 9.6, 9.6],
        "L10m_gain_dB": [None, None, -9.5, -1.82, 2.6, 5.68, 7.73, 9.83, 10.7, 12.16, 13.92, 15.0, 15.85,
                         16.42, 16.99, 17.9, 18.41, 18.81, 19.2, 19.32],
        "L15m_gain_dB": [None, None, None, -9.6, -3.07, 2.56, 6.76, 9.83, 13.12, 15.6, 19.2, 21.48, 22.9,
                         23.98, 24.94, 26.5, 27.39, 28.18, 28.69, 28.98],
        "L20m_gain_dB": [None, None, None, None, None, -8.75, -1.08, 5.3, 10.7, 15.8, 21.82, 25.91, 28.58,
                         30.45, None, 33.6, 35.11, 36.14, 36.99, 37.56],
        "notes": ("None = below -10 dB (off-scale) or hidden by a label. -10 dB onset: L5 ~0.2, L10 ~0.7, "
                  "L15 ~1.0, L20 ~1.45 mW. L20 has a gap at 4.6-5.9 mW (label)."),
    },

    "fig8b_gain_vs_length": {
        "model": True, "digitized": True,
        "cite": "Fig. 8(b), p. 278 (caption mislabels as (c)); Pp = 1..10 mW, 100 nW signal assumed",
        "uncertainty": "+/-0.5 dB, +/-0.5 m",
        "peaks": {  # (optimum length m (flat-top range), max gain dB)
            "Pp1mW": {"Lopt_m_range": (3.7, 6.7), "Gmax_dB": 1.5},
            "Pp2mW": {"Lopt_m_range": (12.0, 14.8), "Gmax_dB": 10.5},
            "Pp3mW": {"Lopt_m_range": (19.7, 22.0), "Gmax_dB": 21.5},
            "Pp4mW": {"Lopt_m_range": (26.4, 29.6), "Gmax_dB": 31.4},
            "Pp5mW": {"Lopt_m_range": (30.7, 34.9), "Gmax_dB": 37.6},
            "Pp6to10mW": "peaks above 40 dB (off-scale)",
        },
        "length_m": [5, 10, 15, 20, 25, 30, 35, 40, 45, 49.5],
        "Pp1mW_gain_dB": [1.47, -0.3, -9.5, None, None, None, None, None, None, None],
        "Pp2mW_gain_dB": [5.44, 9.53, 10.38, 5.44, -4.7, None, None, None, None, None],
        "Pp3mW_gain_dB": [7.09, 13.89, 19.05, 21.55, 19.85, 11.5, -0.06, None, None, None],
        "Pp4mW_gain_dB": [8.0, 15.88, 22.8, 27.96, 30.91, 31.13, 27.45, 17.64, 4.99, -7.0],
        "Pp5mW_gain_dB": [9.07, 18.2, 24.84, 31.36, 35.56, 37.43, 37.49, 35.84, None, 23.0],
        "notes": "Pp1 at 10 m is hidden by a label (approx.). G returns to 0 dB at ~9.5 (1 mW), ~23.5 (2 mW), ~35 (3 mW), ~46 m (4 mW).",
    },

    "fig8c_saturation_L15m": {
        "model": True, "digitized": True,
        "cite": "Fig. 8(c), p. 278 (caption mislabels as (b)); L = 15 m, Pp = 1..10 mW",
        "uncertainty": "+/-0.6 dB gain; +/-0.7 dBm output; curve identity for Pp 7-10 mildly uncertain",
        "small_signal_gain_dB": {"Pp1": -9.2, "Pp2": 10.9, "Pp3": 19.4, "Pp4": 23.2, "Pp5": 25.3,
                                 "Pp6": 26.5, "Pp7": 27.4, "Pp8": 28.3, "Pp9": 28.7, "Pp10": 29.05},
        "output_power_at_3dB_compression_dBm": {"Pp2": -6.7, "Pp3": -4.5, "Pp4": -2.7, "Pp5": -1.6,
                                                "Pp6": -0.3, "Pp7": 0.6, "Pp8": 1.0, "Pp9": 1.7, "Pp10": 2.3},
        "Pp1mW_curve": {"Pout_dBm": [-30, -25, -20, -15, -10, -7.5, -5, -2.5],
                        "gain_dB": [-9.2, -8.9, -7.7, -6.1, -4.0, -3.1, -2.2, -1.5]},
        "notes": "text: saturated output proportional to pump; near-complete pump-to-signal conversion at large compression",
    },

    "fig8d_pump_for_constant_gain_L15m": {
        "model": True, "digitized": True,
        "cite": "Fig. 8(d), p. 278; L = 15 m",
        "uncertainty": "+/-0.1 mW (+/-0.2 mW above 0 dBm where labels and axis ticks overlap)",
        "Pout_dBm": [-19, -15, -10, -7.5, -5, -4, -3, -2, 3, 4, 5],
        "G25dB_pump_mW": [5.14, 5.17, 5.33, 5.51, 5.83, 6.02, 6.26, 6.58, 10.0, None, None],
        "G20dB_pump_mW": [3.23, 3.23, 3.31, 3.44, 3.68, 3.83, 4.00, 4.22, 6.61, 7.38, 8.5],
        "G15dB_pump_mW": [2.43, 2.43, 2.51, 2.61, 2.78, 2.89, 3.02, 3.18, 5.14, 5.87, 6.6],
        "G10dB_pump_mW": [1.93, 1.96, 2.00, 2.11, 2.28, 2.38, 2.49, 2.62, 4.2, 4.74, 5.4],
        "G5dB_pump_mW": [1.64, 1.62, 1.66, 1.72, 1.82, 1.89, 1.97, 2.07, 3.1, 3.62, 4.1],
        "notes": "G=25 dB curve reaches 10 mW at ~+2.9 dBm output.",
    },

    "fig9_conversion_efficiency": {
        "model": True, "digitized": "rough (3-D mesh)",
        "cite": "Fig. 9, p. 279; Al:silicate fibre, eq. (30); signal 1520-1560 nm, pump 1460-1510 nm",
        "notes": ("efficiency lowest near 1530 nm signal. Approx. corner values: 3-dB compression "
                  "~0.15 at (s=1560, p=1460) and up to ~0.6 peak; 10-dB compression ~0.45 at (1560, 1460) and up to ~0.9. "
                  "Read off a 3-D projection: +/-0.1."),
    },

    "fig10_nsp_highly_pumped": {
        "model": True, "digitized": "rough (3-D mesh)",
        "cite": "Fig. 10 and text p. 279; eq. (31), uniform full-pump inversion",
        "text_statement": "pump < 1487 nm gives n_sp < 2 over 1520-1560 nm signal band",
        "notes": ("n_sp ~1.1 at long signal / short pump (1560/1460); rises to ~3 at pump 1510 nm and signal "
                  "~1520-1530 nm. 980 nm pumping (sigma_ep = 0) gives n_sp = 1."),
    },

    "overlap_appendix_figs11_12": {
        "model": True,
        "cite": "Appendix, Figs. 11-12, pp. 281-282",
        "notes": ("Gaussian mode i_k=(1/(pi w^2))exp(-r^2/w^2), uniform Er to radius b. Gamma_k,1(p->0)=1-exp(-b^2/w^2). "
                  "Fig. 11: power-dependent Gamma for 0.2<=b/w<=2 vs normalized power 0-100. Constant-Gamma "
                  "error <8% if b<w (Gamma<0.6). Curves not digitized."),
    },

    # ------------------------------------------------------------------
    # Other scalar statements
    # ------------------------------------------------------------------
    "misc_statements": {
        "cite": "text pp. 271-280",
        "two_level_valid_980nm_pump_below_W": 1.0,
        "noise_bandwidth_typical_GHz": (200.0, 2000.0),
        "fig4_noise_bandwidth_GHz": 1250.0,
        "steady_state_valid_mod_freq_above_kHz": 10.0,
        "gain_corner_frequency_typical_kHz": 5.0,
        "Ge_silicate_inhomog_linewidth_nm": 8.0,
        "room_T_homog_linewidth_nm": (3.0, 4.0),
        "literature_records_quoted": {"gain_dB": 51, "sat_output_mW": 100, "noise_figure_dB": 3,
                                      "fiber_laser_pump_eff_pct": 71},
        "energy_levels_fig1": {"4I13/2_to_4I15/2_nm": 1530, "4I11/2_pump_nm": 980,
                               "4I15/2_stark_spread_cm-1": "~400", "homog_width_cm-1": "~10"},
    },

    # ------------------------------------------------------------------
    # 3. MODEL EQUATIONS / ASSUMPTIONS (transcribed)
    # ------------------------------------------------------------------
    "model_equations": {
        "cite": "Sections III-VIII, eqs. (1)-(32), pp. 272-280",
        "assumptions": [
            "homogeneously broadened two-level system (980 and 1480 nm pumping; no ESA at these wavelengths)",
            "steady state (dn2/dt = 0); CW or modulation > 10 kHz",
            "m = 2 modes (two polarizations of LP01) for spontaneous emission",
            "forward and backward ASE resolved in bins dnu_k (spectral model), or a single effective bandwidth (spatial model)",
            "step-index weakly-guiding LP01 with u = 1.1428V - 0.9960 for 1 < V <= 3 (eq. 9)",
            "effective-overlap models assume Er confined to the mode centre (b < w, Gamma < 0.6)",
            "Saleh model: constant Gamma, no ASE saturation, negligible excess loss (accurate for G <~ 20 dB)",
        ],
        "eq1": "alpha(l) = sigma_a(l) Gamma(l) n_t ; g*(l) = sigma_e(l) Gamma(l) n_t",
        "eq2_LF": "sigma_{a,e}(l) = l_peak^4 I_{a,e}(l) / (8 pi c n^2 tau Int I_{a,e} dl)",
        "eq6": ("dn2/dt = sum_k P_k i_k sigma_ak n1/(h nu_k) - sum_k P_k i_k sigma_ek n2/(h nu_k) - n2/tau ; "
                "eq7: n_t = n1 + n2"),
        "eq13_spatial_propagation": ("dP_k/dz = u_k[(alpha_k+g*_k)/Gamma_k P_k Int(n2/n_t i_k r dr dphi) "
                                     "- (alpha_k + l_k) P_k + g*_k/Gamma_k m h nu_k dnu_k Int(n2/n_t i_k r dr dphi)]"),
        "eq18_19": ("n2bar/n_t = sum_k P_k alpha_k/(h nu_k zeta) / (1 + sum_k P_k (alpha_k+g*_k)/(h nu_k zeta)), "
                    "zeta = pi b_eff^2 n_t / tau"),
        "eq20": ("dP_k/dz = u_k(alpha_k+g*_k)(n2/n_t)P_k + u_k g*_k (n2/n_t) m h nu_k dnu_k "
                 "- u_k(alpha_k + l_k)P_k"),
        "zeta_from_Psat": "zeta = P_sat (alpha_k + g*_k)/(h nu_k)",
        "eq25_saleh": "Q_out = sum_k Q_k_in exp{(alpha_k+g*_k)(Q_in - Q_out)/zeta - alpha_k L}  (photon fluxes)",
        "eq27": "g_k = ln G_k = L[(alpha_k+g*_k)<n2>/n_t - alpha_k]",
        "eq28": "dg_i/dg_j = (alpha_i+g*_i)/(alpha_j+g*_j)",
        "eq30": "eta = 1 - exp{-dg_s (alpha_p+g*_p)/(alpha_s+g*_s)}",
        "eq31": "n_sp = 1/(1 - g*_p alpha_s/(alpha_p g*_s)) at uniform inversion n2/n_t = sigma_ap/(sigma_ap+sigma_ep)",
        "eq32": "P_ASE = 2 n_sp (G - 1) h nu dnu  (per direction)",
        "eq35_37": "P_th = h nu A/(sigma_ak Gamma_k tau)-type threshold/saturation powers; gain responds as a single pole, f_c ~ 5 kHz",
        "transparency_threshold": "P_th = (h nu_p A/((sigma_ap+sigma_ep) tau)) (sigma_as/sigma_es)  (p. 276, verified on scan; printed as-is)",
    },
}
