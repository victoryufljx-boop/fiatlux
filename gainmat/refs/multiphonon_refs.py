"""Multiphonon (nonradiative) relaxation reference data, read from the papers themselves.

Sources (all read page by page; numbers below carry page / table / figure citations):

[LLW77] C. B. Layne, W. H. Lowdermilk, M. J. Weber, "Multiphonon relaxation of rare-earth ions
        in oxide glasses", Phys. Rev. B 16, 10 (1977), pp. 10-20.
[vDS83] J. M. F. van Dijk, M. F. H. Schuurmans, "On the nonradiative and radiative decay rates
        and a modified exponential energy gap law for 4f-4f transitions in rare-earth ions",
        J. Chem. Phys. 78, 5317 (1983), pp. 5317-5323.
[WB04]  B. M. Walsh, N. P. Barnes, "Comparison of Tm:ZBLAN and Tm:silica fiber lasers;
        spectroscopy and tunable pulsed laser operation around 1.9 um", Appl. Phys. B 78,
        325-333 (2004).

Caveats -- read before using these numbers for calibration
---------------------------------------------------------
1. [LLW77] gives NO table of B and alpha.  It plots data and fitted lines only (Fig. 5 for
   silicate, Fig. 6 for all five glasses, p. 17).  The B and alpha in GAP_LAW_LAYNE1977 were
   DIGITIZED here from the Fig. 6 lines: the plot frame and axis ticks were located on a
   300 dpi render, each line was traced pixel by pixel, and ln W was fitted against dE.
   Expect alpha to be good to about +-0.2e-3 cm and B to within a factor of about 2-3.
   B is an extrapolation to dE = 0 and is strongly correlated with alpha, so compare
   models through W at dE = 2000-4000 cm^-1 (the `W_at_3000` field), not through B.
2. The van Dijk Table I values `B_vDS83_tableI` / `alpha_vDS83_tableI` are the numbers
   that paper attributes to [LLW77] (its ref. 17) for tellurite, germanate, phosphate and
   borate.  For silicate and fluoroberyllate it cites Layne & Weber, PRB 16, 3259 (1977)
   (its ref. 18), a DIFFERENT paper on a multicomponent silicate glass.  B = 1.4e12 s^-1,
   alpha = 4.7e-3 cm for "silicate" (the value in gainmat/refs/ln_refs.py GAP_LAW) comes
   from that paper, not from [LLW77].
3. [vDS83] prints alpha = "44.3" (x 1e-3 cm) for borate glass in Table I.  This is a typo
   for 4.3: with 4.3e-3 the table's beta_el column is self-consistent (see note 4), and
   4.25e-3 is what the digitized Fig. 6 of [LLW77] gives.
4. [vDS83] Table I: for every host, tabulated beta_el = 2.0 x beta x exp(-2 alpha hw_max)
   (ratio 1.97-2.03, YAG 2.3).  Eq. (9) alone would give a ratio of 1.  The paper does not
   explain this factor of 2.  We store beta_el exactly as printed.
5. Tm level labels.  [LLW77] uses the old convention: its "Tm 3F4" level, pumped through 3F2
   at 680 nm and monitored at 800 nm, decaying to "3H5" across a gap of about 3800 cm^-1, is
   the level now labelled 3H4 (the ~790 nm pump manifold).  [WB04] p. 327 warns about exactly
   this relabelling.  MULTIPHONON_DATA uses MODERN labels and keeps the original label in
   `label_in_src`.
6. [LLW77] energy gaps are nominal values (p. 15, Sec. III D): the gap from the lowest Stark
   level of the emitting state to the highest Stark level of the terminal state, estimated
   from absorption spectra.  For one ion the gap was taken to be the same in every glass.
   `dE_cm1` holds the gap stated in the text.  `dE_fig_cm1` holds where the point sits in
   Fig. 6, which differs by up to about 120 cm^-1.
7. [WB04] gives no Stark levels and no energy-gap-law fit.  Its implied W_nr values are
   1/tau_meas - 1/tau_rad, measured on a 1.44 wt% Tm2O3 preform (1.12e20 cm^-3).  For 3H4
   this includes Tm-Tm cross-relaxation ("two-for-one"), so it is an UPPER BOUND on the
   multiphonon rate.  The [WB04] gaps given here are barycentre differences computed from
   the Table 3 Judd-Ofelt wavelengths.  They are larger than the Stark-edge gaps that
   [LLW77] uses.
8. Digitized curve values (WB04 Figs. 1-2, LLW77 Figs. 5, 6, 9) are marked `digitized=True`.
"""

# ----------------------------------------------------------------------------------------------
# [LLW77] energy-gap law
# ----------------------------------------------------------------------------------------------
LAYNE1977_FORMULAS = dict(
    gap_law="W = W0 exp(-alpha dE),  alpha = -ln(eps)/(hbar w)   (Eq. 10, p. 13)",
    temperature="W(T) = C [n(T)+1]^p exp(-alpha dE),  n(T) = 1/(exp(hbar w/kT) - 1),  "
                "p = dE/(hbar w)   (Eq. 11, p. 13; n(T) defined p. 13)",
    note="No numerical C, W0 or eps is given anywhere in the paper; all rates are room "
         "temperature measurements except Figs. 8-9.  Multiple excitation of one mode would "
         "give (n+1)(n+2)... instead of (n+1)^p (p. 13).",
    glass_compositions_mol_pct={  # Table I, p. 15
        "borate":    dict(B2O3=67, Na2O=15, BaO=18),
        "phosphate": dict(P2O5=67, Na2O=15, BaO=18),
        "silicate":  dict(SiO2=67, Na2O=15, BaO=18),
        "germanate": dict(GeO2=67, Na2O=15, BaO=18),
        "tellurite": dict(TeO2=82, Al2O3=12, WO3=6),
    },
    RE_doping="Nd, Er, Tm oxides < 0.2 mol% (p. 15)",
    time_resolution_ns=3,
    measurable_W_max_s1=1e8,  # p. 16: "fastest rates that could be measured confidently ~1e8"
    src="[LLW77]",
)

_F6 = "[LLW77] Fig. 6 p. 17 (lines digitized); hw from Fig. 6 caption"
GAP_LAW_LAYNE1977 = {
    #             B (s^-1)  alpha (cm)  hw = Raman high-freq band (Fig. 6 caption)
    "borate": dict(B=1.0e13, alpha=4.21e-3, hw=1350, W_at_3000=3.3e7, fit_range_cm1=(2680, 4960),
                   B_points=1.06e13, alpha_points=4.25e-3,
                   p_Tfit=None, hw_Tfit=None,
                   B_vDS83_tableI=1.4e13, alpha_vDS83_tableI=4.3e-3, hw_vDS83=1400,
                   src=_F6, digitized=True,
                   note="Tm decay varied with T less than the experimental error (p. 19), consistent "
                        "with a low-order process involving 1400 cm^-1 phonons"),
    "phosphate": dict(B=3.5e12, alpha=4.56e-3, hw=1100, W_at_3000=3.9e6, fit_range_cm1=(2320, 4970),
                      B_points=1.4e12, alpha_points=4.30e-3,
                      p_Tfit=4, hw_Tfit=1100,
                      B_vDS83_tableI=6.3e12, alpha_vDS83_tableI=4.7e-3, hw_vDS83=1300,
                      src=_F6 + "; T-fit Fig. 9 caption p. 19", digitized=True),
    "silicate": dict(B=1.9e12, alpha=4.86e-3, hw=1000, W_at_3000=8.6e5, fit_range_cm1=(1980, 4850),
                     B_points=1.2e12, alpha_points=4.69e-3,
                     B_fig5=3.4e12, alpha_fig5=5.00e-3,
                     p_Tfit=4, hw_Tfit=1000,
                     B_vDS83_tableI=1.4e12, alpha_vDS83_tableI=4.7e-3, hw_vDS83=1100,
                     src=_F6 + "; Fig. 5 p. 17 (silicate-only line); T-fit Figs. 8-9 p. 19",
                     digitized=True,
                     note="vDS83 silicate values come from Layne & Weber PRB 16, 3259 (1977), not "
                          "LLW77.  Fig. 8 compares p=4/hw=1000 with p=5/hw=800.  The text (pp. 18-19) "
                          "says 'better fit ... five-phonon decay and ... 1000 cm^-1', while the "
                          "Fig. 9 caption uses 4 x 1000; the paper is internally inconsistent here.  "
                          "Raman peaks at 950, ~1025 and 1100 cm^-1 (p. 19)."),
    "germanate": dict(B=2.6e11, alpha=4.80e-3, hw=900, W_at_3000=1.4e5, fit_range_cm1=(1670, 4510),
                      B_points=7.2e10, alpha_points=4.44e-3,
                      p_Tfit=5, hw_Tfit=825,
                      B_vDS83_tableI=1.2e11, alpha_vDS83_tableI=4.6e-3, hw_vDS83=900,
                      src=_F6 + "; T-fit Fig. 9 caption p. 19", digitized=True),
    "tellurite": dict(B=4.9e10, alpha=4.60e-3, hw=800, W_at_3000=4.9e4, fit_range_cm1=(1330, 4340),
                      B_points=2.8e10, alpha_points=4.42e-3,
                      p_Tfit=5, hw_Tfit=760,
                      B_vDS83_tableI=2.8e10, alpha_vDS83_tableI=4.5e-3, hw_vDS83=800,
                      src=_F6 + "; T-fit Fig. 9 caption p. 19", digitized=True),
}
# B_points / alpha_points: least-squares fit to the digitized data points at the nominal text
# gaps, excluding the 1400 cm^-1 Nd 4G7/2 points, which the paper itself excludes (p. 17).
# Fluoroberyllate is NOT in LLW77 (it is in Layne & Weber PRB 16, 3259); see VANDIJK1983.

# Tm "3F4"(old) = 3H4(modern) -> 3H5 rate vs T: low-T plateau of the Eq. (11) curves in Fig. 9
# (digitized).  Germanate/tellurite curves include the calculated radiative rate (p. 19).
LAYNE1977_FIG9_T0_PLATEAU_s1 = dict(phosphate=1.39e5, silicate=3.0e4, germanate=6.5e3,
                                    tellurite=4.7e3, src="[LLW77] Fig. 9 p. 19", digitized=True)
LAYNE1977_RADIATIVE = {
    # p. 16: Nd 4F3/2 in borate, total decay 1.57e4 s^-1 vs JO radiative 2.45e3 s^-1
    ("Nd3+", "borate", "4F3/2"): dict(W_total=1.57e4, A_rad_JO=2.45e3),
    # p. 17: Tm 3H4 (old "3F4") radiative rates estimated from Reisfeld-Eckstein JO params (+-30%)
    ("Tm3+", "germanate", "3H4"): dict(A_rad_JO=880.0),
    ("Tm3+", "tellurite", "3H4"): dict(A_rad_JO=2.2e3),
    "src": "[LLW77] pp. 16-17",
}

# ----------------------------------------------------------------------------------------------
# Measured multiphonon rates.  Tuples are kept as dicts for clarity.
#   ion, glass, level (MODERN label), dE_cm1, W_s1, T_K, src  (+ extras)
# ----------------------------------------------------------------------------------------------
_L5 = "[LLW77] Fig. 5 p. 17 (digitized)"
_L6 = "[LLW77] Fig. 6 p. 17 (digitized)"
_L9 = "[LLW77] Fig. 9 p. 19 (digitized)"


def _p(ion, glass, level, dE, W, T, src, **kw):
    d = dict(ion=ion, glass=glass, level=level, dE_cm1=dE, W_s1=W, T_K=T, src=src)
    d.update(kw)
    return d


MULTIPHONON_DATA = [
    # ---- LLW77 Fig. 6: five glasses, room temperature.  Gaps from the text (pp. 15-16):
    #      Nd 4G7/2 1400, Nd 2P3/2 2200, Er 4F9/2 2600, Er 4S3/2 2800, Er 4I11/2 3300,
    #      Tm 3H4 (old "3F4") 3800, Nd 4F3/2 4700.
    _p("Er3+", "borate", "4S3/2", 2800, 5.3e7, 295, _L6, dE_fig_cm1=2819, digitized=True),
    _p("Er3+", "borate", "4I11/2", 3300, 7.7e6, 295, _L6, dE_fig_cm1=3302, digitized=True),
    _p("Tm3+", "borate", "3H4", 3800, 2.35e6, 295, _L6, dE_fig_cm1=3791, digitized=True,
       label_in_src="3F4"),
    _p("Nd3+", "borate", "4F3/2", 4700, 1.47e4, 295, _L6, dE_fig_cm1=4767, digitized=True,
       note="W_nr = 1.57e4 - 2.45e3 (JO radiative) = 1.33e4 s^-1 per text p. 16"),
    _p("Er3+", "phosphate", "4F9/2", 2600, 4.7e7, 295, _L6, dE_fig_cm1=2546, digitized=True),
    _p("Er3+", "phosphate", "4S3/2", 2800, 2.76e6, 295, _L6, dE_fig_cm1=2820, digitized=True),
    _p("Er3+", "phosphate", "4I11/2", 3300, 1.03e6, 295, _L6, dE_fig_cm1=3296, digitized=True),
    _p("Tm3+", "phosphate", "3H4", 3800, 1.30e5, 295, _L6, dE_fig_cm1=3801, digitized=True,
       label_in_src="3F4"),
    _p("Nd3+", "silicate", "4G7/2", 1400, 9.0e7, 295, _L6, dE_fig_cm1=1379, digitized=True,
       note="falls BELOW the gap law (1-2 phonon regime, p. 17); near the 1e8 resolution limit"),
    _p("Nd3+", "silicate", "2P3/2", 2200, 2.69e7, 295, _L6, dE_fig_cm1=2134, digitized=True),
    _p("Er3+", "silicate", "4F9/2", 2600, 1.76e7, 295, _L6, dE_fig_cm1=2534, digitized=True),
    _p("Er3+", "silicate", "4S3/2", 2800, 1.61e6, 295, _L6, dE_fig_cm1=2822, digitized=True),
    _p("Er3+", "silicate", "4I11/2", 3300, 1.20e5, 295, _L6, dE_fig_cm1=3284, digitized=True),
    _p("Tm3+", "silicate", "3H4", 3800, 2.92e4, 295, _L6, dE_fig_cm1=3799, digitized=True,
       label_in_src="3F4"),
    _p("Nd3+", "germanate", "4G7/2", 1400, 9.2e7, 295, _L6, dE_fig_cm1=1386, digitized=True,
       note="near 1e8 resolution limit"),
    _p("Nd3+", "germanate", "2P3/2", 2200, 4.87e6, 295, _L6, dE_fig_cm1=2132, digitized=True),
    _p("Er3+", "germanate", "4F9/2", 2600, 1.45e6, 295, _L6, dE_fig_cm1=2554, digitized=True),
    _p("Er3+", "germanate", "4S3/2", 2800, 1.17e5, 295, _L6, dE_fig_cm1=2819, digitized=True),
    _p("Er3+", "germanate", "4I11/2", 3300, 1.52e4, 295, _L6, dE_fig_cm1=3287, digitized=True),
    _p("Tm3+", "germanate", "3H4", 3800, 6.1e3, 295, _L6, dE_fig_cm1=3796, digitized=True,
       label_in_src="3F4", note="radiative 880 s^-1 subtracted by the authors"),
    _p("Nd3+", "tellurite", "4G7/2", 1400, 9.3e7, 295, _L6, dE_fig_cm1=1519, digitized=True,
       note="near 1e8 resolution limit"),
    _p("Er3+", "tellurite", "4F9/2", 2600, 6.65e5, 295, _L6, dE_fig_cm1=2529, digitized=True),
    _p("Er3+", "tellurite", "4S3/2", 2800, 4.6e4, 295, _L6, dE_fig_cm1=2808, digitized=True),
    _p("Er3+", "tellurite", "4I11/2", 3300, 1.01e4, 295, _L6, dE_fig_cm1=3280, digitized=True),
    _p("Tm3+", "tellurite", "3H4", 3800, 1.89e3, 295, _L6, dE_fig_cm1=3789, digitized=True,
       label_in_src="3F4", note="radiative 2.2e3 s^-1 subtracted by the authors"),

    # ---- LLW77 Fig. 9: Tm 3H4 (old "3F4") -> 3H5, dE = 3800, vs temperature (error-bar centres)
    *[_p("Tm3+", "phosphate", "3H4", 3800, W, T, _L9, digitized=True, label_in_src="3F4")
      for T, W in [(300, 1.33e5), (374, 1.55e5), (472, 1.71e5), (523, 1.86e5), (573, 1.83e5),
                   (623, 1.93e5)]],
    *[_p("Tm3+", "silicate", "3H4", 3800, W, T, _L9, digitized=True, label_in_src="3F4")
      for T, W in [(300, 2.95e4), (374, 3.24e4), (423, 3.27e4), (471, 3.54e4), (523, 3.80e4),
                   (573, 3.96e4), (623, 4.48e4), (673, 4.88e4)]],
    *[_p("Tm3+", "germanate", "3H4", 3800, W, T, _L9, digitized=True, label_in_src="3F4",
         note="may include radiative 880 s^-1 (curves include it, p. 19)")
      for T, W in [(298, 7.56e3), (373, 8.3e3), (421, 8.65e3), (471, 8.61e3), (523, 1.01e4),
                   (572, 1.09e4), (623, 1.26e4), (672, 1.48e4), (723, 1.63e4)]],
    *[_p("Tm3+", "tellurite", "3H4", 3800, W, T, _L9, digitized=True, label_in_src="3F4",
         note="may include radiative 2.2e3 s^-1 (curves include it, p. 19)")
      for T, W in [(313, 5.38e3), (373, 5.63e3), (471, 5.71e3), (522, 6.0e3), (572, 6.73e3),
                   (621, 7.56e3)]],

    # ---- WB04: Tm:silica preform, 1.44 wt% Tm2O3 = 1.12e20 cm^-3, room temperature.
    #      W = 1/tau_meas (Table 4) - 1/tau_rad (Table 1).  dE = barycentre difference from
    #      Table 3 JO wavelengths: 3F4->3H6 1847 nm, 3H4->3H5 2428 nm.
    _p("Tm3+", "silica", "3F4", 1e7 / 1847, 1 / 0.42e-3 - 1 / 4.56e-3, 295,
       "[WB04] Tables 1, 3, 4", dE_kind="barycentre (JO lambda)",
       note="W = 2381 - 219 = 2.16e3 s^-1.  Tm conc. 1.12e20 cm^-3; may include concentration "
            "quenching/OH.  The Stark-edge gap is smaller (emission onset ~1.55 um)."),
    _p("Tm3+", "silica", "3H4", 1e7 / 2428, 1 / 21e-6 - 1 / 0.67e-3, 295,
       "[WB04] Tables 1, 3, 4; text p. 331 (21 us main component)", dE_kind="barycentre (JO lambda)",
       upper_bound=True,
       note="W = 4.76e4 - 1.49e3 = 4.6e4 s^-1.  Includes Tm-Tm cross-relaxation; decay is "
            "non-exponential.  Compare LLW77 silicate 3H4: 2.9e4 s^-1 at 0.2 mol% (dE 3800)."),
]

# Silicate-only data of LLW77 Fig. 5, for completeness (same glass and points as Fig. 6).
LAYNE1977_FIG5_SILICATE = dict(
    points=[("Nd3+", "4G7/2", 1402, 8.4e7), ("Nd3+", "2P3/2", 2148, 2.4e7),
            ("Er3+", "4F9/2", 2563, 1.6e7), ("Er3+", "4S3/2", 2850, 1.28e6),
            ("Er3+", "4I11/2", 3290, 1.28e5), ("Tm3+", "3H4 (old 3F4)", 3808, 3.1e4)],
    line=dict(B=3.4e12, alpha=5.0e-3),
    src=_L5, digitized=True,
)

# ----------------------------------------------------------------------------------------------
# [vDS83] modified energy-gap law
# ----------------------------------------------------------------------------------------------
VANDIJK1983 = dict(
    form="k_NR = beta_el exp[-(dE - 2 hw_max) alpha]   (Eq. 9, p. 5319)",
    form_classic="k_NR = beta exp(-alpha dE)   (Eq. 1, p. 5317)",
    theory="F(v) = exp(-D Delta^2/2) (D Delta^2/2)^v / v!,  v = dE/hw_max - 1 (Eqs. 5-6); "
           "the exponential fit at v~3-4 extrapolates to F~1 at v=1, hence Eq. 8-9 (p. 5319)",
    factor_2_range=(1.5, 2.5),       # p. 5319: 'factor of 2 in Eq. (9) ... between 1.5 and 2.5'
    beta_el_spread="beta_el varies only by a factor of 10 over all hosts (p. 5319)",
    universal_alpha=None,            # the paper gives NO universal alpha; alpha is per host
    beta_el_range_s1=(0.3e7, 20.4e7),        # Table I, all hosts
    beta_el_range_glasses_s1=(4.2e7, 20.4e7),
    beta_el_geomean_all_s1=3.75e7,   # computed here from Table I (14 rows)
    beta_el_geomean_glasses_s1=8.9e7,  # computed here (6 glasses)
    beta_el_calc_formula="beta_el ~ (pi c^3 / 4) M dE0 C_sym k_R  (Eq. 28, atomic units)",
    C_sym_range=(1, 10),             # p. 5322 'a posteriori C_sym ~ 1-10'
    T_dependence="(t+1)^N with t = 1/(exp(hw/kT)-1) (p. 5322)",
    # Table I, p. 5318.  beta (s^-1), beta_el (s^-1), beta_el_calc range (s^-1), ratio
    # beta_el/beta_el_calc, dE0 range (cm^-1), kR range (s^-1), hw_max (cm^-1), alpha (cm).
    table_I={
        "tellurite":      dict(beta=2.8e10, beta_el=4.2e7, beta_el_calc=(0.5e7, 2e7), ratio=(10, 2),
                               dE0=(2000, 4000), kR=(4, 30), hw_max=800, alpha=4.5e-3, ref_beta=17, ref_kR=23),
        "germanate":      dict(beta=1.2e11, beta_el=6.1e7, beta_el_calc=(0.5e7, 2e7), ratio=(10, 3),
                               dE0=(2000, 4000), kR=(4, 30), hw_max=900, alpha=4.6e-3, ref_beta=17, ref_kR=23),
        "phosphate":      dict(beta=6.3e12, beta_el=6.2e7, beta_el_calc=(0.5e7, 2e7), ratio=(10, 3),
                               dE0=(2000, 4000), kR=(4, 30), hw_max=1300, alpha=4.7e-3, ref_beta=17, ref_kR=23),
        "borate":         dict(beta=1.4e13, beta_el=16.8e7, beta_el_calc=(0.5e7, 2e7), ratio=(30, 10),
                               dE0=(2000, 4000), kR=(4, 30), hw_max=1400, alpha=4.3e-3, ref_beta=17, ref_kR=23,
                               note="printed as 44.3e-3 (typo); 4.3e-3 is self-consistent"),
        "silicate":       dict(beta=1.4e12, beta_el=9.0e7, beta_el_calc=(0.5e7, 2e7), ratio=(20, 5),
                               dE0=(2000, 4000), kR=(4, 30), hw_max=1100, alpha=4.7e-3, ref_beta=18, ref_kR=18),
        "fluoroberyllate": dict(beta=1.2e12, beta_el=20.4e7, beta_el_calc=(0.4e7, 2e7), ratio=(50, 10),
                                dE0=(2000, 4000), kR=(4, 30), hw_max=900, alpha=5.2e-3, ref_beta=18, ref_kR=18),
        "SrF2":           dict(beta=3.3e8, beta_el=2.6e7, beta_el_calc=(0.4e7, 1e7), ratio=(6, 3),
                               dE0=(1500, 2500), kR=(3, 10), hw_max=360, alpha=4.5e-3, ref_beta=20, ref_kR=20),
        "LaF3":           dict(beta=3.9e8, beta_el=1.9e7, beta_el_calc=(0.3e7, 1e7), ratio=(6, 2),
                               dE0=(1500, 2500), kR=(2, 10), hw_max=350, alpha=5.3e-3, ref_beta=20, ref_kR=20),
        "LaCl3":          dict(beta=5.4e10, beta_el=6.4e7, beta_el_calc=(0.1e7, 0.5e7), ratio=(50, 10),
                               dE0=(800, 1500), kR=(1, 7), hw_max=260, alpha=14.3e-3, ref_beta=19, ref_kR=19),
        "LaBr3":          dict(beta=1.4e10, beta_el=2.9e7, beta_el_calc=(0.1e7, 0.1e7), ratio=(30, 30),
                               dE0=(600, 1000), kR=(0.5, 3), hw_max=175, alpha=19.8e-3, ref_beta=20, ref_kR=20),
        "Y2O3(a)":        dict(beta=1.6e8, beta_el=0.6e7, beta_el_calc=(0.3e7, 1e7), ratio=(2, 1),
                               dE0=(1800, 2800), kR=(3, 13), hw_max=550, alpha=3.6e-3, ref_beta=20, ref_kR=20),
        "Y2O3(b)":        dict(beta=4.6e9, beta_el=3.4e7, beta_el_calc=(0.3e7, 1e7), ratio=(10, 3),
                               dE0=(1800, 2800), kR=(3, 13), hw_max=550, alpha=5.1e-3, ref_beta=26, ref_kR=20),
        "Y3Al5O12":       dict(beta=1.0e8, beta_el=0.3e7, beta_el_calc=(0.2e7, 1e7), ratio=(1, 3),
                               dE0=(1200, 2900), kR=(1, 11), hw_max=700, alpha=3.1e-3, ref_beta=21, ref_kR=(24, 25)),
        "YAlO3":          dict(beta=5e9, beta_el=4.0e7, beta_el_calc=(0.3e7, 2e7), ratio=(10, 2),
                               dE0=(1500, 3200), kR=(2, 24), hw_max=600, alpha=4.6e-3, ref_beta=22, ref_kR=22),
    },
    table_refs={17: "Layne, Lowdermilk, Weber, PRB 16, 10 (1977)",
                18: "Layne & Weber, PRB 16, 3259 (1977)",
                19: "Riseberg, Gandrud, Moos, Phys. Rev. 159, 262 (1967)",
                20: "Riseberg & Moos, Phys. Rev. 174, 429 (1968)",
                21: "Zverev, Kolodnyi, Onishchenko, Sov. Phys. JETP 33, 497 (1971)",
                22: "Weber, PRB 8, 54 (1973)",
                23: "Reisfeld & Eckstein, J. Chem. Phys. 63, 4001 (1975)",
                26: "Weber, Phys. Rev. 171, 283 (1968)"},
    note="Tabulated beta_el = 2.0 * beta * exp(-2 alpha hw_max) for every row (see module note 4).  "
         "Eu3+ 5D1->5D0 in Y2O3: k_NR ~ 5e3 s^-1, ~100x below the gap law (J=1->0 forbidden), p. 5321.",
    src="[vDS83] Table I p. 5318; Eqs. 1, 5-9 pp. 5317-5319; Eq. 28 p. 5322",
)

# ----------------------------------------------------------------------------------------------
# [WB04] Tm3+ in silica (preform / fiber core), and ZBLAN for comparison
# ----------------------------------------------------------------------------------------------
TM_SILICA_WALSH2004 = dict(
    host="silica preform/fiber core (INO); inner cladding SiO2 co-doped F, P2O5 (p. 329)",
    Tm2O3_wt_pct=1.44, N_Tm_cm3=1.12e20,          # Table 4; p. 330
    n=1.44,                                       # p. 329, at 1.8 um
    hw_max_cm1=1100,                              # abstract (cites Reisfeld & Jorgensen 1987)
    Omega_1e20_cm2=dict(O2=6.23, O4=1.91, O6=1.36),   # Table 1
    tau_rad_s={"3F4": 4.56e-3, "3H4": 0.67e-3,        # Table 1
               # Table 3 last column (us), JO radiative lifetimes of every manifold:
               "3H5": 3896e-6, "3F3": 447e-6, "3F2": 631e-6, "1G4": 428e-6, "1D2": 29e-6,
               "1I6": 59e-6},
    tau_rad_3F4_table3_s=4559e-6,
    tau_meas_s={"3F4": 0.42e-3, "3H4": 21e-6},   # Table 4 (3H4 printed 0.02 ms; text p. 331 '21 us')
    tau_meas_conditions=dict(
        T_K=295,  # room temperature (not stated explicitly)
        **{"3F4": "Co:MgF2 excitation at 1.78 um, detection at 1.8 um, InSb (p. 327)",
           "3H4": "alexandrite excitation at 0.78 um, detection at 0.82 um (815 nm filter), "
                  "non-exponential, main component quoted (pp. 327, 331)"}),
    tau_meas_3H5_s=None,                          # not measured in WB04
    W_nr_implied_s1={"3F4": 1 / 0.42e-3 - 1 / 4.56e-3,   # = 2.16e3
                     "3H4": 1 / 21e-6 - 1 / 0.67e-3},    # = 4.6e4 (upper bound; incl. cross-relaxation)
    eta_Q={"3F4": 0.42 / 4.56, "3H4": 0.021 / 0.67},     # = 0.092, 0.031
    # Table 3: A_ED (s^-1) and branching ratio beta from the JO analysis (lambda in nm)
    JO_transitions={
        ("3F4", "3H6"): dict(lam=1847, S_ED=5.060, A=219.36, beta=1.000),
        ("3H5", "3F4"): dict(lam=3563, S_ED=2.072, A=10.22, beta=0.040),
        ("3H5", "3H6"): dict(lam=1216, S_ED=1.977, A=246.42, beta=0.960),
        ("3H4", "3H5"): dict(lam=2428, S_ED=1.009, A=19.22, beta=0.013),
        ("3H4", "3F4"): dict(lam=1444, S_ED=1.333, A=121.17, beta=0.084),
        ("3H4", "3H6"): dict(lam=810, S_ED=2.494, A=1295.24, beta=0.902),
        ("3F3", "3H4"): dict(lam=5200, S_ED=1.574, A=3.92, beta=0.002),
        ("3F3", "3H5"): dict(lam=1655, S_ED=4.579, A=354.79, beta=0.158),
        ("3F3", "3F4"): dict(lam=1130, S_ED=0.243, A=59.39, beta=0.027),
        ("3F3", "3H6"): dict(lam=701, S_ED=1.757, A=1820.65, beta=0.813),
        ("3F2", "3F3"): dict(lam=20449, S_ED=0.166, A=0.01, beta=0.000),
        ("3F2", "3H4"): dict(lam=4145, S_ED=2.486, A=17.12, beta=0.011),
        ("3F2", "3H5"): dict(lam=1531, S_ED=1.354, A=185.58, beta=0.117),
        ("3F2", "3F4"): dict(lam=10722, S_ED=2.036, A=819.48, beta=0.517,
                             note="printed 10722 nm; the 3F2-3F4 gap (~9400 cm^-1) implies ~1072 nm (decimal slip in source)"),
        ("3F2", "3H6"): dict(lam=678, S_ED=0.350, A=562.49, beta=0.355),
        ("1G4", "3F2"): dict(lam=1632, S_ED=0.233, A=14.64, beta=0.006),
        ("1G4", "3F3"): dict(lam=1511, S_ED=0.606, A=47.93, beta=0.020),
        ("1G4", "3H4"): dict(lam=1171, S_ED=1.483, A=253.07, beta=0.108),
        ("1G4", "3H5"): dict(lam=790, S_ED=1.190, A=667.13, beta=0.285),
        ("1G4", "3F4"): dict(lam=647, S_ED=0.157, A=162.01, beta=0.069),
        ("1G4", "3H6"): dict(lam=479, S_ED=0.461, A=1193.17, beta=0.510),
    },
    JO_transitions_units="lam nm; S_ED 1e-20 cm^2; A s^-1 (ED only)",
    JO_note="Same U^(k) and wavelengths are printed for ZBLAN (Table 2); 3F2->3F4 '10722 nm' is "
            "almost certainly 1072.2 nm.  1D2 and 1I6 rows omitted here (see the paper).",
    sigma_em_3F4_3H6=dict(peak_cm2=0.46e-20, peak_um=1.80, FWHM_um=0.31,
                          method="reciprocity/F-L with tau_rad = 4.56 ms (p. 326)",
                          others="Zou & Toratani 6.1e-21 cm^2 (tau_rad 6.3 ms); Jackson & King ~25% higher"),
    sigma_abs_3H6_3H4=dict(peak_cm2=1.00e-20, peak_um=0.789, FWHM_um=0.0154,
                           note="ZBLAN ~3x lower (0.33e-20 at 0.791 um)"),
    sigma_abs_3H6_3F4=None,  # not shown in WB04
    # Digitized curves (Fig. 2 emission, Fig. 1 absorption), um -> 1e-20 cm^2.
    sigma_em_3F4_curve=[(1.550, 0.013), (1.575, 0.020), (1.600, 0.048), (1.625, 0.092), (1.650, 0.151),
                        (1.675, 0.222), (1.700, 0.298), (1.725, 0.367), (1.750, 0.420), (1.775, 0.453),
                        (1.800, 0.464), (1.825, 0.458), (1.850, 0.437), (1.875, 0.414), (1.900, 0.385),
                        (1.925, 0.352), (1.950, 0.316), (1.975, 0.266), (2.000, 0.212), (2.025, 0.160),
                        (2.050, 0.115), (2.075, 0.083), (2.100, 0.058), (2.125, 0.042), (2.150, 0.027),
                        (2.175, 0.018)],
    sigma_abs_3H4_curve=[(0.750, 0.037), (0.755, 0.068), (0.760, 0.110), (0.765, 0.155), (0.770, 0.191),
                         (0.775, 0.238), (0.780, 0.377), (0.785, 0.776), (0.790, 0.996), (0.795, 0.662),
                         (0.800, 0.407), (0.805, 0.284), (0.810, 0.168), (0.815, 0.091), (0.820, 0.043),
                         (0.825, 0.027)],
    curves_digitized=True,
    stark_levels=None,  # WB04 gives no Stark energies
    gap_law_fit=None,   # WB04 fits no energy-gap-law parameters; it only quotes (1+m)^p, p=dE/hw (p. 328)
    ZBLAN=dict(Omega_1e20_cm2=dict(O2=1.96, O4=1.36, O6=1.16), tau_rad_s={"3F4": 11.1e-3, "3H4": 1.52e-3},
               tau_meas_s={"0.32wt%": {"3F4": 11.2e-3, "3H4": 1.51e-3},
                           "2.88wt%": {"3F4": 13.7e-3, "3H4": 0.33e-3}},
               N_Tm_cm3={"0.32wt%": 0.36e20, "2.88wt%": 3.23e20}, hw_max_cm1=500,
               sigma_em_3F4_peak=(0.24e-20, 1.82)),
    src="[WB04] Tables 1-4, Figs. 1-2, pp. 325-331",
)
