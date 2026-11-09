"""Spectroscopic / energy-transfer data for Er3+ and Er3+,Pr3+ ZBLAN.

Source (verified from the PDF masthead, p. 856):
    P. S. Golding, S. D. Jackson, T. A. King, and M. Pollnau,
    "Energy transfer processes in Er3+-doped and Er3+,Pr3+-codoped ZBLAN glasses",
    Phys. Rev. B 62(2), 856-864 (1 July 2000-II).  DOI 10.1103/PhysRevB.62.856
    Received 26 Oct 1999; revised 10 Feb 2000.  Laser Photonics Group,
    Univ. of Manchester.  Bulk glass samples from Le Verre Fluore (France).

Page numbers below are journal page numbers (856-864); PDF page = journal page - 855.

What the paper provides (and what it does NOT):
  * All transfer coefficients are MACROSCOPIC rate-equation parameters W (cm^3/s)
    (Table I, p. 858), NOT Forster/Dexter microparameters.  No C_DA, R0,
    critical concentrations, Inokuti-Hirayama or Burshtein/Yokota-Tanimoto fits.
  * Conventions (important when importing into a rate-equation model):
      - ETU:  dN/dt = -N/tau - 2*W_ETU*N^2   (Eq. 3, p. 859), i.e. the event rate
        is W_ETU*N^2 and each event removes TWO excitations from the level
        (one goes to ground, one is promoted).  Jensen's (Ref. 37, LiYF4) alpha
        must be divided by 2 to compare with W (p. 861).
      - ETU from 4I11/2 was evaluated with back-transfer correction via 4S3/2 CR:
        dN/dt = -N/tau - (2-beta)*W_ETU*N^2, beta = 1 - tau_eff(4S3/2)/tau(4S3/2)
        (Eqs. 6-7, p. 860).  W is still the event-rate coefficient.
      - CR:  R_CR = W_CR * N(4S3/2) * N(4I15/2), with N(4I15/2) ~ N(Er)
        (Eqs. 1-2, p. 858); W_CR obtained from 1/tau_eff = 1/tau + W_CR*N(Er).
      - Er->Pr ET:  R_ET = W_ET * N(Er*) * N(Pr 3H4), N(3H4) ~ N(Pr)
        (Eqs. 8-9, p. 860); W_ET from 1/tau_eff(Er*) = 1/tau(Er*) + W_ET*N(Pr).
  * ETU parameters from the initial (0-~1.6 ms) part of the decay linearised per
    Eq. 5 (Fig. 5, p. 859); intrinsic lifetimes from the late exponential tail.
  * W values depend on Er concentration (roughly linear in N_Er, i.e. a
    migration-accelerated regime, p. 862), and ET depends only weakly on Pr
    concentration.  Do NOT treat any W as a concentration-independent constant.
  * No Judd-Ofelt parameters, no absorption/emission cross sections, no
    branching ratios, no multiphonon rates, no glass composition and no density
    are given.  Table II gives oscillator strengths (taken from Ref. 39,
    Wetenkamp PhD thesis 1991) used only for a qualitative ranking.
  * Pr3+ level lifetimes were NOT measured.  Er 4F9/2 and 4I9/2 lifetimes were
    NOT measured (4I9/2 only "calculated from Ref. 36 to be < 20 us", p. 859).
  * Values marked 'digitized' were read from figures by pixel calibration of
    400-dpi renders (estimated uncertainty ~+-5-10 % for log-axis points in
    Figs. 3/6a, worse (~+-0.5 ms) for the overlapping 4I13/2/4I11/2 points of
    Fig. 3 above the lowest concentration).  Fig. 6(a) digitized lifetimes were
    cross-checked against Eq. 9 with Table I W_ET values and agree to <5 %.
  * Concentration conversion used by the paper: 1 mol% = 1.6e20 cm^-3 (p. 857).
"""

CITATION = ("P. S. Golding, S. D. Jackson, T. A. King, M. Pollnau, "
            "Phys. Rev. B 62, 856 (2000), doi:10.1103/PhysRevB.62.856")

MOL_PCT_TO_CM3 = 1.6e20  # p. 857, Sec. II: "(1 mol% = 1.6x10^20 cm^-3)"

# ---------------------------------------------------------------------------
# Energy-transfer processes (Table I, p. 858; process definitions pp. 858-861;
# Fig. 1, p. 856).  k_cm3_s = macroscopic W at the listed concentrations.
# ---------------------------------------------------------------------------
ENERGY_TRANSFER = [
    # --- ETU from Er 4I13/2 -------------------------------------------------
    {
        "name": "ETU1 (Er-Er upconversion from 4I13/2)",
        "donor_species": "Er3+", "donor_from": "4I13/2", "donor_to": "4I15/2",
        "acceptor_species": "Er3+", "acceptor_from": "4I13/2", "acceptor_to": "4I9/2",
        "rate_law": "R = W*N(4I13/2)^2; dN13/dt = -N13/tau - 2*W*N13^2 (Eq. 3)",
        "k_cm3_s": {2e20: 1.3e-17, 8e20: 2.8e-17, 14e20: 6.7e-17},
        "concentration": "keys = N_Er (cm^-3) = 1.25, 5, 8.75 mol%; Pr = 0",
        "uncertainty": "error bar shown only at 2e20 in Fig. 4, ~1.0-1.7e-17 (digitized)",
        "method": "initial-decay linearisation, Eq. 5 (Fig. 5b), tau = 9.0 ms of lowest-doped "
                  "sample; back-CR (4I9/2,4I15/2)->(4I13/2,4I13/2) neglected (<10 % effect per Ref. 37)",
        "concentration_dependence": "increases with N_Er; log-log slope ~0.84 over 2-14e20 (derived "
                                    "from Table I; paper: 'close to linear', p. 862)",
        "src": "Table I p. 858; Eqs. 3-5 p. 859; Fig. 4 p. 858",
    },
    # --- ETU from Er 4I11/2 -------------------------------------------------
    {
        "name": "ETU2 (Er-Er upconversion from 4I11/2)",
        "donor_species": "Er3+", "donor_from": "4I11/2", "donor_to": "4I15/2",
        "acceptor_species": "Er3+", "acceptor_from": "4I11/2", "acceptor_to": "4F7/2",
        "rate_law": "R = W*N(4I11/2)^2; dN11/dt = -N11/tau - (2-beta)*W*N11^2 (Eq. 7)",
        "k_cm3_s": {2e20: 0.2e-17, 8e20: 1.0e-17, 14e20: 1.9e-17},
        "concentration": "keys = N_Er (cm^-3); Pr = 0",
        "uncertainty": "error bar at 2e20 in Fig. 4 roughly 1.8-2.4e-18 (digitized)",
        "method": "Eq. 5 linearisation with slope [1 + tau_eff(4S3/2)/tau(4S3/2)]*W*N0*tau, "
                  "tau(4S3/2) = 586 us, tau_eff from Fig. 3 (accounts for back transfer "
                  "4F7/2->4S3/2 -CR-> 4I9/2 -> 4I11/2)",
        "concentration_dependence": "log-log slope ~1.16 over 2-14e20 (derived from Table I)",
        "src": "Table I p. 858; Eqs. 6-7 p. 860; Fig. 4 p. 858",
    },
    # --- CR from thermally coupled 4S3/2/2H11/2 (two channels, one lumped W) --
    {
        "name": "CR (4S3/2/2H11/2 self-quenching), channel a",
        "donor_species": "Er3+", "donor_from": "2H11/2 (thermalised with 4S3/2)", "donor_to": "4I9/2",
        "acceptor_species": "Er3+", "acceptor_from": "4I15/2", "acceptor_to": "4I13/2",
        "rate_law": "R_CR = W_CR*N(4S3/2+2H11/2)*N(4I15/2) (Eq. 1), lumped over channels a+b",
        "k_cm3_s": {0.4e20: 0.6e-17, 2e20: 0.6e-17, 8e20: 2.4e-17, 14e20: 3.3e-17},
        "concentration": "keys = N_Er (cm^-3) = 0.25, 1.25, 5, 8.75 mol%; Pr = 0",
        "method": "1/tau_eff(4S3/2) = 1/586us + W_CR*N_Er (Eq. 2) using measured 4S3/2 lifetimes "
                  "(532-nm excitation, 550-nm detection)",
        "note": "W_CR is ONE lumped coefficient for channels a and b (and endothermic CR from 4S3/2, "
                "p. 861). Paper writes channel a as (2H11/2,4I15/2)->(4I9/2,4I13/2).",
        "src": "Table I p. 858; Eqs. 1-2 p. 858; text p. 858 col. 1",
    },
    {
        "name": "CR (4S3/2/2H11/2 self-quenching), channel b",
        "donor_species": "Er3+", "donor_from": "2H11/2 (thermalised with 4S3/2)", "donor_to": "4I13/2",
        "acceptor_species": "Er3+", "acceptor_from": "4I15/2", "acceptor_to": "4I9/2",
        "rate_law": "same lumped W_CR as channel a (do not double count)",
        "k_cm3_s": None,
        "concentration": "see channel a",
        "note": "Paper writes (2H11/2,4I15/2)->(4 13/2,4I9/2) [sic, typo for 4I13/2]. Not separately quantified.",
        "src": "p. 858 col. 1; p. 861 (oscillator-strength discussion)",
    },
    # --- Er -> Pr ET ---------------------------------------------------------
    {
        "name": "ET Er 4I13/2 -> Pr 3F3 (lower-laser-level quenching)",
        "donor_species": "Er3+", "donor_from": "4I13/2", "donor_to": "4I15/2",
        "acceptor_species": "Pr3+", "acceptor_from": "3H4", "acceptor_to": "3F3",
        "rate_law": "R_ET = W_ET*N(Er 4I13/2)*N(Pr 3H4) (Eq. 8); 1/tau_eff = 1/tau + W_ET*N_Pr (Eq. 9)",
        "k_cm3_s": {  # (N_Er, N_Pr) in cm^-3
            (8e20, 0.4e20): 9.0e-17, (8e20, 0.8e20): 10.6e-17, (8e20, 2.5e20): 13.2e-17,
            (14e20, 0.4e20): 16.8e-17, (14e20, 0.8e20): 17.2e-17, (14e20, 2.0e20): 20.5e-17,
        },
        "concentration": "keys = (N_Er, N_Pr) cm^-3; Er 5 mol% with Pr 0.25/0.5/1.55 mol%, "
                         "Er 8.75 mol% with Pr 0.25/0.5/1.25 mol%",
        "method": "effective-lifetime quenching (Eq. 9) of the 1560-nm decay after 1510-nm excitation",
        "concentration_dependence": "~linear in N_Er (donor), nearly independent of N_Pr (acceptor) -> "
                                    "fast migration within Er 4I13/2 (abstract; p. 862)",
        "note": "Back transfer inhibited by fast multiphonon relaxation from Pr 3F3 (abstract). "
                "Acceptor level named '3F3' by the paper (p. 861, p. 862 conclusions); the Pr 3F3/3F4 "
                "pair is near-degenerate - a model may lump into 3F3+3F4 then relax to 3F2/3H6.",
        "src": "Table I p. 858; Eqs. 8-9 p. 860; Fig. 6b p. 860; Fig. 4 open squares (N_Pr = 8e19)",
    },
    {
        "name": "ET Er 4I11/2 -> Pr 1G4",
        "donor_species": "Er3+", "donor_from": "4I11/2", "donor_to": "4I15/2",
        "acceptor_species": "Pr3+", "acceptor_from": "3H4", "acceptor_to": "1G4",
        "rate_law": "R_ET = W_ET*N(Er 4I11/2)*N(Pr 3H4) (Eq. 8)",
        "k_cm3_s": {
            (8e20, 0.4e20): 0.9e-17, (8e20, 0.8e20): 1.0e-17, (8e20, 2.5e20): 1.0e-17,
            (14e20, 0.4e20): 2.1e-17, (14e20, 0.8e20): 2.0e-17, (14e20, 2.0e20): 2.2e-17,
        },
        "concentration": "keys = (N_Er, N_Pr) cm^-3",
        "method": "effective-lifetime quenching (Eq. 9) of the 1010-nm decay after 979-nm excitation",
        "concentration_dependence": "increases with N_Er, ~independent of N_Pr",
        "src": "Table I p. 858; Eqs. 8-9 p. 860; Fig. 6b p. 860",
    },
    # --- mentioned, not quantified --------------------------------------------
    {
        "name": "Back-CR from 4I9/2 (inverse of ETU1) - mentioned only",
        "donor_species": "Er3+", "donor_from": "4I9/2", "donor_to": "4I13/2",
        "acceptor_species": "Er3+", "acceptor_from": "4I15/2", "acceptor_to": "4I13/2",
        "k_cm3_s": None,
        "concentration": None,
        "note": "Neglected: 4I9/2 lifetime < 20 us (calc. from Ref. 36); in LiYF4 changes ETU1 W by <10 %.",
        "src": "p. 859 col. 2",
    },
]

# ---------------------------------------------------------------------------
GOLDING2000 = {
    "citation": CITATION,
    "doi": "10.1103/PhysRevB.62.856",

    # ---------------- samples ----------------
    "samples": {
        "src": "Sec. II p. 857; Table I p. 858",
        "host": "ZBLAN bulk glass (Le Verre Fluore); composition NOT given",
        "density_g_cm3": None,  # not given
        "dimensions_mm": (10, 10, 5),
        "mol_pct_to_cm3": MOL_PCT_TO_CM3,
        "Er_only": [  # (Er mol%, N_Er cm^-3)
            {"Er_mol_pct": 0.25, "N_Er_cm3": 0.4e20},
            {"Er_mol_pct": 1.25, "N_Er_cm3": 2e20},
            {"Er_mol_pct": 5.0, "N_Er_cm3": 8e20},
            {"Er_mol_pct": 8.75, "N_Er_cm3": 14e20},
        ],
        "Er_Pr": [
            {"Er_mol_pct": 8.75, "N_Er_cm3": 14e20, "Pr_mol_pct": 0.25, "N_Pr_cm3": 0.4e20},
            {"Er_mol_pct": 8.75, "N_Er_cm3": 14e20, "Pr_mol_pct": 0.5, "N_Pr_cm3": 0.8e20},
            {"Er_mol_pct": 8.75, "N_Er_cm3": 14e20, "Pr_mol_pct": 1.25, "N_Pr_cm3": 2.0e20},
            {"Er_mol_pct": 5.0, "N_Er_cm3": 8e20, "Pr_mol_pct": 0.25, "N_Pr_cm3": 0.4e20},
            {"Er_mol_pct": 5.0, "N_Er_cm3": 8e20, "Pr_mol_pct": 0.5, "N_Pr_cm3": 0.8e20},
            {"Er_mol_pct": 5.0, "N_Er_cm3": 8e20, "Pr_mol_pct": 1.55, "N_Pr_cm3": 2.5e20},
        ],
    },

    "experiment": {
        "src": "Sec. II p. 857; Fig. 5 p. 859",
        "excitation_nm": {"4I13/2": 1510, "4I11/2": 979, "4S3/2/2H11/2": 532},
        "detection_nm": {"4I13/2": 1560, "4I11/2": 1010, "4S3/2": 550},
        "pump_pulse_ns": 6, "pump_bandwidth_nm": 5, "detection_bandwidth_nm": 4,
        "initial_excitation_N0_cm3_range_Fig5": (8.6e18, 2.2e19),  # 4I13/2, 8.75 mol% sample
        "Fig5_windows_ms": {"ETU_fit": (0, 1.6), "exponential_tail_fit": (4.8, 8.0)},  # digitized
    },

    # ---------------- lifetimes, Er-only samples ----------------
    "lifetimes_Er_only": {
        "tau_intrinsic_4I13_2_ms": 9.0, "tau_intrinsic_4I13_2_err_ms": 0.2,
        "tau_intrinsic_4I11_2_ms": 6.9, "tau_intrinsic_4I11_2_err_ms": 0.1,
        "intrinsic_src": "p. 858 col. 2 (lowest-doped 0.25 mol% sample, exponential tail); "
                         "abstract: 4I13/2 intrinsic 9 ms",
        "tau_intrinsic_4S3_2_us": 586,
        "tau_intrinsic_4S3_2_src": "p. 858 col. 1: extrapolation to N_Er -> 0 (Ref. 34 gives 570 us)",
        "tau_4I9_2_us_upper_bound": 20,
        "tau_4I9_2_src": "p. 859 col. 2: 'calculated from Ref. 36 to be < 20 us' (not measured)",
        "tau_4F9_2": None,  # not measured / not given
        # Fig. 3 (p. 858): lifetimes vs N_Er
        "tau_eff_4S3_2_us_vs_NEr": {
            0.4e20: 518,   # +-6, stated p. 858
            2e20: 345,     # digitized Fig. 3 (Eq. 2 with Table I W_CR gives 344)
            8e20: 47,      # digitized Fig. 3 (Eq. 2 gives 47.8)
            14e20: 21,     # +-1, stated p. 858
        },
        "tau_4I13_2_ms_vs_NEr": {0.4e20: 9.0, 2e20: 8.6, 8e20: 8.5, 14e20: 8.3},
        "tau_4I11_2_ms_vs_NEr": {0.4e20: 6.9, 2e20: 7.6, 8e20: 7.6, 14e20: 7.7},
        "tau_vs_NEr_note": "Values at N_Er >= 2e20 digitized from Fig. 3 (overlapping markers, "
                           "~+-0.5 ms). Text p. 858: 4I11/2 lifetime increases and 4I13/2 decreases "
                           "with N_Er. These are exponential-TAIL lifetimes, not 1/e times.",
        "decay_shape": "4I13/2 and 4I11/2: fast non-exponential initial part (ETU, pump-energy "
                       "dependent) followed by exponential tail independent of pump (p. 858, Fig. 5a). "
                       "4S3/2: single effective lifetime reported (shape not discussed).",
        "src": "Fig. 3 p. 858; text p. 858",
    },

    # ---------------- lifetimes, Er,Pr co-doped (Fig. 6a, digitized) ----------
    "lifetimes_Er_Pr_codoped": {
        "src": "Fig. 6(a) p. 860 (digitized); cross-checked with Eq. 9 + Table I",
        "note": "Effective lifetimes; how much non-exponentiality remained is not stated.",
        "tau_eff_4I13_2_us": {  # (N_Er, N_Pr) cm^-3
            (8e20, 0.4e20): 267, (8e20, 0.8e20): 114, (8e20, 2.5e20): 31,
            (14e20, 0.4e20): 147, (14e20, 0.8e20): 70, (14e20, 2.0e20): 25,
        },
        "tau_eff_4I11_2_us": {
            (8e20, 0.4e20): 2030, (8e20, 0.8e20): 1110, (8e20, 2.5e20): 394,
            (14e20, 0.4e20): 1060, (14e20, 0.8e20): 566, (14e20, 2.0e20): 230,
        },
        "eq9_check_4I13_2_us": {  # 1/(1/9ms + W_ET*N_Pr)
            (8e20, 0.4e20): 269, (8e20, 0.8e20): 116, (8e20, 2.5e20): 30,
            (14e20, 0.4e20): 146, (14e20, 0.8e20): 72, (14e20, 2.0e20): 24,
        },
        "eq9_check_4I11_2_us": {  # 1/(1/6.9ms + W_ET*N_Pr)
            (8e20, 0.4e20): 1980, (8e20, 0.8e20): 1060, (8e20, 2.5e20): 378,
            (14e20, 0.4e20): 1015, (14e20, 0.8e20): 573, (14e20, 2.0e20): 220,
        },
    },

    # ---------------- Pr effect on Er 4I13/2 (laser lower level) --------------
    "Pr_effect_on_4I13_2": {
        "statement_abstract": "ET effectively quenches the 4I13/2 intrinsic lifetime of 9 ms down to "
                              "20 us for the highest Er3+ and Pr3+ concentrations investigated",
        "tau_4I13_2_quenched_us_abstract": 20,  # Eq. 9 gives ~24 us for Er 14e20 / Pr 2e20
        "W_ET_4I13_2_over_W_ET_4I11_2": "~9-13 (Table I)",
        "statements": [
            "ET from 4I13/2 more efficient than from 4I11/2: large Pr 3H4->3F3 oscillator strength "
            "and back transfer inhibited by fast multiphonon relaxation from Pr level (abstract)",
            "Small amounts of Pr3+ rapidly deplete the lower laser level while having much lesser "
            "effect on the upper level (p. 863)",
            "Pr codoping preferred when high Er concentrations for efficient ETU recycling are "
            "unavailable (p. 863)",
        ],
        "src": "abstract p. 856; Fig. 6 p. 860; Sec. V pp. 862-863",
    },

    # ---------------- comparison values ----------------
    "LiYF4_comparison_at_8e20": {
        "W_ETU_4I13_2_cm3_s": {"LiYF4": 1.7e-17, "ZBLAN": 2.8e-17},
        "W_ETU_4I11_2_cm3_s": {"LiYF4": 1.6e-17, "ZBLAN": 1.0e-17},
        "note": "LiYF4 from Ref. 37 (Jensen thesis) alpha/2; 8e20 cm^-3 = 5.8 at.% LiYF4 = 5 mol% ZBLAN",
        "src": "Sec. IV.B p. 861",
    },

    # ---------------- oscillator strengths (Table II) ----------------
    "oscillator_strengths_x1e8": {
        "src": "Table II p. 861 (from Ref. 39, Wetenkamp thesis 1991); caption '[10^8]' "
               "= values are f x 1e8. 'YLF' column reproduced as printed (meaning not explained).",
        "Er3+": {
            "4I15/2-4I13/2": {"ZBLAN": 100, "YLF": 135.8},
            "4I15/2-4I11/2": {"ZBLAN": 44, "YLF": 126.7},
            "4I15/2-4I9/2": {"ZBLAN": 33, "YLF": 54.3},
            "4I13/2-4I9/2": {"ZBLAN": 45, "YLF": 48.5},
            "4I13/2-4S3/2": {"ZBLAN": 42, "YLF": 474.2},
            "4I13/2-2H11/2": {"ZBLAN": 30, "YLF": 58.2},
            "4I11/2-4F7/2": {"ZBLAN": 44, "YLF": 129.6, "note": "assumed = 4I15/2-4I11/2 (asterisk)"},
            "4I9/2-4S3/2": {"ZBLAN": 34, "YLF": 52.1},
            "4I9/2-2H11/2": {"ZBLAN": 95, "YLF": 43.8},
        },
        "Pr3+": {"3H4-3F3": 840, "3H4-1G4": 21},
        "products_x1e16_as_quoted": {  # p. 861; printed as 'x10^16' (should be 1e-16 in f units)
            "ETU1": 4500, "ETU2": 1963, "ET_4I13_2_to_3F3": 84000, "ET_4I11_2_to_1G4": 924,
            "CR_combined": 2133,
        },
        "products_note": "ETU2 quoted 1963 whereas 44*44 = 1936 (likely typo).",
        "ranking_Table_III": {
            "measured": ["ET 4I13/2", "ETU 4I13/2", "CR 4S3/2", "ET 4I11/2", "ETU 4I11/2"],
            "calculated": ["ET 4I13/2", "ETU 4I13/2", "CR 4S3/2", "ETU 4I11/2", "ET 4I11/2"],
            "src": "Table III p. 861",
        },
    },

    "judd_ofelt": None,       # not given
    "cross_sections": None,   # not given (ESA at 979 nm only described as 'small', p. 857)
    "branching_ratios": None,  # not given (only CR fraction beta, Eq. 6)
    "multiphonon_rates": None,  # not given

    "energy_transfer": ENERGY_TRANSFER,
}
