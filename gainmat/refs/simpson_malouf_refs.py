"""Extracted numbers (no text or figures reproduced).

SIMPSON2008: D. A. Simpson, "Spectroscopy of thulium doped silica glass", PhD thesis, Victoria
University (Melbourne).  MCVD/solution-doped Tm fibres made with LPMC Nice.  Tm2O3 in ppm (taken
as weight ppm when converting to cm^-3), Al2O3 in mol%.
MALOUF2016: A. Malouf, O. Henderson-Sapir, M. Gorjan, D. J. Ottaway, "Numerical modeling of
3.5 um dual-wavelength pumped erbium doped mid-infrared fiber lasers", arXiv:1606.07842 (IEEE JQE
2016), Tables I-IV: parameters of their model of the Henderson-Sapir 2014 laser (H2014 column).
"""

SIMPSON2008 = dict(
    src="Simpson, PhD thesis, Victoria University, Tables 4.1-4.7",
    # Table 4.2 / 4.3: peak absorption wavelength (nm), FWHM (nm) from 3H6
    abs_peaks_aluminosilicate={  # sample: {upper: (peak_nm, fwhm_nm)}, Al2O3 mol%
        "Al-2 (3.4 % Al2O3)": {"3F2,3": (681.0, 27.5), "3H4": (787.0, 19.0), "3H5": (1209.5, 10.5), "3F4": (1632.5, 181.5)},
        "Al-5 (4.8 % Al2O3)": {"3F2,3": (682.0, 26.5), "3H4": (787.0, 18.0), "3H5": (1208.5, 10.0), "3F4": (1641.5, 182.0)},
        "Al-7 (9.0 % Al2O3)": {"3F2,3": (682.0, 24.5), "3H4": (788.0, 19.0), "3H5": (1208.5, 10.0), "3F4": (1660.5, 179.5)},
    },
    abs_peaks_germanosilicate={
        "Ge-1": {"3F4": (1576.5, 160.5)}, "Ge-2": {"3H5": (1209.5, 7.5), "3F4": (1580.5, 168.0)},
    },
    note_3H5_fwhm="FWHM of the narrow peak near 1208 nm only",
    # Tables 4.5-4.7: 1/e lifetimes (us)
    lifetimes_us=[  # (sample, Tm2O3 ppm, Al2O3 mol%, tau_3H4, tau_3F4)
        ("Al-1", 100, 1.9, 16.8, 384), ("Al-2", 240, 3.4, 20.0, 417), ("Al-3", 450, 4.4, 20.7, 428),
        ("Al-4", 1850, 4.6, 21.0, 441), ("Al-5", 220, 4.8, 22.6, 448), ("Al-6", 1950, 5.4, 21.2, 447),
        ("Al-7", 250, 9.0, 32.5, 573), ("Al-9 (+Sb)", 330, 3.6, 23.2, 477), ("Al-10 (+Sn)", 260, 2.7, 20.7, 445),
        ("Ge-1 (+Sb)", 300, None, 21.0, 479), ("Ge-2 (+Sn)", 1400, None, 20.1, 467),
    ],
)

MALOUF2016 = dict(
    src="Malouf et al., arXiv:1606.07842, Tables I-IV, column H2014 (Henderson-Sapir 2014 laser)",
    Er_mol_pct=1.7, core_radius_um=5.25, NA=0.150, L_m=0.18, loss_per_m_3p5um=0.060, R_HR=0.99, R_OC=0.87,
    sigma_cm2={  # 1e-26 m^2 = 1e-22 cm^2
        ("GSA", "4I15/2", "4I11/2", 985): 9.30e-22,
        ("ESA", "4I11/2", "4F7/2", 985): 2.00e-22,
        ("ESA", "4I11/2", "4F9/2", 1973): 3.00e-21,
        ("EM", "4F9/2", "4I11/2", 1973): 3.61e-21,
        ("EM", "4F9/2", "4I9/2", 3470): 1.20e-21,
    },
    lifetimes_s={"4I13/2": 9.9e-3, "4I11/2": 7.9e-3, "4I9/2": 8.0e-6, "4F9/2": 177e-6, "4S3/2": 530e-6, "4F7/2": 5.0e-6},
    energy_transfer_m3_s={"W2214": 0.40e-24, "W3317": 0.08e-24, "W6142": 0.10e-24, "W5362": 17.0e-24},
)

NEWELL2007 = dict(
    src="Newell, Peterson, Gavrielides & Sharma, Opt. Commun. 273, 256 (2007), Fig. 1 and text",
    host="Yb-doped silica fibre core (composition not stated)", T_K=77,
    levels_2F7_2_cm1=[0, 425, 530, 865],             # from 77 K emission lines 977, 1020, 1032, 1069 nm
    gap_2F5_2_ab_cm1=630,                            # lifetime-vs-T fit and emission peaks
    zpl_nm=977, absorption_lines_nm=[977, 960, 920],
    tau_a_ms=0.884, tau_b_ms=0.547,
)

WANG2014 = dict(
    src="Wang et al., J. Mater. Chem. C (2014), Yb/Al/P co-doped silica (Fig. 7 inset; values not tabulated)",
    zpl_nm={"series A, B": 975, "series C": 976}, abs_1to7_nm=915,
)

LINCOLN1991 = dict(
    src="Lincoln, Brocklesby, Cusso, Townsend, Tropper & Pearson, J. Lumin. 50, 297 (1991)",
    label_note="uses OLD labels: their '3F4' = modern 3H4 (~12600 cm^-1), their '3H4' = modern 3F4",
    fln_9K_GeSiO2_3H4_to_3H6_peaks_nm=[787, 803, 815],  # first split into a doublet; >=2 weaker peaks beyond 820 nm
    note="ground 3H6 Stark components ~0, ~250, ~430 cm^-1 and more above ~500 (Ge-silica); Al-silica much less resolved",
    homogeneous_fwhm_3H4_9K_cm1_upper=10,
    tau_3H4_us=20, tau_3H4_concentration_dependence="none (20-12000 ppm)",
)
