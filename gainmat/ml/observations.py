"""Observation table for the host-microstructure layer (L2).

Each observation is something the site-ensemble spectra must reproduce and that depends on
the host's microstructure (crystal-field strength, site disorder, nephelauxetic scatter,
homogeneous width) rather than on the Judd-Ofelt intensities:

  eff_width_nm   effective emission/absorption width  int(sigma d lambda) / sigma_max
                 (Wetenkamp 1992 Table 2 uses exactly this width in its Fuchtbauer-Ladenburg
                 cross sections; for digitised curves it is computed the same way)
  curve          a digitised cross-section curve: peak-normalised shape, effective width, peak
  split_cm1      total Stark splitting of a multiplet
  zpl_fwhm_cm1   inhomogeneous width of the lowest-to-lowest Stark line (low temperature)
  ratio_log10    log10(sigma_e / sigma_a) at given wavelengths (McCumber energy)

Every entry carries the source.  Only public, repository data are used here.
"""
from __future__ import annotations

import numpy as np


def _arr(x):
    return np.array([np.nan if v is None else v for v in x], float)


def build(include=("ZBLAN", "aluminosilicate", "germanosilicate", "silica")):
    from ..refs.zblan_yb_refs import WETENKAMP1992, PASK1995
    from ..refs.tm_silica_refs import AGGER2006
    from ..refs.midir_ho_refs import LEGOUET2020, JACKSON2012
    from ..refs.edfa_refs import GILES1991
    obs = []

    def add(**kw):
        if kw["family"] in include:
            obs.append(kw)

    # --- Wetenkamp 1992 Table 2: effective emission widths in ZBLAN (4 ions, 15 transitions)
    T2 = WETENKAMP1992["table2_emission_ZBLAN"]
    src = "Wetenkamp, West & Toebben, JNCS 140, 35 (1992) Table 2"
    for ion_el in ("Nd3+", "Ho3+", "Er3+", "Tm3+"):
        for tr, row in T2.get(ion_el, {}).items():
            lam, A, beta, dl = row[0], row[1], row[2], row[3]
            if dl is None:
                continue
            up, lo = [x.strip() for x in tr.split("->")]
            # '5S2,5F4' stays a thermally coupled group (Boltzmann-weighted in host_calib._band)
            add(family="ZBLAN", ion=ion_el, kind="eff_width_nm", pair=(up, lo), band_nm=(lam * 0.85, lam * 1.18),
                center_nm=lam, value=float(dl), sigma_rel=0.15, src=src, which="em")
    # --- Jackson 2012 Fig. 4 emission curves (ZBLAN and aluminosilicate)
    J = JACKSON2012["FIG4_EMISSION_XS_CURVES_um_1e-25m2"]
    jm = {"Ho_5I7-5I8_aluminosilicate": ("Ho3+", "aluminosilicate", ("5I7", "5I8")),
          "Tm_3F4-3H6_aluminosilicate": ("Tm3+", "aluminosilicate", ("3F4", "3H6")),
          "Er_4I13/2-4I15/2_aluminosilicate": ("Er3+", "aluminosilicate", ("4I13/2", "4I15/2")),
          "Er_4I11/2-4I13/2_ZBLAN": ("Er3+", "ZBLAN", ("4I11/2", "4I13/2")),
          "Ho_5I6-5I7_ZBLAN": ("Ho3+", "ZBLAN", ("5I6", "5I7")),
          "Dy_6H13/2-6H15/2_ZBLAN": ("Dy3+", "ZBLAN", ("6H13/2", "6H15/2"))}
    for key, (ion, fam, pair) in jm.items():
        d = np.array(J[key], float)
        add(family=fam, ion=ion, kind="curve", pair=pair, which="em", lam_nm=d[:, 0] * 1e3, sigma=d[:, 1] * 1e-21,
            sigma_rel=0.15, src="Jackson, Nat. Photon. 6, 423 (2012) Fig. 4 (digitised)")
    # --- Agger 2006 Tm aluminosilicate (fibre Tm1)
    ag = AGGER2006["cross_sections"]["spectra_digitized_1e-25_m2"]
    for k, which in (("Tm1_sigma_a", "abs"), ("Tm1_sigma_e", "em")):
        add(family="aluminosilicate", ion="Tm3+", kind="curve", pair=("3F4", "3H6"), which=which,
            lam_nm=np.array(ag["lambda_nm"], float), sigma=_arr(ag[k]) * 1e-21, sigma_rel=0.1,
            src="Agger & Povlsen, OE 14, 50 (2006) Fig. 4")
    # --- Le Gouet 2020 Ho aluminosilicate
    lg = np.array(LEGOUET2020["XS_CURVE_DIGITIZED_nm_abs_em_1e-25m2"], float)
    for col, which in ((1, "abs"), (2, "em")):
        add(family="aluminosilicate", ion="Ho3+", kind="curve", pair=("5I7", "5I8"), which=which,
            lam_nm=lg[:, 0], sigma=lg[:, col] * 1e-21, sigma_rel=0.1, src="Le Gouet et al., arXiv:2003.14188")
    # --- Giles 1991 Fig. 2: Er alpha and g* (shape only; scale unknown) and their ratio
    sp2 = GILES1991["spectra_fig2"]; lam = np.array(sp2["wavelength_nm"], float)
    for key, fam in (("a_Ge_silicate", "germanosilicate"), ("b_Al_Ge_silicate", "aluminosilicate")):
        a = _arr(sp2[key]["alpha_dB_per_m"]); g = _arr(sp2[key]["gstar_dB_per_m"])
        src = "Giles & Desurvire, JLT 9, 271 (1991) Fig. 2"
        add(family=fam, ion="Er3+", kind="curve", pair=("4I13/2", "4I15/2"), which="abs", lam_nm=lam, sigma=a,
            shape_only=True, sigma_rel=0.1, src=src)
        add(family=fam, ion="Er3+", kind="curve", pair=("4I13/2", "4I15/2"), which="em", lam_nm=lam, sigma=g,
            shape_only=True, sigma_rel=0.1, src=src)
        for lw in (1500.0, 1531.0, 1550.0, 1565.0):
            i = int(np.argmin(abs(lam - lw)))
            add(family=fam, ion="Er3+", kind="ratio_log10", pair=("4I13/2", "4I15/2"), lam_point=float(lam[i]),
                value=float(np.log10(g[i] / a[i])), sigma_abs=0.05, src=src)
    # --- Pask 1995 Yb germanosilicate: curve and Stark splittings
    pk = PASK1995["cross_sections"]["curve_DIGITIZED"]
    for k, which in (("sigma_abs_1e-24m2", "abs"), ("sigma_em_1e-24m2", "em")):
        add(family="germanosilicate", ion="Yb3+", kind="curve", pair=("2F5/2", "2F7/2"), which=which,
            lam_nm=np.array(pk["wavelength_nm"], float), sigma=_arr(pk[k]) * 1e-20, sigma_rel=0.15,
            src="Pask et al., JSTQE 1, 2 (1995) Fig. 1(b)")
    for lab, v in (("2F7/2", 1490.0), ("2F5/2", 1370.0)):
        add(family="germanosilicate", ion="Yb3+", kind="split_cm1", label=lab, value=v, sigma_abs=150.0,
            src="Pask et al. 1995 Fig. 1(a) ('approximate')")
    # --- Desurvire & Simpson 1990, Er aluminosilicate: ZPL inhomogeneous width; other silicates' splittings
    add(family="aluminosilicate", ion="Er3+", kind="zpl_fwhm_cm1", upper="4I13/2", lower="4I15/2", value=49.0,
        sigma_abs=8.0, src="Desurvire & Simpson, OL 15, 547 (1990), 4.2 K")
    add(family="silica", ion="Er3+", kind="split_cm1", label="4I15/2", value=445.0, sigma_abs=60.0,
        src="Desurvire 1990 Table 1: silicate glasses 395-500")
    add(family="silica", ion="Er3+", kind="split_cm1", label="4I13/2", value=330.0, sigma_abs=40.0,
        src="Desurvire 1990 Table 1: silicate glasses 305-357")
    # --- Agger: sigma_e/sigma_a anchors for Tm aluminosilicate
    for lw, r in ((1740.0, 1.20), (1600.0, 0.10)):
        add(family="aluminosilicate", ion="Tm3+", kind="ratio_log10", pair=("3F4", "3H6"), lam_point=lw,
            value=float(np.log10(r)), sigma_abs=0.08, src="Agger & Povlsen 2006 p. 53")
    return obs
