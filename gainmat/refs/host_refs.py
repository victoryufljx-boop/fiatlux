"""Experimental reference values for the host tier, with citations.

Each quantity is a dict:
    value    : number (eV, dimensionless n, or cm^-1)
    kind     : what exactly is measured (e.g. "exciton peak", "optical gap",
               "n_D at 589 nm", "max phonon (multiphonon-law fit)", "highest LO")
    source   : citation (authors, journal, year / handbook)
    verified : True  -> the number was checked against the cited source (or a
                        page quoting it) during this build;
               False -> the number is a widely-quoted literature value that was
                        NOT re-checked against the primary source in this
                        build; treat as provisional and confirm before using
                        it for calibration.
Do not add values without a source.  Do not flip `verified` to True without
reading the source.

Band-gap convention: for wide-gap ionic crystals the first absorption peak is
an exciton, and the fundamental (transport/quasiparticle) gap lies ~1 eV above
it.  GLLB-SC targets the fundamental gap, so `gap_fund_eV` is the comparison
value; `gap_exciton_eV` is the optical edge and sets the practical UV edge.
"""
from __future__ import annotations


def _r(value, kind, source, verified=False, note=None):
    d = dict(value=value, kind=kind, source=source, verified=verified)
    if note:
        d["note"] = note
    return d


HOST_REFS = {
    # ------------------------------------------------------------------ CaF2
    "CaF2": dict(
        gap_exciton_eV=_r(
            11.11, "first exciton absorption peak (room T)",
            "M. Letz, L. Parthier, A. Gottwald, M. Richter, 'Spatial anisotropy "
            "of the exciton level in CaF2 at 11.1 eV ...', arXiv:cond-mat/0306510 "
            "(2003); fitted 11.111 eV (100) / 11.131 eV (111)",
            verified=True),
        gap_fund_eV=_r(
            12.1, "fundamental (band-to-band) gap from VUV reflectance",
            "G.W. Rubloff, Phys. Rev. B 5, 662 (1972); T. Tomiki, T. Miyata, "
            "J. Phys. Soc. Jpn. 27, 658 (1969)",
            note="commonly quoted 11.8-12.1 eV"),
        n=_r(1.4338, "n_D at 589.3 nm, room T",
             "I.H. Malitson, Appl. Opt. 2, 1103 (1963)"),
        phonon_max_cm1=_r(
            466, "effective max phonon used in energy-gap-law literature",
            "CaF2 Raman T2g mode 322 cm^-1 and LO ~ 463-474 cm^-1; "
            "466 cm^-1 as tabulated e.g. in M.J. Weber (ed.), Handbook of Laser "
            "Science and Technology / CRC Handbook of Optical Materials (2003)"),
        phonon_LO_cm1=_r(
            472, "highest LO phonon (IR reflectivity / neutron)",
            "W. Kaiser, W.G. Spitzer, R.H. Kaiser, L.E. Howarth, Phys. Rev. 127, "
            "1950 (1962): TO 257, LO ~ 463-474 cm^-1",
            note="literature LO range 463-482 cm^-1"),
    ),
    # ----------------------------------------------------------------- LiYF4
    "LiYF4": dict(
        gap_exciton_eV=_r(
            10.5, "VUV absorption edge / exciton (~118 nm)",
            "Commonly quoted for undoped YLF (e.g. P. Dorenbos, J. Lumin. 91, 155 "
            "(2000), Table of host gaps; M. Kirm et al. VUV studies)",
            note="literature values 10.5-11.8 eV depending on definition"),
        gap_fund_eV=_r(
            11.8, "fundamental gap (exciton + binding energy)",
            "P. Dorenbos, J. Lumin. 91, 155 (2000) (E_ex x 1.08 rule)"),
        n=_r(1.4536, "n_o at 589 nm (n_e = 1.4758)",
             "N.P. Barnes, D.J. Gettemy, J. Opt. Soc. Am. 70, 1244 (1980)"),
        phonon_max_cm1=_r(
            460, "effective phonon in energy-gap law",
            "D.N. Patel et al. / A. Ellens et al. ; YLF often quoted 400-490 cm^-1",
            note="highest Raman mode 425-446 cm^-1; IR LO up to ~560 cm^-1 "
                 "(S.A. Miller, H.E. Rast, H.H. Caspers, J. Chem. Phys. 52, 4172 "
                 "(1970))"),
        phonon_LO_cm1=_r(
            560, "highest LO phonon (IR reflectance)",
            "S.A. Miller, H.E. Rast, H.H. Caspers, J. Chem. Phys. 52, 4172 (1970)"),
    ),
    # ----------------------------------------------------------- alpha-Al2O3
    "Al2O3": dict(
        gap_exciton_eV=_r(
            9.0, "exciton absorption peak (~9.0-9.1 eV)",
            "R.H. French, J. Am. Ceram. Soc. 73, 477 (1990)"),
        gap_fund_eV=_r(
            8.8, "optical band gap from VUV spectroscopy (8.8 eV, 300 K)",
            "R.H. French, J. Am. Ceram. Soc. 73, 477 (1990)",
            note="many papers quote 8.8-9.4 eV"),
        n=_r(1.768, "n_o at 589 nm (n_e = 1.760)",
             "I.H. Malitson, J. Opt. Soc. Am. 52, 1377 (1962)"),
        phonon_max_cm1=_r(
            870, "effective max phonon (energy-gap law)",
            "Widely quoted for sapphire; e.g. M.J. Weber, Phys. Rev. 157, 262 "
            "(1967) multiphonon fits",
            note="range 750-900 cm^-1 used in literature"),
        phonon_LO_cm1=_r(
            907, "highest LO phonon (Eu, E||c ~ 881; E_u LO 907 cm^-1)",
            "A.S. Barker, Phys. Rev. 132, 1474 (1963)"),
    ),
    # --------------------------------------------------------------- Y3Al5O12
    "Y3Al5O12": dict(
        gap_exciton_eV=_r(
            7.0, "exciton / absorption edge (~6.5-7.0 eV, ~180 nm)",
            "P. Dorenbos, J. Lumin. 91, 155 (2000); Y.N. Xu, W.Y. Ching, "
            "Phys. Rev. B 59, 10530 (1999) (expt. comparison)"),
        gap_fund_eV=_r(
            7.5, "fundamental gap estimate",
            "P. Dorenbos, J. Lumin. 91, 155 (2000)",
            note="Not confirmed. Zhunusbekov et al., Crystals 14, 897 (2024) "
                 "(https://mdpi.com/2073-4352/14/10/897) quote gap >6.9-7.0 eV "
                 "(Zych 2000; Tanner 2007) and 7.85-7.95 eV (Zorenko 2018), "
                 "exciton 6.77 eV"),
        n=_r(1.8233, "n at 589 nm",
             "D.E. Zelmon, D.L. Small, R. Page, Appl. Opt. 37, 4933 (1998)"),
        phonon_max_cm1=_r(
            700, "effective max phonon (energy-gap law)",
            "Standard YAG value (~700 cm^-1), e.g. C.B. Layne, W.H. Lowdermilk, "
            "M.J. Weber, Phys. Rev. B 16, 10 (1977) fits; Weber handbook"),
        phonon_LO_cm1=_r(
            857, "highest Raman/IR optical mode",
            "J.P. Hurrell, S.P.S. Porto, I.F. Chang, S.S. Mitra, R.P. Bauman, "
            "Phys. Rev. 173, 851 (1968): highest Raman A1g 783, T2g 857 cm^-1"),
    ),
    # ---------------------------------------------------------- alpha-quartz
    "SiO2-alpha-quartz": dict(
        gap_exciton_eV=_r(
            10.3, "first exciton peak (~10.3 eV)",
            "H.R. Philipp, Solid State Commun. 4, 73 (1966)"),
        gap_fund_eV=_r(
            9.65, "optical gap (T.H. DiStefano, D.E. Eastman 9.0 eV "
            "photoconductivity edge; 9.65 eV commonly used)",
            "R.B. Laughlin, Phys. Rev. B 22, 3021 (1980)",
            note="8.9-9.7 eV depending on definition"),
        n=_r(1.5443, "n_o at 589 nm (n_e = 1.5534)",
             "G. Ghosh, Opt. Commun. 163, 95 (1999)"),
        phonon_max_cm1=_r(
            1100, "effective max phonon (Si-O stretch)",
            "Generic silicate value used in energy-gap law (Layne et al. 1977)"),
        phonon_LO_cm1=_r(
            1230, "highest LO phonon (E mode TO 1162 / LO 1230 cm^-1)",
            "F. Gervais, B. Piriou, Phys. Rev. B 11, 3944 (1975)"),
    ),
    # ----------------------------------------------------------- glasses
    "SiO2-glass": dict(
        gap_exciton_eV=_r(
            10.2, "exciton absorption peak of fused silica",
            "H.R. Philipp, J. Phys. Chem. Solids 32, 1935 (1971)"),
        gap_fund_eV=_r(
            9.0, "optical gap / band edge (~9 eV); practical UV edge 160 nm "
                 "(7.75 eV) limited by defects/OH",
            "H.R. Philipp, J. Phys. Chem. Solids 32, 1935 (1971)"),
        n=_r(1.4585, "n_D at 589.3 nm, 20 C",
             "I.H. Malitson, J. Opt. Soc. Am. 55, 1205 (1965)"),
        phonon_max_cm1=_r(
            1100, "maximum phonon energy (energy-gap-law value)",
            "X. Zhu, N. Peyghambarian, Adv. OptoElectron. 2010, 501956, Table 1 "
            "(silica: 1100 cm^-1); also arXiv:2110.06515 ('nominal maximum phonon "
            "energy of 1100 cm-1 associated with silicate glass'); originally "
            "attributed to Layne, Lowdermilk, Weber, PRB 16, 10 (1977) (not re-read)",
            verified=True,
            note="Raman: 1060 / 1200 cm^-1 Si-O stretches"),
        phonon_LO_cm1=_r(
            1250, "highest LO (IR, TO4/LO4 1065/1250 cm^-1)",
            "F.L. Galeener, G. Lucovsky, Phys. Rev. Lett. 37, 1474 (1976)"),
    ),
    "ZBLAN": dict(
        gap_exciton_eV=_r(
            None, "fundamental absorption edge not well-defined",
            "unverified"),
        gap_fund_eV=_r(
            6.5, "intrinsic UV absorption edge ~190-200 nm (6.2-6.5 eV); "
                 "practical edge ~ 250 nm (impurity limited)",
            "M. Poulain, in 'Fluoride glass fiber optics' (Aggarwal & Lu eds, "
            "Academic Press 1991); X. Zhu, N. Peyghambarian, Adv. OptoElectron. "
            "2010, 501956 (review)"),
        n=_r(1.499, "n at 0.589 um (Table 1)",
             "X. Zhu, N. Peyghambarian, Adv. OptoElectron. 2010, 501956, Table 1",
             verified=True),
        # CORRECTED 580 -> 600: Zhu & Peyghambarian Table 1 lists "Maximum phonon
        # energy (cm-1)" = 600 for ZBLAN (the earlier "565-580" quote was not in it).
        phonon_max_cm1=_r(
            600, "maximum phonon energy (Table 1)",
            "X. Zhu, N. Peyghambarian, Adv. OptoElectron. 2010, 501956, Table 1 "
            "(silica 1100, ZBLAN 600)", verified=True,
            note="Raman Zr-F symmetric stretch ~ 580 cm^-1 (not re-checked); "
                 "was 580 before verification"),
        phonon_LO_cm1=_r(
            580, "highest Raman band (Zr-F stretch)",
            "same as above"),
    ),
}


def get_ref(host: str, quantity: str):
    """Return the reference dict (or None)."""
    return HOST_REFS.get(host, {}).get(quantity)


def ref_value(host: str, quantity: str):
    r = get_ref(host, quantity)
    return None if r is None else r["value"]
