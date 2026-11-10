"""External benchmark from open-access papers (round 4).  Numbers stated in the text or tables of
each paper (no values read off figures), with citation and DOI; nothing here was used to build or
fit the model.  `family` = the model's closest glass family; proxy=True marks a mapping that is
only approximate (e.g. a bismuth-germanate or AlF3-fluoride glass run as 'germanate' / 'ZBLAN').

Record: dict(ion, family, proxy, pair=(upper, lower), q, value, src, note)
  q: em_peak_nm | em_fwhm_nm | sigma_em_cm2 | sigma_abs_cm2 | abs_peak_nm | tau_rad_ms | tau_ms | omegas
"""

CANKAYA2018 = "Cankaya et al., Appl. Sci. 8, 333 (2018), doi:10.3390/app8030333"
HUANG2014A = "Huang, Liu, Hu & Chen, Sci. Rep. 4, 5053 (2014), doi:10.1038/srep05053"
HUANG2014B = "Huang et al., Sci. Rep. 4, 3607 (2014), doi:10.1038/srep03607"
WANG2017 = "Wang et al., Materials 10, 241 (2017), doi:10.3390/ma10030241"
CAI2015 = "Cai et al., Opt. Mater. Express 5, 1431 (2015), doi:10.1364/OME.5.001431"
KOCH2019 = "Kochanowicz et al., Opt. Mater. Express 9, 1450 (2019), doi:10.1364/OME.9.001450"
CHOUD = "Choudhary, Kannan, Mackenzie, Feng & Shepherd, 'An ion-exchanged Tm3+:glass channel waveguide laser' (Southampton eprint 350893)"
TURRI2008 = "Turri et al., J. Appl. Phys. 103, 093104 (2008), doi:10.1063/1.2912952 (open-access copy, core.ac.uk 236316474)"
HAN2016 = "Han et al., Sci. Rep. 6, 31207 (2016), doi:10.1038/srep31207"
SEGURA = "Segura et al., 'Fabrication and spectroscopy of high-quality Tm3+-doped germanate glass for 2 um laser emission' (INRIM post-print)"
ZMOJDA2010 = "Zmojda, Dorosz, Kochanowicz & Dorosz, Photonics Lett. Pol. 2, 2 (2010), doi:10.4302/plp.2010.2.08"
CHANTHIMA2019 = "Chanthima et al., J. Met. Mater. Miner. 29 (2019), doi:10.14456/jmmm.2019.19"
BELL2014 = "Bell et al., JOSA B 31, 1590 (2014), doi:10.1364/JOSAB.31.001590"
MOREIRA2016 = "Moreira et al., Opt. Mater. (2016), doi:10.1016/j.optmat.2016.05.024 (open-access copy, UNESP)"
DANTAS2013 = "Dantas et al., Braz. J. Phys. (2013), doi:10.1007/s13538-013-0145-y"


def _r(ion, family, pair, q, value, src, proxy=False, note=""):
    return dict(ion=ion, family=family, proxy=proxy, pair=pair, q=q, value=value, src=src, note=note)


TM_F4 = ("3F4", "3H6"); TM_H4 = ("3H4", "3H6")
ER_13 = ("4I13/2", "4I15/2"); ER_11_13 = ("4I11/2", "4I13/2")
YB = ("2F5/2", "2F7/2"); HO = ("5I7", "5I8"); ND_11 = ("4F3/2", "4I11/2")
_CK = "1 mol% Tm2O3 sample (Table 1 Omegas). Table columns for 3F4/3H4 were swapped in extraction; the " \
      "1.8 um-band values are assigned to 3F4 by physics (tau_rad 2.6 ms is 3F4, not 3H4)"

RECORDS = [
    # ---- Tm tellurite (Cankaya 2018): three glasses, 1 mol% Tm
    _r("Tm3+", "tellurite", TM_F4, "omegas", (4.51, 0.76, 1.13), CANKAYA2018, note="TeO2-K2O-Nb2O5; " + _CK),
    _r("Tm3+", "tellurite", TM_F4, "tau_rad_ms", 2.57, CANKAYA2018, note="TeO2-K2O-Nb2O5; " + _CK),
    _r("Tm3+", "tellurite", TM_F4, "sigma_em_cm2", 5.61e-21, CANKAYA2018, note="TeO2-K2O-Nb2O5; " + _CK),
    _r("Tm3+", "tellurite", TM_F4, "em_peak_nm", 1860, CANKAYA2018, note="TeO2-K2O-Nb2O5"),
    _r("Tm3+", "tellurite", TM_F4, "em_fwhm_nm", 170, CANKAYA2018, note="TeO2-K2O-Nb2O5; " + _CK),
    _r("Tm3+", "tellurite", TM_F4, "omegas", (4.02, 0.93, 1.12), CANKAYA2018, note="TeO2-ZnO; " + _CK),
    _r("Tm3+", "tellurite", TM_F4, "tau_rad_ms", 2.55, CANKAYA2018, note="TeO2-ZnO; " + _CK),
    _r("Tm3+", "tellurite", TM_F4, "sigma_em_cm2", 5.48e-21, CANKAYA2018, note="TeO2-ZnO; " + _CK),
    _r("Tm3+", "tellurite", TM_F4, "em_fwhm_nm", 192, CANKAYA2018, note="TeO2-ZnO; " + _CK),
    _r("Tm3+", "tellurite", TM_F4, "omegas", (4.09, 0.69, 1.11), CANKAYA2018, note="TeO2-Nb2O5; " + _CK),
    _r("Tm3+", "tellurite", TM_F4, "tau_rad_ms", 2.76, CANKAYA2018, note="TeO2-Nb2O5; " + _CK),
    _r("Tm3+", "tellurite", TM_F4, "sigma_em_cm2", 4.99e-21, CANKAYA2018, note="TeO2-Nb2O5; " + _CK),
    _r("Tm3+", "tellurite", TM_F4, "em_peak_nm", 1795, CANKAYA2018, note="TeO2-Nb2O5"),
    _r("Tm3+", "tellurite", TM_F4, "em_fwhm_nm", 237, CANKAYA2018, note="TeO2-Nb2O5; " + _CK),
    _r("Tm3+", "tellurite", TM_H4, "abs_peak_nm", 794, CANKAYA2018, note="all TeO2 glasses"),
    # ---- Tm germanate family
    _r("Tm3+", "germanate", TM_F4, "omegas", (4.4, 1.2, 0.9), TURRI2008, note="Tm germanate (composition not disclosed)"),
    _r("Tm3+", "germanate", TM_F4, "tau_rad_ms", 5.3, TURRI2008),
    _r("Tm3+", "germanate", TM_F4, "em_peak_nm", 1870, TURRI2008),
    _r("Tm3+", "germanate", TM_F4, "em_peak_nm", 1865, CHOUD, note="65.5GeO2-12Al2O3-4.5BaO-18Na2O-1Tm2O3"),
    _r("Tm3+", "germanate", TM_F4, "sigma_em_cm2", 2.3e-21, CHOUD, note="McCumber; 1 mol% Tm2O3"),
    _r("Tm3+", "germanate", TM_H4, "abs_peak_nm", 790, CHOUD),
    _r("Tm3+", "germanate", TM_F4, "em_peak_nm", 1800, KOCH2019, note="60GeO2-15Ga2O3-10BaO-15Na2O, 0.7 mol% Tm2O3"),
    _r("Tm3+", "germanate", TM_H4, "abs_peak_nm", 792, SEGURA),
    _r("Tm3+", "germanate", TM_H4, "sigma_abs_cm2", 1.34e-20, SEGURA, note="2.6e20 cm^-3"),
    _r("Tm3+", "germanate", TM_F4, "omegas", (6.26, 0.24, 1.6), HAN2016, proxy=True, note="36Bi2O3-29GeO2-25Ga2O3-10Na2O (bismuth germanate)"),
    _r("Tm3+", "germanate", TM_F4, "tau_rad_ms", 2.25, HAN2016, proxy=True, note="bismuth germanate"),
    _r("Tm3+", "germanate", TM_F4, "em_peak_nm", 1865, HAN2016, proxy=True, note="bismuth germanate"),
    _r("Tm3+", "germanate", TM_F4, "sigma_em_cm2", 7.01e-21, HAN2016, proxy=True, note="bismuth germanate"),
    _r("Tm3+", "germanate", TM_H4, "abs_peak_nm", 790, HAN2016, proxy=True, note="bismuth germanate"),
    # ---- Tm silica / phosphate (Turri 2008)
    _r("Tm3+", "aluminosilicate", TM_F4, "omegas", (3.7, 2.3, 0.6), TURRI2008, note="10.3 wt% Al2O3 silica"),
    _r("Tm3+", "aluminosilicate", TM_F4, "tau_rad_ms", 7.0, TURRI2008),
    _r("Tm3+", "phosphate", TM_F4, "omegas", (3.6, 0.7, 1.0), TURRI2008, note="Kigre borophosphate"),
    _r("Tm3+", "phosphate", TM_F4, "tau_rad_ms", 8.3, TURRI2008),
    # ---- Er fluorides (Huang 2014 a, b)
    _r("Er3+", "ZBLAN", ER_13, "omegas", (3.27, 1.3, 1.75), HUANG2014A, note="ZBLAN, 1 mol% ErF3"),
    _r("Er3+", "ZBLAN", ER_13, "tau_rad_ms", 6.44, HUANG2014A),
    _r("Er3+", "ZBLAN", ER_11_13, "tau_rad_ms", 5.50, HUANG2014A, note="4I11/2 total radiative lifetime"),
    _r("Er3+", "ZBLAN", ER_13, "em_peak_nm", 1535, HUANG2014A),
    _r("Er3+", "ZBLAN", ER_13, "sigma_em_cm2", 7.26e-21, HUANG2014A, note="Table 3 alignment uncertain (extraction note)"),
    _r("Er3+", "ZBLAN", ER_13, "sigma_abs_cm2", 5.79e-21, HUANG2014A, note="Table 3 alignment uncertain (extraction note)"),
    _r("Er3+", "ZBLAN", ER_11_13, "em_peak_nm", 2708, HUANG2014A),
    _r("Er3+", "ZBLAN", ER_11_13, "sigma_em_cm2", 1.003e-20, HUANG2014A),
    _r("Er3+", "ZBLAN", ER_13, "omegas", (3.01, 1.03, 1.64), HUANG2014A, note="ZBYA fluorozirconate"),
    _r("Er3+", "ZBLAN", ER_13, "tau_rad_ms", 6.82, HUANG2014A, note="ZBYA"),
    _r("Er3+", "ZBLAN", ER_11_13, "sigma_em_cm2", 8.87e-21, HUANG2014A, note="ZBYA"),
    _r("Er3+", "ZBLAN", ER_13, "omegas", (2.70, 1.68, 1.18), HUANG2014B, note="53ZrF4-20BaF2-4LaF3-3AlF3-20NaF"),
    _r("Er3+", "ZBLAN", ER_13, "tau_rad_ms", 8.29, HUANG2014B),
    _r("Er3+", "ZBLAN", ER_11_13, "tau_rad_ms", 7.74, HUANG2014B),
    _r("Er3+", "ZBLAN", ER_11_13, "em_peak_nm", 2710, HUANG2014B),
    _r("Er3+", "ZBLAN", ER_11_13, "sigma_em_cm2", 9.5e-21, HUANG2014B),
    _r("Er3+", "ZBLAN", ER_13, "omegas", (1.56, 1.10, 1.05), HUANG2014B, proxy=True, note="AlF3-YF3 fluoroaluminate AYF3"),
    _r("Er3+", "ZBLAN", ER_13, "tau_rad_ms", 9.6, HUANG2014B, proxy=True, note="AYF3"),
    _r("Er3+", "ZBLAN", ER_11_13, "sigma_em_cm2", 8.8e-21, HUANG2014B, proxy=True, note="AYF3"),
    # ---- Er phosphates
    _r("Er3+", "phosphate", ER_13, "sigma_abs_cm2", 6.8e-21, ZMOJDA2010, proxy=True, note="antimony-phosphate"),
    _r("Er3+", "phosphate", ER_13, "sigma_em_cm2", 7.74e-21, ZMOJDA2010, proxy=True, note="antimony-phosphate"),
    _r("Er3+", "phosphate", ER_13, "em_peak_nm", 1536, ZMOJDA2010, proxy=True),
    _r("Er3+", "phosphate", ER_13, "abs_peak_nm", 1537, CHANTHIMA2019, note="15Na2O-20Al2O3-63P2O5, 2 mol% Er2O3"),
    _r("Er3+", "phosphate", ER_13, "em_peak_nm", 1536, CHANTHIMA2019),
    # ---- Yb phosphates (Wang 2017): BaO-P2O5-SiO2 with P2O5 >= 53 mol%
    _r("Yb3+", "phosphate", YB, "abs_peak_nm", 975, WANG2017, note="20BaO-71P2O5-9SiO2 (BPS9)"),
    _r("Yb3+", "phosphate", YB, "sigma_abs_cm2", 1.09e-20, WANG2017, note="BPS9"),
    _r("Yb3+", "phosphate", YB, "sigma_em_cm2", 1.46e-20, WANG2017, note="BPS9"),
    _r("Yb3+", "phosphate", YB, "tau_ms", 1.17, WANG2017, note="BPS9, 1 mol% Yb2O3"),
    _r("Yb3+", "phosphate", YB, "sigma_abs_cm2", 0.98e-20, WANG2017, note="BPS16"),
    _r("Yb3+", "phosphate", YB, "sigma_em_cm2", 1.30e-20, WANG2017, note="BPS16"),
    _r("Yb3+", "phosphate", YB, "sigma_abs_cm2", 0.85e-20, WANG2017, note="BPS26.7"),
    _r("Yb3+", "phosphate", YB, "sigma_em_cm2", 1.13e-20, WANG2017, note="BPS26.7"),
    # ---- Ho germanate (Cai 2015, Yb co-doped; fluoro-germanate-gallate)
    _r("Ho3+", "germanate", HO, "em_peak_nm", 2000, CAI2015, proxy=True, note="'centred at about 2 um'"),
    _r("Ho3+", "germanate", HO, "sigma_em_cm2", 4.9e-21, CAI2015, proxy=True),
    _r("Ho3+", "germanate", HO, "em_fwhm_nm", 171, CAI2015, proxy=True),
    _r("Ho3+", "germanate", HO, "tau_rad_ms", 1e3 / 103.38, CAI2015, proxy=True, note="from A_rad = 103.38 s^-1"),
    # ---- Nd tellurite (Bell 2014, Moreira 2016) and phosphate (Dantas 2013)
    _r("Nd3+", "tellurite", ND_11, "omegas", (3.09, 3.11, 3.25), BELL2014, note="85TeO2-15ZnO wt%, n = 2.0"),
    _r("Nd3+", "tellurite", ND_11, "tau_rad_ms", 0.2175, BELL2014),
    _r("Nd3+", "tellurite", ND_11, "em_peak_nm", 1062.5, BELL2014),
    _r("Nd3+", "tellurite", ND_11, "sigma_em_cm2", 3.14e-20, BELL2014),
    _r("Nd3+", "tellurite", ND_11, "em_fwhm_nm", 29.3, BELL2014),
    _r("Nd3+", "tellurite", ND_11, "omegas", (3.865, 4.018, 3.792), MOREIRA2016, note="TZO 0.5 wt% Nd2O3"),
    _r("Nd3+", "tellurite", ND_11, "tau_rad_ms", 0.198, MOREIRA2016, note="paper's '198 ms' read as us"),
    _r("Nd3+", "tellurite", ND_11, "sigma_em_cm2", 4.2e-20, MOREIRA2016),
    _r("Nd3+", "tellurite", ND_11, "em_fwhm_nm", 24, MOREIRA2016),
    _r("Nd3+", "tellurite", ND_11, "em_peak_nm", 1062, MOREIRA2016),
    _r("Nd3+", "phosphate", ND_11, "omegas", (2.62, 2.87, 2.92), DANTAS2013, note="40P2O5-20Al2O3-35Na2O-5K2O"),
    _r("Nd3+", "phosphate", ND_11, "tau_ms", 0.385, DANTAS2013, note="dilute-limit fit tau_0"),
]
