"""Measured total Stark splittings of rare-earth ground (and some excited) multiplets in glasses.

Only extracted numbers with their sources; no figures or text are reproduced.  'family' maps
each glass to the forward model's closest glass family (see notes).  Values are total
splittings (highest minus lowest Stark level of the multiplet) in cm^-1.
"""

STARK_SPLITTINGS = [
    # ---------------- Yb3+ 2F7/2, room temperature, Lorentzian decomposition of abs + emission
    # Yang, Liu, Wang, Zhang, Hu & Zhang, Opt. Lett. 39, 1772 (2014), Table 2 ("Ed", highest
    # Stark level of 2F7/2).  Compositions (mol%) from their Table 1.
    dict(ion="Yb3+", level="2F7/2", value=901, sigma=50, family="germanate", glass="55GeO2-20PbO-15BaO-10ZnO",
         src="Yang et al., OL 39, 1772 (2014) Table 2", T="RT"),
    dict(ion="Yb3+", level="2F7/2", value=874, sigma=50, family="aluminosilicate", glass="55SiO2-30PbO-15Al2O3",
         src="Yang et al., OL 39, 1772 (2014) Table 2", T="RT",
         note="lead aluminosilicate; mapped to the aluminosilicate family"),
    dict(ion="Yb3+", level="2F7/2", value=722, sigma=50, family="tellurite", glass="70TeO2-17WO3-8La2O3-5Na2O",
         src="Yang et al., OL 39, 1772 (2014) Table 2", T="RT"),
    dict(ion="Yb3+", level="2F7/2", value=596, sigma=50, family="phosphate",
         glass="70P2O5-10H3BO3-10BaO-3Al2O3-3Nb2O5-4K2O", src="Yang et al., OL 39, 1772 (2014) Table 2", T="RT"),
    dict(ion="Yb3+", level="2F7/2", value=624, sigma=50, family="ZBLAN", glass="54ZrF4-20BaF2-20NaF-3LaF3-3AlF3",
         src="Yang et al., OL 39, 1772 (2014) Table 2", T="RT", note="fluorozirconate (ZBNLA), mapped to ZBLAN"),
    # borate 830, bismuthate 758, fluorophosphate 734 (Yang Table 2): no matching model family
    # Cao et al., J. Appl. Phys. 123, 215106 (2018), Table II: scalar field N_J; Delta E(7/2) = 0.245 N_J (their Eq. 2)
    dict(ion="Yb3+", level="2F7/2", value=0.245 * 3657, sigma=50, family="phosphosilicate",
         glass="93.98SiO2-3.4P2O5-2.5Al2O3 (SPA)", src="Cao et al., JAP 123, 215106 (2018) Table II (N_J x 0.245)", T="RT",
         note="Yb3+ in a local phosphorus environment in silica"),
    dict(ion="Yb3+", level="2F7/2", value=0.245 * 3257, sigma=50, family="phosphate",
         glass="71P2O5-27Al2O3 (PA)", src="Cao et al., JAP 123, 215106 (2018) Table II (N_J x 0.245)", T="RT",
         note="aluminophosphate"),
    # Zhang, Xue, He, Guzik & Boulon, Opt. Express 23, 1505 (2015): low-temperature levels
    dict(ion="Yb3+", level="2F7/2", value=589, sigma=30, family="phosphate", glass="Yb:phosphate laser glass",
         src="Zhang et al., OE 23, 1505 (2015) p. 1507", T="LT"),
    # Pask et al., JSTQE 1, 2 (1995) Fig. 1(a), germanosilicate, labelled approximate.  Larger than every
    # modern silicate value above; kept for reference with a wide uncertainty.
    dict(ion="Yb3+", level="2F7/2", value=1490, sigma=300, family="germanosilicate", glass="Ge-doped silica fibre",
         src="Pask et al., JSTQE 1, 2 (1995) Fig. 1(a) (approximate)", T="?"),
    dict(ion="Yb3+", level="2F5/2", value=1370, sigma=300, family="germanosilicate", glass="Ge-doped silica fibre",
         src="Pask et al., JSTQE 1, 2 (1995) Fig. 1(a) (approximate)", T="?"),
    # ---------------- Er3+ 4I15/2, 4.2 K fluorescence line narrowing
    # Zemon et al., J. Appl. Phys. 69, 6799 (1991), Sec. III.A (total spread A..H)
    dict(ion="Er3+", level="4I15/2", value=335, sigma=25, family="ZBLAN", glass="fluorozirconate F114",
         src="Zemon et al., JAP 69, 6799 (1991) p. 6806", T="4.2 K"),
    dict(ion="Er3+", level="4I15/2", value=325, sigma=30, family="phosphate", glass="phosphate L28 (50P2O5-31K2O-15.5BaO)",
         src="Zemon et al., JAP 69, 6799 (1991) p. 6807 ('300-350')", T="4.2 K"),
    dict(ion="Er3+", level="4I15/2", value=400, sigma=30, family="silica", glass="silicate L22 (70.5SiO2-12CaO-10.5Na2O-4K2O)",
         src="Zemon et al., JAP 69, 6799 (1991) p. 6807", T="4.2 K", note="multicomponent silicate"),
    # Desurvire & Simpson, OL 15, 547 (1990) Table 1 summary for silicate glasses (395-500)
    dict(ion="Er3+", level="4I15/2", value=445, sigma=50, family="silica", glass="silicate glasses (range 395-500)",
         src="Desurvire & Simpson, OL 15, 547 (1990) Table 1", T="4.2 K"),
    dict(ion="Er3+", level="4I13/2", value=330, sigma=40, family="silica", glass="silicate glasses (range 305-357)",
         src="Desurvire & Simpson, OL 15, 547 (1990) Table 1", T="4.2 K"),
]
