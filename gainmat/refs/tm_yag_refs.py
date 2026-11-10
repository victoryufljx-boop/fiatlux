"""Tm3+:Y3Al5O12 (D2 site): measured Stark levels and fitted crystal-field parameters.
Gruber, Hills, Macfarlane, Morrison, Turner, Quarles, Kintz & Esterowitz, Phys. Rev. B 40, 9464 (1989),
Table V (experimental levels; '...' = not observed -> None) and Table V footnote b (fitted B_nm, Wybourne
normalisation C_nm = sqrt(4pi/(2n+1)) Y_nm, cm^-1; rms 11 cm^-1 over 66 levels).  Numbers only."""
GRUBER1989 = dict(
    src="Gruber et al., Phys. Rev. B 40, 9464 (1989), Table V",
    bkq={(2, 0): 474.0, (2, 2): 47.0, (4, 0): -213.0, (4, 2): -1571.0, (4, 4): -824.0,
         (6, 0): -984.0, (6, 2): -310.0, (6, 4): 591.0, (6, 6): -193.0},
    levels={  # absolute cm^-1 above the ground level; None = not observed
        "3H6": [0, 27, 216, 241, 247, 252, 588, 610, None, 690, None, 730, None],
        "3F4": [5556, 5736, 5832, 5901, 6041, 6108, 6170, 6224, 6233],
        "3H5": [8339, 8345, 8516, 8530, None, 8556, 8711, 8773, 8800, None, 8882],
        "3H4": [12607, 12679, 12747, None, 12824, None, 13072, 13139, 13159],
    },
    calc_levels={  # Gruber's own calculated levels (Table V), for checking our Hamiltonian
        "3H6": [-3, 24, 215, 225, 262, 253, 518, 610, 650, 684, 698, 751, 765],
        "3F4": [5536, 5757, 5810, 5912, 6040, 6111, 6164, 6228, 6243],
        "3H5": [8343, 8353, 8507, 8517, 8524, 8558, 8712, 8774, 8804, 8873, 8885],
        "3H4": [12614, 12677, 12749, 12818, 12832, 12954, 13070, 13127, 13159],
    },
    note="levels 588/610/690/730 of 3H6 from 80 K emission (tentative symmetries); several sites present",
)
