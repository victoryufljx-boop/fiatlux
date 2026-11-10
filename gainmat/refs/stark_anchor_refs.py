"""Measured Stark levels used by the measured-level mode (spectra.ensemble_spectra(level_anchor=...)).
Numbers only, with their sources.  Levels are cm^-1 above the multiplet's lowest level, one entry per
Kramers doublet for Er/Yb (None = not resolved); E0 = lowest level of the multiplet above the lowest
ground level.  Each anchor comes from ONE source and is tested on data from OTHER sources
(scripts/validate_level_anchor.py)."""
ANCHORS = {
    "Yb3+": dict(
        levels={"2F7/2": dict(levels=[0, 425, 530, 865]), "2F5/2": dict(levels=[0, 630], E0=1e7 / 977.0)},
        families=("silica", "aluminosilicate", "germanosilicate"),
        src="Newell, Peterson, Gavrielides & Sharma, Opt. Commun. 273, 256 (2007): 77 K emission lines "
            "977, 1020, 1032, 1069 nm; 2F5/2 gap 630 cm^-1 (Yb silica fibre, composition not stated)"),
    # Zero-phonon line only (no level pattern: the 2F7/2 splitting differs by host, e.g. phosphate
    # ~590 vs silica ~870 cm^-1, Yang 2014 / Zhang 2015): Yb 2F5/2 lowest level at 977 nm (Newell),
    # applied to the other oxide families.  Held-out check: Wang et al. 2017 phosphate, 975 nm.
    "Yb3+:zpl": dict(
        ion="Yb3+", levels={"2F5/2": dict(levels=[0], E0=1e7 / 977.0)},
        families=("phosphosilicate", "phosphate", "germanate", "tellurite"),
        src="Newell et al., Opt. Commun. 273, 256 (2007): Yb zero-phonon line 977 nm (level pattern not "
            "transferred to non-silicate hosts)"),
    "Er3+": dict(
        levels={"4I15/2": dict(levels=[0, 30, 55, 129, 234, None, None, 420]),
                "4I13/2": dict(levels=[0, None, None, None, None, None, 330])},
        families=("silica", "aluminosilicate", "germanosilicate"),
        src="Desurvire & Simpson, Opt. Lett. 15, 547 (1990) (4.2 K FLN, silicate glasses): 4I15/2 lowest "
            "five levels and total splitting 395-500 (420 used); 4I13/2 total 305-357 (330 used)"),
    "Tm3+": dict(
        levels={"3H6": dict(levels=[0, 250, 430])},
        families=(), auto=False,     # partial (3 of 13 levels): tested, makes spectra worse, not used
        src="Lincoln et al., J. Lumin. 50, 297 (1991): 9 K FLN in Ge-silica, ground 3H6 components "
            "~0, ~250, ~430 cm^-1 (more above ~500, not resolved)"),
}

# Yb3+ full level patterns per glass family from Yang et al., Opt. Lett. 39, 1772 (2014), Fig. 2 (room-
# temperature Lorentzian decomposition; one glass per family).  Held-out checks: Wang et al. 2017 phosphate
# cross sections at 975 and 1005 nm, and the open-access benchmark.
from .yang2014_yb_refs import YANG2014 as _Y, SRC as _YSRC
for _k, _fam in (("Te", "tellurite"), ("Ge", "germanate"), ("P", "phosphate"), ("F", "ZBLAN")):
    _v = _Y[_k]
    ANCHORS[f"Yb3+:{_fam}"] = dict(
        ion="Yb3+", families=(_fam,), auto=False,   # held-out test (Wang 2017 phosphate 1005/975 nm) worse than ZPL-only

        levels={"2F7/2": dict(levels=list(_v["f72"])),
                "2F5/2": dict(levels=[0, _v["f52"][1] - _v["f52"][0], _v["f52"][2] - _v["f52"][0]], E0=float(_v["f52"][0]))},
        src=f"{_YSRC}; glass {_v['comp']}")
