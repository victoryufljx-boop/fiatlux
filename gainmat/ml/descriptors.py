"""Physically meaningful composition descriptors for a host glass.

From a composition {compound formula: mol%} (dopant compound removed):
  x_<anion>          anion fractions (O, F, S, Se, Cl, Br, I)
  basicity           optical basicity (Duffy): charge-weighted mean of 1 / (1.36 (chi - 0.26))
                     over cations, chi = Pauling electronegativity
  field_mean/max     cation field strength z / r^2 (Shannon radius, A) -- cation-count weighted
                     mean and maximum
  chi_mean           mean cation electronegativity
  mass_cation        mean cation atomic mass (heavier lattices -> lower phonon energies)
  polarizability     mean anion polarizability (approximate free-anion values, A^3)
  f_<group>          cation fractions of network formers / intermediates / modifiers
  former             label of the dominant network former (for leave-one-family-out splits)
All descriptors are intensive and composition-weighted, so they are defined for any glass.
"""
from __future__ import annotations

import functools

import numpy as np

ANIONS = {"O": -2, "F": -1, "S": -2, "Se": -2, "Cl": -1, "Br": -1, "I": -1, "N": -3}
# approximate anion polarizabilities (A^3): ionic-refraction scale values, for trends only
POLARIZABILITY = {"O": 1.8, "F": 0.8, "S": 5.5, "Se": 7.0, "Cl": 3.0, "Br": 4.2, "I": 6.2, "N": 2.5}
GROUPS = {
    "Si": ["Si"], "B": ["B"], "P": ["P"], "Ge": ["Ge"], "Te": ["Te"], "Al": ["Al"], "Ga": ["Ga"],
    "heavy": ["Pb", "Bi", "Tl", "Sb", "W", "Nb", "Ta", "Mo"], "Zr": ["Zr", "Hf"], "In": ["In"],
    "alkali": ["Li", "Na", "K", "Rb", "Cs"], "alkearth": ["Mg", "Ca", "Sr", "Ba"], "Zn": ["Zn", "Cd"],
    "La": ["La", "Y", "Gd", "Lu", "Sc"], "Ti": ["Ti", "V"],
}
RE_DOPANTS = {"Pr", "Nd", "Pm", "Sm", "Eu", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Ce"}
FORMERS = ["Si", "B", "P", "Ge", "Te", "heavy", "Ga", "Al"]
# Approximate optical basicities of oxides per cation (Duffy; Dimitrov & Komatsu compilations),
# used where the electronegativity formula is known to fail (lone-pair and heavy cations).
LAMBDA_OXIDE = {"Si": 0.48, "B": 0.42, "P": 0.40, "Ge": 0.61, "Te": 0.93, "Al": 0.60, "Ga": 0.72,
                "Li": 1.00, "Na": 1.15, "K": 1.40, "Rb": 1.52, "Cs": 1.70, "Mg": 0.78, "Ca": 1.00,
                "Sr": 1.07, "Ba": 1.15, "Zn": 1.03, "Cd": 1.12, "Pb": 1.18, "Bi": 1.19, "Sb": 1.18,
                "La": 1.18, "Y": 1.00, "Gd": 1.18, "Lu": 1.05, "Zr": 0.87, "Ti": 0.95, "Nb": 1.05,
                "W": 1.04, "Ta": 1.00, "Mo": 0.95, "Tl": 1.20, "Sn": 0.85, "In": 0.85, "V": 0.85}
_SUB = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")
NAMES = (["x_" + a for a in ["O", "F", "S", "Se", "Cl", "Br", "I"]] +
         ["basicity", "field_mean", "field_max", "chi_mean", "mass_cation", "polarizability"] +
         ["f_" + g for g in GROUPS])


@functools.lru_cache(maxsize=None)
def _el(sym):
    from pymatgen.core import Element
    return Element(sym)


@functools.lru_cache(maxsize=None)
def _oxi(formula):
    """Oxidation states of a compound: charge balance with fixed anion charges, else pymatgen."""
    from pymatgen.core import Composition
    c = Composition(formula).get_el_amt_dict()
    an = {k: v for k, v in c.items() if k in ANIONS}
    ca = {k: v for k, v in c.items() if k not in ANIONS}
    if len(ca) == 1 and an:
        el, n = next(iter(ca.items()))
        q = -sum(ANIONS[a] * m for a, m in an.items()) / n
        return {el: q, **{a: ANIONS[a] for a in an}}, c
    try:
        g = Composition(formula).oxi_state_guesses(max_sites=-1)
        if g:
            return dict(g[0]), c
    except Exception:
        pass
    return {k: (ANIONS[k] if k in ANIONS else 3.0) for k in c}, c


def _radius(sym, q):
    e = _el(sym)
    r = e.ionic_radii
    if r:
        qq = min(r, key=lambda k: abs(k - q))
        return float(r[qq])
    return float(e.average_ionic_radius or 1.0)


def describe(components: dict, dopant: str | None = None):
    """Descriptor vector (np.ndarray, order NAMES) and dominant-former label."""
    cat, an, charge = {}, {}, {}
    for f, x in components.items():
        f = f.translate(_SUB).replace(" ", "").replace("·", "")
        try:
            ox, c = _oxi(f)
        except Exception:
            continue
        els = set(c)
        if dopant and dopant in els and len(els - set(ANIONS)) == 1:
            continue                                   # the dopant compound itself
        for el, n in c.items():
            if el in ANIONS:
                an[el] = an.get(el, 0.0) + n * x
            elif el in RE_DOPANTS and el == dopant:
                continue
            else:
                cat[el] = cat.get(el, 0.0) + n * x
                charge[el] = abs(ox.get(el, 3.0)) or 3.0
    sa, sc = sum(an.values()), sum(cat.values())
    if sc <= 0 or sa <= 0:
        return None, None
    xa = {k: v / sa for k, v in an.items()}
    xc = {k: v / sc for k, v in cat.items()}
    # optical basicity: weights = cation count x charge (oxygen-neutralising share)
    wz = {k: xc[k] * charge[k] for k in xc}; Z = sum(wz.values())
    oxide = xa.get("O", 0.0) > 0.5
    lam = sum(wz[k] / Z * (LAMBDA_OXIDE[k] if (oxide and k in LAMBDA_OXIDE)
                           else 1.0 / (1.36 * (max(_el(k).X or 1.5, 0.5) - 0.26))) for k in xc)
    fs = {k: charge[k] / _radius(k, charge[k]) ** 2 for k in xc}
    field_mean = sum(xc[k] * fs[k] for k in xc); field_max = max(fs[k] for k in xc if xc[k] > 0.02) \
        if any(xc[k] > 0.02 for k in xc) else max(fs.values())
    chi = sum(xc[k] * (_el(k).X or 1.5) for k in xc)
    mass = sum(xc[k] * float(_el(k).atomic_mass) for k in xc)
    pol = sum(xa[k] * POLARIZABILITY.get(k, 2.0) for k in xa)
    fg = {g: sum(xc.get(e, 0.0) for e in els) for g, els in GROUPS.items()}
    v = ([xa.get(a, 0.0) for a in ["O", "F", "S", "Se", "Cl", "Br", "I"]] +
         [lam, field_mean, field_max, chi, mass, pol] + [fg[g] for g in GROUPS])
    # dominant former label
    if xa.get("F", 0) + xa.get("Cl", 0) > 0.5:
        label = "fluoride/halide" if xa.get("O", 0) < 0.2 else "oxyfluoride"
    elif xa.get("S", 0) + xa.get("Se", 0) > 0.5:
        label = "chalcogenide"
    else:
        fr = {g: fg[g] for g in FORMERS}
        top = max(fr, key=fr.get)
        label = {"Si": "silicate", "B": "borate", "P": "phosphate", "Ge": "germanate", "Te": "tellurite",
                 "heavy": "heavy-metal oxide", "Ga": "gallate", "Al": "aluminate"}[top] if fr[top] > 0.1 else "other oxide"
        if xa.get("F", 0) > 0.2:
            label = "oxyfluoride"
    return np.array(v, float), label
