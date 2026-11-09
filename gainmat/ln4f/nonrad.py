"""Multiphonon (nonradiative) relaxation via the energy-gap law.

  W_nr = B exp(-alpha dE)                       (family form; Layne-Lowdermilk-Weber / Reisfeld)
  W_nr = beta_el exp[-alpha (dE - 2 hw_max)]    (van Dijk-Schuurmans modified law)
  W_nr(host-agnostic, hw) : log B and alpha*hw interpolated over the family table as a
       function of hw_max, i.e. W = B(hw) exp(-[alpha hw](hw) * p),  p = dE/hw_max
       (phonon-order form; consumes the host tier's computed phonon cutoff).
"""
from __future__ import annotations
import numpy as np
from ..refs.ln_refs import GAP_LAW, VAN_DIJK

FAMILIES = {k: v for k, v in GAP_LAW.items() if isinstance(v, dict)}


def w_family(dE, family):
    p = FAMILIES[family]
    return p["B"] * np.exp(-p["alpha"] * dE)


def w_van_dijk(dE, hw_max, beta_el=VAN_DIJK["beta_el"], alpha=VAN_DIJK["alpha"]):
    return min(beta_el * np.exp(-alpha * (dE - 2 * hw_max)), 1e13)


def _fit_phonon_order():
    """Least-squares lines ln B(hw) and alpha*hw (hw) over the oxide+fluoride families.
    Chalcogenide is excluded from the fit (its B is an order-of-magnitude value)."""
    tab = [v for k, v in FAMILIES.items() if k != "chalcogenide" and v.get("fit", True)]
    hw = np.array([t["hw"] for t in tab], float)
    cb = np.polyfit(hw, np.log([t["B"] for t in tab]), 1)
    ca = np.polyfit(hw, [t["alpha"] * t["hw"] for t in tab], 1)
    return cb, ca


_CB, _CA = _fit_phonon_order()


def w_phonon_order(dE, hw_max):
    """Host-agnostic energy-gap law parameterised by the max phonon energy:
    W = B(hw) exp(-a(hw) p), p = dE/hw, with ln B and a = alpha*hw linear in hw
    (fitted to the family table). Smooth and monotone in hw for p >~ 2."""
    lb = np.polyval(_CB, hw_max)
    a = max(np.polyval(_CA, hw_max), 0.5)
    return float(np.exp(lb - a * dE / hw_max))


def w_nr(dE, hw_max=None, family=None, model="auto"):
    """Nonradiative rate (s^-1) for a gap dE (cm^-1)."""
    if dE <= 0:
        return 1e13
    if model == "family" or (model == "auto" and family in FAMILIES):
        return float(w_family(dE, family))
    if model == "van_dijk":
        return float(w_van_dijk(dE, hw_max))
    if hw_max is None:
        raise ValueError("need hw_max or a known host family")
    return w_phonon_order(dE, hw_max)


def quantum_eff(A, W):
    return A / (A + W) if A + W > 0 else 0.0
