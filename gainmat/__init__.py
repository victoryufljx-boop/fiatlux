"""gainmat: first-principles optical properties of laser / fiber-amplifier gain media.

Tiers (independent; each returns pieces of gainmat.schema.Result):
  gainmat.host     - band gap, refractive index, phonon cutoff (GPAW + MACE)
  gainmat.ln4f     - lanthanide 4f^N levels, Judd-Ofelt rates, multiphonon quenching
  gainmat.cluster  - embedded-cluster CASSCF/NEVPT2 for 5d (Ce3+) and 3d (Ti3+) centers
"""
from .schema import *  # noqa: F401,F403
