"""Host (dopant-independent) tier: band gap, refractive index, phonon cutoff.

Submodules import their heavy dependencies (gpaw, torch/mace, phonopy) lazily.
"""
from .run_host import run_host, compare_to_experiment, CRYSTAL_KEYS, GLASS_KEYS  # noqa: F401
