# MD_SITES — Er³⁺ site geometry in fiber glasses from MLIP melt-quench MD

Branch `forward/md-sites`. Driver: `md_sites/run_md_sites.py`. Outputs:
`md_sites/<glass>.json` (+ `<glass>_final.extxyz`, `<glass>.log`).

## What was run (scheduled run, 2026-11-09, 55-min hard budget)

Hardware for this run: **2 CPU cores, 7 GB RAM, no GPU**. ~17 of the 55 min
went to installing torch (the CPU-only PyTorch index was blocked by the proxy,
so the full CUDA wheel set from PyPI had to be installed). The two glasses ran
in parallel, one core each.

| | Er-silica | Er-ZBLAN |
|---|---|---|
| composition | 64 SiO₂ + 2 Er₂O₃ (202 atoms, 4 Er) | 27 ZrF₄ + 10 BaF₂ + 2 ErF₃ + 2 AlF₃ + 10 NaF (201 atoms, 2 Er; Er replaces La of 53-20-4-3-20 ZBLAN) |
| density (fixed, NVT) | 2.30 g/cm³ | 4.33 g/cm³ |
| model | MACE-OMAT-0 medium (float32, contains Er/Yb/Tm/Pr) | same |
| protocol | 4000 K melt 0.104 ps → linear quench to 300 K over 0.288 ps (~1.3×10¹⁶ K/s) → 0.132 ps production, dt 1 fs; wall 1656 s | 1800 K melt 0.077 ps → quench over 0.213 ps (~7×10¹⁵ K/s) → 0.099 ps production, dt 1 fs; wall 1542 s |
| step time (1 core) | 2.74 s/step | 3.61 s/step |

Protocol: random packing at the experimental density (min. distance 1.7 Å),
40 FIRE steps, NVT Langevin (friction 0.02 fs⁻¹) melt, linear quench to 300 K
in 40 stages, then weakly-damped (0.002 fs⁻¹) 300 K production; 20 frames
saved evenly through production. **Volume is held fixed** at the experimental
density throughout (no NPT — universal MLIPs give poor glass densities under
constant pressure). Step counts were set at run time from the measured step
time so the run fit its wall-clock budget. Model: the repo default
`DEFAULT_MODEL` from `gainmat/host/phonons.py` (MACE-OMAT-0 medium; the
driver checks Z = 59, 68, 69, 70 are in its element table and would fall back
to MACE-MP-0 medium otherwise — no fallback was needed).

First shell: cutoff at the first minimum of the (smoothed, r²-normalised)
Er–ligand pair distribution, searched between peak+0.2 Å and 3.6 Å.

## Results

| | Er-silica (4 Er × 20 frames) | Er-ZBLAN (2 Er × 20 frames) | EXAFS (lit.) |
|---|---|---|---|
| ligand | O | F | |
| Er–L RDF first peak | **2.13 Å** | **2.28 Å** | Er–O 2.2–2.3 Å (silica); Er–F 2.2–2.3 Å (ZBLAN) |
| first-min cutoff r_c | 3.38 Å | 3.28 Å | |
| CN (r < r_c) mean | 6.56 | 8.75 | ~6 (silica); ~8 (ZBLAN) |
| CN histogram | 6:47, 7:22, 8:10, 9:1 | 8:10, 9:30 | |
| R_mean ± R_std (all ligands < r_c) | 2.49 ± 0.40 Å | 2.41 ± 0.31 Å | |
| **actual T of production** | **~1700 K** (final frame) | **~1230 K** (final frame) | 300 K |
| Er–Er nearest | 4.9 Å | 5.9 Å | |
| VDOS 99 % cutoff / top peak (by-product, unreliable) | 1516 / 1011 cm⁻¹ | 1361 / 1021 cm⁻¹ | |

Reading: the **RDF first-peak positions** (2.13 Å Er–O, 2.28 Å Er–F) and the
**CNs** (≈6–7 in silica, ≈8–9 in ZBLAN) land in or near the EXAFS ranges —
encouraging for the model's Er chemistry. The mean distance R_mean and
R_std are inflated because the first-minimum cutoff (~3.3 Å) is wide in hot,
poorly-ordered structures and catches a long tail of second-shell-ish
ligands; for SPM input, the ligands at r < ~2.8 Å (the core shell) are the
physically meaningful set and can be re-cut from the JSON coordinates.
Er-silica peak at 2.13 Å is slightly short of 2.2–2.3 Å (possibly Er–NBO
bonding in a very under-relaxed network, and MLIP error).

**Big problem found in this run:** the target thermostat temperature was
300 K, but the systems were *not* at 300 K in production — the final frame
temperatures are ~1700 K (silica) and ~1230 K (ZBLAN). Relaxation of the
random packing released potential energy faster than the Langevin
thermostat (0.02 fs⁻¹ during quench, 0.002 fs⁻¹ in production) could remove
it within a few hundred fs. So the "glass" is really a hot, partially formed
network and the site geometries are thermally broadened accordingly. Fix:
much longer melt/quench (≥10 ps), stronger friction during quench, and a
check of ⟨T⟩ before production.

Literature EXAFS for comparison:

* Er in silica: Er–O ≈ 2.2–2.3 Å, CN ≈ 6 (Marcus & Polman, "Local structure
  around Er in silica and sodium silicate glasses", J. Non-Cryst. Solids 136,
  260 (1991); MeV-implanted silica gives similar values, Polman et al., MRS
  Proc.). 
* Er in fluorozirconate (ZBLAN-type) glass: Er–F ≈ 2.2–2.3 Å, CN ≈ 8 — this is
  the target range given in the task brief; I did not verify a specific ZBLAN
  EXAFS paper during this run, so treat the citation as pending.

## Honest caveats (read before using these sites)

1. **The quenches are absurdly fast.** Sub-picosecond to ~1 ps quenches
   (10¹⁵–10¹⁶ K/s) vs ~10³ K/s in a fiber draw — 12+ orders of magnitude.
   Even "normal" GPU MLIP glass work (10–100 ps) is 10¹²–10¹³ K/s. These
   structures are quenched liquids with essentially the melt's medium-range
   order. The *first coordination shell* of a strongly-bound ion forms in
   ~0.1 ps and is the least quench-sensitive quantity, which is why it is the
   only thing reported, but CN and the distance spread should be considered
   biased towards liquid-like disorder (larger R_std, possibly CN off by ~1).
2. **Production is ~0.1 ps**, so the 20 "snapshots" per Er are strongly
   correlated thermal fluctuations of the *same* site, not independent sites.
   The number of genuinely independent sites is the number of Er atoms
   (4 in silica, 2 in ZBLAN). That is far too few for an inhomogeneous
   site ensemble; it is a pipeline proof + first estimate only.
3. Er concentration is very high (Er₂O₃ ~3 mol% in silica, needed to get >1
   Er in 200 atoms); Er–Er clustering at this level is plausible and would not
   be representative of ~100-ppm fibers. No Al/P co-dopant was included.
4. Melt temperatures are short and from a random packing; equilibration of the
   melt is not demonstrated.
5. VDOS cutoff from <0.1 ps of 300 K velocities has ~150 cm⁻¹ resolution —
   reported in the JSON as a by-product only, do not use it.

Tellurite and phosphate were **not run** (no time left on 2 cores).

## To do this properly (one A100 or a 32-core node)

`python md_sites/run_md_sites.py Er-silica --budget-min 240` — the driver is
budget-driven, so the same command scales: on an A100 MACE-medium runs ~600
atoms at ~20–50 ms/step, i.e. ~100 ps quench + 50 ps production in a few
hours. Recommended: 4–8 independent melts per glass (different seeds),
≤1 Er per 100 cations, Al co-doping variant for silica, and production frames
spaced ≥1 ps apart. `GLASSES["Er-tellurite"]` is already defined.
