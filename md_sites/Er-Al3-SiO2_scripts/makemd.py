import json
J=json.load(open('/home/claude/fiatlux/md_sites/Er-Al3-SiO2.json'));S=J['summary'];s=S['sanity'];P=J['protocol']
t=P['wall_time_s_per_replica']
nrep=S['n_replicas']
cv=S['CN_vs_cutoff']
md=f"""# Er-Al3-SiO2 — classical-MD rare-earth sites (Pedone potential, slow quench to 300 K)

Branch `md2/Er-Al3-SiO2`. Data: `Er-Al3-SiO2.json` (keys glass, composition, potential, protocol, summary, er_sites),
final 300 K configurations `Er-Al3-SiO2_rep*_300K.extxyz`, inputs/scripts in `Er-Al3-SiO2_scripts/`.

## What was run
* Glass: 966 SiO2 + 34 Al2O3 + 4 Er2O3 = 3088 atoms (Si 966, Al 68, Er 8, O 2046); Al2O3 = 3.386 mol%, Er2O3 = 0.398 mol%.
  Cubic box 36.127 A, fixed density 2.22 g/cm3 (NVT).
* Potential: Pedone et al. J. Phys. Chem. B 110, 11780 (2006), rigid ion, q(O)=-1.2, q(Si)=2.4, q(Al)=q(Er)=1.8;
  Morse + C/r^12 for cation-O and O-O (LAMMPS `pair_style pedone`, cutoff 5.5 A); cation-cation Coulomb only.

  | pair | D (eV) | a (1/A) | r0 (A) | C (eV A^12) |
  |---|---|---|---|---|
  | Si-O | 0.340554 | 2.006700 | 2.100000 | 1.0 |
  | Al-O | 0.361581 | 1.900442 | 2.164818 | 0.9 |
  | Er-O | 0.040448 | 2.294078 | 2.837722 | 3.0 |
  | O-O  | 0.042395 | 1.379316 | 3.618701 | 22.0 |

  **Source of every number:** the PMMCS (Pedone 2006) parameter dictionary in `amorphouspy` v0.8.1 (PyPI wheel,
  `amorphouspy/lammps/potentials/pmmcs_potential.py`, BAM / A. Atila), which presents them as the Pedone 2006 table.
  I could not fetch the original ACS paper or its SI in this run (blocked), and found no second independent reproduction
  of the **Er-O** row, so Er-O is single-sourced (not invented, but not cross-checked against the paper). Si-O, Al-O and
  O-O match the commonly reproduced Pedone values. Ge-O is not used (no Ge in this glass).
* Electrostatics: **DSF** (`coul/dsf`, alpha 0.25 1/A, cutoff 8 A — the amorphouspy PMMCS default), not PPPM with a 10-11 A
  real-space cutoff. Reason: only 2 CPU cores; PPPM(1e-5)/10 A cost 19 ms/step/core vs 6.5 ms/step for DSF/8 A
  (PPPM 1e-4 was 10 ms). erfc(0.25*8)=0.005, so the damped sum is well converged, but this is a deviation.
* Protocol per replica (1 fs, Nose-Hoover NVT, Tdamp 0.1 ps): random placement with minimum distances -> CG minimisation ->
  3 ps Langevin 3000->6000 K with `nve/limit` (contact relaxation) -> **6000 K 20 ps** -> 2 ps ramp -> **3500 K 50 ps** ->
  **linear quench 3500->300 K at 10 K/ps (320 ps)** -> **300 K 50 ps**; 10 frames every 2 ps over the last 20 ps are analysed.
  6000 K was stable (no lost atoms, no dangerous neighbour builds).
* **Quench rate 10 K/ps, not 5 K/ps** — the machine had `nproc` = 2. Measured ~7 ms/step per single-core replica,
  two replicas concurrently (faster in aggregate than one 2-core run: 14.4 ms/step with 2 OpenMP threads + PPPM).
  A 5 K/ps protocol (765k steps, ~90 min/replica) would have allowed only 4 replicas; 10 K/ps gave {nrep}.
* Replicas: **{nrep}** fully independent (different random placement seed and velocity seed each) -> **{S['n_distinct_er']} distinct Er ions**,
  {S['n_site_frames']} Er site-frames in `er_sites`. The target of >=60 distinct Er sites was **not** reached (8 Er per box,
  composition fixed by the task). Wall time per replica: {', '.join(f'rep{k} {v/60:.0f} min' for k,v in t.items())}.
  The cloud container was restarted twice during the session; each time the two replicas then running (3+4, then 5+6), which were only minutes in, were killed and re-run from scratch
  with the same seeds (no partial trajectories were stitched).

## Sanity checks (all frames, all replicas)
| quantity | result | target |
|---|---|---|
| Si-O RDF peak / mean first shell | {s['SiO_rdf_peak']} A / {s['SiO_mean_first_shell']['mean']} +- {s['SiO_mean_first_shell']['std']} A | ~1.61 A |
| Si CN (r<{s['rmin_SiO']} A) | {s['Si_CN_mean']} ({s['Si_CN_hist']}) | 4 |
| O-Si-O | {s['O_Si_O']['mean']} +- {s['O_Si_O']['std']} deg | ~109 deg |
| Si-O-Si | mean {s['Si_O_Si']['mean']} +- {s['Si_O_Si']['std']}, mode ~{s['Si_O_Si']['peak']} deg | ~145-150 deg |
| Al-O RDF peak / mean first shell | {s['AlO_rdf_peak']} A / {s['AlO_mean_first_shell']['mean']} +- {s['AlO_mean_first_shell']['std']} A | ~1.75 A |
| Al CN (r<{s['rmin_AlO']} A) | {s['Al_CN_mean']} ({s['Al_CN_hist']}) | mostly 4 |
| O-Al-O | {s['O_Al_O']['mean']} +- {s['O_Al_O']['std']} deg | ~109 deg |
| O speciation (whole glass) | {s['oxygen_speciation_all_frames']} | |
| Er-O RDF peak | {S['ErO_rdf_peak']} A | ~2.25 A (EXAFS) |
| Er-O first minimum r_min | {S['r_min']} A | ~2.9-3.2 A |
| Er CN (r<r_min) | {S['CN_mean']} +- {S['CN_std']}, hist {S['CN_hist']} | ~6-7 |
| Er-O distance (r<r_min) | {S['R_mean']} +- {S['R_std']} A | ~2.2-2.35 A |

Er CN vs cutoff: {', '.join(f"r<{k} A: {v['mean']}" for k,v in cv.items())}.

Mean pressure in the 300 K stage is about -2 GPa (tensile) at the imposed 2.22 g/cm3, i.e. this potential with DSF/8 A
would prefer a denser glass; the task specified fixed experimental density, so the box was not relaxed.

## Er site statistics
* Ligand O types: {S['ligand_O_types']} -> NBO fraction in the first shell **{S['nbo_fraction']}** (NBO = O bonded to one former;
  Si within 2.0 A, Al within 2.1 A). Almost all NBO in the glass sit on Er (whole-glass NBO count ~= Er first-shell NBO count).
* Er-O-former linkages: {S['er_o_former_linkages']['Er-O-Si']} Er-O-Si, {S['er_o_former_linkages']['Er-O-Al']} Er-O-Al -> Al fraction
  **{S['er_o_former_linkages']['frac_Al']}** vs {S['er_o_former_linkages']['random_expectation_frac_Al_by_cation_count']} for random
  former choice: Er is enriched near Al ({S['er_o_former_linkages']['frac_Al']/S['er_o_former_linkages']['random_expectation_frac_Al_by_cation_count']:.1f}x), but most of its ligands still bond to Si.
* Er-Er nearest neighbour: {S['er_er_nn']['mean']} +- {S['er_er_nn']['std']} A, min {S['er_er_nn']['min']} A; fraction of Er with another
  Er within 5 A: {S['er_er_nn']['frac_below_5A']} (some clustering/O-sharing; 0.4 mol% Er2O3 is far above fibre doping levels).
* O-Er-O angle: {S['O_Er_O_angle']['mean']} +- {S['O_Er_O_angle']['std']} deg (histogram in JSON); the first shell is strongly distorted.
* Er-O RDF (r, g) sampled to 8 A is in `summary.ErO_rdf`.

## Assessment / what did not work
* Compared with forward/md-sites (ML potential, >1000 K, Er-O 2.49 +- 0.40 A, CN 6.6) these sites are properly quenched to 300 K
  and have a much tighter, shorter first shell (peak {S['ErO_rdf_peak']} A, mean {S['R_mean']} +- {S['R_std']} A).
* **Er is under-coordinated relative to EXAFS**: CN {S['CN_mean']} within r_min={S['r_min']} A and only ~{cv['2.6']['mean']} within 2.6 A,
  versus CN ~6 at 2.25 A measured. This is a property of the Pedone Er-O pair (very shallow Morse well, D=0.04 eV), possibly
  aggravated by DSF/8 A and the -2 GPa tensile state; it is not a quench artefact. Downstream crystal-field use should keep this in mind
  (sites with CN 4-5 will give larger, more axial crystal fields than real 6-coordinated Er).
* Al: {s['Al_CN_hist']} — a sizeable AlO3 fraction plus O triclusters, expected because only 8 Er (charge compensators) are available for 68 Al
  (peraluminous, no alkali/alkaline-earth), but it may be over-represented by the potential.
* Quench at 10 K/ps rather than 5 K/ps; DSF rather than PPPM; {S['n_distinct_er']} instead of >=60 distinct Er ions; Er-O parameters single-sourced.
"""
open('/home/claude/fiatlux/md_sites/Er-Al3-SiO2.md','w').write(md)
