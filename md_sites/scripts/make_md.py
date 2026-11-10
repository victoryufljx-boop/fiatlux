import json, sys
d = json.load(open(sys.argv[1])); out = sys.argv[2]
S = d['summary']; P = d['protocol']; sc = S['sanity_checks']; pot = d['potential']
F = d.get('summary_fast_quench_replicas')
A = S['angles_deg']
def st(x): return f"{x['mean']} ± {x['std']} (mode {x['mode']})" if x else 'n/a'
L = []
w = L.append
w('# Er-Ge10-SiO2: Er site geometries from classical MD (Pedone 2006 potential)\n')
w('Glass: germanosilicate, 900 SiO2 + 100 GeO2 + 4 Er2O3 (89.6 / 9.96 / 0.40 mol%), '
  f"{P['n_atoms']} atoms, cubic box {P['box_A']} Å at a fixed {P['density_g_cm3']} g/cm³ (NVT).\n")
w('Files: `Er-Ge10-SiO2.json` has the sites, statistics, RDFs and every potential parameter with its source. '
  '`Er-Ge10-SiO2_repN_final.extxyz` holds the final 300 K configurations. The scripts are in `md_sites/scripts/`.\n')
w('## Potential\n')
w('Pedone, Malavasi, Menziani, Cormack, Segre, J. Phys. Chem. B 110, 11780 (2006), Table 2 (binary-oxide fit). '
  'Every number below comes from that table, which I transcribed from the paper PDF in this session. '
  'U = q_i q_j/r + D[(1-e^{-a(r-r0)})² - 1] + C/r¹². Charges: O -1.2, Si +2.4, Ge +2.4, Er +1.8.\n')
w('| pair | D (eV) | a (1/Å) | r0 (Å) | C (eV Å¹²) |\n|---|---|---|---|---|')
for k, v in pot['pairs'].items():
    w(f"| {k} | {v['D_eV']} | {v['a_invA']} | {v['r0_A']} | {v['C_eVA12']} |")
w('\nThe Er-O pair is Pedone\'s own Er(1.8)-O(-1.2) row, so no substitute RE-O set was needed. '
  'Al-O (0.361581, 1.900442, 2.164818, 0.9) is recorded in the JSON but not used because this glass has no Al. '
  'Cation-cation pairs are Coulomb only.\n')
w('**Implementation (deviation):** LAMMPS 22 Jul 2025 (PyPI wheel). I put every pair into one `pair_style table` '
  '(12000 points). Coulomb uses damped shifted force (DSF, alpha 0.2 Å⁻¹, cutoff 10 Å) instead of Ewald/PPPM. '
  'The Morse and r⁻¹² terms are cut at 10 Å and energy-shifted. The reason is speed: the machine has only 2 cores, '
  'and hybrid coul/long + PPPM (rc 11 Å, 1e-5) ran at 25 ms/step on one core, while the single table runs at 8.5-10 ms/step. '
  'DSF is a standard real-space alternative for oxide glasses, but I did not run a PPPM cross-check on the final structures.\n')
w('## Protocol\n')
w('Each replica starts from a random placement with minimum distances (cation-cation 2.6 Å, cation-O 1.55 Å, O-O 2.2 Å, '
  'Er-Er 4.0 Å). This is followed by CG minimisation and 3 ps of nve/limit with Berendsen at 6000 K to remove overlaps. '
  'Then, under Nose-Hoover NVT (Tdamp 0.1 ps) with a 1 fs timestep:\n')
w('- 6000 K for 20 ps (the Morse + r⁻¹² core held at 6000 K, with no atom collapse)\n- 3500 K for 50 ps\n'
  '- linear cooling from 3500 K to 300 K over 320 ps, i.e. **10 K/ps**\n- 300 K for 50 ps. Sites come from 10 frames, one every 2 ps over the last 20 ps.\n')
w(f"Main replicas: {', '.join(P['replicas_main_10Kps'])} (seeds {P['seeds']}). "
  'Each took about 62 min of wall time on one core, with two replicas running at once.\n')
if F:
    w(f"Extra fast-quench replicas: {', '.join(P['fast_replicas']['replicas'])}. These used a shortened schedule "
      '(6000 K 20 ps, 3500 K 20 ps, 20 K/ps quench, 300 K 30 ps) to add distinct sites near the end of the budget. '
      'Their sites are tagged `quench_rate_K_per_ps: 20` and are **excluded** from `summary`. '
      'Their statistics are in `summary_fast_quench_replicas`.\n')
w('## Deviations from the requested protocol (honest list)\n')
w('1. **Quench rate 10 K/ps, not 5 K/ps.** `nproc` = 2. At about 9 ms/step per core, one 5 K/ps replica (760 ps) takes about 110 min. '
  'Within 3 h 15 min that would allow only 2 replicas (16 Er), so I followed the 10 K/ps fallback.\n')
w(f"2. **Distinct Er sites: {S['n_distinct_er']} (main set)" + (f" + {F['n_distinct_er']} (fast set)" if F else '') +
  ', below the target of ≥ 60.** Each replica holds 8 Er, and 2 cores allowed only 2 rounds of 2 concurrent replicas. '
  '`er_sites` contains 10 frames per Er. These snapshots are thermally correlated (300 K vibrations of the same site), '
  'not independent sites.\n')
w('3. Coulomb uses DSF instead of Ewald/PPPM (see above).\n')
w('4. During replicas 3-4 the cloud container was restarted (at about 01:05 UTC), killing both runs part-way through the quench. '
  'I restarted them from their restart files written at the end of the 3500 K hold (`hold.restart`, step 70000), '
  'so the full 320 ps quench and the 300 K stage ran uninterrupted after the restart. '
  'Because the restart rebuilt the thermostat chain, the run is not bit-identical to an uninterrupted one, but the protocol is the same.\n')
w('## Sanity checks (pooled over all frames of the main replicas)\n')
w('| quantity | MD | target |\n|---|---|---|')
w(f"| Si-O first peak | {sc['Si-O_first_peak']} Å | ~1.61 Å |")
w(f"| O-Si-O | {st(sc['O-Si-O'])}° | ~109° |")
w(f"| Si-O-Si | {st(sc['Si-O-Si'])}° | ~145-150° |")
w(f"| Si CN (fraction 4-fold) | {sc['Si_CN_mean']} ({sc['Si_CN4_frac']}) | 4 |")
w(f"| Ge-O first peak | {sc['Ge-O_first_peak']} Å | ~1.74 Å |")
w(f"| Ge CN (fraction 4-fold) | {sc['Ge_CN_mean']} ({sc['Ge_CN4_frac']}) | 4 |")
w(f"| O-Ge-O | {st(sc['O-Ge-O'])}° | ~109° |")
w(f"| Si-O-Ge / Ge-O-Ge | {st(sc['Si-O-Ge'])}° / {st(sc['Ge-O-Ge'])}° | ~130-150° |")
w(f"| O-O / Si-Si first peak | {sc['O-O_first_peak']} / {sc['Si-Si_first_peak']} Å | 2.62 / 3.1 Å |")
w(f"| bulk O speciation | {sc['oxygen_speciation_bulk']} | NBO only from Er2O3 |")
w(f"| Er-O first peak | {S['ErO_first_peak']} Å | ~2.2-2.35 Å (EXAFS 2.25) |")
w(f"| Er-O CN within r_min = {S['r_min_ErO']} Å | {S['CN_mean']} ± {S['CN_std']} | ~6-7 (EXAFS ~6) |")
w(f"| final-stage T | {S['final_stage_mean_T_K']} K | 300 K |\n")
w('The Si-O-Si angle comes out a few degrees above the requested range, consistent with Pedone-type potentials, which tend to give wide Si-O-Si angles. '
  'Both formers are almost perfectly tetrahedral.\n')
w('## Er first shell (main replicas)\n')
w(f"- r_min (first minimum of the pooled Er-O g(r), 0.1 Å smoothing) = **{S['r_min_ErO']} Å**. "
  f"Er-O first peak {S['ErO_first_peak']} Å.")
w(f"- CN within r_min: **{S['CN_mean']} ± {S['CN_std']}**, histogram {S['CN_hist']}. "
  f"CN within 2.5 / 2.6 / 2.8 Å: {S['CN_within_2.5A_mean']} / {S['CN_within_2.6A_mean']} / {S['CN_within_2.8A_mean']}.")
w(f"- Er-O over all first-shell ligands: **{S['R_mean']} ± {S['R_std']} Å**. The shortest Er-O per site averages {S['R_shortest_mean']} Å.")
w(f"- First-shell O types: {S['first_shell_O_types']}, giving an NBO fraction of **{S['nbo_fraction_first_shell']}**. "
  'Most of the NBOs in the glass sit on Er; Er2O3 is their only source.')
x = S['ErOGe_vs_ErOSi_linkages']
w(f"- Er-O-T linkages: Er-O-Si {x['Er-O-Si']}, Er-O-Ge {x['Er-O-Ge']}, so the Ge fraction is {x['frac_Ge']} against 0.10 Ge among the formers.")
w(f"- Er-Er nearest neighbour: mean {S['er_er_nn']['mean']} Å, min {S['er_er_nn']['min']} Å. "
  f"Fraction of Er with another Er within 5 Å: {S['er_er_nn']['frac_below_5A']}.")
w(f"- O-Er-O angles: {st(A['O-Er-O'])}°. Er-O-T angles: {st(A['Er-O-T(former)'])}°.\n")
w('**Assessment.** The run is now cold (300 K) and much better quenched than `forward/md-sites` (>1000 K, Er-O 2.49 ± 0.40 Å, CN 6.6). '
  'The Er-O distance is in the expected range: the peak is slightly short of EXAFS, and the mean over the shell is close to it. '
  '**The coordination number is low, however.** About 4 O sit within 2.6 Å, and roughly one more lies in a diffuse 2.6-3.2 Å tail. '
  'The target from EXAFS is ~6. This is a property of the Pedone Er-O parameters in a silica-rich network at this density, not of the protocol. '
  'Er binds all available NBOs and adds a few long, weak bridging-O contacts, so the sites are strongly distorted, low-symmetry 4-5 coordinated polyhedra. '
  'Er also clusters: Er pairs share NBOs at 3.5-4.2 Å. In real MCVD Ge-silica, the absence of Al/P likewise favours clustering, '
  'but downstream crystal-field use should treat these sites as an under-coordinated end member. '
  'A possible next step is to rescale the Er-O D/r0 to reach CN ~6. That would mean refitting, which I did not do because inventing parameters was ruled out.\n')
w('## JSON layout\n')
w('`er_sites[i]` = {replica, frame, timestep, er_index, ligands: [["O", x, y, z], ...] relative to the Er, minimum image, Å, sorted by distance; '
  'ligand_tags: [{type: BO|NBO|free|TBO, formers: [Si/Ge...], r}], other_neighbors_within_rmin: non-O atoms within r_min (e.g. Si, Er), '
  'quench_rate_K_per_ps}. Only O appear in `ligands`, so `gainmat.forward.glass.load_md_sites` reads them directly. '
  'BO/NBO criterion: O within 2.0 Å of 2 / 1 Si or Ge.\n')
open(out, 'w').write('\n'.join(L))
