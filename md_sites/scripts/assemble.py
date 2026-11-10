"""Assemble md_sites/Er-Ge10-SiO2.json from analysis output + metadata.
usage: python assemble.py analysis.json OUT.json rep1 rep2 ...
"""
import sys, json, os, re, numpy as np
sys.path.insert(0, '/home/claude/md')
from build import PED, CHG, N, DENS, box_length

ana = json.load(open(sys.argv[1])); out = sys.argv[2]; reps = sys.argv[3:]
SRC = ("A. Pedone, G. Malavasi, M. C. Menziani, A. N. Cormack, A. Segre, J. Phys. Chem. B 110, 11780 (2006), "
       "Table 2 'Potential Parameters of Eq 4 Derived from Binary Oxides' (pp. 11782-11783); "
       "transcribed from the paper PDF in this session")
pot = {
    'name': 'Pedone 2006 partial-charge rigid-ion potential',
    'form': 'U(r) = q_i q_j e^2/r + D_ij[(1-exp(-a_ij(r-r0)))^2 - 1] + C_ij/r^12',
    'charges_e': {k: v for k, v in CHG.items() if k != 'Al'},
    'charges_source': SRC + ' (q_O = -1.2 fixed; cation charges in the Table 2 row labels, e.g. Er1.8, Ge2.4, Si2.4)',
    'pairs': {f'{a}-{b}': {'D_eV': p[0], 'a_invA': p[1], 'r0_A': p[2], 'C_eVA12': p[3], 'source': SRC}
              for (a, b), p in PED.items()},
    'not_used_but_available_Al-O': {'D_eV': 0.361581, 'a_invA': 1.900442, 'r0_A': 2.164818, 'C_eVA12': 0.9,
                                    'source': SRC, 'note': 'no Al in this glass'},
    'cation_cation': 'Coulomb only (as in Pedone 2006)',
    'implementation': ('LAMMPS 22 Jul 2025 (PyPI wheel), pair_style table linear 12000 for ALL pairs: '
                       'Coulomb by damped shifted force (DSF, Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006)), '
                       'alpha = 0.2 1/A, cutoff 10.0 A; Morse + C/r^12 tabulated to 10.0 A and energy-shifted to zero '
                       'at the cutoff. Table spacing 0.0008 A from 0.5 A.'),
    'note_on_coulomb': ('PPPM (rc 11 A, 1e-5) was benchmarked at 25 ms/step on one core; with only 2 cores available '
                        'the tabulated DSF form (8.5-10 ms/step) was used to fit the time budget.'),
}
FAST = os.environ.get('FAST', '').split()
reps = reps + FAST
timing = {}
for r in reps:
    try:
        t0 = int(open(os.path.join(r, 't_start')).read()); t1 = int(open(os.path.join(r, 't_end')).read())
        timing[r] = round((t1 - t0) / 60, 1)
    except Exception:
        pass
# final-stage temperatures from logs
Tfin = {}
for r in reps:
    rows = []
    for line in open(os.path.join(r, 'log.lammps')):
        p = line.split()
        if len(p) == 5 and p[0].isdigit():
            rows.append([float(x) for x in p])
    rows = np.array(rows)
    # last 25 thermo lines = last 50 ps (300 K stage)
    nl = 15 if r in FAST else 25
    Tfin[r] = round(float(rows[-nl:, 1].mean()), 1)
summ = ana['summary']
summ['final_stage_mean_T_K'] = Tfin
sites = ana['sites']
import numpy as _np
for rc in (2.5, 2.6, 2.8):
    summ[f'CN_within_{rc}A_mean'] = round(float(_np.mean([sum(t['r'] < rc for t in s['ligand_tags']) for s in sites])), 3)
FSITES, FSUMM = [], None
if FAST:
    fa = json.load(open(os.environ['FASTANA']))
    FSITES = fa['sites']; FSUMM = fa['summary']
    FSUMM['note'] = 'replicas ' + ', '.join(FAST) + ': shortened protocol 6000 K 20 ps, 3500 K 20 ps, 3500->300 K at 20 K/ps (160 ps), 300 K 30 ps'
    FSUMM['final_stage_mean_T_K'] = {r: Tfin[r] for r in FAST}
    for rc in (2.5, 2.6, 2.8):
        FSUMM[f'CN_within_{rc}A_mean'] = round(float(_np.mean([sum(t['r'] < rc for t in s['ligand_tags']) for s in FSITES])), 3)
    summ['final_stage_mean_T_K'] = {r: Tfin[r] for r in reps if r not in FAST}
SEEDS = {'rep1': 11, 'rep2': 22, 'rep3': 33, 'rep4': 44, 'rep5': 55, 'rep6': 66}
doc = {
    'glass': 'Er-Ge10-SiO2',
    'description': 'Er-doped germanosilicate (10 mol% GeO2), MCVD-type fibre core; classical MD melt-quench',
    'composition': {'formula_units': {'SiO2': 900, 'GeO2': 100, 'Er2O3': 4},
                    'mol_percent': {'SiO2': round(900 / 1004 * 100, 3), 'GeO2': round(100 / 1004 * 100, 3),
                                    'Er2O3': round(4 / 1004 * 100, 3)},
                    'atoms': N, 'n_atoms': sum(N.values())},
    'potential': pot,
    'protocol': {
        'ensemble': 'NVT, Nose-Hoover (Tdamp 0.1 ps); fixed density',
        'density_g_cm3': DENS, 'box_A': round(box_length(), 4), 'n_atoms': sum(N.values()),
        'timestep_fs': 1.0,
        'start': 'random placement with minimum distances (cation-cation 2.6, cation-O 1.55, O-O 2.2, Er-Er 4.0 A); '
                 'CG minimisation; 3 ps nve/limit + Berendsen at 6000 K to relax overlaps',
        'stages': [
            {'T_K': 6000, 'time_ps': 20}, {'T_K': 3500, 'time_ps': 50},
            {'T_K': '3500->300 linear', 'time_ps': 320, 'quench_rate_K_per_ps': 10},
            {'T_K': 300, 'time_ps': 50}],
        'quench_rate_K_per_ps': 10,
        'quench_rate_note': '10 K/ps instead of 5 K/ps: only 2 CPU cores available (see .md)',
        'replicas': len(reps), 'replicas_main_10Kps': [r for r in reps if r not in FAST], 'seeds': {r: SEEDS.get(os.path.basename(r)) for r in reps}, 'replica_wall_minutes': timing,
        'fast_replicas': {'replicas': FAST, 'stages': [{'T_K': 6000, 'time_ps': 20}, {'T_K': 3500, 'time_ps': 20}, {'T_K': '3500->300 linear', 'time_ps': 160, 'quench_rate_K_per_ps': 20}, {'T_K': 300, 'time_ps': 30}], 'note': 'added at the end of the budget for more distinct sites; sites tagged quench_rate_K_per_ps=20, statistics in summary_fast_quench_replicas'} if FAST else None,
        'sampling': '10 frames every 2 ps over the last 20 ps (300 K), per replica',
    },
    'er_sites': [dict(x, quench_rate_K_per_ps=10) for x in ana['sites']] + [dict(x, quench_rate_K_per_ps=20) for x in FSITES],
    'summary': summ,
    'summary_fast_quench_replicas': FSUMM,
    'rdf': {'r_A': ana['rdf_r'], **ana['rdfs']},
}
json.dump(doc, open(out, 'w'), separators=(',', ':'))
print('wrote', out, os.path.getsize(out) / 1e6, 'MB')
