import os, json, sys
from pyscf import lib; lib.num_threads(2)
from gainmat.cluster import ti_3d, ce_5d
o = {}
o['ti'] = ti_3d.run('fast'); print('TI', o['ti']['abs_eg_mean_cm1'], o['ti']['timings'], flush=True)
o['ce'] = ce_5d.run('fast'); print('CE', o['ce']['e_pt2']['so_5d_cm1'], o['ce']['e_pt2']['em_to_2F52_cm1'], o['ce']['timings'], flush=True)
json.dump(o, open('runs/baseline.json','w'), default=float)
