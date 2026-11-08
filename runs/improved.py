import json, sys
from pyscf import lib; lib.num_threads(1)
from gainmat.cluster import ti_3d, ce_5d
which = sys.argv[1]
if which == 'ti':
    o = ti_3d.run('fast+')
    print('TI abs', o['abs_eg_mean_cm1'], 'em', o.get('em_cm1'), 'stokes', o.get('stokes_cm1'), o['timings'], flush=True)
else:
    o = ce_5d.run('fast+', relax_grid=(0.0, 0.05, 0.10))
    r = o['e_pt2']
    print('CE abs', r['so_5d_cm1'], 'em_vert', r['em_to_2F52_cm1'], r['em_to_2F72_cm1'], 'em_rel', r.get('em_relaxed_to_2F52_cm1'), r.get('em_relaxed_to_2F72_cm1'), 'stokes', r.get('stokes_cm1'), o['timings'], flush=True)
json.dump(o, open(f'runs/improved_{which}.json','w'), default=lambda x: x.tolist() if hasattr(x,'tolist') else str(x))
