import json
from pyscf import lib; lib.num_threads(1)
from gainmat.cluster import ti_3d
o = ti_3d.run('fast+', relax=False)
print('TI fast+ vertical', o['abs_eg_mean_cm1'], o['e_pt2_cm1'], o['timings'], flush=True)
json.dump(o, open('runs/ti_vert.json','w'), default=lambda x: x.tolist() if hasattr(x,'tolist') else str(x))
