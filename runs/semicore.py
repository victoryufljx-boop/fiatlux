import json
from pyscf import lib; lib.num_threads(1)
from gainmat.cluster import ce_5d
o = ce_5d.run('fast+', relax=False)
print('SEMI', o['semicore'], o['e_pt2']['so_5d_cm1'], flush=True)
json.dump(o, open('runs/semicore_ce.json','w'), default=lambda x: x.tolist() if hasattr(x,'tolist') else str(x))
