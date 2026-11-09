
import json, os, sys, time, traceback, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, '/home/claude/fiatlux')
from ase.io import read
from gpaw.mpi import world
from gainmat.host.gap import compute_gap, compute_eps_inf
atoms = read('/home/claude/fiatlux/bench/yag/caf2_mpi/CaF2.traj')
gap_kw = {'ecut': 500, 'kdens': 3.5, 'band_edges': 'path', 'method': 'gllbsc'}
eps_kw = None
out = dict(nproc=world.size, errors={})
t0 = time.time()
if gap_kw is not None:
    try:
        out["gap"] = compute_gap(atoms, txt='/home/claude/fiatlux/bench/yag/caf2_mpi/CaF2_gap.txt', **gap_kw)
    except Exception as e:
        out["errors"]["gap"] = "%s: %s" % (type(e).__name__, e)
        out["traceback_gap"] = traceback.format_exc(limit=4)
if eps_kw is not None:
    try:
        g = out.get("gap", {})
        if eps_kw.pop("scissor_from_gap", False) and g.get("gap_eV") and g.get("gap_pbe_eV"):
            eps_kw["scissor_eV"] = g["gap_eV"] - g["gap_pbe_eV"]
        out["eps"] = compute_eps_inf(atoms, workdir='/home/claude/fiatlux/bench/yag/caf2_mpi', txt='/home/claude/fiatlux/bench/yag/caf2_mpi/CaF2_eps.txt', **eps_kw)
    except Exception as e:
        out["errors"]["eps"] = "%s: %s" % (type(e).__name__, e)
        out["traceback_eps"] = traceback.format_exc(limit=4)
out["t_wall_s"] = round(time.time() - t0, 1)
if world.rank == 0:
    with open('/home/claude/fiatlux/bench/yag/caf2_mpi/CaF2_gpaw.json', "w") as f:
        json.dump(out, f, default=float)
