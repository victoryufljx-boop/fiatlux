"""Build a random Er2O3-SiO2 configuration + LAMMPS input for the Pedone (2006) potential.

usage: python build_er_sio2.py OUTDIR SEED [cool_rate_K_per_ps]
"""
import sys, os
import numpy as np

outdir, seed = sys.argv[1], int(sys.argv[2])
rate = float(sys.argv[3]) if len(sys.argv) > 3 else 10.0
os.makedirs(outdir, exist_ok=True)
rng = np.random.default_rng(seed)

N_SIO2, N_ER2O3 = 1000, 4
n = {"Si": N_SIO2, "Er": 2 * N_ER2O3, "O": 2 * N_SIO2 + 3 * N_ER2O3}
mass = {"O": 15.999, "Si": 28.0855, "Er": 167.259}
q = {"O": -1.2, "Si": 2.4, "Er": 1.8}
types = {"O": 1, "Si": 2, "Er": 3}
NA = 6.02214076e23
rho = 2.20  # g/cm3
M = sum(n[e] * mass[e] for e in n)  # g/mol
V = M / NA / rho * 1e24  # A^3
L = V ** (1 / 3)

# random placement with minimum-distance check (cation-cation 3.0, cation-O 1.5, O-O 2.2)
dmin = {("Er", "Er"): 4.5, ("Si", "Si"): 2.9, ("Er", "Si"): 3.0,
        ("O", "O"): 2.2, ("O", "Si"): 1.5, ("Er", "O"): 2.0}
order = ["Er"] * n["Er"] + ["Si"] * n["Si"] + ["O"] * n["O"]
pos, el = [], []
for e in order:
    for tries in range(200000):
        x = rng.random(3) * L
        if pos:
            d = np.array(pos) - x
            d -= L * np.round(d / L)
            r = np.sqrt((d * d).sum(1))
            lim = np.array([dmin[tuple(sorted((e, f)))] for f in el])
            if np.any(r < lim):
                continue
        pos.append(x); el.append(e)
        break
    else:
        raise RuntimeError("placement failed")
pos = np.array(pos)
with open(f"{outdir}/init.data", "w") as f:
    f.write(f"Er-SiO2 seed {seed}\n\n{len(el)} atoms\n3 atom types\n\n")
    f.write(f"0 {L:.6f} xlo xhi\n0 {L:.6f} ylo yhi\n0 {L:.6f} zlo zhi\n\nMasses\n\n")
    for e in ["O", "Si", "Er"]:
        f.write(f"{types[e]} {mass[e]}\n")
    f.write("\nAtoms # charge\n\n")
    for i, (e, x) in enumerate(zip(el, pos)):
        f.write(f"{i+1} {types[e]} {q[e]} {x[0]:.5f} {x[1]:.5f} {x[2]:.5f}\n")

# ---- Pedone et al. JPCB 110, 11780 (2006), Table 2: D (eV), a (1/A), r0 (A), C (eV A^12)
# Everything (DSF-damped Coulomb + Morse + C/r^12) is tabulated per pair: one pair loop, no k-space.
from scipy.special import erfc
KE = 14.399645  # e^2/(4 pi eps0) in eV*A
ALPHA = 0.2     # DSF damping (1/A), Fennell & Gezelter, JCP 124, 234104 (2006)
RC = 10.0
P = {("O", "O"): (0.042395, 1.379316, 3.618701, 22.0),
     ("O", "Si"): (0.340554, 2.006700, 2.100000, 1.0),
     ("Er", "O"): (0.040448, 2.294078, 2.837722, 3.0)}
pairs = [("O", "O"), ("O", "Si"), ("Er", "O"), ("Si", "Si"), ("Er", "Si"), ("Er", "Er")]
rs = np.linspace(0.4, RC, 19201)
ec = erfc(ALPHA * RC) / RC
fc = erfc(ALPHA * RC) / RC**2 + 2 * ALPHA / np.sqrt(np.pi) * np.exp(-(ALPHA * RC) ** 2) / RC
with open(f"{outdir}/pedone.table", "w") as f:
    for p in pairs:
        qq = KE * q[p[0]] * q[p[1]]
        u = qq * (erfc(ALPHA * rs) / rs - ec + fc * (rs - RC))
        fr = qq * (erfc(ALPHA * rs) / rs**2 + 2 * ALPHA / np.sqrt(np.pi) * np.exp(-(ALPHA * rs) ** 2) / rs - fc)
        if p in P:
            D, a, r0, C = P[p]
            ex = np.exp(-a * (rs - r0))
            u = u + D * ((1 - ex) ** 2 - 1) + C / rs ** 12
            fr = fr - (D * 2 * (1 - ex) * a * ex - 12 * C / rs ** 13)
            exc = np.exp(-a * (RC - r0))
            u = u - (D * ((1 - exc) ** 2 - 1) + C / RC ** 12)  # shift (< 2e-5 eV)
        f.write(f"{p[0]}_{p[1]}\nN {len(rs)} R {rs[0]} {rs[-1]}\n\n")
        for k, (r, uu, ff) in enumerate(zip(rs, u, fr)):
            f.write(f"{k+1} {r:.6f} {uu:.12e} {ff:.12e}\n")
        f.write("\n")
pc = "\n".join(f"pair_coeff {min(types[a],types[b])} {max(types[a],types[b])} pedone.table {a}_{b} {RC}" for a, b in pairs)
t_cool = (3500 - 300) / rate  # ps
nsteps_cool = int(round(t_cool * 1000))
inp = f"""units metal
atom_style charge
boundary p p p
read_data init.data
pair_style table linear 19201
{pc}
neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
timestep 0.001
thermo 2000
thermo_style custom step time temp pe etotal press cpu
# soft start: minimise
min_style cg
minimize 1e-6 1e-8 2000 20000
reset_timestep 0
velocity all create 6000 {seed} mom yes rot yes dist gaussian
fix th all nvt temp 6000 6000 0.1
run 20000
unfix th
fix th all nvt temp 3500 3500 0.1
run 50000
write_restart melt.restart
unfix th
fix th all nvt temp 3500 300 0.1
run {nsteps_cool}
write_restart quenched.restart
unfix th
fix th all nvt temp 300 300 0.1
run 30000
dump d all custom 2000 final20ps.dump id type q x y z
dump_modify d sort id
run 20000
write_data final.data
"""
open(f"{outdir}/in.lmp", "w").write(inp)
print(f"L={L:.4f} A, N={len(el)}, cool {t_cool:.0f} ps at {rate} K/ps")
