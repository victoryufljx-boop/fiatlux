"""Build random initial config + Pedone tables for Er-Ge10-SiO2.
usage: python build.py SEED OUTDIR
"""
import sys, os, numpy as np

# Pedone et al. J. Phys. Chem. B 110, 11780 (2006), Table 2 (verbatim)
# pair: (D_ij eV, a_ij 1/A, r0 A, C_ij eV A^12)
PED = {
    ('Si', 'O'): (0.340554, 2.006700, 2.100000, 1.0),
    ('Al', 'O'): (0.361581, 1.900442, 2.164818, 0.9),
    ('Ge', 'O'): (0.158118, 2.294230, 2.261313, 5.0),
    ('Er', 'O'): (0.040448, 2.294078, 2.837722, 3.0),
    ('O', 'O'): (0.042395, 1.379316, 3.618701, 22.0),
}
CHG = {'O': -1.2, 'Si': 2.4, 'Ge': 2.4, 'Er': 1.8, 'Al': 1.8}
MASS = {'O': 15.9994, 'Si': 28.0855, 'Ge': 72.630, 'Er': 167.259}
TYPES = ['O', 'Si', 'Ge', 'Er']
N = {'Si': 900, 'Ge': 100, 'Er': 8}
N['O'] = 2 * 900 + 2 * 100 + 3 * 4
DENS = 2.35  # g/cm3
RC_SR = 6.0


def box_length():
    m = sum(N[e] * MASS[e] for e in N) / 6.02214076e23
    return (m / DENS * 1e24) ** (1 / 3)


def write_tables(path):
    r = np.linspace(0.6, RC_SR, 5401)
    with open(path, 'w') as f:
        f.write('# Pedone 2006 Morse + C/r^12 (no Coulomb), energy-shifted to 0 at rc\n')
        for (a, b), (D, al, r0, C) in PED.items():
            def U(x):
                return D * ((1 - np.exp(-al * (x - r0))) ** 2 - 1) + C / x ** 12
            def F(x):
                ex = np.exp(-al * (x - r0))
                return -(2 * D * al * (1 - ex) * ex) + 12 * C / x ** 13
            e = U(r) - U(RC_SR)
            f.write(f'\n{a}_{b}\nN {len(r)} R {r[0]} {r[-1]}\n\n')
            for i, (x, ee, ff) in enumerate(zip(r, e, F(r))):
                f.write(f'{i+1} {x:.6f} {ee:.10e} {ff:.10e}\n')


def place(seed, L):
    rng = np.random.default_rng(seed)
    order = ['Er'] * N['Er'] + ['Ge'] * N['Ge'] + ['Si'] * N['Si'] + ['O'] * N['O']
    dmin = {('cat', 'cat'): 2.6, ('cat', 'O'): 1.55, ('O', 'O'): 2.2}
    dmin_er = 4.0  # keep RE ions apart initially (no forced clustering)
    pos, el = [], []
    for e in order:
        for tries in range(200000):
            x = rng.random(3) * L
            if pos:
                P = np.array(pos)
                d = P - x
                d -= L * np.round(d / L)
                d = np.sqrt((d ** 2).sum(1))
                ok = True
                for e2 in set(el):
                    m = np.array([q == e2 for q in el])
                    k = ('O', 'O') if (e == 'O' and e2 == 'O') else (('cat', 'O') if 'O' in (e, e2) else ('cat', 'cat'))
                    lim = dmin[k]
                    if e == 'Er' and e2 == 'Er':
                        lim = dmin_er
                    if d[m].min() < lim:
                        ok = False
                        break
                if not ok:
                    continue
            pos.append(x); el.append(e)
            break
        else:
            raise RuntimeError('placement failed')
    return np.array(pos), el


def write_data(path, pos, el, L):
    with open(path, 'w') as f:
        f.write(f'Er-Ge10-SiO2 random start\n\n{len(el)} atoms\n{len(TYPES)} atom types\n\n')
        f.write(f'0 {L:.6f} xlo xhi\n0 {L:.6f} ylo yhi\n0 {L:.6f} zlo zhi\n\nMasses\n\n')
        for i, t in enumerate(TYPES):
            f.write(f'{i+1} {MASS[t]}\n')
        f.write('\nAtoms # charge\n\n')
        for i, (p, e) in enumerate(zip(pos, el)):
            f.write(f'{i+1} {TYPES.index(e)+1} {CHG[e]} {p[0]:.5f} {p[1]:.5f} {p[2]:.5f}\n')


if __name__ == '__main__':
    seed, out = int(sys.argv[1]), sys.argv[2]
    os.makedirs(out, exist_ok=True)
    L = box_length()
    pos, el = place(seed, L)
    q = sum(CHG[e] for e in el)
    assert abs(q) < 1e-9, q
    write_data(os.path.join(out, 'start.data'), pos, el, L)
    write_tables(os.path.join(out, 'pedone.table'))
    print('L =', L, 'natoms', len(el), 'charge', q)
