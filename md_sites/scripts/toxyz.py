import sys
from ase.io import read, write
for r in sys.argv[1:]:
    a = read(f'{r}/final.data', format='lammps-data', atom_style='charge', Z_of_type={1: 8, 2: 14, 3: 32, 4: 68})
    from ase import Atoms
    a = Atoms(a.get_chemical_symbols(), positions=a.positions, cell=a.cell, pbc=True)
    a.info['comment'] = f'Er-Ge10-SiO2 {r} final 300 K, Pedone 2006 (DSF Coulomb)'
    write(f'/home/claude/fiatlux/md_sites/Er-Ge10-SiO2_{r}_final.extxyz', a, format='extxyz')
