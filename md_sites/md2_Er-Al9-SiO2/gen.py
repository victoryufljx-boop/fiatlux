import numpy as np, sys
seed=int(sys.argv[1]); out=sys.argv[2]
rng=np.random.default_rng(seed)
nSi,nAl2,nEr2=870,86,4
species=['O','Si','Al','Er']; q={'O':-1.2,'Si':2.4,'Al':1.8,'Er':1.8}
mass={'O':15.9994,'Si':28.0855,'Al':26.981539,'Er':167.259}
counts={'Si':nSi,'Al':2*nAl2,'Er':2*nEr2,'O':2*nSi+3*nAl2+3*nEr2}
M=sum(mass[s]*n for s,n in counts.items())
rho=2.28; V=M/6.02214076e23/rho*1e24; L=V**(1/3)
types=[]
for s in ['Er','Al','Si','O']: types+= [s]*counts[s]
pos=np.zeros((len(types),3)); dmin=1.7
for i,s in enumerate(types):
    for t in range(100000):
        p=rng.random(3)*L
        if i==0: break
        d=pos[:i]-p; d-=L*np.round(d/L)
        if (d*d).sum(1).min()>dmin**2: break
    else: raise SystemExit('placement failed')
    pos[i]=p
with open(out,'w') as f:
    f.write(f'Er-Al9-SiO2 seed {seed}\n\n{len(types)} atoms\n4 atom types\n\n0 {L:.6f} xlo xhi\n0 {L:.6f} ylo yhi\n0 {L:.6f} zlo zhi\n\nMasses\n\n')
    for k,s in enumerate(species): f.write(f'{k+1} {mass[s]}\n')
    f.write('\nAtoms # charge\n\n')
    for i,(s,p) in enumerate(zip(types,pos)):
        f.write(f'{i+1} {species.index(s)+1} {q[s]} {p[0]:.5f} {p[1]:.5f} {p[2]:.5f}\n')
print(len(types),counts,L, 'mol% Al2O3', nAl2/(nSi+nAl2+nEr2)*100, 'net q', sum(q[s] for s in types))
