import numpy as np, sys
seed=int(sys.argv[1]); out=sys.argv[2]
nSi,nAl,nEr=966,68,8
nO=966*2+34*3+4*3
rng=np.random.default_rng(seed)
mass={'Si':28.0855,'Al':26.981538,'Er':167.259,'O':15.9994}
M=nSi*mass['Si']+nAl*mass['Al']+nEr*mass['Er']+nO*mass['O']
L=(M/6.02214076e23/2.22*1e24)**(1/3)
species=['Er']*nEr+['Al']*nAl+['Si']*nSi+['O']*nO
dmin={('O','O'):2.2}; 
def dm(a,b):
    if a=='O' and b=='O': return 2.2
    if a=='O' or b=='O': return 1.5
    return 2.8
pos=[];sp=[]
for s in species:
    for t in range(100000):
        x=rng.random(3)*L
        if pos:
            P=np.array(pos); d=P-x; d-=L*np.round(d/L); r=np.sqrt((d*d).sum(1))
            ok=all(r[i]>=dm(s,sp[i]) for i in np.nonzero(r<3.0)[0])
        else: ok=True
        if ok: pos.append(x); sp.append(s); break
    else: raise SystemExit('placement failed')
tid={'Si':1,'Al':2,'Er':3,'O':4}; q={'Si':2.4,'Al':1.8,'Er':1.8,'O':-1.2}
with open(out,'w') as f:
    f.write(f'Er-Al3-SiO2 seed {seed}\n\n{len(pos)} atoms\n4 atom types\n0 {L:.6f} xlo xhi\n0 {L:.6f} ylo yhi\n0 {L:.6f} zlo zhi\n\nMasses\n\n')
    for s in ['Si','Al','Er','O']: f.write(f'{tid[s]} {mass[s]}\n')
    f.write('\nAtoms # charge\n\n')
    for i,(s,x) in enumerate(zip(sp,pos)): f.write(f'{i+1} {tid[s]} {q[s]} {x[0]:.5f} {x[1]:.5f} {x[2]:.5f}\n')
print(len(pos),L, sum(q[s] for s in sp))
