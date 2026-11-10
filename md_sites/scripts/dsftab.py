import numpy as np, sys
from scipy.special import erfc
sys.path.insert(0,'/home/claude/md')
from build import PED, CHG, TYPES
KE=14.399645
def write(path, rc=10.0, alpha=0.2, rin=0.5, n=12000):
    r=np.linspace(rin, rc, n)
    erc=erfc(alpha*rc)/rc; frc=erfc(alpha*rc)/rc**2+2*alpha/np.sqrt(np.pi)*np.exp(-alpha**2*rc**2)/rc
    with open(path,'w') as f:
        f.write('# DSF Coulomb (Fennell-Gezelter, alpha=%g, rc=%g) + Pedone 2006 Morse + C/r^12\n'%(alpha,rc))
        for i,a in enumerate(TYPES):
            for b in TYPES[i:]:
                qq=KE*CHG[a]*CHG[b]
                E=qq*(erfc(alpha*r)/r-erc+frc*(r-rc))
                F=qq*(erfc(alpha*r)/r**2+2*alpha/np.sqrt(np.pi)*np.exp(-alpha**2*r**2)/r-frc)
                key=tuple(sorted((a,b),key=lambda s:(s!='O',s)))
                k2=(b,a) if (b,a) in PED else (a,b)
                if k2 in PED:
                    D,al,r0,C=PED[k2]
                    ex=np.exp(-al*(r-r0)); exc=np.exp(-al*(rc-r0))
                    E=E+D*((1-ex)**2-1)+C/r**12-(D*((1-exc)**2-1)+C/rc**12)
                    F=F-2*D*al*(1-ex)*ex+12*C/r**13
                f.write(f'\n{a}_{b}\nN {n} R {rin} {rc}\n\n')
                for j in range(n): f.write(f'{j+1} {r[j]:.6f} {E[j]:.10e} {F[j]:.10e}\n')
if __name__=='__main__': write(sys.argv[1])
