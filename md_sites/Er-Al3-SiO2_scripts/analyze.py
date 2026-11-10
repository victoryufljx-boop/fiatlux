import numpy as np, json, sys, glob, os
EL={1:'Si',2:'Al',3:'Er',4:'O'}
def read_dump(fn):
    frames=[]
    with open(fn) as f:
        while True:
            l=f.readline()
            if not l: break
            if l.startswith('ITEM: TIMESTEP'):
                step=int(f.readline()); f.readline(); n=int(f.readline()); f.readline()
                b=[list(map(float,f.readline().split())) for _ in range(3)]; f.readline()
                a=np.array([f.readline().split() for _ in range(n)],float)
                L=np.array([x[1]-x[0] for x in b]); lo=np.array([x[0] for x in b])
                frames.append((step,a[:,1].astype(int),a[:,2:5]-lo,L))
    return frames
def mic(d,L): return d-L*np.round(d/L)
def dist_to(pos,i,idx,L):
    d=mic(pos[idx]-pos[i],L); return d,np.sqrt((d*d).sum(1))
def rdf(pos,typ,L,ta,tb,rmax=6.0,nb=300):
    A=np.where(typ==ta)[0];B=np.where(typ==tb)[0];h=np.zeros(nb)
    for i in A:
        _,r=dist_to(pos,i,B,L); r=r[(r>1e-6)&(r<rmax)]; h+=np.histogram(r,nb,(0,rmax))[0]
    e=np.linspace(0,rmax,nb+1);c=0.5*(e[1:]+e[:-1]);sh=4/3*np.pi*(e[1:]**3-e[:-1]**3)
    rho=len(B)/np.prod(L); return c,h/(len(A)*rho*sh),h
def first_min(c,g,start,stop):
    m=(c>start)&(c<stop); i=np.argmin(g[m]); return float(c[m][i])
def angles(pos,typ,L,center,lig,rc):
    out=[]
    C=np.where(typ==center)[0];Lg=np.where(typ==lig)[0]
    for i in C:
        d,r=dist_to(pos,i,Lg,L); d=d[(r<rc)&(r>1e-6)]
        for a in range(len(d)):
            for b in range(a+1,len(d)):
                out.append(np.degrees(np.arccos(np.clip(d[a]@d[b]/np.linalg.norm(d[a])/np.linalg.norm(d[b]),-1,1))))
    return np.array(out)
def hist(x,lo,hi,n): h,e=np.histogram(x,n,(lo,hi)); return {'bin_lo':lo,'bin_hi':hi,'n_bins':n,'counts':h.tolist()}
def stats(x): return {'mean':round(float(np.mean(x)),3),'std':round(float(np.std(x)),3),'n':int(len(x))} if len(x) else None

reps=sorted(glob.glob('rep*_final.dump'),key=lambda s:int(s[3:].split('_')[0]))
reps=[r for r in reps if os.path.exists(r.replace('_final.dump','.end'))]
allF={}
for fn in reps:
    rid=int(fn[3:].split('_')[0]); allF[rid]=read_dump(fn)[-10:]
# pooled RDFs for cutoffs
acc={}
for rid,F in allF.items():
    for step,typ,pos,L in F:
        for key,(a,b) in {'SiO':(1,4),'AlO':(2,4),'ErO':(3,4),'OO':(4,4),'SiSi':(1,1)}.items():
            c,g,h=rdf(pos,typ,L,a,b,rmax=6.0 if key!='ErO' else 8.0)
            acc.setdefault(key,[c,0,0]); acc[key][1]=acc[key][1]+g/1.0; acc[key][2]=acc[key][2]+1
rd={k:(v[0],v[1]/v[2]) for k,v in acc.items()}
peak=lambda k,lo,hi:float(rd[k][0][(rd[k][0]>lo)&(rd[k][0]<hi)][np.argmax(rd[k][1][(rd[k][0]>lo)&(rd[k][0]<hi)])])
rmin_ErO=first_min(*rd['ErO'],2.6,3.6)
rmin_SiO=first_min(*rd['SiO'],1.8,2.6); rmin_AlO=first_min(*rd['AlO'],1.9,2.8)
sites=[];SiO_d=[];AlO_d=[];ErO_d=[];AlCN=[];SiCN=[];OSiO=[];SiOSi=[];OAlO=[];OErO=[];ErErnn=[];CN=[];nbo=0;nO=0;linkSi=0;linkAl=0
ObondTot={'BO':0,'NBO':0,'free':0,'tri':0}
for rid,F in allF.items():
    for fi,(step,typ,pos,L) in enumerate(F):
        Si=np.where(typ==1)[0];Al=np.where(typ==2)[0];Er=np.where(typ==3)[0];O=np.where(typ==4)[0]
        # former bonds per O
        formers={}
        for o in O:
            _,rs=dist_to(pos,o,Si,L);_,ra=dist_to(pos,o,Al,L)
            formers[o]=(['Si']*int((rs<2.0).sum())+['Al']*int((ra<2.1).sum()))
        for o in O:
            n=len(formers[o]); ObondTot['free' if n==0 else 'NBO' if n==1 else 'BO' if n==2 else 'tri']+=1
        for i in Si:
            _,r=dist_to(pos,i,O,L); SiCN.append(int((r<rmin_SiO).sum())); SiO_d+=r[r<rmin_SiO].tolist()
        for i in Al:
            _,r=dist_to(pos,i,O,L); AlCN.append(int((r<rmin_AlO).sum())); AlO_d+=r[r<rmin_AlO].tolist()
        OSiO+=angles(pos,typ,L,1,4,2.0).tolist(); OAlO+=angles(pos,typ,L,2,4,2.1).tolist()
        for o in O:
            d,r=dist_to(pos,o,Si,L); d=d[r<2.0]
            if len(d)==2: SiOSi.append(np.degrees(np.arccos(np.clip(d[0]@d[1]/np.linalg.norm(d[0])/np.linalg.norm(d[1]),-1,1))))
        for e in Er:
            others=Er[Er!=e]; _,re_=dist_to(pos,e,others,L); ErErnn.append(float(re_.min()))
            allidx=np.arange(len(typ)); allidx=allidx[allidx!=e]
            d,r=dist_to(pos,e,allidx,L); m=r<rmin_ErO
            ligs=[];tags=[]
            for j,dv,rv in sorted(zip(allidx[m],d[m],r[m]),key=lambda t:t[2]):
                el=EL[typ[j]]; ligs.append([el,round(float(dv[0]),4),round(float(dv[1]),4),round(float(dv[2]),4)])
                if el=='O':
                    fm=formers[j]; nf=len(fm); kind='free' if nf==0 else 'NBO' if nf==1 else 'BO' if nf==2 else 'TBO'
                    tags.append({'type':kind,'formers':fm,'r':round(float(rv),4)}); ErO_d.append(float(rv)); nO+=1
                    nbo+= kind in ('NBO','free'); linkSi+=fm.count('Si'); linkAl+=fm.count('Al')
                else: tags.append({'type':'cation','r':round(float(rv),4)})
            CN.append(sum(1 for l in ligs if l[0]=='O'))
            sites.append({'replica':rid,'frame':fi,'step':int(step),'er_index':int(e),'ligands':ligs,'ligand_tags':tags})
            dd=np.array([l[1:] for l in ligs if l[0]=='O'])
            for a in range(len(dd)):
                for b in range(a+1,len(dd)): OErO.append(np.degrees(np.arccos(np.clip(dd[a]@dd[b]/np.linalg.norm(dd[a])/np.linalg.norm(dd[b]),-1,1))))
CN=np.array(CN);ErO_d=np.array(ErO_d)
nAl_total=sum(len(np.where(F[0][1]==2)[0]) for F in allF.values()); nSi_total=sum(len(np.where(F[0][1]==1)[0]) for F in allF.values())
summary={'n_replicas':len(allF),'replicas':sorted(allF),'n_distinct_er':len(allF)*8,'n_site_frames':len(sites),
 'r_min':round(rmin_ErO,3),'ErO_rdf_peak':round(peak('ErO',1.8,3.0),3),
 'CN_mean':round(float(CN.mean()),3),'CN_std':round(float(CN.std()),3),'CN_hist':{int(k):int((CN==k).sum()) for k in range(CN.min(),CN.max()+1)},
 'R_mean':round(float(ErO_d.mean()),4),'R_std':round(float(ErO_d.std()),4),
 'nbo_fraction':round(nbo/nO,4),'nbo_definition':'O bonded to <=1 network former (Si<2.0 A, Al<2.1 A); free O (0 formers) counted as NBO',
 'ligand_O_types':{k:int(sum(1 for s in sites for t in s['ligand_tags'] if t.get('type')==k)) for k in ['BO','NBO','free','TBO']},
 'er_o_former_linkages':{'Er-O-Si':linkSi,'Er-O-Al':linkAl,'frac_Al':round(linkAl/max(1,linkSi+linkAl),4),
     'random_expectation_frac_Al_by_cation_count':round(nAl_total/(nAl_total+nSi_total),4)},
 'er_er_nn':{**stats(ErErnn),'min':round(min(ErErnn),3),'frac_below_5A':round(float(np.mean(np.array(ErErnn)<5.0)),4)},
 'O_Er_O_angle':{**stats(OErO),'hist':hist(OErO,40,180,28)},
 'sanity':{'SiO_rdf_peak':round(peak('SiO',1.4,2.0),3),'SiO_mean_first_shell':stats(SiO_d),'Si_CN_mean':round(float(np.mean(SiCN)),3),
   'Si_CN_hist':{int(k):int((np.array(SiCN)==k).sum()) for k in sorted(set(SiCN))},'rmin_SiO':round(rmin_SiO,3),
   'O_Si_O':{**stats(OSiO),'hist':hist(OSiO,60,180,24)},'Si_O_Si':{**stats(SiOSi),'peak':None,'hist':hist(SiOSi,90,180,18)},
   'AlO_rdf_peak':round(peak('AlO',1.5,2.2),3),'AlO_mean_first_shell':stats(AlO_d),'rmin_AlO':round(rmin_AlO,3),'Al_CN_mean':round(float(np.mean(AlCN)),3),
   'Al_CN_hist':{int(k):int((np.array(AlCN)==k).sum()) for k in sorted(set(AlCN))},'O_Al_O':stats(OAlO),
   'OO_rdf_peak':round(peak('OO',2.2,3.0),3),'SiSi_rdf_peak':round(peak('SiSi',2.6,3.6),3),
   'oxygen_speciation_all_frames':ObondTot}}
h=summary['sanity']['Si_O_Si']['hist']; e=np.linspace(90,180,19); summary['sanity']['Si_O_Si']['peak']=float(0.5*(e[np.argmax(h['counts'])]+e[np.argmax(h['counts'])+1]))
ErCN_alt={}
for rc in (2.5,2.6,2.8,3.0):
    c=[sum(1 for t in st['ligand_tags'] if t.get('type') in ('BO','NBO','free','TBO') and t['r']<rc) for st in sites]
    ErCN_alt[str(rc)]={'mean':round(float(np.mean(c)),3),'hist':{int(k):int(np.sum(np.array(c)==k)) for k in sorted(set(c))}}
summary['CN_vs_cutoff']=ErCN_alt
summary['ErO_rdf']={'r':[round(x,3) for x in rd['ErO'][0][::2].tolist()],'g':[round(x,3) for x in rd['ErO'][1][::2].tolist()]}
json.dump({'summary':summary,'er_sites':sites},open('analysis.json','w'))
print(json.dumps({k:v for k,v in summary.items() if k!='ErO_rdf'},indent=1))
