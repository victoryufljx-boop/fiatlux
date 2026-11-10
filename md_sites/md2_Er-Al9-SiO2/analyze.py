import numpy as np, json, sys, glob
EL={1:'O',2:'Si',3:'Al',4:'Er'}
def read_dump(fn):
    frames=[]
    with open(fn) as f:
        lines=f.read().split('\n')
    i=0
    while i<len(lines):
        if lines[i].startswith('ITEM: TIMESTEP'):
            step=int(lines[i+1]); n=int(lines[i+3])
            b=np.array([list(map(float,lines[i+5+k].split()[:2])) for k in range(3)])
            L=b[:,1]-b[:,0]
            a=np.array([lines[i+9+k].split() for k in range(n)],float)
            frames.append(dict(step=step,L=L,ids=a[:,0].astype(int),types=a[:,1].astype(int),pos=a[:,2:5]-b[:,0]))
            i+=9+n
        else: i+=1
    return frames
def mic(d,L): return d-L*np.round(d/L)
def dists(pos,i_idx,j_idx,L):
    d=mic(pos[j_idx][None,:,:]-pos[i_idx][:,None,:],L)
    return d,np.sqrt((d*d).sum(-1))
def rdf(fr_list,ti,tj,rmax=6.0,nb=300):
    h=np.zeros(nb); norm=0
    for fr in fr_list:
        ii=np.where(fr['types']==ti)[0]; jj=np.where(fr['types']==tj)[0]
        _,r=dists(fr['pos'],ii,jj,fr['L'])
        if ti==tj: r=r[r>1e-6]
        h+=np.histogram(r.ravel(),bins=nb,range=(0,rmax))[0]
        V=np.prod(fr['L']); norm+=len(ii)*len(jj)/V
    e=np.linspace(0,rmax,nb+1); rc=0.5*(e[1:]+e[:-1]); shell=4/3*np.pi*(e[1:]**3-e[:-1]**3)
    return rc,h/(norm*shell)
def first_min(rc,g,rlo,rhi):
    m=(rc>rlo)&(rc<rhi); gs=np.convolve(g,np.ones(5)/5,'same')
    return float(rc[m][np.argmin(gs[m])])
def angles(fr,tc,tl,rcut):
    """angles l-c-l for centre type tc with ligands of type tl within rcut"""
    out=[]
    cc=np.where(fr['types']==tc)[0]; ll=np.where(fr['types']==tl)[0]
    d,r=dists(fr['pos'],cc,ll,fr['L'])
    for k in range(len(cc)):
        nb=np.where(r[k]<rcut)[0]
        v=d[k,nb]/r[k,nb][:,None]
        for a in range(len(nb)):
            for b in range(a+1,len(nb)):
                out.append(np.degrees(np.arccos(np.clip(v[a]@v[b],-1,1))))
    return out
def hist_stats(x):
    x=np.asarray(x); 
    if len(x)==0: return None
    h,e=np.histogram(x,bins=36,range=(0,180))
    return dict(mean=float(x.mean()),std=float(x.std()),peak=float(0.5*(e[np.argmax(h)]+e[np.argmax(h)+1])),n=int(len(x)),hist_5deg=h.tolist())

def analyze(dumps, out_prefix):
    reps={}
    for fn in dumps:
        rep=int(fn.split('/')[-1].split('.')[0][1:])
        fr=read_dump(fn)[-10:]
        reps[rep]=fr
    allfr=[f for r in sorted(reps) for f in reps[r]]
    res={}
    rc,gEr=rdf(allfr,4,1); rmin=first_min(rc,gEr,2.5,3.6)
    res['er_o_rdf_peak']=float(rc[np.argmax(gEr)]); res['r_min']=rmin
    rc,gSi=rdf(allfr,2,1); rc,gAl=rdf(allfr,3,1); rc,gOO=rdf(allfr,1,1)
    res['si_o_peak']=float(rc[np.argmax(gSi*(rc<2.5))]); res['al_o_peak']=float(rc[np.argmax(gAl*(rc<2.5))])
    res['o_o_peak']=float(rc[np.argmax(gOO*(rc<3.5))])
    res['si_o_rmin']=first_min(rc,gSi,1.8,2.6); res['al_o_rmin']=first_min(rc,gAl,1.9,2.8)
    res['er_o_first_min_note']='first minimum of 5-bin-smoothed pooled Er-O g(r) searched in 2.5-3.6 A'
    # coordination of Si, Al
    def cn(fr,tc,rc_):
        cc=np.where(fr['types']==tc)[0]; oo=np.where(fr['types']==1)[0]
        _,r=dists(fr['pos'],cc,oo,fr['L']); return (r<rc_).sum(1)
    siCN=np.concatenate([cn(f,2,res['si_o_rmin']) for f in allfr]); alCN=np.concatenate([cn(f,3,res['al_o_rmin']) for f in allfr])
    alCN21=np.concatenate([cn(f,3,2.1) for f in allfr])
    res['si_cn_hist']={int(k):float((siCN==k).mean()) for k in np.unique(siCN)}
    res['al_cn_hist']={int(k):float((alCN==k).mean()) for k in np.unique(alCN)}
    res['al_cn_hist_2.1A']={int(k):float((alCN21==k).mean()) for k in np.unique(alCN21)}
    # angles
    A={k:[] for k in ['O-Si-O','O-Al-O','Si-O-Si','Si-O-Al','Al-O-Al','O-Er-O']}
    oxy_stats=dict(BO=0,NBO=0,free=0,tri=0)
    sites=[]; erer=[]; CNs=[]; Rs=[]; link=dict(Si=0,Al=0); bulk_link=dict(Si=0,Al=0); nbo_lig=0; nlig=0; tagcount={}
    for rep in sorted(reps):
        for fi,fr in enumerate(reps[rep]):
            A['O-Si-O']+=angles(fr,2,1,2.0); A['O-Al-O']+=angles(fr,3,1,2.1)
            A['O-Er-O']+=angles(fr,4,1,rmin)
            t=fr['types']; L=fr['L']; pos=fr['pos']
            oo=np.where(t==1)[0]; ff=np.where((t==2)|(t==3))[0]
            d,r=dists(pos,oo,ff,L)
            cut=np.where(t[ff]==2,2.0,2.1)[None,:]
            bond=r<cut
            ofm={}
            for k,o in enumerate(oo):
                nb=np.where(bond[k])[0]
                ofm[o]=[('Si' if t[ff[j]]==2 else 'Al') for j in nb]
                for x in ofm[o]: bulk_link[x]+=1
                n=len(nb); oxy_stats['BO' if n==2 else 'NBO' if n==1 else 'free' if n==0 else 'tri']+=1
                if n==2:
                    v=d[k,nb]/r[k,nb][:,None]; ang=np.degrees(np.arccos(np.clip(v[0]@v[1],-1,1)))
                    key='-O-'.join(sorted(ofm[o],key=lambda s:{'Si':0,'Al':1}[s])); A[key].append(ang)
            ee=np.where(t==4)[0]
            de,re=dists(pos,ee,oo,L)
            _,rr=dists(pos,ee,ee,L); np.fill_diagonal(rr,np.inf); erer+=rr.min(1).tolist()
            for k,e in enumerate(ee):
                nb=np.where(re[k]<rmin)[0]; nb=nb[np.argsort(re[k,nb])]
                ligs=[['O']+[round(float(x),4) for x in de[k,j]] for j in nb]
                tags=[]
                for j in nb:
                    fm=ofm[oo[j]]; n=len(fm)
                    kind='BO' if n==2 else 'NBO' if n==1 else 'free_O' if n==0 else 'TBO'
                    tg=kind+('-'+'-'.join(sorted(fm)) if fm else '')
                    tags.append(tg); tagcount[tg]=tagcount.get(tg,0)+1
                    for s in fm: link[s]+=1
                    nlig+=1; nbo_lig+= (n==1)
                CNs.append(len(nb)); Rs+= re[k,nb].tolist()
                sites.append(dict(replica=rep,frame=fi,timestep=fr['step'],er_index=int(fr['ids'][e]),ligands=ligs,ligand_tags=tags,
                                  ligand_dist=[round(float(x),4) for x in re[k,nb]]))
    CNs=np.array(CNs); Rs=np.array(Rs); erer=np.array(erer)
    tot=sum(oxy_stats.values())
    summ=dict(CN_mean=float(CNs.mean()),CN_std=float(CNs.std()),CN_hist={int(k):int((CNs==k).sum()) for k in np.unique(CNs)},
        R_mean=float(Rs.mean()),R_std=float(Rs.std()),R_median=float(np.median(Rs)),r_min=rmin,er_o_rdf_peak=res['er_o_rdf_peak'],
        nbo_fraction_first_shell=nbo_lig/nlig, ligand_tag_counts=tagcount,
        er_o_former_linkages=dict(Er_O_Si=link['Si'],Er_O_Al=link['Al'],frac_Al=link['Al']/max(1,link['Si']+link['Al'])),
        bulk_O_species={k:v/tot for k,v in oxy_stats.items()},
        er_er_nn=dict(mean=float(erer.mean()),min=float(erer.min()),max=float(erer.max()),frac_below_4A=float((erer<4).mean()),frac_below_5A=float((erer<5).mean())),
        n_sites=len(sites),n_distinct_er=len(reps)*8,
        sanity=dict(si_o_peak=res['si_o_peak'],al_o_peak=res['al_o_peak'],o_o_peak=res['o_o_peak'],si_o_rmin=res['si_o_rmin'],al_o_rmin=res['al_o_rmin'],
            si_cn_hist=res['si_cn_hist'],al_cn_hist=res['al_cn_hist'],al_cn_hist_2p1A=res['al_cn_hist_2.1A'],
            angles={k:hist_stats(v) for k,v in A.items()}))
    # composition of bulk Al fraction for reference: random expectation of Er-O-Al linkage
    summ['er_o_former_linkages']['bulk_frac_Al_of_all_O-former_bonds']=bulk_link['Al']/(bulk_link['Al']+bulk_link['Si'])
    np.savetxt(out_prefix+'_rdf.txt',np.c_[rc,gEr,gSi,gAl,gOO],fmt='%.4f',header='r gErO gSiO gAlO gOO (pooled final-20ps frames)')
    return sites,summ
if __name__=='__main__':
    import pickle
    sites,summ=analyze(sys.argv[2:],sys.argv[1])
    json.dump(dict(er_sites=sites,summary=summ),open(sys.argv[1]+'_analysis.json','w'))
    s=dict(summ); s['sanity']=dict(s['sanity']); s['sanity']['angles']={k:(v and {kk:vv for kk,vv in v.items() if kk!='hist_5deg'}) for k,v in s['sanity']['angles'].items()}
    print(json.dumps(s,indent=1))
