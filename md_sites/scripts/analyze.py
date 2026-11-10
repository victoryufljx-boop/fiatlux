"""Analyse final 300 K stage of all replicas -> md_sites JSON + stats.
usage: python analyze.py OUT.json rep1 rep2 ...
"""
import sys, json, os, numpy as np

TYPES = {1: 'O', 2: 'Si', 3: 'Ge', 4: 'Er'}
FORMER_CUT = {'Si': 2.0, 'Ge': 2.0, 'Al': 2.1}


def read_dump(path):
    frames = []
    with open(path) as f:
        lines = f.read().split('\n')
    i = 0
    while i < len(lines):
        if lines[i].startswith('ITEM: TIMESTEP'):
            step = int(lines[i + 1]); n = int(lines[i + 3])
            b = [list(map(float, lines[i + 5 + k].split()[:2])) for k in range(3)]
            L = np.array([hi - lo for lo, hi in b])
            arr = np.array([l.split() for l in lines[i + 9:i + 9 + n]], float)
            arr = arr[np.argsort(arr[:, 0])]
            frames.append((step, L, arr[:, 1].astype(int), arr[:, 2:5] - np.array([lo for lo, hi in b])))
            i += 9 + n
        else:
            i += 1
    return frames


def mic(d, L):
    return d - L * np.round(d / L)


def dists(P, x, L):
    d = mic(P - x, L)
    return d, np.sqrt((d ** 2).sum(-1))


def rdf(P, el, a, b, L, rmax=8.0, dr=0.02):
    ia = np.where(el == a)[0]; ib = np.where(el == b)[0]
    h = np.zeros(int(rmax / dr))
    for i in ia:
        _, r = dists(P[ib], P[i], L)
        r = r[(r > 1e-6) & (r < rmax)]
        h += np.histogram(r, bins=len(h), range=(0, rmax))[0]
    rr = (np.arange(len(h)) + 0.5) * dr
    nb = len(ib) - (1 if a == b else 0)
    rho = nb / np.prod(L)
    g = h / (len(ia) * rho * 4 * np.pi * rr ** 2 * dr)
    return rr, g


def first_peak_min(rr, g, rlo=1.2):
    gs = np.convolve(g, np.ones(5) / 5, mode='same')
    m = rr > rlo
    ip = np.argmax(np.where(m, gs, -1))
    seg = gs[ip:ip + int(1.2 / (rr[1] - rr[0]))]
    mn = seg.min(); tol = 0.001 * max(gs[ip], 1e-9)
    idx = np.where(seg <= mn + tol)[0]
    # first contiguous block of near-minimal values: take its midpoint
    blk = [idx[0]]
    for k in idx[1:]:
        if k == blk[-1] + 1: blk.append(k)
        else: break
    jm = ip + (blk[0] + blk[-1]) // 2
    return float(rr[ip]), float(rr[jm])


def angle(u, v):
    c = (u * v).sum(-1) / np.linalg.norm(u, axis=-1) / np.linalg.norm(v, axis=-1)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def stats(a):
    a = np.asarray(a, float)
    if len(a) == 0:
        return None
    h, e = np.histogram(a, bins=90, range=(0, 180))
    hs = np.convolve(h, np.ones(3) / 3, mode='same')
    return {'mean': round(float(a.mean()), 2), 'std': round(float(a.std()), 2),
            'mode': round(float(e[np.argmax(hs)] + 1.0), 1), 'n': int(len(a))}


def main():
    out = sys.argv[1]; reps = sys.argv[2:]
    allframes = {}
    for rp in reps:
        fr = read_dump(os.path.join(rp, 'final.dump'))
        allframes[rp] = fr[-10:]
    # pooled RDFs over all frames
    pairs = [('Er', 'O'), ('Si', 'O'), ('Ge', 'O'), ('O', 'O'), ('Si', 'Si'), ('Er', 'Er')]
    G = {p: None for p in pairs}
    for rp, fr in allframes.items():
        for step, L, t, P in fr:
            el = np.array([TYPES[k] for k in t])
            for p in pairs:
                rr, g = rdf(P, el, p[0], p[1], L)
                G[p] = g if G[p] is None else G[p] + g
    nfr = sum(len(f) for f in allframes.values())
    peaks = {}
    for p in pairs:
        G[p] /= nfr
        if p != ('Er', 'Er'):
            peaks['-'.join(p)] = first_peak_min(rr, G[p])
    r_min = peaks['Er-O'][1]
    ge_cut = peaks['Ge-O'][1]; si_cut = peaks['Si-O'][1]

    sites = []; cn_all = []; r_all = []; tags_all = []; link = {'Si': 0, 'Ge': 0}
    oero = []; ereenn = []; er_o_former = []
    cn_si = []; cn_ge = []; osio = []; sios = []; ogeo = []; tox = {'Si-O-Si': [], 'Si-O-Ge': [], 'Ge-O-Ge': []}
    oxy_spec = {'NBO': 0, 'BO': 0, 'free': 0, 'TBO': 0}
    for rp, fr in allframes.items():
        for fi, (step, L, t, P) in enumerate(fr):
            el = np.array([TYPES[k] for k in t])
            iO = np.where(el == 'O')[0]
            iF = np.where((el == 'Si') | (el == 'Ge'))[0]
            # O -> formers bonded
            obond = {}
            for o in iO:
                d, r = dists(P[iF], P[o], L)
                cut = np.array([FORMER_CUT[el[j]] for j in iF])
                obond[o] = [iF[k] for k in np.where(r < cut)[0]]
            for o in iO:
                nf = len(obond[o])
                oxy_spec[{0: 'free', 1: 'NBO', 2: 'BO'}.get(nf, 'TBO')] += 1
            # former geometry (first frame of each replica is enough: use all, cheap)
            for fct, lst, ang in (('Si', cn_si, osio), ('Ge', cn_ge, ogeo)):
                cut = si_cut if fct == 'Si' else ge_cut
                for j in np.where(el == fct)[0]:
                    d, r = dists(P[iO], P[j], L)
                    nb = np.where(r < cut)[0]
                    lst.append(len(nb))
                    v = d[nb]
                    for a in range(len(v)):
                        for b in range(a + 1, len(v)):
                            ang.append(angle(v[a], v[b]))
            for o in iO:
                bs = obond[o]
                if len(bs) == 2:
                    u = mic(P[bs[0]] - P[o], L); v = mic(P[bs[1]] - P[o], L)
                    k = '-O-'.join(sorted([el[bs[0]], el[bs[1]]], key=lambda s: s != 'Si'))
                    tox[k].append(angle(u, v))
            iEr = np.where(el == 'Er')[0]
            for e in iEr:
                d, r = dists(P, P[e], L)
                r[e] = 1e9
                nb = np.where(r < r_min)[0]
                nb = nb[np.argsort(r[nb])]
                ligs = []; tags = []; vO = []
                others = []
                for j in nb:
                    if el[j] != 'O':
                        others.append([str(el[j])] + [round(float(x), 4) for x in d[j]])
                        continue
                    ligs.append([str(el[j])] + [round(float(x), 4) for x in d[j]])
                    if True:
                        f = obond[j]
                        nf = len(f)
                        kind = {0: 'free', 1: 'NBO', 2: 'BO'}.get(nf, 'TBO')
                        fs = sorted(str(el[k]) for k in f)
                        tags.append({'type': kind, 'formers': fs, 'r': round(float(r[j]), 4)})
                        for k in f:
                            link[str(el[k])] += 1
                            er_o_former.append(angle(-d[j], mic(P[k] - P[j], L)))
                        vO.append(d[j]); r_all.append(r[j]); tags_all.append(kind)
                    else:
                        tags.append({'type': 'non-O', 'r': round(float(r[j]), 4)})
                for a in range(len(vO)):
                    for b in range(a + 1, len(vO)):
                        oero.append(angle(vO[a], vO[b]))
                cn_all.append(len(vO))
                rEr = r[iEr]; rEr = rEr[rEr < 1e8]
                ereenn.append(float(rEr.min()) if len(rEr) else None)
                sites.append({'replica': rp, 'frame': fi, 'timestep': int(step), 'er_index': int(e),
                              'ligands': ligs, 'ligand_tags': tags, 'other_neighbors_within_rmin': others})
    cn_all = np.array(cn_all); r_all = np.array(r_all)
    hist = {int(k): int((cn_all == k).sum()) for k in range(cn_all.min(), cn_all.max() + 1)}
    nO = sum(oxy_spec.values())
    summary = {
        'n_sites': len(sites), 'n_distinct_er': len(sites) // 10,
        'r_min_ErO': round(r_min, 3),
        'ErO_first_peak': round(peaks['Er-O'][0], 3),
        'CN_mean': round(float(cn_all.mean()), 3), 'CN_std': round(float(cn_all.std()), 3), 'CN_hist': hist,
        'R_mean': round(float(r_all.mean()), 4), 'R_std': round(float(r_all.std()), 4),
        'R_shortest_mean': round(float(np.nanmean([min([t['r'] for t in s['ligand_tags']] or [np.nan]) for s in sites])), 4),
        'nbo_fraction_first_shell': round(tags_all.count('NBO') / len(tags_all), 4),
        'first_shell_O_types': {k: tags_all.count(k) for k in ('BO', 'NBO', 'free', 'TBO')},
        'ErOGe_vs_ErOSi_linkages': {'Er-O-Si': link['Si'], 'Er-O-Ge': link['Ge'],
                                    'frac_Ge': round(link['Ge'] / max(1, link['Si'] + link['Ge']), 4),
                                    'glass_Ge_fraction_of_formers': 0.1},
        'er_er_nn': {'mean': round(float(np.mean(ereenn)), 3), 'min': round(float(np.min(ereenn)), 3),
                     'frac_below_5A': round(float(np.mean(np.array(ereenn) < 5.0)), 3),
                     'per_er_nn_last_frame': [round(x, 3) for x in ereenn[-8:]]},
        'angles_deg': {'O-Er-O': stats(oero), 'Er-O-T(former)': stats(er_o_former)},
        'sanity_checks': {
            'Si-O_first_peak': round(peaks['Si-O'][0], 3), 'Si-O_rmin': round(si_cut, 3),
            'Ge-O_first_peak': round(peaks['Ge-O'][0], 3), 'Ge-O_rmin': round(ge_cut, 3),
            'O-O_first_peak': round(peaks['O-O'][0], 3), 'Si-Si_first_peak': round(peaks['Si-Si'][0], 3),
            'Si_CN_mean': round(float(np.mean(cn_si)), 3), 'Si_CN4_frac': round(float(np.mean(np.array(cn_si) == 4)), 4),
            'Ge_CN_mean': round(float(np.mean(cn_ge)), 3), 'Ge_CN4_frac': round(float(np.mean(np.array(cn_ge) == 4)), 4),
            'O-Si-O': stats(osio), 'O-Ge-O': stats(ogeo),
            'Si-O-Si': stats(tox['Si-O-Si']), 'Si-O-Ge': stats(tox['Si-O-Ge']), 'Ge-O-Ge': stats(tox['Ge-O-Ge']),
            'oxygen_speciation_bulk': {k: round(v / nO, 5) for k, v in oxy_spec.items()},
        },
    }
    rdfs = {'-'.join(p): [round(float(x), 3) for x in G[p][::2]] for p in pairs}
    json.dump({'summary': summary, 'sites': sites, 'rdf_r': [round(float(x), 3) for x in rr[::2]], 'rdfs': rdfs},
              open(out, 'w'))
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
