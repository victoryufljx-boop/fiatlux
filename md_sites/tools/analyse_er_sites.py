"""Extract Er first-shell sites + structural sanity checks from finished replicas.

usage: python analyse_er_sites.py WORKDIR OUT_JSON_PREFIX rep1 rep2 ...
Each replica dir must contain final20ps.dump (LAMMPS custom dump: id type q x y z) and run.log.
"""
import sys, json, os, re
import numpy as np

work, outpref, reps = sys.argv[1], sys.argv[2], sys.argv[3:]
EL = {1: "O", 2: "Si", 3: "Er"}
SIO_CUT = 2.0  # Si-O bond cutoff for BO/NBO tagging and angles


def read_dump(fn):
    frames = []
    with open(fn) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if line.startswith("ITEM: TIMESTEP"):
                step = int(f.readline())
                f.readline(); n = int(f.readline())
                f.readline()
                lo_hi = [list(map(float, f.readline().split()[:2])) for _ in range(3)]
                f.readline()
                a = np.array([f.readline().split() for _ in range(n)], float)
                L = np.array([h - l for l, h in lo_hi])
                frames.append((step, a[:, 0].astype(int), a[:, 1].astype(int), a[:, 3:6] - np.array([l for l, h in lo_hi]), L))
    return frames


def mic(d, L):
    return d - L * np.round(d / L)


def dists(A, B, L):
    d = mic(B[None, :, :] - A[:, None, :], L)
    return d, np.sqrt((d ** 2).sum(-1))


def angles(center, nbr_lists, X, L):
    out = []
    for c, nb in zip(center, nbr_lists):
        if len(nb) < 2:
            continue
        v = mic(X[nb] - X[c], L)
        v /= np.linalg.norm(v, axis=1)[:, None]
        for i in range(len(nb)):
            for j in range(i + 1, len(nb)):
                out.append(np.degrees(np.arccos(np.clip(v[i] @ v[j], -1, 1))))
    return out


bins = np.linspace(0, 8, 401)
rdf_acc = {k: np.zeros(400) for k in ["Si-O", "O-O", "Si-Si", "Er-O", "Er-Si", "Er-Er"]}
nfr_tot = 0
OSiO, SiOSi, SiCN = [], [], []
frames_by_rep = {}
for r in reps:
    fr = read_dump(f"{work}/{r}/final20ps.dump")
    fr = fr[-10:]
    frames_by_rep[r] = fr
    for step, ids, typ, X, L in fr:
        nfr_tot += 1
        V = np.prod(L)
        idx = {e: np.where(typ == t)[0] for t, e in EL.items()}
        for k in rdf_acc:
            a, b = k.split("-")
            _, R = dists(X[idx[a]], X[idx[b]], L)
            if a == b:
                R = R[np.triu_indices(len(idx[a]), 1)]
                norm = len(idx[a]) * (len(idx[a]) - 1) / 2
            else:
                norm = len(idx[a]) * len(idx[b])
            h, _ = np.histogram(R.ravel(), bins)
            shell = 4 / 3 * np.pi * (bins[1:] ** 3 - bins[:-1] ** 3)
            rdf_acc[k] += h / (norm / V * shell)
        # Si-O network
        Si, O = idx["Si"], idx["O"]
        _, R = dists(X[Si], X[O], L)
        bond = R < SIO_CUT
        SiCN += list(bond.sum(1))
        OSiO += angles(Si, [O[np.where(b)[0]] for b in bond], X, L)
        SiOSi += angles(O, [Si[np.where(b)[0]] for b in bond.T], X, L)
rdf = {k: (v / nfr_tot) for k, v in rdf_acc.items()}
rc = 0.5 * (bins[1:] + bins[:-1])


def first_peak_min(g, lo=1.0, hi=4.0):
    m = (rc > lo) & (rc < hi)
    gs = np.convolve(g, np.ones(5) / 5, mode="same")
    ip = np.argmax(np.where(m, gs, -1))
    j = ip
    while j < len(gs) - 3 and not (gs[j + 1] >= gs[j] and gs[j + 2] >= gs[j]) :
        j += 1
    # take min within region after peak up to hi
    mm = (rc > rc[ip]) & (rc < hi)
    jmin = np.argmin(np.where(mm, gs, 1e9))
    return rc[ip], rc[jmin]

p_sio, m_sio = first_peak_min(rdf["Si-O"], 1.2, 2.6)
p_ero, r_min = first_peak_min(rdf["Er-O"], 1.8, 3.6)
r_min = float(os.environ.get("RMIN", round(r_min, 3)))

sites, CN, Rall, nbo_n, o_n, ernn, tagcount = [], [], [], 0, 0, [], {}
for r in reps:
    for fi, (step, ids, typ, X, L) in enumerate(frames_by_rep[r]):
        idx = {e: np.where(typ == t)[0] for t, e in EL.items()}
        Er, Si, O = idx["Er"], idx["Si"], idx["O"]
        _, RSiO = dists(X[Si], X[O], L)
        formers_of_O = [Si[np.where(RSiO[:, k] < SIO_CUT)[0]] for k in range(len(O))]
        Opos = {o: k for k, o in enumerate(O)}
        _, REE = dists(X[Er], X[Er], L)
        np.fill_diagonal(REE, np.inf)
        for ie, e in enumerate(Er):
            d = mic(X - X[e], L)
            R = np.sqrt((d ** 2).sum(1))
            R[e] = np.inf
            nb = np.where(R < r_min)[0]
            nb = nb[np.argsort(R[nb])]
            lig, tags = [], []
            cn = 0
            for j in nb:
                elj = EL[typ[j]]
                lig.append([elj] + [round(float(x), 4) for x in d[j]])
                if elj == "O":
                    cn += 1
                    Rall.append(R[j])
                    f = formers_of_O[Opos[j]]
                    nf = len(f)
                    t = {0: "free_O", 1: "NBO", 2: "BO"}.get(nf, "O_3c")
                    tagcount[t] = tagcount.get(t, 0) + 1
                    o_n += 1
                    nbo_n += (nf < 2)
                    tags.append({"tag": t, "formers": [f"Si{int(ids[s])}" for s in f], "atom_id": int(ids[j])})
                else:
                    tags.append({"tag": "cation", "atom_id": int(ids[j])})
            CN.append(cn)
            ernn.append(float(REE[ie].min()))
            sites.append({"replica": r, "frame": fi, "timestep": int(step), "er_index": int(ids[e]),
                          "ligands": lig, "ligand_tags": tags})

CN = np.array(CN); Rall = np.array(Rall)
hist = {int(k): int((CN == k).sum()) for k in range(CN.min(), CN.max() + 1)}
OSiO = np.array(OSiO); SiOSi = np.array(SiOSi); SiCN = np.array(SiCN)


def peak(a, lo, hi, w=1.0):
    h, e = np.histogram(a, np.arange(lo, hi + w, w))
    return float(0.5 * (e[np.argmax(h)] + e[np.argmax(h) + 1]))

summary = {
    "n_sites": len(sites), "n_distinct_er": len(sites) // 10 if sites else 0,
    "r_min": r_min, "ErO_rdf_peak": float(round(p_ero, 3)),
    "CN_mean": float(CN.mean()), "CN_std": float(CN.std()), "CN_hist": hist,
    "R_mean": float(Rall.mean()), "R_std": float(Rall.std()),
    "nbo_fraction": nbo_n / o_n, "first_shell_O_tag_counts": tagcount,
    "er_er_nn": {"mean": float(np.mean(ernn)), "min": float(np.min(ernn)),
                  "per_er_last_frame": {f"{st['replica']}:Er{st['er_index']}": round(x, 3) for st, x in zip(sites, ernn) if st["frame"] == 9},
                  "fraction_below_5A": float(np.mean(np.array(ernn) < 5.0))},
    "sanity": {
        "SiO_rdf_peak": float(round(p_sio, 3)), "SiO_rdf_min": float(round(m_sio, 3)),
        "Si_CN_mean": float(SiCN.mean()), "Si_CN4_fraction": float(np.mean(SiCN == 4)),
        "OSiO_mean": float(OSiO.mean()), "OSiO_std": float(OSiO.std()), "OSiO_mode": peak(OSiO, 60, 180),
        "SiOSi_mean": float(SiOSi.mean()), "SiOSi_std": float(SiOSi.std()), "SiOSi_mode": peak(SiOSi, 90, 180, 2.0),
        "ErO_peak_target": "EXAFS Er-doped silica ~2.25 A, CN~6",
    },
    "angle_hist": {
        "O-Si-O": dict(zip(["edges", "counts"], [list(map(float, np.arange(60, 181, 2))), np.histogram(OSiO, np.arange(60, 181, 2))[0].tolist()])),
        "Si-O-Si": dict(zip(["edges", "counts"], [list(map(float, np.arange(90, 181, 2))), np.histogram(SiOSi, np.arange(90, 181, 2))[0].tolist()])),
    },
    "rdf": {"r": [round(x, 3) for x in rc[::2]],
            **{k: [round(float(x), 3) for x in v[::2]] for k, v in rdf.items()}},
}
# CN sensitivity to the shell cutoff (Er-O RDF has no deep first minimum)
sens = {}
for rcut in (2.6, 2.8, 3.0, 3.2):
    cns = []
    for r in reps:
        for step, ids, typ, X, L in frames_by_rep[r]:
            _, R = dists(X[typ == 3], X[typ == 1], L)
            cns += list((R < rcut).sum(1))
    sens[str(rcut)] = {"CN_mean": float(np.mean(cns)), "CN_std": float(np.std(cns))}
summary["CN_vs_cutoff"] = sens
# O-Er-O angles in the first shell (angular distortion)
oeo = []
for s in sites:
    v = np.array([l[1:] for l in s["ligands"] if l[0] == "O"])
    if len(v) > 1:
        v = v / np.linalg.norm(v, axis=1)[:, None]
        c = v @ v.T
        oeo += list(np.degrees(np.arccos(np.clip(c[np.triu_indices(len(v), 1)], -1, 1))))
summary["O-Er-O"] = {"edges": list(map(float, np.arange(40, 181, 5))),
                     "counts": np.histogram(oeo, np.arange(40, 181, 5))[0].tolist(),
                     "mean": float(np.mean(oeo))}
json.dump({"sites": sites, "summary": summary}, open(outpref + ".analysis.json", "w"))
print(json.dumps({k: v for k, v in summary.items() if k not in ("rdf", "angle_hist", "O-Er-O")}, indent=1))
