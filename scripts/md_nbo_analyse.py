"""Does a separate strength for non-bridging oxygen explain the Tm composition series?
Reads the grids from md_nbo_grid.py (bridging strength a, non-bridging b, both x handbook).

  Test 1 (leave one glass out): fit (a, b) on the Tm 3F4 peak of two glasses, predict the third.
  Test 2 (leave Tm out): fit (a, b) on Er + Yb Stark levels in pure silica, predict all three Tm peaks.
  Baseline: one O strength (a = b), same protocol.
Every (a, b) within loss 1 of the best (i.e. within ~1 sigma) is kept; predictions are given as
the range over that set, so non-unique fits show up as wide ranges.

    python scripts/md_nbo_analyse.py GRID_DIR OUT.json
"""
import sys, json, numpy as np
from scipy.interpolate import RegularGridInterpolator as RGI

D, OUT = sys.argv[1:3]
MEAS = {"Er-Ge10-SiO2": 1578.5, "Er-Al3-SiO2": 1632.5, "Er-Al9-SiO2": 1660.5}
SIG = 10.0
a_f = np.linspace(0.6, 2.6, 201); b_f = np.linspace(1.0, 4.6, 361)
AA, BB = np.meshgrid(a_f, b_f, indexing="ij")
peaks = {}
for g in MEAS:
    d = json.load(open(f"{D}/{g}.json")); A, B = d["A"], d["B"]
    Z = np.array([[next(r["tm_abs"] for r in d["rows"] if r["a"] == a and r["b"] == b) for b in B] for a in A])
    peaks[g] = RGI((A, B), Z)(np.stack([AA, BB], -1))
st = json.load(open(f"{D}/stark_SiO2.json")); A, B = st["A"], st["B"]


def grid(key, i):
    return RGI((A, B), np.array([[next(r[key][i] for r in st["rows"] if r["a"] == a and r["b"] == b) for b in B] for a in A]))(np.stack([AA, BB], -1))


er_low = [0, 30, 55, 129, 234]; yb = [0, 425, 530, 865]
L_stark = sum(((grid("er15", i) - er_low[i]) / 30) ** 2 for i in range(1, 5)) + ((grid("er15", 7) - 420) / 30) ** 2 \
    + ((grid("er13", 6) - 330) / 30) ** 2 + sum(((grid("yb7", i) - yb[i]) / 40) ** 2 for i in range(1, 4)) \
    + ((grid("yb5", 1) - 630) / 50) ** 2
diag = np.isclose(AA, BB, atol=0.006)            # one-O-strength baseline: a == b


def predict(L, target, mask=None):
    L = np.where(mask, L, np.inf) if mask is not None else L
    ok = L <= L.min() + 1.0
    i = np.unravel_index(np.argmin(L), L.shape)
    p = peaks[target]
    return dict(best_a=float(AA[i]), best_b=float(BB[i]), fit_loss=float(L[i]), pred_best=float(p[i]),
                pred_range=[float(p[ok].min()), float(p[ok].max())], measured=MEAS[target],
                err_best=float(p[i] - MEAS[target]))


res = {"leave_one_glass_out": {}, "baseline_one_O": {}, "from_Er_Yb_silica": {}, "stark_fit": {}}
for g in MEAS:
    L = sum(((peaks[h] - MEAS[h]) / SIG) ** 2 for h in MEAS if h != g)
    res["leave_one_glass_out"][g] = predict(L, g)
    res["baseline_one_O"][g] = predict(L, g, diag)
i = np.unravel_index(np.argmin(L_stark), L_stark.shape)
res["stark_fit"] = dict(a=float(AA[i]), b=float(BB[i]), loss=float(L_stark[i]),
                        loss_one_O_best=float(np.min(np.where(diag, L_stark, np.inf))),
                        er15=[float(grid("er15", k)[i]) for k in range(8)], yb7=[float(grid("yb7", k)[i]) for k in range(4)],
                        yb5_gap=float(grid("yb5", 1)[i]))
for g in MEAS:
    res["from_Er_Yb_silica"][g] = predict(L_stark, g)
# measured shifts vs model shifts at the Er/Yb-fitted point
res["shifts_from_Er_Yb_fit"] = {g: float(peaks[g][i] - peaks["Er-Al3-SiO2"][i]) for g in MEAS}
res["shifts_measured"] = {g: MEAS[g] - MEAS["Er-Al3-SiO2"] for g in MEAS}
json.dump(res, open(OUT, "w"), indent=1)
for k, v in res.items():
    print(k); print(json.dumps(v, indent=None, default=lambda x: round(x, 1))[:900])
