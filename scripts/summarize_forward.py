"""Markdown tables from results/forward/<run>/{validation,grid}.json."""
import json
import sys

run = sys.argv[1] if len(sys.argv) > 1 else "results/forward/final"
V = json.load(open(f"{run}/validation.json"))
G = json.load(open(f"{run}/grid.json"))


def f(x, fmt="{:.2g}", scale=1.0):
    if x is None:
        return "—"
    if isinstance(x, (list, tuple)):
        xs = [v for v in x if isinstance(v, (int, float))]
        return "–".join(fmt.format(v * scale) for v in xs) if xs else "—"
    if isinstance(x, (int, float)):
        return fmt.format(x * scale)
    return str(x)


print("### Validation against known fibre gain media\n")
print("| Medium | Measured band (nm) | Model window(s) (nm) | sigma_em model / exp (1e-21 cm2) | FWHM model / exp (nm) | tau model / exp (ms) | Omega source |")
print("|---|---|---|---|---|---|---|")
for r in V:
    if "error" in r or "skipped" in r:
        print(f"| {r.get('entry')} | | {r.get('error') or r.get('skipped')} | | | | |")
        continue
    band = r.get("band_exp_nm")
    win = r.get("window_model_nm")
    if win is None:
        win_s = "—"
    else:
        win_s = ", ".join(f"{a:.0f}–{b:.0f}" for a, b in win) or "none"
    se_m, se_e = r.get("model_sigma_em"), r.get("exp_sigma_em")
    sig = f"{f(se_m, '{:.1f}', 1e21)} / {f(se_e, '{:.1f}', 1e21)}"
    fw = f"{f(r.get('model_fwhm_em_nm'), '{:.0f}')} / {f(r.get('exp_fwhm_nm'), '{:.0f}')}"
    te = r.get("exp_tau_ms")
    tau = f"{f(r.get('model_tau_ms'), '{:.3g}')} / {f(te, '{:.3g}')}"
    band_s = f"{band[0]:.0f}–{band[1]:.0f}" if band else "—"
    print(f"| {r['entry']} | {band_s} | {win_s} | {sig} | {fw} | {tau} | {r.get('omega_src', '')} |")

print("\n### Predicted gain windows, all ions x glass families\n")
print("| Ion | " + " | ".join(sorted({g['family'] for g in G}, key=lambda s: [x['family'] for x in G].index(s))) + " |")
fams = []
for g in G:
    if g["family"] not in fams:
        fams.append(g["family"])
ions = []
for g in G:
    if g["ion"] not in ions:
        ions.append(g["ion"])
print("|---|" + "---|" * len(fams))
for ion in ions:
    cells = []
    for fam in fams:
        g = next((x for x in G if x["ion"] == ion and x["family"] == fam), None)
        if not g or "error" in g:
            cells.append("err")
            continue
        ws = [w for w in g.get("windows", []) if w["hi_nm"] - w["lo_nm"] >= 10]
        cells.append("<br>".join(f"{w['lo_nm']/1000:.2f}–{w['hi_nm']/1000:.2f}" for w in ws) or "—")
    print(f"| {ion.replace('3+', '')} | " + " | ".join(cells) + " |")
