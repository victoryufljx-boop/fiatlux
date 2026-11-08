"""Entry point: run_center(name, level) -> schema.DopantProps.

  run_center("Ce:LiYF4", "fast")   # minutes on 2 CPU cores (ECP Ce, SVP F, 12 A PC sphere)
  run_center("Ti:Al2O3", "full")   # Colab: def2-TZVP / ANO-RCC + sfX2C-1e, 20 A PC sphere

Threads: PySCF uses OpenMP; set OMP_NUM_THREADS before import (Colab ~12).
gpu4pyscf is NOT used: its CASSCF/NEVPT2 support is limited and QM/MM point
charges + ECP-only caps + ROHF are not all supported on GPU; CPU PySCF is used.
"""
from __future__ import annotations

import json
import numpy as np

from ..schema import DopantProps, Level, cm1_to_nm
from ..refs import cluster_refs as R

CENTERS = {"Ce:LiYF4": "ce", "Ce3+:LiYF4": "ce", "Ti:Al2O3": "ti", "Ti3+:Al2O3": "ti"}


def run_center(name: str, level: str = "fast", verbose: int = 0, log=print,
               raw_out: dict | None = None, **flags) -> DopantProps:
    """Never raises on a calculation failure: returns DopantProps with method['error'].
    ``flags`` are passed to ce_5d.run / ti_3d.run (cap_mode, r_cap, relax, semicore,
    basis, relax_grid ...); None values keep the level default."""
    try:
        return _run_center(name, level, verbose, log, raw_out, **flags)
    except Exception as e:  # noqa: BLE001
        import traceback
        log(f"[cluster] {name} {level} FAILED: {e!r}")
        return DopantProps(method=dict(tier="cluster", level=level, error=repr(e),
                                       traceback=traceback.format_exc()))


def _run_center(name, level, verbose, log, raw_out, **flags):
    kind = CENTERS[name]
    if kind == "ti":
        from . import ti_3d
        o = ti_3d.run(level, verbose=verbose, log=log,
                      **{k: v for k, v in flags.items() if k in ("relax", "relax_grid", "basis")})
        e = np.array(o["e_pt2_so_cm1"])
        abs_cm1 = [o["abs_eg_mean_cm1"]]
        dp = DopantProps(
            levels=[Level(label=f"3d-{i}", energy_cm1=float(x), J=0.5) for i, x in enumerate(e)],
            abs_peaks_nm=[cm1_to_nm(x) for x in abs_cm1],
            em_peak_nm=cm1_to_nm(o["em_cm1"]) if o.get("em_cm1") else None,
            stokes_cm1=o.get("stokes_cm1"),
            method=dict(tier="cluster", level=level, relax=o.get("relax"),
                        scheme="embedded [TiO6]9- ROHF -> SA-CASSCF(1,5) -> SC-NEVPT2; "
                               "vertical 2T2g->2Eg (mean of eg pair); SO zeta_3d=154 cm-1",
                        basis=o["basis"], x2c=o["x2c"], nao=o["nao"],
                        casscf_cm1=o["e_cas_cm1"], nevpt2_cm1=o["e_pt2_cm1"],
                        t2g_split_cm1=o["t2g_split_cm1"], timings=o["timings"],
                        cluster=o["cluster"],
                        not_computed="2Eg Jahn-Teller relaxation -> emission/Stokes shift, "
                                     "JT shoulder (~560 nm)",
                        exp_ref=R.TI_AL2O3))
    else:
        from . import ce_5d
        o = ce_5d.run(level, verbose=verbose, log=log, **flags)
        r = o.get("e_pt2", o["e_cas"])
        em = r.get("em_relaxed_to_2F52_cm1", r["em_to_2F52_cm1"])
        levels = [Level(label=f"4f-{i}", energy_cm1=float(x), J=0.5) for i, x in enumerate(r["so_4f_cm1"])]
        levels += [Level(label=f"5d-{i}", energy_cm1=float(x), J=0.5) for i, x in enumerate(r["so_5d_cm1"])]
        dp = DopantProps(
            levels=levels,
            abs_peaks_nm=[cm1_to_nm(x) for x in r["so_5d_cm1"]],
            em_peak_nm=cm1_to_nm(em),
            stokes_cm1=r.get("stokes_cm1"),
            method=dict(tier="cluster", level=level, flags=o.get("flags"), relax=o.get("relax"),
                        semicore=o.get("semicore"),
                        em_relaxed_nm=([cm1_to_nm(r["em_relaxed_to_2F52_cm1"]),
                                        cm1_to_nm(r["em_relaxed_to_2F72_cm1"])]
                                       if "em_relaxed_to_2F52_cm1" in r else None),
                        scheme="embedded [CeF8]5- ROHF -> SA-CASSCF(1,12) 4f+5d -> SC-NEVPT2; "
                               "12x12 ab initio LF matrix + SO (zeta_4f=615, zeta_5d=1082 cm-1)",
                        basis=o["basis"], x2c=o["x2c"], nao=o["nao"],
                        casscf=o["e_cas"], nevpt2=o.get("e_pt2"),
                        em_vertical_nm=[cm1_to_nm(r["em_to_2F52_cm1"]),
                                        cm1_to_nm(r["em_to_2F72_cm1"]) if r["em_to_2F72_cm1"] else None],
                        em_note="vertical at ground-state geometry (no 5d relaxation): "
                                "upper bound on emission energy, Stokes shift not computed",
                        timings=o["timings"], cluster=o["cluster"], exp_ref=R.CE_LIYF4))
    if raw_out is not None:
        raw_out.update(o)
    return dp


def compare_table(name: str, dp: DopantProps) -> list:
    """Rows (quantity, computed, experiment, unit) for printing."""
    rows = []
    if not dp.abs_peaks_nm:
        return rows
    if "Ce" in name:
        exp = R.CE_LIYF4["abs_5d_cm1"]
        for i, (lam, ex) in enumerate(zip(dp.abs_peaks_nm, exp)):
            rows.append((f"4f->5d{i+1}", 1e7 / lam, ex, "cm-1"))
        lab = "5d1->2F5/2 (relaxed)" if dp.stokes_cm1 else "5d1->2F5/2 (vertical)"
        rows.append((lab, 1e7 / dp.em_peak_nm, 1e7 / R.CE_LIYF4["em_nm"][0], "cm-1"))
        rows.append(("Stokes shift", dp.stokes_cm1, R.CE_LIYF4["stokes_cm1"], "cm-1"))
    else:
        rows.append(("2T2g->2Eg", 1e7 / dp.abs_peaks_nm[0], 1e7 / R.TI_AL2O3["abs_nm"][0], "cm-1"))
        rows.append(("emission (2Eg JT min)", 1e7 / dp.em_peak_nm if dp.em_peak_nm else None,
                     1e7 / R.TI_AL2O3["em_nm"], "cm-1"))
        rows.append(("Stokes shift", dp.stokes_cm1, R.TI_AL2O3["stokes_cm1"], "cm-1"))
    return rows


def print_table(name, dp):
    if not dp.abs_peaks_nm:
        print(f"\n{name}: no result ({dp.method.get('error')})"); return
    print(f"\n{name}  ({dp.method.get('level')})")
    print(f"{'quantity':26s} {'calc cm-1':>10s} {'calc nm':>8s} {'exp cm-1':>9s} {'exp nm':>7s} {'diff':>7s}")
    for q, c, e, _ in compare_table(name, dp):
        if q == "Stokes shift":
            print(f"{q:26s} {(f'{c:10.0f}' if c else f'{chr(45):>10s}')} {'':8s} {e:9.0f}"); continue
        cs = f"{c:10.0f} {1e7/c:8.1f}" if c else f"{'-':>10s} {'-':>8s}"
        ds = f"{c-e:7.0f}" if c else f"{'-':>7s}"
        print(f"{q:26s} {cs} {e:9.0f} {1e7/e:7.1f} {ds}")


if __name__ == "__main__":
    import sys
    from dataclasses import asdict
    lvl = sys.argv[1] if len(sys.argv) > 1 else "fast"
    out = {}
    for n in ("Ti:Al2O3", "Ce:LiYF4"):
        dp = run_center(n, lvl)
        print_table(n, dp)
        out[n] = asdict(dp)
    print(json.dumps(out, default=float)[:2000])
