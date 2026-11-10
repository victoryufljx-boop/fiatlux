"""Assemble md_sites/Er-SiO2.json from analysed replicas.
usage: python finalize_er_sio2.py WORKDIR rep1 rep2 ...
"""
import sys, json, os, re, subprocess
import numpy as np

work, reps = sys.argv[1], sys.argv[2:]
here = os.path.dirname(os.path.abspath(__file__))
outdir = os.path.dirname(here)
subprocess.run([sys.executable, f"{here}/analyse_er_sites.py", work, f"{work}/all", *reps], check=True)
A = json.load(open(f"{work}/all.analysis.json"))

timing = {}
for r in reps:
    log = open(f"{work}/{r}/run.log").read()
    if os.path.exists(f"{work}/{r}/run_resume.log"):  # resumed after a container restart
        log += open(f"{work}/{r}/run_resume.log").read()
    loops = [float(x) for x in re.findall(r"Loop time of ([0-9.]+)", log)]
    timing[r] = {"resumed_after_container_restart": os.path.exists(f"{work}/{r}/run_resume.log"), "wall_s_total": round(sum(loops), 1), "stages_wall_s": [round(x, 1) for x in loops],
                 "seed": int(re.search(r"velocity all create 6000 (\d+)", log).group(1))}
    # final temperature of production
    th = re.findall(r"^\s+(\d+)\s+([0-9.]+)\s+([0-9.]+)\s+(-?[0-9.]+)", log, re.M)
    timing[r]["T_last_thermo"] = float(th[-1][2])
    prod = [float(t[2]) for t in th if int(t[0]) > 440000 - 50000]
    timing[r]["T_mean_300K_stage"] = round(float(np.mean(prod)), 1) if prod else None
    # extxyz of final structure
    try:
        from ase.io import read, write
        at = read(f"{work}/{r}/final.data", format="lammps-data", atom_style="charge",
                  Z_of_type={1: 8, 2: 14, 3: 68})
        fn = f"{outdir}/Er-SiO2_{r}_final.extxyz"
        write(fn, at, format="extxyz")
    except Exception as e:
        print("extxyz failed", r, e)

src = "A. Pedone, G. Malavasi, M. C. Menziani, A. N. Cormack, A. Segre, J. Phys. Chem. B 110, 11780 (2006), Table 2 (p. 11783) 'Potential Parameters of Eq 4 Derived from Binary Oxides'; read from the PDF hosted at https://www.lptmc.jussieu.fr/user/mmi/IMI_lecture/lecture_7/Pedone_JPCB_2006.pdf"
out = {
    "glass": "Er-SiO2",
    "composition": {"formula": "1000 SiO2 + 4 Er2O3", "mol_pct_Er2O3": round(400 / 1004, 3),
                    "atoms": {"Si": 1000, "O": 2012, "Er": 8}, "n_atoms": 3020},
    "potential": {
        "name": "Pedone et al. 2006 partial-charge rigid-ion (Coulomb + Morse + C/r^12)",
        "form": "U(r) = qi qj e^2/(4 pi eps0 r) + D[(1 - exp(-a(r - r0)))^2 - 1] + C/r^12",
        "charges": {"O": -1.2, "Si": 2.4, "Er": 1.8, "source": src},
        "pairs": {
            "O-O": {"D_eV": 0.042395, "a_invA": 1.379316, "r0_A": 3.618701, "C_eVA12": 22.0, "source": src},
            "Si-O": {"D_eV": 0.340554, "a_invA": 2.006700, "r0_A": 2.100000, "C_eVA12": 1.0, "source": src},
            "Er-O": {"D_eV": 0.040448, "a_invA": 2.294078, "r0_A": 2.837722, "C_eVA12": 3.0, "source": src},
            "Al-O (not used here)": {"D_eV": 0.361581, "a_invA": 1.900442, "r0_A": 2.164818, "C_eVA12": 0.9, "source": src},
            "Ge-O (not used here)": {"D_eV": 0.158118, "a_invA": 2.294230, "r0_A": 2.261313, "C_eVA12": 5.0, "source": src},
            "cation-cation": "Coulomb only (as in Pedone 2006)",
        },
        "notes": ["Table 2 column header for a_ij prints A^-2 in the extracted text; a appears inside exp(-a(r-r0)) so the unit is 1/A (standard reading).",
                  "Parameters were taken from a text extraction of the PDF; values were not cross-checked against a second source in this run.",
                  "Coulomb: damped shifted force (DSF, Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006)), alpha = 0.2 1/A, cutoff 10 A, instead of Ewald/PPPM, for speed (3.5x on this 2-core machine). Short-range cutoff 10 A (paper: 15 A; Morse terms are < 2e-5 eV at 10 A).",
                  "Implemented as one LAMMPS pair_style table (linear, 19201 points, 0.4-10 A) per pair containing DSF Coulomb + Morse + C/r^12."],
    },
    "protocol": {
        "engine": "LAMMPS 22 Jul 2025 (PyPI wheel lammps 2025.7.22.4.0), 1 MPI rank per replica",
        "ensemble": "NVT, Nose-Hoover (Tdamp 0.1 ps)", "density_g_cm3": 2.20, "box_A": 35.9612,
        "timestep_fs": 1.0,
        "stages": ["random placement (min-distance check) + CG minimisation", "6000 K 20 ps", "3500 K 50 ps",
                   "linear cool 3500 -> 300 K at 10 K/ps (320 ps)", "300 K 50 ps (sites from last 20 ps, 10 frames every 2 ps)"],
        "quench_rate_K_per_ps": 10.0, "replicas": len(reps), "n_atoms": 3020,
        "replica_timing": timing,
    },
    "dsf_vs_ewald_check": {
        "what": "Replica rep3 final 300 K DSF structure: single-point comparison with PPPM/Ewald (1e-6), then 19 ps NVT 300 K re-equilibration with PPPM 1e-5 (10 frames over the last 9 ps) and the same site analysis",
        "single_point": {"pressure_bar_DSF": -32166, "pressure_bar_Ewald": -30622,
                         "rms_force_Ewald_eV_A": {"O": 1.345, "Si": 1.616, "Er": 1.336},
                         "rms_force_difference_eV_A": {"O": 0.247, "Si": 0.496, "Er": 0.356}},
        "Er_shell_at_fixed_cutoff_3.0A": {"DSF": {"CN_mean": 4.94, "R_mean": 2.341, "R_std": 0.232, "nbo_fraction": 0.608},
                                          "PPPM_relaxed": {"CN_mean": 4.99, "R_mean": 2.353, "R_std": 0.248, "nbo_fraction": 0.602}},
        "Er_shell_at_fixed_cutoff_2.8A": {"DSF": {"CN_mean": 4.64, "R_mean": 2.306}, "PPPM_relaxed": {"CN_mean": 4.61, "R_mean": 2.309}},
        "conclusion": "instantaneous DSF forces differ noticeably from Ewald, but the relaxed Er first shell is statistically unchanged",
    },
    "er_sites": A["sites"],
    "summary": A["summary"],
}
json.dump(out, open(f"{outdir}/Er-SiO2.json", "w"), separators=(",", ":"))
print("wrote", f"{outdir}/Er-SiO2.json", os.path.getsize(f"{outdir}/Er-SiO2.json") // 1024, "kB")
