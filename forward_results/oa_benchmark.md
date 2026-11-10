# External benchmark: open-access literature

`scripts/oa_benchmark.py`, data in `gainmat/refs/oa_benchmark_refs.py` (never used to build the model).

| quantity | unit | n | model (own Omega) median abs err | model + paper's Omega | lookup (other papers) | model beats lookup |
|---|---|---|---|---|---|---|
| abs_peak_nm | nm | 8 | 0.877 (n=8) | 0.523 (n=2) | 1.5 (n=4) | 4/4 |
| em_fwhm_nm | nm | 7 | 27.3 (n=7) | 15 (n=5) | 5.3 (n=2) | 0/2 |
| em_peak_nm | nm | 18 | 12.7 (n=18) | 19.3 (n=9) | 1.5 (n=18) | 3/18 |
| omegas | dex | 16 | 0.0866 (n=16) | — | 0.0918 (n=15) | 7/15 |
| sigma_abs_cm2 | dex | 7 | 0.0794 (n=7) | 0.00674 (n=1) | 0.0698 (n=2) | 1/2 |
| sigma_em_cm2 | dex | 21 | 0.185 (n=21) | 0.137 (n=10) | 0.124 (n=17) | 5/17 |
| tau_ms | dex | 2 | 0.0669 (n=2) | — | — | 0/0 |
| tau_rad_ms | dex | 20 | 0.112 (n=20) | 0.0181 (n=14) | 0.159 (n=20) | 14/20 |
