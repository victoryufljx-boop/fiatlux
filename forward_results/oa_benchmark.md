# External benchmark: open-access literature

`scripts/oa_benchmark.py`, data in `gainmat/refs/oa_benchmark_refs.py` (never used to build the model).

| quantity | unit | n | model (own Omega) median abs err | model + paper's Omega | lookup (other papers) | model beats lookup |
|---|---|---|---|---|---|---|
| abs_peak_nm | nm | 6 | 1.11 (n=6) | 0.764 (n=1) | 1.5 (n=4) | 4/4 |
| em_fwhm_nm | nm | 6 | 20.1 (n=6) | 15 (n=5) | 5.3 (n=2) | 0/2 |
| em_peak_nm | nm | 14 | 19.5 (n=14) | 19.3 (n=9) | 2 (n=13) | 2/13 |
| omegas | dex | 14 | 0.0866 (n=14) | — | 0.0816 (n=14) | 6/14 |
| sigma_abs_cm2 | dex | 6 | 0.11 (n=6) | — | 0.0698 (n=2) | 1/2 |
| sigma_em_cm2 | dex | 17 | 0.185 (n=17) | 0.137 (n=10) | 0.0944 (n=13) | 3/13 |
| tau_ms | dex | 2 | 0.0669 (n=2) | — | — | 0/0 |
| tau_rad_ms | dex | 16 | 0.118 (n=16) | 0.0153 (n=13) | 0.161 (n=15) | 11/15 |
