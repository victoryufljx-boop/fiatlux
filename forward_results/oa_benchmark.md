# External benchmark: open-access literature

`scripts/oa_benchmark.py`, data in `gainmat/refs/oa_benchmark_refs.py` (never used to build the model).

| quantity | unit | n | model (own Omega) median abs err | model + paper's Omega | lookup (other papers) | model beats lookup |
|---|---|---|---|---|---|---|
| abs_peak_nm | nm | 46 | 2.53 (n=46) | 4.57 (n=14) | 2.48 (n=37) | 17/37 |
| em_fwhm_nm | nm | 19 | 28.7 (n=19) | 18.9 (n=10) | 23.6 (n=14) | 4/14 |
| em_peak_nm | nm | 68 | 12.3 (n=68) | 12.4 (n=29) | 13.2 (n=62) | 27/62 |
| n_refr | - | 18 | 0.073 (n=18) | 0.095 (n=10) | 0.151 (n=12) | 8/12 |
| omegas | dex | 38 | 0.127 (n=38) | — | 0.0735 (n=38) | 14/38 |
| sigma_abs_cm2 | dex | 19 | 0.141 (n=19) | 0.251 (n=7) | 0.217 (n=16) | 10/16 |
| sigma_em_cm2 | dex | 36 | 0.189 (n=36) | 0.15 (n=20) | 0.0941 (n=33) | 14/33 |
| tau_ms | dex | 70 | 0.322 (n=70) | 0.2 (n=18) | 0.164 (n=63) | 24/63 |
| tau_rad_ms | dex | 59 | 0.147 (n=59) | 0.0576 (n=47) | 0.157 (n=51) | 27/51 |

Hybrid = model band shifted/stretched to the emission peak and FWHM reported by OTHER papers on the same ion and family (library), integrated strength from the model.  Rows where both hybrid and lookup exist:

| quantity | n | model | hybrid | hybrid + paper's Omega | lookup | hybrid beats lookup |
|---|---|---|---|---|---|---|
| abs_peak_nm | 37 | 3.24 | 4.28 | 8.04 | 2.48 | 14/37 |
| em_fwhm_nm | 14 | 27.3 | 27.2 | 22.2 | 23.6 | 3/14 |
| em_peak_nm | 62 | 12.3 | 11.4 | 11.6 | 13.2 | 33/62 |
| sigma_abs_cm2 | 16 | 0.144 | 0.144 | 0.271 | 0.217 | 10/16 |
| sigma_em_cm2 | 33 | 0.187 | 0.187 | 0.129 | 0.0941 | 11/33 |

Without Tm (median abs error; paired = rows where a lookup exists):

| quantity | model (all) | model (paired) | lookup (paired) | model beats lookup |
|---|---|---|---|---|
| abs_peak_nm | 2.53 (n=28) | 5.37 (n=20) | 2.03 (n=20) | 8/20 |
| em_fwhm_nm | 28.7 (n=15) | 23.6 (n=10) | 13.9 (n=10) | 1/10 |
| em_peak_nm | 8.47 (n=43) | 7.17 (n=37) | 6.5 (n=37) | 16/37 |
| n_refr | 0.083 (n=14) | 0.131 (n=8) | 0.216 (n=8) | 6/8 |
| omegas | 0.127 (n=30) | 0.127 (n=30) | 0.0664 (n=30) | 12/30 |
| sigma_abs_cm2 | 0.141 (n=17) | 0.155 (n=14) | 0.246 (n=14) | 8/14 |
| sigma_em_cm2 | 0.173 (n=26) | 0.127 (n=24) | 0.102 (n=24) | 11/24 |
| tau_ms | 0.393 (n=52) | 0.392 (n=47) | 0.164 (n=47) | 16/47 |
| tau_rad_ms | 0.21 (n=47) | 0.168 (n=39) | 0.148 (n=39) | 19/39 |

Records: 99 hand-collected + 282 collected by extraction agents (rejected by checks: {'abstract only': 4}).
