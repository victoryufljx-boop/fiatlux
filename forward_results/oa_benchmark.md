# External benchmark: open-access literature

`scripts/oa_benchmark.py`, data in `gainmat/refs/oa_benchmark_refs.py` (never used to build the model).

| quantity | unit | n | model (own Omega) median abs err | model + paper's Omega | lookup (other papers) | model beats lookup | model + library Omega (paired) | beats lookup |
|---|---|---|---|---|---|---|---|---|
| abs_peak_nm | nm | 59 | 2.88 (n=59) | 4.57 (n=14) | 1.95 (n=50) | 17/50 | 3.74 (n=50) | 16/50 |
| em_fwhm_nm | nm | 19 | 28.7 (n=19) | 18.9 (n=10) | 23.6 (n=14) | 4/14 | 27.8 (n=14) | 4/14 |
| em_peak_nm | nm | 72 | 12 (n=72) | 12.4 (n=31) | 13.5 (n=65) | 30/65 | 12.2 (n=65) | 30/65 |
| n_refr | - | 19 | 0.07 (n=19) | 0.08 (n=11) | 0.131 (n=13) | 9/13 | 0.07 (n=13) | 9/13 |
| omegas | dex | 41 | 0.119 (n=41) | — | 0.0736 (n=41) | 16/41 | 0.0743 (n=41) | 3/41 |
| sigma_abs_cm2 | dex | 24 | 0.112 (n=24) | 0.251 (n=7) | 0.173 (n=21) | 10/21 | 0.115 (n=21) | 10/21 |
| sigma_em_cm2 | dex | 41 | 0.188 (n=41) | 0.15 (n=20) | 0.0717 (n=38) | 15/38 | 0.144 (n=38) | 15/38 |
| tau_ms | dex | 93 | 0.318 (n=93) | 0.2 (n=18) | 0.167 (n=87) | 30/87 | 0.266 (n=87) | 34/87 |
| tau_rad_ms | dex | 64 | 0.135 (n=64) | 0.047 (n=51) | 0.122 (n=54) | 29/54 | 0.0845 (n=54) | 33/54 |

model + library Omega = the model given the geometric-mean Omegas of the same ion and family from OTHER papers in this data set (physics + data library, blind to the record; falls back to the model's own Omega).

Hybrid = model band shifted/stretched to the emission peak and FWHM reported by OTHER papers on the same ion and family (library), integrated strength from the model.  Rows where both hybrid and lookup exist:

| quantity | n | model | hybrid | hybrid + paper's Omega | lookup | hybrid beats lookup |
|---|---|---|---|---|---|---|
| abs_peak_nm | 50 | 3.74 | 4.47 | 8.04 | 1.95 | 14/50 |
| em_fwhm_nm | 14 | 27.3 | 27.2 | 22.2 | 23.6 | 3/14 |
| em_peak_nm | 65 | 12.4 | 12.4 | 11.7 | 13.5 | 35/65 |
| sigma_abs_cm2 | 21 | 0.115 | 0.115 | 0.271 | 0.173 | 10/21 |
| sigma_em_cm2 | 38 | 0.186 | 0.184 | 0.129 | 0.0717 | 12/38 |

Without Tm (median abs error; paired = rows where a lookup exists):

| quantity | model (all) | model (paired) | lookup (paired) | model beats lookup |
|---|---|---|---|---|
| abs_peak_nm | 2.18 (n=33) | 4.25 (n=25) | 1.9 (n=25) | 10/25 |
| em_fwhm_nm | 28.7 (n=15) | 23.6 (n=10) | 13.9 (n=10) | 1/10 |
| em_peak_nm | 8.47 (n=47) | 9.68 (n=40) | 10.1 (n=40) | 19/40 |
| n_refr | 0.08 (n=15) | 0.11 (n=9) | 0.138 (n=9) | 7/9 |
| omegas | 0.119 (n=33) | 0.119 (n=33) | 0.0644 (n=33) | 14/33 |
| sigma_abs_cm2 | 0.112 (n=22) | 0.115 (n=19) | 0.145 (n=19) | 8/19 |
| sigma_em_cm2 | 0.16 (n=31) | 0.158 (n=29) | 0.078 (n=29) | 12/29 |
| tau_ms | 0.374 (n=58) | 0.364 (n=53) | 0.167 (n=53) | 19/53 |
| tau_rad_ms | 0.158 (n=52) | 0.147 (n=42) | 0.0992 (n=42) | 21/42 |

Records: 99 hand-collected + 341 collected by extraction agents (rejected by checks: {'abstract only': 4, 'duplicate paper': 10}).
