# Explicit model vs explicit + ML layer

Residual = log10(model / measured).  ML layer = learned ZBLAN host microstructure (L2, adopted because it transfers across ions) + extrinsic quenching rule for silica-family 2 um levels (L3).  Rows whose measurement trained an ML component are shown separately.

| check | subset | n | changed by ML | rms dex: explicit | rms dex: +ML | within 30 %: explicit | +ML | within 2x: explicit | +ML |
|---|---|---|---|---|---|---|---|---|---|
| Stark splitting | independent | 9 | 0 | 0.177 | 0.177 | 67% | 67% | 78% | 78% |
| integrated cross section | independent | 14 | 0 | 0.223 | 0.223 | 29% | 29% | 79% | 79% |
| integrated cross section | same source as L2 widths | 3 | 2 | 0.083 | 0.085 | 100% | 67% | 100% | 100% |
| lifetime | independent | 24 | 11 | 0.389 | 0.387 | 62% | 67% | 88% | 88% |
| lifetime | in-sample (trained L3) | 4 | 2 | 0.610 | 0.400 | 25% | 75% | 25% | 75% |
| peak cross section | independent | 11 | 1 | 0.197 | 0.179 | 64% | 64% | 82% | 91% |
| peak cross section | same source as L2 widths | 7 | 7 | 0.140 | 0.092 | 43% | 57% | 100% | 100% |
| radiative lifetime | independent | 12 | 1 | 0.065 | 0.065 | 92% | 92% | 100% | 100% |

## Rows the ML layer changed

| check | ion | host | item | measured | x explicit | x +ML | subset |
|---|---|---|---|---|---|---|---|
| integrated cross section | Dy3+ | ZBLAN | Dy:ZBLAN 6H13/2 em (Jackson 2012) | 1.0 | 1.18 | 1.15 | same source as L2 widths |
| integrated cross section | Er3+ | ZBLAN | Er:ZBLAN 4I11/2 em (Jackson 2012) | 1.0 | 0.79 | 0.77 | same source as L2 widths |
| lifetime | Dy3+ | ZBLAN | ['6H13/2', '6H15/2'] | 0.65 | 6.03 | 4.03 | independent |
| lifetime | Er3+ | ZBLAN | ['4I11/2', '4I13/2'] | 6.8999999999999995 | 0.96 | 0.92 | independent |
| lifetime | Er3+ | ZBLAN | ['4F9/2', '4I9/2'] | 0.177 | 1.32 | 1.14 | independent |
| lifetime | Er3+ | ZBLAN | 4I11/2 | 6.71 | 0.99 | 0.94 | independent |
| lifetime | Er3+ | ZBLAN | 4S3/2 | 0.57 | 1.00 | 0.96 | independent |
| lifetime | Ho3+ | ZBLAN | ['5I6', '5I7'] | 3.5 | 1.28 | 1.22 | independent |
| lifetime | Ho3+ | ZBLAN | 5I6 | 3.5 | 1.28 | 1.22 | independent |
| lifetime | Ho3+ | ZBLAN | 5F5 | 0.29 | 0.03 | 0.02 | independent |
| lifetime | Ho3+ | ZBLAN | 5S2 | 0.3 | 0.98 | 0.93 | independent |
| lifetime | Ho3+ | silica | ['5I7', '5I8'] | 0.7 | 3.05 | 0.83 | in-sample (trained L3) |
| lifetime | Pr3+ | ZBLAN | ['1G4', '3H5'] | 0.11 | 1.95 | 1.31 | independent |
| lifetime | Pr3+ | ZBLAN | 1G4 | 0.11 | 1.95 | 1.31 | independent |
| lifetime | Tm3+ | aluminosilicate | ['3F4', '3H6'] | 0.65 | 6.15 | 1.02 | in-sample (trained L3) |
| peak cross section | Dy3+ | ZBLAN | 6H13/2->6H15/2 em | 2.95e-21 | 1.38 | 1.06 | same source as L2 widths |
| peak cross section | Er3+ | ZBLAN | 4I11/2->4I13/2 em | 5.6e-21 | 1.01 | 0.76 | same source as L2 widths |
| peak cross section | Er3+ | ZBLAN | 4I13/2->4I15/2 em | 4.2e-21 | 1.40 | 1.02 | same source as L2 widths |
| peak cross section | Ho3+ | ZBLAN | 5I6->5I7 em | 5.2e-21 | 1.83 | 1.22 | same source as L2 widths |
| peak cross section | Ho3+ | ZBLAN | 5I7->5I8 em | 5.6e-21 | 1.01 | 0.76 | same source as L2 widths |
| peak cross section | Nd3+ | ZBLAN | 4F3/2->4I11/2 em | 2.64e-20 | 1.45 | 0.87 | same source as L2 widths |
| peak cross section | Pr3+ | ZBLAN | 1G4->3H5 em | 3.48e-21 | 2.47 | 1.91 | independent |
| peak cross section | Tm3+ | ZBLAN | 3F4->3H6 em | 5.7e-21 | 0.94 | 0.73 | same source as L2 widths |
| radiative lifetime | Pr3+ | ZBLAN | 1G4 | 3.24 | 0.72 | 0.71 | independent |
