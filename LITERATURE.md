# Literature review: where gainmat sits, and what to upgrade

Compiled 2026-11-08 from web sources. Every claim cites a source in the [Sources](#sources) list as [Sn].
Where a source could not be opened (paywall, bot check), the claim is marked **(unverified)** or left out.
The numbers come from abstracts or full texts as they were fetched.

---

## A. Prior art: screening and designing luminescent and gain materials

### A1. What exists

| Area | Representative work | Scale and method | Outcome |
|---|---|---|---|
| Ce³⁺ 4f→5d by ΔSCF-DFT | Jia, Poncé, Miglio, Mikami, Gonze [S1] | 13 Ce³⁺ compounds. Constrained DFT/ΔSCF with relaxed ground and excited states (ABINIT) | Absorption/emission within 0.3 eV (two cases at 0.5 eV). Stokes shift MARE 19% vs 37% for the Dorenbos semi-empirical model |
| Eu²⁺ 5d by ΔSCF-DFT | Jia et al. [S2] | 15 Eu²⁺ phosphors, PBE+U (U = 7.5 eV) | Energies within 0.3 eV. Stokes shift error about 40%. Fails when 5d is not in the gap (CaO:Eu) |
| Beyond 1D configuration-coordinate model | Jia et al. [S3]; Bouquiaux et al. [S4] | ΔSCF plus multi-D phonon coupling (DFPT, generating function) | 1D-CCM overestimates the 4f–5d crossing barrier [S3]. Multi-D model fits the Eu²⁺ line shape with 6.5% error vs 12% for 1D, Huang–Rhys S = 3.05 [S4] |
| Ce scintillator screening | Canning et al., PRB 83, 125115 (2011) [S5] | PBE+U with constrained occupancy, about 100 hosts. Criteria: gap, Ce 4f–VBM offset, 5d localisation | Ba₂YCl₇:Ce predicted, then made and bright |
| Physics-informed ML for scintillators | Pilania et al., J. Chem. Phys. 148 (2018) [S6] | ML on Dorenbos-type 4f/5d binding-energy trends plus DFT vacuum-referred band edges | First-pass down-selection tool for perovskites and elpasolites |
| ML phosphor host (Brgoch) | Zhuo et al., Nat. Commun. 9 (2018) [S7] | SVM-predicted Debye temperature as a quantum-yield proxy, plus DFT gaps | NaBaB₉O₁₅:Eu²⁺, QY 95% |
| ML Ce³⁺ excitation (Brgoch) | Lee …, Scanlon, Brgoch (c. 2025) [S8] | XGBoost on 357 measured Ce³⁺ 5d₁ energies. Screened 54,885 sites | MAE 0.15 eV. Ca₂SrSc₆O₁₂:Ce³⁺ validated (0.12 eV error) |
| ML Eu²⁺ emission | Koyama et al. (2022) [S9]; Zhang et al. GPT+CGCNN (2023) [S10]; Li, Amachraa, Ong, Xie (NIR nitrides) [S11] | 129 / 264 / 223-compound sets | MAE 25 nm [S9]. MAE 0.20 eV [S10]. 3 of 19 synthesized candidates showed Eu²⁺ emission [S9]. (Sr,Ba)₃Li₄Si₂N₆:Eu NIR validated [S11] |
| DFT-guided phosphor host (Ong) | Wang, Ong et al. Sr₂LiAlO₄ [S12] | Structure prediction plus DFT screening | New quaternary phase. Eu²⁺/Ce³⁺ phosphor with CRI > 90 |
| Embedded-cluster wavefunction | Seijo & Barandiarán (AIMP + CASSCF/CASPT2/RASSI-SO) [S13–S16]; Ning, Duan, Seijo [S17]; Oleynichenko et al. (FS-RCC) [S18] | One system per paper; 21 garnet hosts in [S14] | YAG:Ce absorption overestimated by 2300–3300 cm⁻¹, emission by about 4600 cm⁻¹, Stokes shift strongly underestimated [S15]. FS-RCC gives Ce³⁺ 4f→5d within about 0.3 eV and τ within ×2 [S18] |
| ML Judd–Ofelt Ω (glasses) | Ahmadi et al., Sci. Rep. 14, 15505 (2024), GEP [S19]; ANN for Dy³⁺ lead borosilicate (MDPI Electronics 2022) [S20]; RF for Dy³⁺ glass (2025) [S21]; statistical learning for Er³⁺ tellurite [S22] | 60–70 samples, one glass family each | R² of 0.93–0.95 on the test set, valid only inside the training composition range [S19] |
| JO database and tools | LOMS.cz [S23]; JOYSpectra [S24] | LOMS has more than 1,200 JO records for 12 Ln³⁺ ions in more than 550 hosts (literature-curated, not independently verified) | Useful as an ML label source. LOMS has no ab initio part |
| Glass property ML | GlassNet (Cassar, 2023) [S25]; SciGlass predictor [S26] | MLP on 218k compositions, 85 properties incl. n_D and Abbe number | Composition→n is a solved problem for oxide glasses at this level |
| Laser-glass design (no ML) | Lun et al., Materials Futures 2025, "neighboring glassy compounds" [S27] | Statistical-ensemble model plus MD for germanate laser glasses | Multi-property luminescence charts. The authors say current ML models are inadequate for this problem |
| High-throughput optical data | Petousis et al., Sci. Data 4, 160134 (2017): DFPT ε∞/n for 1,056 compounds in Materials Project [S28]; Kim et al., Sci. Data 7, 387 (2020): HSE gaps for 10,481 materials [S29]; JARVIS TB-mBJ, 7,358 gaps [S30]; Trinquet et al.: about 2,200 DFPT SHG tensors via active learning [S31] | — | Ready-made host-tier labels and benchmarks |
| Multi-fidelity ML | Chen, Ong et al., multi-fidelity MEGNet [S32] | PBE (52k) + GLLB-SC (2.3k) + HSE (6k) + experiment (3k) | Experimental-gap MAE falls 22–45% vs single fidelity. This is the template for our ML layer |
| Autonomous glass labs | GlasDigital (Jena/BAM, 2021–24) [S33] | Aims for an automated, robot-assisted glass development workflow | The project page reports no published closed-loop optical-glass results. I found no autonomous lab for doped gain glasses |

### A2. What has NOT been done (the gap gainmat fills)

- **I found no paper that screens laser or fiber-amplifier gain media (4f–4f emitters such as Er, Nd, Tm, Ho, Yb, or Ti/Cr broadband centres) end to end from first principles or with ML.** The searches covered laser crystals, gain fibers and AI-designed lasers from 2023 to 2026. All high-throughput and ML work targets **5d→4f phosphors and scintillators** (Ce³⁺/Eu²⁺: [S1, S5–S12]) or single **glass-family JO regressions** on fewer than 100 samples [S19–S22].
- **No pipeline joins** a host band edge, ħω_max (multiphonon), 4f^N levels, JO intensities and an amplifier figure of merit (σ_e·τ, η, gain window) per (host, dopant) pair. Phosphor screens use Debye temperature [S7] or 5d position [S8] as the only proxy.
- **Ab initio Judd–Ofelt Ω for solids is still open.** Reid's 2019 notes state that ab initio methods "are not yet reliable" for 4f^N electric-dipole intensities [S34]. Molecular CASSCF work on Eu³⁺ (Janicki, Kędziorski, Mondry 2016) computed only energies and left intensities for future work [S35]. The Dijon "extended JO" models still fit three crystal-field parameters to experiment [S36, S37]. Hovhannesyan et al. cite a Wen et al. (2014) "purely ab initio" intensity calculation [S36]; **I could not open or confirm that paper.**
- **No multi-fidelity dataset** spans PBE → GLLB-SC/MLIP → CASSCF/NEVPT2 → experiment for dopant centres. The MEGNet multi-fidelity work covers host gaps only [S32].
- **No validated uMLIP workflow for fluoride or oxide laser-glass phonon spectra.** uMLIP melt-quench has known failures [S38–S40].

---

## B. Standard practice and accuracy-per-cost upgrades, by tier

### B1. Band gap and refractive index

- **Benchmarks.** On 472 solids, mBJ is the best density functional, with HSE06 close behind [S41]. Großmann et al. (2025) give MAEs of 0.50 eV for mBJ, 0.53 eV for HSE06, 0.54–0.60 eV for G₀W₀-PPA, 0.37 eV for full-frequency QP G₀W₀, and 0.30 eV for QSGŴ [S42]. Plain G₀W₀-PPA "gives only a marginal gain over the best DFT functionals" [S42].
- **GLLB-SC.** On 76 solids, GLLB-SC has an MAE of 0.64 eV (MARE 24%), vs 0.47 eV for mBJ and 0.82 eV for HSE06. It tends to **overestimate gaps above about 4 eV** [S43]. Every validation host here (fluorides, Al₂O₃, SiO₂, YAG) is in that regime, so GLLB-SC may push the UV edge too far into the UV. Calibrate it, or cross-check with mBJ.
- **Recommendation.** Keep GLLB-SC as the cheap tier. Add **TB-mBJ** as a second cheap fidelity (GPAW supports it via libxc; **unverified for GPAW specifically**). Use the HSE [S29] and TB-mBJ [S30] databases as free fidelity labels. Reserve QP G₀W₀ for a handful of anchors.
- **Refractive index.** DFPT-PBE ε∞ is systematically **overestimated**, mainly from the XC functional. Even so, n has MARD of 5.7% with local-field effects (7.4% without), Spearman 0.96 [S44]. That is good enough for the local-field factor in JO rates. For the 1,056 Materials Project compounds, use the MP DFPT values directly [S28]. Our code notes that the RPA+LFE ε∞ overestimates and the scissored value underestimates (`host/gap.py`). Petousis supports reporting the unscissored RPA/DFPT value as the default, with a systematic correction calibrated on the validation set.

### B2. Phonons and glass MD

- **uMLIP phonon benchmarks (MDR set, about 10k DFT phonon calculations; Loew et al.)** [S45]. MAE of ω_max in K: MatterSim 17, SevenNet 40, MACE(-MP-0) 61, CHGNet 89, ORB(v1/v2) 291, eqV2-OMat24 780. The MACE/CHGNet/M3GNet family **systematically softens** phonons. Direct-force models (early ORB, eqV2) give unreliable small-displacement phonons [S45].
- **Updated MDR table (eSEN paper)** [S46]. ω_max MAE: eSEN-30M-OAM 15, SevenNet-MF-ompa 15, MatterSim 17, GRACE-2L-OAM 19, eSEN-30M-MP 21. Conservative-force models clearly beat direct-force ones [S46].
- **ORNL INS benchmark (12 uMLIPs, 4,869 crystals)** [S47]. Best: ORB v3, SevenNet-MF-ompa, GRACE-2L-OAM. Then MatterSim, MACE-MPA-0, eSEN-OAM. "Most models underestimate high-frequency optical bands." **ω_max is exactly the quantity that sets our multiphonon cutoff.**
- **MACE docs** recommend **MACE-OMAT-0 as "most accurate for phonons"** [S48]. Matbench-Discovery's phonon (κ_SRME) leaderboard is led by newer OAM-trained models (TECE-OAM, EquFlashV2, EquiformerV3, PET-OAM-XL, GRACE-3L-OAM) [S49].
- **Cheap fix for softening.** Phonon fine-tuning (PFT) on DFT Hessians cuts phonon MAEs by about 50–55%, at about 30 A100-h for MACE-MP-0 [S50]. Cheaper still: compute a few DFT finite-difference ω_max per chemistry and fit a linear correction.
- **LO–TO splitting.** Local MLIPs cannot produce LO–TO splitting. Add the non-analytic correction from DFPT Born charges and ε∞ (phonopy NAC). Monacelli & Marzari show this needs no retraining [S51]. **In ionic fluorides and oxides the highest LO mode often *is* ħω_max**, so omitting NAC biases ħω_max low.
- **Glass MD with uMLIPs.** Eight uMLIPs (MACE, MatterSim, ORB, SevenNet, GRACE, UMA, …) give **far too low densities in NPT melt-quench**; 19 of 30 materials failed. Fix: quench under NVT at the reference density, then equilibrate under NPT at low temperature. This cut density MAE from 2.46 to 0.35 g/cm³ [S38]. On amorphous materials (41 uMLIPs incl. a-SiO₂), good force errors do not guarantee good structures, and fine-tuning helps numerically but not reliably structurally [S39]. MACE short-range models over-structure silica's medium-range order; long-range terms are "necessary but not sufficient" [S40].

### B3. Lanthanide 4f^N levels

- **Full Carnall Hamiltonian.** Carnall et al. (LaF₃, ANL-88-8) use F², F⁴, F⁶, ζ, α, β, γ, T² T³ T⁴ T⁶ T⁷ T⁸, and M⁰ (with M²=0.56 M⁰, M⁴=0.38 M⁰) and P² (with P⁴=0.5 P², P⁶=0.1 P²) [S52]. Reid: the higher-order terms "are important for accurate calculations, but the effects are subtle" [S34]. Our `fn_ci.py` omits T^i, M^k and P^k. Adding them is cheap (they are fixed tables of matrix elements) and is the standard way to reach Carnall-level RMS. The size of the gain for our levels is not quantified here.
- **Ab initio crystal-field parameters.** Novák, Kuneš, Knížek use DFT (4f in core) → Wannier projection → B^k_q. For 13 Ln³⁺ ions in LaF₃ they get RMS < 30 cm⁻¹ against experimental Stark splittings, with one adjustable hybridisation parameter Δ [S53, S54]. **This is the most direct route to host-specific B^k_q for screening.**
- **Superposition model** (Newman): intrinsic parameters per ligand with a power law in distance [S34, S55]. It is the cheap, transferable B^k_q generator for glasses and MD snapshots.
- **Wavefunction routes.** CASSCF/NEVPT2 with AILFT on Cs₂NaLnCl₆ across the series captures dynamic correlation and nephelauxetic effects adequately but **underestimates ligand-field splittings** [S56]. Minimal-active-space CASSCF-SO overestimates F^k by 12–30% and ζ by 2–10% [S57]. Its errors are systematic, so a 0.7–0.9 scaling is usable [S57]. Corrected-Racah CASSCF gives near-NEVPT2 Nd³⁺ emission at about 1/100 of the cost (2 h vs 237 h) [S58].
- **Implication.** For 4f–4f, a fitted-parameter Hamiltonian (Carnall free-ion plus ab initio or superposition B^k_q) is both more accurate and cheaper than direct CASSCF/NEVPT2. Keep wavefunction methods for 5d and charge-transfer states.

### B4. Ce³⁺ 5d and Ti³⁺ d–d

- **Embedded cluster.** AIMP embedding + SA-CASSCF/MS-CASPT2 + RASSI-SO (MOLCAS) is the standard [S13–S17]. Known errors for YAG:Ce: absorption +2300–3300 cm⁻¹, emission +4600 cm⁻¹, **Stokes shift strongly underestimated** because cluster relaxation is too small [S15]. Using the undistorted host gives "reasonable" trends across 21 garnets [S14]. Supercell DFT geometry followed by embedded-cluster CASPT2 gives good agreement for CaF₂:Ce centres [S17].
- **ΔSCF-DFT.** It is cheaper, gets 0.3 eV on energies, and gets the Stokes shift better (19% MARE) because the excited-state geometry relaxes in a periodic supercell [S1]. **Best hybrid:** take the relaxed geometries from ΔSCF-DFT (or MLIP-relaxed supercells) and run the cluster CASSCF/NEVPT2 on both the ground- and excited-state geometries. This follows [S1, S15, S17].
- **Higher level.** Compound-tunable embedding + FS-RCC reaches about 0.3 eV for Ce:YPO₄ and gives radiative τ within a factor of 2 [S18].
- **ML shortcut.** XGBoost reaches 0.15 eV MAE for Ce³⁺ 5d₁ excitation [S8], better than any ab initio route at negligible cost. Use it as a baseline and an extra fidelity.
- **Ti:sapphire.** Experimental targets: ²E JT components at 17,760 and 20,500 cm⁻¹; emission components at 12,680 and 14,210 cm⁻¹; Stokes shift about 5430 cm⁻¹ [S59]. A PBE supercell overestimates 10Dq (2.79 vs 2.37 eV) [S60]. An AIMD-averaged oscillator-strength approach reproduces the absorption band [S61]. **I found no recent embedded CASSCF/NEVPT2 benchmark for Ti:Al₂O₃**; our calculation would be new.

### B5. Nonradiative rates

- **Standard.** Energy-gap law, in either the Layne–Weber or the van Dijk–Schuurmans form. Our code implements both.
- **Beyond it.** Use a first-principles configuration-coordinate or multi-D Huang–Rhys treatment: ΔSCF geometries plus DFPT phonons [S3, S4]. The Nonrad code (Alkauskas/Turiansky/Van de Walle) implements the static-coupling nonradiative rate [S62]. It suits **strong-coupling 5d→4f and Ti³⁺** quenching. For weak-coupling 4f–4f multiphonon decay, **I found no first-principles rate calculation for a lanthanide in a solid**. A calculated phonon DOS of NaYF₄ matched INS after a 1.08 energy scaling [S63]; that supports using the computed PDOS rather than only ħω_max.
- **Thermal quenching of 5d.** The 4f–5d crossing barrier from ΔSCF is much larger than the measured quenching barriers, which points to auto-ionisation into the CB instead [S3]. So **the 5d–CBM gap (our "5d vs conduction band" output) is the right quenching descriptor**, consistent with [S5, S6].

### B6. Judd–Ofelt Ω from first principles

- **State of the art.** Not reliable [S34]. Partial routes exist. Hovhannesyan/Boudon/Lepers compute free-ion quantities with Cowan codes but fit three crystal-field parameters [S36, S37]. JOYSpectra computes Ω from structure with covalency and polarisability models; its exact models could not be confirmed from the abstract [S24]. Hu, Reid, Duan et al. extract crystal-field parameters (including odd-rank ones in principle) from DV-Xα cluster calculations [S64, S65].
- **Practical path.**
  1. Data-driven: train on LOMS's 1,200+ curated Ω's [S23] with host descriptors (n, ħω, coordination, covalency). This is better than glass-family regressions [S19–S22].
  2. Physics-based: odd-rank B^k_q from the superposition or point-charge model on DFT/MLIP geometries, plus dynamic coupling, giving A^λ_tp and then Ω. Calibrate on the validation set.

---

## C. Ranked recommendations: (accuracy gain) / (effort + compute)

| # | Change | Why | Effort / compute | Source |
|---|---|---|---|---|
| 1 | **Switch phonons from MACE-MP-0 to a conservative OAM/OMat model** (MACE-OMAT-0, SevenNet-MF-ompa, MatterSim, or eSEN-30M-OAM). Report a model ensemble spread | MACE-MP-0 ω_max MAE is about 61 K vs 15–17 K, and it softens optical modes. ħω_max drives W_nr exponentially | Hours; drop-in ASE calculators | [S45–S48] |
| 2 | **Add LO–TO NAC (Born charges + ε∞ from DFPT/GPAW, or MP)** to the phonopy runs | The top LO mode sets ħω_max in fluorides and oxides; local MLIPs miss it | One DFPT per host, or MP lookup | [S51, S28] |
| 3 | **Fix the glass melt-quench protocol**: NVT quench at experimental density, then low-T NPT. Validate the glass VDOS cutoff vs crystal analogues | uMLIP NPT quenches give drastically wrong densities | Protocol change only | [S38–S40] |
| 4 | **Add the three-body T^i, M^k, P^k terms** with Carnall LaF₃ values to `fn_ci.py` | The standard Hamiltonian; needed for Carnall-level RMS (gain size not quantified) | Days of coding; seconds to run | [S52, S34] |
| 5 | **Host-specific B^k_q** via the superposition model on relaxed geometries (glass snapshots included), with a DFT+Wannier (Novák) check on 2–3 crystals | Replaces "LaF₃ for everything". RMS < 30 cm⁻¹ demonstrated | SPM is cheap; Wannier needs one DFT per host | [S53–S55] |
| 6 | **Calibrate the GLLB-SC gap** (> 4 eV overestimation) by adding TB-mBJ plus the HSE/TB-mBJ databases as fidelities in a multi-fidelity graph net | Wide-gap hosts are exactly where GLLB-SC errs; the multi-fidelity method cuts experimental MAE by 22–45% | Low; data is public | [S43, S41, S42, S29, S30, S32] |
| 7 | **Ce³⁺ 5d: use ΔSCF-DFT (or MLIP)-relaxed ground and excited geometries for the cluster CASSCF/NEVPT2**, and add the Brgoch-type ML 5d₁ model as a baseline fidelity | Cluster-only relaxation underestimates the Stokes shift and overestimates emission by about 0.5 eV; ΔSCF gets 0.3 eV / 19% | Moderate (2 supercell DFT runs per centre) | [S1, S15, S17, S8] |
| 8 | **Learn Ω_λ from LOMS (1,200+ records, 550 hosts)** with physical host descriptors, as the JO tier until an ab initio route exists | Literature-Ω only covers known hosts; ab initio Ω is not reliable | Low | [S23, S34, S19] |
| 9 | Replace the single-ħω energy-gap law with a **PDOS-weighted multiphonon model**; use Huang–Rhys/Nonrad for 5d and Ti³⁺ quenching | Computed PDOS matches INS; 5d quenching goes via the CB, not crossing | Moderate | [S63, S4, S62, S3] |
| 10 | For 4f-only CASSCF checks, use minimal-CAS CASSCF with Racah/F^k scaling (0.7–0.9) instead of NEVPT2 | NEVPT2 costs about 100× more for similar band placement | Low | [S57, S58] |

---

## Sources

- [S1] Jia et al., Ce³⁺ first-principles vs semi-empirical — https://arxiv.org/pdf/1708.06225
- [S2] Jia et al., Eu²⁺ phosphors ΔSCF — https://arxiv.org/abs/1709.00311
- [S3] Jia et al., beyond 1D configuration-coordinate model — https://arxiv.org/abs/1909.09713
- [S4] Bouquiaux et al., Sr[Li₂Al₂O₂N₂]:Eu²⁺ multi-D coupling — https://www.arxiv.org/abs/2010.00423v2
- [S5] Canning et al., PRB 83, 125115 (2011) — https://link.aps.org/accepted/10.1103/PhysRevB.83.125115
- [S6] Pilania et al., physics-informed ML for scintillators (JCP 2018) — https://www.osti.gov/pages/servlets/purl/1467213
- [S7] Zhuo et al., Nat. Commun. 9 (2018) — https://ideas.repec.org/a/nat/natcom/v9y2018i1d10.1038_s41467-018-06625-z.html
- [S8] Lee et al. (Brgoch), ML of Ce³⁺ excitation band — https://arxiv.org/pdf/2502.18859
- [S9] Koyama et al., data-driven Eu²⁺ phosphors (2022) — https://pubs.rsc.org/en/content/articlehtml/2023/ma/d2ma00881e
- [S10] Zhang et al., GPT-assisted CGCNN for Eu²⁺ — https://arxiv.org/pdf/2306.14238
- [S11] Li, Amachraa, Wang, Ong, Xie, NIR Eu²⁺ nitrides (ECS 2022 abstract) — https://ecs.confex.com/ecs/242/webprogram/Paper164270.html
- [S12] Ong et al., Sr₂LiAlO₄ phosphor (APS March Meeting 2018 abstract) — https://meetings-archive.aps.org/mar/2018/s12/4
- [S13] Seijo & Barandiarán, Ce³⁺ in D₂ 8-fold O coordination — https://arxiv.org/abs/1203.5391
- [S14] Seijo & Barandiarán, host effects on Ce³⁺ in garnets, PCCP 2013 — https://dissem.in/p/29363011/host-effects-on-the-optically-active-4f-and-5d-levels-of-ce3-in-garnets
- [S15] Gracia, Seijo, Barandiarán et al., YAG:Ce ab initio, J. Lumin. 2007 — https://data.zhenzhi.com/details/1/875990020
- [S16] Seijo & Barandiarán, YAG:Ce red shifts (RSC Adv. 2016) — https://pubs.rsc.org/en/content/articlehtml/2016/ra/c6ra02611g
- [S17] Ning, Duan, Seijo et al., CaF₂:Ce charge-compensated centres (JPCC 2012) — https://core.ac.uk/works/33701467
- [S18] Oleynichenko et al., CTEP + FS-RCC for Ce/Th in YPO₄ — https://arxiv.org/abs/2310.09240v2
- [S19] Ahmadi et al., Sci. Rep. 14, 15505 (2024) — https://pmc.ncbi.nlm.nih.gov/articles/PMC11226621/
- [S20] ANN Judd–Ofelt Dy³⁺ lead borosilicate — https://mdpi.com/2079-9292/11/7/1045
- [S21] Singh et al., RF for JO parameters (2025) — https://update.library.uitm.edu.my/Record/2-s2.0-85215939814/Description
- [S22] Statistical learning for JO in Er³⁺ tellurite — https://portaldelaciencia.uva.es/documentos/607e9d2c9f431e6cf77704bb
- [S23] LOMS.cz JO platform and database — https://arxiv.org/pdf/2507.01218
- [S24] JOYSpectra (Optical Materials: X 2021) — https://zenodo.org/records/6201849 ; https://doaj.org/article/0d09e04027734ba08a329df7263f611b
- [S25] GlassNet — https://arxiv.org/abs/2303.15538
- [S26] SciGlass property predictor — https://docs.sciglass.uni-jena.de/guide/sidebar/property-predictor/
- [S27] Lun et al., Materials Futures 2025 — https://openaccess.library.uitm.edu.my/find/Record/doaj-art-01475b78020d41cb8ef4f93a4ee33049/Description
- [S28] Petousis et al., Sci. Data 4, 160134 (2017) — https://escholarship.org/content/qt0v43v5d6/qt0v43v5d6.pdf?t=oszjpj
- [S29] Kim et al., HSE band-gap database, Sci. Data 7, 387 (2020) — https://www.nature.com/articles/s41597-020-00723-8.pdf
- [S30] Choudhary et al., JARVIS OptB88vdW/TB-mBJ — https://www.arxiv.org/abs/1804.01032
- [S31] Trinquet, Evans, Rignanese, NLO active learning — https://web3.arxiv.org/abs/2504.01526
- [S32] Chen, Ong et al., multi-fidelity graph networks — https://arxiv.org/abs/2005.04338
- [S33] GlasDigital project — https://www.glas.uni-jena.de/202/glasdigital
- [S34] M. F. Reid, lecture notes on rare-earth electronic structure and intensities (2019) — https://spcs.canterbury.ac.nz/~mfr24/electronicstructure/00electronic.pdf
- [S35] Janicki, Kędziorski, Mondry, PCCP 2016 (ab initio Eu³⁺ energies) — https://pubs.rsc.org/en/content/articlehtml/2016/cp/c6cp05284c
- [S36] Hovhannesyan, Boudon, Lepers, revisiting JO — https://ar5iv.labs.arxiv.org/html/2106.02502
- [S37] Hovhannesyan et al., JO extension to Eu, Nd, Er — https://arxiv.org/pdf/2307.14176
- [S38] Li et al., melt-quench failures of uMLIPs — https://arxiv.org/pdf/2606.16385
- [S39] Fragapane & Deringer, amorphous materials vs uMLIPs — https://arxiv.org/pdf/2607.11384
- [S40] Balantrapu et al., silica MRO limits of MACE (2026) — https://www.emergentmind.com/papers/2604.21222
- [S41] Borlido et al., JCTC 15, 5069 (2019) — https://cris.fau.de/publications/296862709
- [S42] Großmann et al., MBPT vs DFT band-gap benchmark — https://arxiv.org/pdf/2508.05247
- [S43] Tran, Ehsan, Blaha, GLLB-SC assessment — https://arxiv.org/pdf/1712.06508
- [S44] Petousis et al., PRB 93, 115151 (2016), DFPT benchmark — https://link.aps.org/accepted/10.1103/PhysRevB.93.115151
- [S45] Loew et al., uMLIPs are ready for phonons — https://arxiv.org/html/2412.16551v1
- [S46] eSEN, smooth conservative potentials (Meta FAIR) — https://arxiv.org/html/2502.12147v2
- [S47] Han & Cheng (ORNL), uMLIPs for INS — https://arxiv.org/html/2506.01860
- [S48] MACE foundation models documentation — https://mace-docs.readthedocs.io/en/latest/guide/foundation_models.html
- [S49] Matbench Discovery phonon task — https://matbench-discovery.materialsproject.org/tasks/phonons
- [S50] Koker et al., Phonon Fine-Tuning (PFT) — https://arxiv.org/pdf/2601.07742
- [S51] Monacelli & Marzari, electrostatics in MLIPs for polar materials — https://arxiv.org/html/2412.01642v1
- [S52] Carnall, Goodman, Rajnak, Rana, LaF₃ systematic analysis (ANL-88-8) — https://www.osti.gov/servlets/purl/6995771
- [S53] Novák, Kuneš, Knížek, CF of RE impurities in LaF₃ — https://arxiv.org/abs/1407.7669
- [S54] Novák et al., CF with Wannier functions, RE aluminates — https://arxiv.org/pdf/1303.1281
- [S55] Superposition model chapter, Crystal Field Handbook — https://www.cambridge.org/core/books/crystal-field-handbook/superposition-model/7A282CCB6C7B01466FA6539BE02159D3
- [S56] Aravena, Atanasov, Neese, Inorg. Chem. 55, 4457 (2016) — https://ersearch2.cvtisr.sk/vufind/SummonRecord/FETCH-LOGICAL-a417t-ddf03484cbfe730a09967253e7b8f47e4259c2215a0a2592b5850eae31f239e12?lng=de
- [S57] Walisinghe & Chilton, minimal-CAS CASSCF-SO parameters (Dalton Trans. 2021) — https://pubs.rsc.org/en/content/articlehtml/2021/dt/d1dt02346b
- [S58] Gil et al., PCCP 27, 9847 (2025), corrected CASSCF for Nd³⁺ — https://pubs.rsc.org/en/content/articlepdf/2025/cp/d5cp00508f
- [S59] García-Revilla et al., Al₂O₃:Ti³⁺ under pressure, JPCM 14, 447 (2002) — https://personales.unican.es/rodriguf/TiAl2O3-JPCM.pdf
- [S60] Kravchenko & Fil, defect complexes in Ti:sapphire (DFT) — https://pdf.arxiv.org/pdf/1708.08850
- [S61] Bao & Ruan, Ti:sapphire absorption from first principles (J. Heat Transfer 2016) — https://ouci.dntb.gov.ua/works/4YepRw64
- [S62] Turiansky et al., Nonrad code — https://arxiv.org/abs/2011.07433
- [S63] Li et al., phonon DOS of NaYF₄/NaGdF₄ nanocrystals (INS/IXS + DFPT) — https://arxiv.org/abs/2012.01056
- [S64] Hu, Reid, Duan et al., CF parameters from DV-Xα cluster calculations — https://arxiv.org/pdf/1011.1709
- [S65] Reid et al., spectroscopy of high-energy lanthanide states — https://arxiv.org/abs/1002.3179

### Items searched but not confirmed
- Wen et al. (2014) "purely ab initio" 4f–4f intensities: cited in [S36]; full reference not opened.
- Pustovarov and Smentek ab initio intensity work: not located in accessible sources.
- Corning/Schott ML glass-design papers: not located in accessible sources (GlassNet and SciGlass cover the public state of the art).
- Dispersion-in-MLIP silicate-glass paper (PMC12079790): behind a bot check; not used.
- Loew et al. journal venue: arXiv version used.
