# Does GENIE have a default list of dials, like T2KReWeight?

Build examined: GENIE Generator R-3_06_02 (`$GENIE=/home/xtao/work/projects/GENIE`), GENIE-Reweight 1.04.02 (`VERSION`), tune G18_10b_02_11b.
Date: 2026-09-29. This was a read-only source review. No builds were run and grwght1p was not run on events. The NEUT σ values in §6 come from one read-only ROOT read of `TweakMetaData`.

Path abbreviations used in citations:

- `RW/` = `/home/xtao/work/projects/GENIE-Reweight/`
- `GC/` = `/home/xtao/work/projects/GENIE/config/`
- `TUNE/` = `GC/G18_10b/G18_10b_02_11b/`
- `MLI/` = `/storage/shared/xtao/mock_lbl_inputs/`

---

## 1. Direct answer

**No. GENIE ships no curated default set (reading (b)). It ships only a menu (reading (a)) plus a table of per-dial 1σ values.**

- The framework starts empty. `GReWeight` adopts no calculators by default (`RW/src/RwFramework/GReWeight.cxx:33-37`), and `GSystSet` starts with no dials (`RW/src/RwFramework/GSystSet.cxx:21-24`). Every calculator and every dial has to be added explicitly by the caller (`GReWeight.cxx:44`, `GSystSet.cxx:36`).
- The one data file that looks like a default, `GC/GSystUncertaintyTable.xml`, only gives ±1σ for 90 of the 111 named dials. It says nothing about which dials to fit.
- The two apps that are built disagree on which calculators to use (§4). Neither is described as a fit configuration.
- The README states the policy directly. GENIE Reweight "does *not* provide the full systematic error for any GENIE comprehensive model or tune", and "the GENIE tuning procedure makes no use of the ReWeight product" (`RW/README.md:6-7`).
- Choosing the dials is left to each experiment. For your tune, many dials on the menu are unusable or have no effect (§3).

---

## 2. The full dial menu in this build

### Counts, reported separately without reconciling

| Source | Count | Citation |
|---|---|---|
| `GSyst_t` enum slots between `kNullSystematic` and `kNTwkDials` | **113** (so `kNTwkDials` = 114) | `RW/src/RwFramework/GSyst.h:46-241` |
| …of which have no name string (`kINVALID_INukeTwkDial_FrElas_pi`, `kINVALID_INukeTwkDial_FrElas_N`) | 2 | `GSyst.h:126,131`; string lines commented out at `GSyst.cxx:89,94` |
| Named dials (selectable with `-s`) | **111** | `RW/src/RwFramework/GSyst.cxx:28-140` (111 `temp_map[...]` assignments) |
| Extra enum value after `kNTwkDials`: `kProfRew` | 1, with no string | `GSyst.h:249` |
| Dials with an active entry in `GSystUncertaintyTable.xml` | **90** (plus 2 commented-out FrElas entries and 1 `ZExpELFF@CovarianceMatrix` param set) | `GC/GSystUncertaintyTable.xml:14-295, 191-193, 207-209, 298-300` |
| Named dials with **no** XML entry | **21**, which get σ = 0 | loader defaults to 0: `RW/src/RwFramework/GSystUncertaintyTable.cxx:67-68`; lookup returns 0 when missing: `GSystUncertainty.cxx:64` |

The 21 dials with no XML entry are AxFFCCQEshape, CCQEMomDistroFGtoSF, DecayAngMEC, the nine EmpMEC_* dials, NC, NormCCQEenu, NormDISCC, RnubarnuCC, RPA_CCQE, Theta_Delta2Npi, ThetaDelta2NRad, VecFFCCQEshape and XSecShape_CCMEC. Most are 0→1 morphing dials, where σ is not used. The **EmpMEC_*** dials are the exception: they *do* multiply by σ (`GReWeightXSecEmpiricalMEC.cxx:214-269`), so with σ = 0 they have no effect.

The XML is registered at `GC/master_config.xml:271` and loaded as param set "Default" (`GSystUncertainty.cxx:96-100`). It is the only copy I found. The install tree has no `config/` directory, and `$GENIE` points to the source tree.

**About the σ values.** Almost every calculator treats σ as a fractional change of the tune's nominal parameter, for example `Ma = Ma_nom·(1 + dial·σ)` (`GReWeightNuXSecCCQE.cxx:242`). A dial value of 1 means +1σ.

### Menu table

In the "Applies to G18_10b_02_11b?" column, "yes" means a nonzero weight can come out for this tune. `RW/src/RwCalculators/` is written `C/` below.

#### GReWeightNuXSecNCEL (`xsec_ncel`) — `IsHandled` at `C/GReWeightNuXSecNCEL.cxx:64-75`
| dial | enum | 1σ (XML line) | applies? | notes |
|---|---|---|---|---|
| MaNCEL | kXSecTwkDial_MaNCEL | 0.25 (14) | yes | uses the tune's QEL-NC model (`NCEL.cxx:207`), which is AhrensNCELPXSec (`GC/G18_10b/ModelConfiguration.xml:74`) |
| EtaNCEL | kXSecTwkDial_EtaNCEL | 0.30 (17) | yes | |

#### GReWeightNuXSecCCQE (`xsec_ccqe`) — `IsHandled` at `C/GReWeightNuXSecCCQE.cxx:90-158`
The tune uses `NievesQELCCPXSec/Dipole` (`ModelConfiguration.xml:71`), which uses `DipoleAxialFormFactorModel` (`GC/LwlynSmithFFCC.xml:37-38`). So `fModelIsDipole` is true (`CCQE.cxx:44,434`), and the default mode is NormAndMaShape (`CCQE.cxx:448-450`).

| dial | enum | 1σ | applies? | notes |
|---|---|---|---|---|
| NormCCQE | kXSecTwkDial_NormCCQE | +0.20/−0.15 (20-21) | yes | default mode |
| MaCCQEshape | kXSecTwkDial_MaCCQEshape | 0.025 (25) | yes | default mode; XML comment says "changed according to best fit" (23-24) |
| MaCCQE | kXSecTwkDial_MaCCQE | 0.03 (28) | yes | needs kModeMa; grwght1p switches to it automatically (`gRwght1Param.cxx:270-278`) |
| E0CCQE, E0CCQEshape | kXSecTwkDial_E0CCQE(shape) | +0.18/−0.16 (31-35) | **no** | needs MArunAxial FF (`CCQE.cxx:111-129`) |
| ZNormCCQE | kXSecTwkDial_ZNormCCQE | +0.20/−0.15 (37) | **no** | needs z-expansion FF (`CCQE.cxx:131-138`) |
| ZExpA1–A4CCQE | kXSecTwkDial_ZExpA{1..4}CCQE | 0.14, 0.67, 1.00, 0.75 (40-50) | **no** | needs z-expansion FF (`CCQE.cxx:140-150`) |
| RPA_CCQE | kXSecTwkDial_RPA_CCQE | none (morph) | yes | nonunit only for Nieves (`GSyst.h:205-207`); tune is Nieves |
| CoulombCCQE | kXSecTwkDial_CoulombCCQE | 0.20 (257) | yes | Nieves Coulomb on (`GC/NievesQELCCPXSec.xml:46`) |

#### GReWeightNuXSecCCQEELFF (`xsec_ccqe_elff`) — `C/GReWeightNuXSecCCQEELFF.cxx:74-116`
| dial | enum | 1σ | applies? | notes |
|---|---|---|---|---|
| ZExpELFFCCQE + 16 coefficients ZExpELFF{A,B}{P,N}{1..4}CCQE | kXSecTwkDial_ZExpELFF, _AP1…_BN4 | 1.0, then 0.0065…41.45 (262-295) | **no** | needs `ZExpELFormFactorModel` (`ELFF.cxx:45,97,415`); tune uses BBA07 (`TUNE/CommonParam.xml:253`) |

#### GReWeightNuXSecCCQEaxial (`xsec_ccqe_axial`) — `C/GReWeightNuXSecCCQEaxial.cxx:61-72`
| AxFFCCQEshape | kXSecTwkDial_AxFFCCQEshape | none (morph) | runs | morphs Dipole→ZExp using `<QEL-CC model>/ZExp` (`axial.cxx:212-213`); `NievesQELCCPXSec` has a ZExp set (`GC/NievesQELCCPXSec.xml:81`) |
|---|---|---|---|---|

#### GReWeightNuXSecCCQEvec (`xsec_ccqe_vec`) — `C/GReWeightNuXSecCCQEvec.cxx:61-72`
| VecFFCCQEshape | kXSecTwkDial_VecFFCCQEshape | none (morph) | runs, with a caveat | morphs the tune's EL FF → dipole (`vec.cxx:228-229`); the app itself warns it is "problematic for various tunes" (`gRwght1Param.cxx:239`) |
|---|---|---|---|---|

#### GReWeightNuXSecCCRES (`xsec_ccres`) — `C/GReWeightNuXSecCCRES.cxx:71-102`; GReWeightNuXSecNCRES (`xsec_ncres`) — `C/GReWeightNuXSecNCRES.cxx:64-95`
The tune uses BergerSehgalRESPXSec2014 (`ModelConfiguration.xml:78-79`). The calculator reads it from the GPL (`CCRES.cxx:211`, `NCRES.cxx:204`). The default mode is NormAndMaMvShape (`CCRES.cxx:233`, `NCRES.cxx:226`).

| dial | 1σ | applies? | notes |
|---|---|---|---|
| NormCCRES / NormNCRES | 0.20 (52, 67) | yes | default mode |
| MaCCRESshape / MaNCRESshape | 0.10 (55, 70) | yes | default mode |
| MvCCRESshape / MvNCRESshape | 0.05 (58, 73) | yes | default mode |
| MaCCRES / MaNCRES | 0.20 (61, 76) | yes | kModeMaMv; grwght1p switches automatically (`gRwght1Param.cxx:280-294`) |
| MvCCRES / MvNCRES | 0.10 (64, 79) | yes | same |

#### GReWeightNuXSecCOH (`xsec_coh`) — `C/GReWeightNuXSecCOH.cxx:59-75`
| dial | 1σ | applies? | notes |
|---|---|---|---|
| MaCOHpi | 0.40 (82) | yes | tune uses BergerSehgalCOHPiPXSec2015 (`ModelConfiguration.xml:114`), which has COH-Ma/COH-Ro (`GC/BergerSehgalCOHPiPXSec2015.xml:10-11`; paths at `COH.cxx:252-253`) |
| R0COHpi | 0.10 (85) | yes | |
| NormCCCOH | 1.0 (248) | yes | "placeholder" per XML comment (247) |
| NormNCCOH | 1.0 (252) | yes | "placeholder" (251) |

#### GReWeightNonResonanceBkg (`xsec_nonresbkg`) — `C/GReWeightNonResonanceBkg.cxx:47-78`
| NonRESBG{vp,vn,vbarp,vbarn}{CC,NC}{1pi,2pi} (16 dials; kXSecTwkDial_Rv…) | 0.50 each (88-134) | yes | scales the tune's DIS-HMultWgt-* values (`NonResBkg.cxx:202-209`). The calculator hard-codes Wmin = 2.0 GeV (`:198`), while the tune's Wcut is 1.809 (`TUNE/CommonParam.xml:99`); see §7 |
|---|---|---|---|

#### GReWeightNuXSecDIS (`xsec_dis`) — `C/GReWeightNuXSecDIS.cxx:71-100`
| dial | 1σ | applies? | notes |
|---|---|---|---|
| AhtBY, BhtBY, CV1uBY, CV2uBY | 0.25, 0.25, 0.30, 0.40 (142-152) | runs, with a caveat | default mode kModeABCV12u (`DIS.cxx:321`). The calculator hard-codes `genie::QPMDISPXSec/Default` (`DIS.cxx:306`), but the tune uses KNOTunedQPMDISPXSec (`ModelConfiguration.xml:97`) |
| AhtBYshape … CV2uBYshape | same values (154-164) | same caveat | grwght1p switches mode automatically (`gRwght1Param.cxx:296-307`) |

#### GReWeightNuXSecNC (`xsec_nc`) — `C/GReWeightNuXSecNC.cxx:45-56`
| NC | kXSecTwkDial_NC | none | yes | **not in σ units**: the weight equals the dial value itself (`NC.cxx:88,115`), and dial = 0 means "untweaked", giving weight 1 (`:93-94`) |
|---|---|---|---|---|

#### GReWeightDISNuclMod (`nuclear_dis`) — `C/GReWeightDISNuclMod.cxx:35-46`
| DISNuclMod | 1.00 (166) | **no — calls exit(-1)** | `CalcWeight` calls `exit(-1)` whenever dial ≠ 0 ("Not implemented", `DISNuclMod.cxx:82-86`) |
|---|---|---|---|

#### GReWeightFGM (`nuclear_qe`) — `C/GReWeightFGM.cxx:50-62`
| dial | 1σ | applies? | notes |
|---|---|---|---|
| CCQEPauliSupViaKF | 0.30 (169) | runs | uses `genie::PauliBlocker/Default` (`FGM.cxx:149`) |
| CCQEMomDistroFGtoSF | none (morph) | runs, but compares against the wrong model | hard-codes FGMBodekRitchie vs SpectralFunc (`FGM.cxx:315-318`), while the tune uses LocalFGM (`ModelConfiguration.xml:37`) |

#### GReWeightAGKY (`hadro_agky`) — `C/GReWeightAGKY.cxx:58-71`
| AGKYxF1pi | 0.20 (172) | yes (my assumption: I did not verify that G18_10b uses AGKY) | |
|---|---|---|---|
| AGKYpT1pi | 0.03 (175) | same | |

#### GReWeightResonanceDecay (`hadro_res_decay`) — `C/GReWeightResonanceDecay.cxx:50-63`; GReWeightDeltaradAngle (`delta_rad`) — `C/GReWeightDeltaradAngle.cxx:46-57`
| RDecBR1gamma | 0.50 (220) | yes | |
|---|---|---|---|
| RDecBR1eta | 0.50 (223) | yes | |
| Theta_Delta2Npi | none (morph) | yes | |
| ThetaDelta2NRad | none (morph) | yes | handled by DeltaradAngle |

#### GReWeightXSecMEC (`xsec_mec`) — `C/GReWeightXSecMEC.cxx:99-113`
The tune uses CC MEC = NievesSimoVacasMECPXSec2016, and NC/EM MEC = EmpiricalMECPXSec2015 (`ModelConfiguration.xml:135-138`).

| dial | 1σ | applies? | notes |
|---|---|---|---|
| NormCCMEC, NormNCMEC, NormEMMEC | 0.50 each (230-237) | yes | "placeholder" (226-229); map at `MEC.cxx:46-48` |
| FracPN_CCMEC | 0.20 (240) | yes | Valencia branch at `MEC.cxx:452-510` |
| FracDelta_CCMEC | 0.20 (244) | yes | Valencia only (`GSyst.h:197-199`) |
| DecayAngMEC | none (morph) | yes | |
| XSecShape_CCMEC | none (morph) | yes | morphs toward `EmpiricalMECPXSec2015/Reweight` (`MEC.cxx:261`) |

#### GReWeightXSecEmpiricalMEC (`xsec_empmec`) — `C/GReWeightXSecEmpiricalMEC.cxx:109-123`
| EmpMEC_Mq2d, _Mass, _Width, _FracPN_NC, _FracPN_CC, _FracCCQE, _FracNCQE, _FracPN_EM, _FracEMQE | **none → σ = 0** | **no effect** | the weight is 1 unless that channel's MEC model is Empirical (`EmpMEC.cxx:298-300`), which rules out CC for this tune. For NC/EM, σ = 0 makes the tweak a no-op (`EmpMEC.cxx:214-215`) |
|---|---|---|---|

#### GReWeightFZone (`hadro_fzone`) and GReWeightINuke (`hadro_intranuke`) — disabled in your build
| FormZone | 0.50 (178) | **no** | §3 |
|---|---|---|---|
| MFP_pi, MFP_N | 0.20, 0.20 (182-186) | **no** | §3 |
| FrCEx_pi, FrInel_pi, FrAbs_pi, FrPiProd_pi | 0.50, 0.40, 0.30, 0.20 (188-202) | **no** | hA fate fractions |
| FrCEx_N, FrInel_N, FrAbs_N, FrPiProd_N | 0.50, 0.40, 0.20, 0.20 (204-218) | **no** | hA fate fractions |

#### Enum entries that no calculator handles
| dial | enum | 1σ | status |
|---|---|---|---|
| **NormCCQEenu** | kXSecTwkDial_NormCCQEenu (`GSyst.h:50`) | none | **no calculator references it** (grep of `RwCalculators/*.cxx` and `Apps/*.cxx` found nothing). `-s NormCCQEenu` parses (`gRwght1Param.cxx:513`) and silently gives all weights = 1 |
| **NormDISCC** | kXSecTwkDial_NormDISCC (`GSyst.h:94`) | none | same: silent weight = 1 |
| **RnubarnuCC** | kXSecTwkDial_RnubarnuCC (`GSyst.h:95`) | none | same: silent weight = 1 |
| (FrElas_pi, FrElas_N) | kINVALID_… (`GSyst.h:126,131`) | commented out (XML 191-193, 207-209) | no string, so not selectable; "hA no longer has elastic fate" (`GSyst.h:125`) |
| (kProfRew) | `GSyst.h:249` | none | beyond `kNTwkDials`, no string; only GReWeightProfessor handles it (§3) |

---

## 3. Dials unusable in this setup

**3a. GReWeightFZone and GReWeightINuke. Your belief is confirmed.**

- GReWeightINuke constructor, `C/GReWeightINuke.cxx:59-81`:
  1. It reads `HadronTransp-Enable` and calls `std::exit(1)` if FSI is off (`:62-68`).
  2. It then runs `dynamic_cast<HAIntranuke2018*>` on the configured `HadronTransp-Model` (`:75`), and calls `std::exit(1)` if the cast returns null (`:77-81`).
- GReWeightFZone has the same guard, in `Init()`: `C/GReWeightFZone.cxx:225-245` (cast at `:239`, exit at `:233` and `:244`).
- Your tune sets `HadronTransp-Model = genie::HNIntranuke2018/Default` (`GC/G18_10b/ModelConfiguration.xml:60`). `HNIntranuke2018` and `HAIntranuke2018` are sibling classes; both derive from `Intranuke2018` (`$GENIE/src/Physics/HadronTransport/HNIntranuke2018.h:51`, `HAIntranuke2018.h:51`). The cast therefore returns null, and the program exits.
- Your edit is real: `diff gRwght1Param.cxx.orig gRwght1Param.cxx` shows only lines 248-249 changed to `// DISABLED (hN2018): …`.
- The FSI fate dials themselves are defined against the hA fate enum `INukeFateHA_t` (`GSyst.h:406-434`). So even apart from the guard, they parametrise hA physics.
- This removes 11 dials: FormZone, MFP_pi, MFP_N, and the 8 Fr* dials.

**3b. GReWeightProfessor. This is not a return-1 stub; it aborts.**

- `AppliesTo` returns true for every event (`C/GReWeightProfessor.cxx:18-20`).
- `SetSystematic`, `Reset` and `Reconfigure` are empty. `SetSystematic` is marked "TODO: figure out how to do it" (`:28-40`).
- `CalcWeight`, without `__GENIE_PROFESSOR2_ENABLED__`, logs FATAL and calls `std::exit(1)` (`:52-56`).
- In this build the macro is off: `//#define __GENIE_PROFESSOR2_ENABLED__` (`GENIE-install/R-3_06_02/include/GENIE/Framework/Conventions/GBuild.h:22`), and `GOPT_ENABLE_PROFESSOR2=NO` (`$GENIE/src/make/Make.config:36`).
- Its app `grwproftest` is built only when Professor2 is enabled (`RW/src/Apps/Makefile:25-27`), so it is absent from `RW/bin`. Neither built app adopts this calculator.

**3c. Other dials that are gated by model or broken**

- DISNuclMod: hard `exit(-1)` for any nonzero value (`DISNuclMod.cxx:82-86`). grwght1p still adopts it, with a comment warning not to use it (`gRwght1Param.cxx:252-254`).
- E0CCQE, E0CCQEshape, ZNormCCQE, ZExpA1-4CCQE: gated on the axial-FF model (`CCQE.cxx:111-150`), so they are not handled for Dipole.
- The 17 ZExpELFF* dials: gated on the z-exp EL FF (`ELFF.cxx:97-108`); the tune uses BBA07.
- The 9 EmpMEC_* dials: σ = 0 and CC MEC is not Empirical, so they have no effect (§2).
- NormCCQEenu, NormDISCC, RnubarnuCC: no calculator exists, so the weight is silently 1 (§2).
- Mode-exclusive groups: MaCCQE vs {NormCCQE, MaCCQEshape}, and MaCCRES/MvCCRES vs {NormCCRES, Ma/MvCCRESshape}, cannot be varied together in one calculator (`CCQE.cxx:102-129`, `CCRES.cxx:77-94`). grwghtnp enforces this (`gRwghtNCorrelatedParams.cxx:661-781`).
- Dials that run but with a model mismatch, with correctness unverified: the DIS BY dials (hard-coded QPMDISPXSec, `DIS.cxx:306`), CCQEMomDistroFGtoSF (hard-coded RFG, `FGM.cxx:315-318`), and VecFFCCQEshape (warned at `gRwght1Param.cxx:239`).

**Net effect for G18_10b_02_11b:** of the 111 named dials, 11 (FSI/FZone) + 1 (DISNuclMod) + 7 (E0CCQE, E0CCQEshape, ZNormCCQE, ZExpA1-4CCQE) + 17 (ZExpELFF) + 9 (EmpMEC) + 3 (orphans) = **48 do nothing useful or abort**. That leaves **63**, 10 of which (8 DIS BY dials, CCQEMomDistroFGtoSF, VecFFCCQEshape) carry model-mismatch caveats. This count is my own tally of the tables above.

---

## 4. Which calculators each app adopts

| | grwght1p (`gRwght1Param.cxx`) | grwghtnp (`gRwghtNCorrelatedParams.cxx`) | grwproftest |
|---|---|---|---|
| How calculators are chosen | Always the same fixed list (`:235-261`) | Only calculators needed by the requested `-s` dials (`AdoptWeightCalcs`, `:783-923`) | Professor only (`gRwProfTest.cxx:40,240`) |
| Built here? | yes (`RW/bin/grwght1p`) | yes (`RW/bin/grwghtnp`) | **no** (`Apps/Makefile:25-27`) |
| Calculators | NCEL, CCQE, CCQEELFF, CCQEaxial, CCQEvec, CCRES, NCRES, NonResBkg, COH, DIS, FGM, ResonanceDecay, ~~FZone~~, ~~INuke~~ (upstream: adopted, `.orig:248-249`), AGKY, DISNuclMod, NC, EmpiricalMEC, XSecMEC, DeltaradAngle: 19 upstream, 17 after your edit | CCQE, CCRES, DIS, NCRES, XSecMEC, CCQEELFF only | GReWeightProfessor |
| Dials per run | exactly one (`syst.Init(gOptSyst)`, `:266`) | several correlated dials, given a covariance matrix | — |
| CLI flags that change the calculator set | **none**. Options are only f, n, t, s, o, p, seed, min/max-tweak, cross-sections (`:440-590`). Modes are switched automatically from `-s` (`:270-307`) | none; the set follows from `-s` | — |

How grwghtnp handles dials outside its whitelist:

- NCEL, COH, NonResBkg, FSI, AGKY, ResDecay, FGM, RPA/Coulomb, the MEC FracPN/FracDelta/Shape/DecayAng dials, EmpMEC, NC and VecFF/AxFF all fall to `default: // no fine-tuning needed` (`:919-920`). **No calculator is adopted for them, so their weights are silently 1.**
- Upstream bug: for NormCCRES, MaCCRESshape and MvCCRESshape, grwghtnp adopts a **`GReWeightNuXSecCCQE`** under the name `"xsec_ccres"` (`:837-840`). The CCQE calculator does not handle those dials (`CCQE.cxx:95-155`), so they also have no effect in grwghtnp.

Separate upstream bug in grwght1p: the default `--max-tweak` is **−5**, not +5 (`gRwght1Param.cxx:584`; also present in `.orig`). The doc comment says +5 (`:45-47`). If you omit `--max-tweak`, every tweak point sits at −5σ. Always pass both limits.

**Conclusion for the supervisor:** the two built apps adopt different calculator sets (17 vs at most 6), and neither is presented as a set of dials to fit. In GENIE, "default" means only "the calculators this utility happens to instantiate". It is not a blessed analysis set.

---

## 5. Mapping: the 15 NEUT dials (`MLI/src/MakeCovariances.C:179-193`) vs GENIE

I could not read the NEUT-side semantics from source: T2KReWeight is not installed here, and `neut_dial_tweaks.toml` (`MLI/scripts/weight.sh:15`, `$T2KReWeight_ROOT/...`) was not found. The NEUT meanings in this table are inferred from the dial names and are labelled **[guess]**.

| NEUT dial | GENIE dial | Same parameter? |
|---|---|---|
| MaCCQE | MaCCQE (or MaCCQEshape + NormCCQE) | **Same parameter** (dipole axial mass). The CCQE model underneath differs: in GENIE it is Nieves+LFG (`ModelConfiguration.xml:37,71`) |
| MaRES | MaCCRES (+ MaNCRES) | **Different parametrisation.** Same name, but GENIE's is the Berger-Sehgal RES-Ma (`CCRES.cxx:240`). NEUT's is the Graczyk-Sobczyk FF axial mass **[guess]** |
| CA5RES | none | **Absent.** GENIE exposes no C5A. NormCCRES is only a normalisation proxy |
| BgSclRES | NonRESBG* (16 dials) | **Different parametrisation.** GENIE scales the DIS-extrapolated non-res bkg per channel and multiplicity (`NonResBkg.cxx:202-209`). NEUT scales an I=1/2 background term inside the RES model **[guess]** |
| NucleonFSI_Total | MFP_N | Different model (NEUT cascade vs GENIE hA MFP); **unusable here** (§3a) |
| NucleonFSI_Elastic | none (FrElas_N removed, `GSyst.h:125,131`) | **Absent** |
| NucleonFSI_Single | none clean (FrPiProd_N?) | **Absent / no clean match [guess]** |
| NucleonFSI_Double | none clean | **Absent [guess]** |
| PionFSI_AbsProb | FrAbs_pi | Different parametrisation (NEUT cascade-probability scale vs hA fate fraction); **unusable here** |
| PionFSI_QELowMomProb | FrInel_pi (closest) | Different; GENIE has no low/high-momentum split; **unusable** |
| PionFSI_QEHighMomProb | FrInel_pi (closest) | Different; **unusable** |
| PionFSI_CExLowMomProb | FrCEx_pi | Different; no momentum split; **unusable** |
| PionFSI_CExHighMomProb | FrCEx_pi | Different; **unusable** |
| PionFSI_InelProb | FrPiProd_pi (closest) **[guess]** | Different; **unusable** |
| SRCFrac_O | none (CCQEMomDistroFGtoSF is a different idea and hard-codes RFG) | **Absent** |

The FSI point is the important one. NEUT tunes microscopic cascade probabilities, and several have low/high-momentum variants. GENIE's FSI dials are hA2018 *fate fractions* plus a *mean free path*, which is a different set of physics quantities. In your hN2018 tune even those are unavailable. None of the 10 NEUT FSI dials has a usable GENIE counterpart in this setup.

---

## 6. Default 1σ comparison

**NEUT σ source.** `MakeCovariances.C` does not hard-code the σ values. It recovers them at run time as σ = (tweak − from)/mult from `TweakMetaData` (`MakeCovariances.C:169-172, 202-242`). I applied the same formula to the knots read (read-only) from `/storage/shared/haowei/neut_output/Mock_FD_Input_Data/numu/numu_x_nue_50K_FHC_Weight.flat.root`. I read one sample, and I assume the others are identical.

Both codes' σ are compared below as fractions of nominal. GENIE uses fractional σ natively (`CCQE.cxx:242`); T2KReWeight's are absolute, so I divided by `_from`.

### Same or same-named parameter

| NEUT dial | nominal, σ (abs) | NEUT σ (frac) | GENIE dial | GENIE σ (frac) | ratio NEUT/GENIE |
|---|---|---|---|---|---|
| MaCCQE | 1.21, 0.196 GeV | **16.2 %** | MaCCQE | 3 % (XML:28) | **5.4** |
| MaRES | 0.95, 0.15 GeV | **15.8 %** | MaCCRES | 20 % (XML:61) | **0.79** |

In absolute GeV: for MaCCQE, GENIE's 3 % of the tune's QEL-Ma = 0.994989 (`TUNE/CommonParam.xml:223`) is ±0.030 GeV, against NEUT's ±0.196 GeV (**6.6×**). For MaRES, GENIE's 20 % of RES-Ma = 1.088962 (`TUNE/CommonParam.xml:27`) is ±0.218 GeV, against NEUT's ±0.15 GeV (0.69×). I did not trace the config-inheritance chain to confirm that the dipole FF reads this exact QEL-Ma key; the calculator reads the path `FormFactorsAlg/AxialFormFactorModel/QEL-Ma` (`CCQE.cxx:444,451`).

### Analogues only (different physics; shown for scale, not equivalence)

| NEUT dial | NEUT σ (frac) | GENIE analogue | GENIE σ | ratio |
|---|---|---|---|---|
| CA5RES | 14.9 % | NormCCRES (proxy) | 20 % | 0.74 |
| BgSclRES | 11.5 % | NonRESBG* | 50 % | 0.23 |
| NucleonFSI_Total | 50 % | MFP_N | 20 % | 2.5 |
| NucleonFSI_Elastic/Single/Double | 50 % each | — | — | — |
| PionFSI_AbsProb | 30.8 % | FrAbs_pi | 30 % | 1.03 |
| PionFSI_QELowMomProb | 29.3 % | FrInel_pi | 40 % | 0.73 |
| PionFSI_QEHighMomProb | 47.1 % | FrInel_pi | 40 % | 1.18 |
| PionFSI_CExLowMomProb | 43.8 % | FrCEx_pi | 50 % | 0.88 |
| PionFSI_CExHighMomProb | 16.0 % | FrCEx_pi | 50 % | 0.32 |
| PionFSI_InelProb | 109.9 % | FrPiProd_pi | 20 % | 5.5 |
| SRCFrac_O | 97 % (0.34 ± 0.33) | — | — | — |

**The factor-of-5 pattern is specific to MaCCQE.** For the only other same-named dial, MaRES, GENIE's σ is *larger* than NEUT's (ratio 0.79). The analogue ratios range from 0.23 to 5.5 with no pattern. GENIE's small MaCCQE σ comes from a deliberate "best fit" update noted in the XML (`GSystUncertaintyTable.xml:23-24`). It is not a general GENIE convention of tighter priors.

---

## 7. Open questions and things I could not determine from source

1. **NEUT dial semantics.** T2KReWeight source and `neut_dial_tweaks.toml` are not on this machine (`weight.sh:15` points to `$T2KReWeight_ROOT`, which is unset, and `$HOME/setup_neut.sh` does not exist). The NEUT meanings in §5 are inferred from names.
2. **NEUT σ from a single file.** I read only `numu_x_nue_50K_FHC_Weight.flat.root`. `MakeCovariances.C` also reads only one file (per `MLI/README.md`), but I did not check that the other samples carry identical knots.
3. **DIS calculator correctness for this tune.** It builds `QPMDISPXSec/Default` (`DIS.cxx:306`) while events come from `KNOTunedQPMDISPXSec`. Whether that gives correct BY-dial weights is untested.
4. **NonResBkg Wmin = 2.0 GeV (`NonResBkg.cxx:198`) vs tune Wcut = 1.809 GeV (`TUNE/CommonParam.xml:99`).** I did not trace how the two cuts interact in `CalcWeight` (`:149`).
5. **CCQEMomDistroFGtoSF** morphs between hard-coded RFG and SF (`FGM.cxx:315-318`), while the tune is LocalFGM. The physical meaning for this tune is unclear.
6. **CoulombCCQE vs `ApplyCoulombCorrection=false`** (`TUNE/CommonParam.xml:284`). I did not determine whether the dial acts on the Nieves-internal Coulomb term only.
7. **AGKY hadronization in G18_10b.** Assumed, not verified in the tune config.
8. **Other copies of `GSystUncertaintyTable.xml`.** Only `GC/GSystUncertaintyTable.xml` was found. A filesystem-wide search did not finish in the time available, so an experiment-specific override elsewhere cannot be ruled out. GENIE would load whichever copy `GXMLPATH`/`$GENIE/config` resolves first.
9. **GENIE-provided covariances.** The README mentions only *plans* to release Professor response functions and global-fit covariance matrices (`README.md:10-11`). Nothing like that is shipped in this build.
