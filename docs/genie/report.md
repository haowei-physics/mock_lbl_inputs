# GENIE dial scan: which dials get splines (G18_10b_02_11b, 12 HK samples)

> Status note (2026-10-08): this report was written before the spline production. The pipeline item below (open question 4) has since been done; see [README.md](README.md). The scan is rerun with `validation/genie/run_checks.sh` (scan, summarize, aux), which writes its outputs under `$data_dir/validation/dial_scan/`.

Updated 2026-10-06. The scan reweighted each of the 53 GENIE Reweight (R-1_04_02) dials one at a time to −2, −1, +1 and +2σ. It then compared oscillated rates in the five analysis samples: FHC 1Rμ, 1Re0de, 1Re1de; RHC 1Rμ, 1Re0de. Full per-dial tables are in the scan outputs `summary.log` and `ratio_plots.pdf` (not in the repository; `$data_dir/validation/dial_scan/` when rerun). The auxiliary array `dialaux` is not needed for any conclusion here.

**Data checks.**
- All 12 scan tasks and the summary job finished cleanly: no errors, no NaN weights, 0 GHEP/flat mismatches.
- Nominal rates agree with `mode_stats` to 0.01 events: in total (104580.1), per flavour, per mode, and per sample. The scan's CCDIS is `mode_stats`' CCnonRES + CCDIS.
- No dial has all weights equal to 1.

## Rules used

1. **Broken for this tune → excluded for now.**
   - **BY dials:** plain QPM is divided by the stored KNOTuned dσ, which includes the DIS/RES join factor R < 1 below Wcut. This gives a nominal offset of about +0.2% total for an infinitesimal tweak. *Fix:* rescale each weight by dσ_stored/dσ_QPM,def.
   - **MaCOHpi/R0COHpi:** the generator stores dσ/dQ²dy × e^{−bt}/tsum as the event dσ, so every weight carries a factor tsum·e^{bt} and COH drops by 85% at any tweak. *Fix:* `UseOldWeightFromFile(false)` on the COH calculator.
2. **No effect possible on this chain → excluded.**
   - `FakeSelect.C` assigns samples from true Eν and RespCat only (flavour, CC/NC, CCQE vs non-QE). The fit variable comes from lepton kinematics.
   - So dials that change only the hadronic final state cannot act: AGKY, RDecBR and the MEC pn/Δ fractions.
   - Their small nonzero scan values are reweighting artifacts.
3. **Pure normalisation** (constant weight on a fixed event class) → GUNDAM normalisation parameter.
4. **Everything else → per-event splines,** however small.

The model choice follows GENIE's own models (Nieves CCQE, Valencia MEC, Berger-Sehgal RES/COH, KNO/BY DIS, hN FSI), not the T2KReWeight parameter set.

## Mode by mode

Max |Δ| is the largest change at ±1σ over the five analysis samples, not in the total rate. The total is dominated by unselected νμ and NC: MaNCEL, for example, is ±7% of the total but ≤0.7% in any selected sample.

| mode | GENIE dials (max \|Δ\| at ±1σ, sample) | verdict |
|---|---|---|
| CCQE | NormCCQE (18.2%, F 1Re0de), CCQEPauliSupViaKF (4.9%, R 1Re0de; one-sided), MaCCQE (2.6%, F 1Re0de), CoulombCCQE (0.4%), MaCCQEshape (0.05%) | NormCCQE norm; MaCCQE, Pauli k_F, Coulomb splines. MaCCQEshape is an alternative parametrisation (decision 2) |
| 2p2h / MEC | NormCCMEC (10.8%, F 1Re1de), NormNCMEC (0.14%), FracPN_CCMEC, FracDelta_CCMEC | two norms. pn/Δ fractions excluded (final state only). No dial for the Valencia MEC energy/shape dependence |
| CC1π (RES) | NormCCRES (11.6%, F 1Re1de), MaCCRES (8.9%), MvCCRES (6.8%), MaCCRESshape (2.3%), MvCCRESshape (1.8%) | MaCCRES, MvCCRES splines. NormCCRES vs the shape pair is decision 1 |
| NC1π⁰ (RES) | MaNCRES (0.3%, R 1Re0de), NormNCRES (0.2%), MvNCRES (0.06%), MaNCRESshape, MvNCRESshape (≤0.07%) | MaNCRES, MvNCRES splines; NormNCRES / shape pair per decision 1 |
| RES decays | RDecBR1gamma, RDecBR1eta | excluded (final state only) |
| non-resonant 1π/2π | 16 × NonRESBG{vp,vn,vbarp,vbarn}{CC,NC}{1pi,2pi}: largest vnCC2pi 1.3% (F 1Rμ), vpCC2pi 0.9%, vnCC1pi 0.6%; NC ones ≤0.04% | all 16 norms (weight 1 + 0.5n on DIS, W < 2 GeV, 1 nucleon + 1 or 2 π) |
| multi-π / DIS | AhtBY, BhtBY, CV1uBY, CV2uBY (+4 *shape twins; ~0.1–0.3% in F 1Rμ, partly the offset), AGKYxF1pi, AGKYpT1pi | **no usable dial now**: BY broken (rule 1), AGKY excluded (rule 2) |
| coherent π | NormCCCOH (1.0%, F 1Re1de), NormNCCOH (0.05%), MaCOHpi, R0COHpi | two norms; MaCOHpi/R0COHpi broken (rule 1) |
| NC elastic | MaNCEL (0.7%, R 1Re0de), EtaNCEL (0.05%) | both splines |
| FSI | none: the tune uses hN and the GENIE FSI calculators are hA-only (they `exit()`) | **gap**; the NEUT fit has 10 FSI parameters |
| other | CCother (hyperon QE, CC diffractive, charm; ≤3.9% of R 1Re0de), NCother (NC diffractive + ν-e) | no dial covers them |

**Response shape.**
- Most dials are close to linear over ±2σ: the ratio of the ±2σ change to twice the ±1σ change is within 0.83–1.17 for the RES and NC-elastic axial masses.
- CCQEPauliSupViaKF is one-sided (asymmetry ≈ −0.9).
- NormCCCOH/NormNCCOH clip at 0 for −1σ and −2σ.
- Splines (5 knots) handle all of these.

**Weight pathologies.**
- No NaN.
- MaCCQE gives negative weights for 12 events at +1σ and 42 at +2σ, all ν̄ CCQE on oxygen near the kinematic edge → clamp to 0.
- w = 0 occurs at −2σ for the σ = 50% norms (the parameter reaches the boundary 0). NormCOH reaches −1 at −2σ (unphysical) and is clipped.
- w > 10 appears in the high-Q² tails of the axial masses at +2σ, which is physical.
- The other pathologies belong to excluded dials.
- The five shape dials were run on 5000 events per sample. That is enough to rank them; they are not proposed anyway.

## Checks: normalisation dials and RES/CCQE σ (read-only: source + existing scan output)

### Normalisation dials

**How the weight is formed (from the source):**
- **NormCCQE:** `GReWeightNuXSecCCQE`, in its default mode `kModeNormAndMaShape`, returns `fNormCurr = max(0, 1 + n·σ±)` for every CC QE event that is not charm or strange. It uses the signed σ, so the response is asymmetric.
- **NormCCMEC / NormNCMEC:** `GReWeightXSecMEC::CalcWeightNorm` gives `max(0, 1 + n·σ)` per MEC interaction type, using the mean σ.
- **NormCCCOH / NormNCCOH:** `GReWeightNuXSecCOH` gives `max(0, 1 + n·σ)` per CC/NC COH event. This norm-only path does not touch the stored dσ, so it is unaffected by the COH bug.
- **NonRESBG (all 16):** `GReWeightNonResonanceBkg` gives `R_curr/R_def` with `R_curr = max(0, R_def·(1 + 0.5n))`.
  - It applies to DIS events with W < 2 GeV, 2–3 primary hadrons and exactly one nucleon, split by ν/ν̄, p/n, CC/NC and 1π/2π.
  - So the weight is `max(0, 1 + 0.5n)`, independent of R_def.

At each tweak point every affected event therefore gets one and the same weight.

**Scan evidence:**
- The weight extremes per point are exactly {1, 1 + nσ}. For example NonRESBG: −2σ [0, 1], −1σ [0.5, 1], +1σ [1, 1.5], +2σ [1, 2]; NormCCQE [0.70, 1], [0.85, 1], [1, 1.20], [1, 1.40].
- The change of the affected mode is exactly nσ:
  - NormCCQE on CCQE: −30.000, −15.000, +20.000, +40.000%;
  - MEC: ±50, ±100%;
  - COH: −100, −100, +100, +200%.
- The ratio of the ±2σ change to twice the ±1σ change is 1.0000 (to 4 decimals) on the total and on the main mode for all 21 dials. The only exception is COH at −2σ (0.5000, from clipping).
- Caveat: the stored statistics keep min/max over *all* events, including unaffected w = 1, not over w ≠ 1 only. So "one weight per point" rests on the source plus these exact numbers, not on a direct per-event histogram.

| dial | physical value at −2, −1, +1, +2σ | reaches 0 / negative | proposed GUNDAM lower bound |
|---|---|---|---|
| NormCCQE | 0.70, 0.85, 1.20, 1.40 | no | 0 |
| NormCCMEC, NormNCMEC | 0, 0.5, 1.5, 2.0 | 0 at −2σ | 0 |
| NormCCCOH, NormNCCOH | −1 (clipped to 0), 0, 2, 3 | 0 at −1σ; **−1 (unphysical) at −2σ** | 0 |
| 16 × NonRESBG | scale 0, 0.5, 1.5, 2.0 of R_def | 0 at −2σ | 0 |

The NonRESBG R_def values are the tune's `DIS-HMultWgt` values:

| ν | CC 1π / 2π | NC 1π / 2π |
|---|---|---|
| vp, vbarn | 0.008 / 0.943 | 0.1 / 1.0 |
| vn, vbarp | 0.0308 / 2.338 | 0.3 / 1.0 |

All of these are norms, so a lower bound of 0 is the physical limit; the knots used for the scan don't apply.

- **Asymmetric ±1σ** in `GSystUncertaintyTable.xml`: among these dials, **only NormCCQE** (+0.20 / −0.15). The other asymmetric entries (E0CCQE, E0CCQEshape, ZNormCCQE) are not applicable to this tune. The MEC and COH calculators use the mean σ anyway.
- **NonRESBG W range:** the calculator reweights DIS with W < 2.0 GeV, but the tune applies R only below Wcut = 1.809 GeV. Events at 1.809–2.0 GeV, which in the model carry no R, are scaled as well, so the dials slightly over-cover.

### MaCCQE and RES masses: raw σ → physical knots

All of these calculators set `M_curr = M_def·(1 + n·σ)`, so the table entry is a *fractional* error.
- CCQE `kModeMa` uses the signed σ; CCRES/NCRES `kModeMaMv` use the mean of ±. All five entries are symmetric anyway.
- `M_def` is read at run time from the model registry. The values below are the tune's: QEL-Ma from the `QuasiElastic` set; RES-Ma from the tune's `Tunable` set, which overrides the model default of 1.120; RES-Mv is the model default.
- The scan did not print them; the producer must log them.

| dial | XML +1σ / −1σ | nominal | 1σ (phys.) | knots −2, −1, 0, +1, +2σ |
|---|---|---|---|---|
| MaCCQE | 0.03 / 0.03 | 0.994989 GeV | 0.0298 GeV | 0.935, 0.965, 0.995, 1.025, 1.055 |
| MaCCRES | 0.20 / 0.20 | 1.088962 GeV | 0.218 GeV | 0.653, 0.871, 1.089, 1.307, 1.525 |
| MvCCRES | 0.10 / 0.10 | 0.840 GeV | 0.084 GeV | 0.672, 0.756, 0.840, 0.924, 1.008 |
| MaNCRES | 0.20 / 0.20 | 1.088962 GeV | 0.218 GeV | as MaCCRES |
| MvNCRES | 0.10 / 0.10 | 0.840 GeV | 0.084 GeV | as MvCCRES |

**GENIE's MaCCQE 1σ is indeed 3%** (0.0298 GeV).
- The entry is 0.03 in the default table, next to the comment "changed according to best fit" (which precedes MaCCQEshape = 0.025).
- The scan agrees: ±1σ moved the CCQE rate by −2.66% / +2.70%, as expected for a ~3% M_A shift.
- For comparison, NEUT's MaCCQE is 0.196 / 1.21 = **16.2%**, so GENIE's prior is ~5× narrower in relative terms.
- For the RES axial mass the order is reversed: GENIE 20% vs NEUT MaRES 0.15 / 0.95 = 15.8%.

## Proposed spline list (per-event TGraph)

MaCCQE, CCQEPauliSupViaKF, CoulombCCQE, MaCCRES, MvCCRES, MaNCRES, MvNCRES, MaNCEL, EtaNCEL.

These may be added after the open decisions: the 4 BY dials (after the fix), MaCOHpi/R0COHpi (after the fix), and MaCCRESshape/MvCCRESshape/MaNCRESshape/MvNCRESshape (option B of decision 1).

## Normalisation list

- NormCCQE, NormCCMEC, NormNCMEC, NormCCCOH, NormNCCOH
- the 16 NonRESBG dials
- NormCCRES and NormNCRES only under option B of decision 1.

The FakeSelect tree has no GENIE mode branch, so GUNDAM `Norm` conditions need a class flag carried through. See [spline_plan.md](spline_plan.md) §2.

## Excluded

| dial(s) | reason |
|---|---|
| AhtBY, BhtBY, CV1uBY, CV2uBY, AhtBYshape, BhtBYshape, CV1uBYshape, CV2uBYshape | broken: nominal offset from plain QPM vs stored KNOTuned dσ (fix: rescale by dσ_stored/dσ_QPM,def; shape twins redundant) |
| MaCOHpi, R0COHpi | broken: t factor in stored COH dσ (fix: `UseOldWeightFromFile(false)`) |
| AGKYxF1pi, AGKYpT1pi, RDecBR1gamma, RDecBR1eta, FracPN_CCMEC, FracDelta_CCMEC | final-state only; no effect possible with FakeSelect |
| MaCCQEshape | alternative to MaCCQE + NormCCQE (decision 2) |
| MaCCRESshape, MvCCRESshape, MaNCRESshape, MvNCRESshape | alternative RES parametrisation (decision 1) |

The GENIE dials that were not scanned are inapplicable to this tune: z-expansion and running-M_A form factors, BBA07-incompatible vector form-factor dials, hA FSI, DISNuclMod (not implemented), and NormEMMEC.

## Open decisions (proposals only; nothing acted on)

1. **RES parametrisation.** Proposed (A): MaCCRES, MvCCRES, MaNCRES, MvNCRES splines with no RES norms.
   - MaCCRES ±1σ already moves CCRES by −22/+29%, more than NormCCRES's ±20%, so adding the norm would double-count.
   - Option A is fast.
   - The alternative (B) is NormCCRES/NormNCRES plus the shape-mode Ma/Mv. It costs ~7 s per RES event for 4 points (~490 CPU-h per CC dial), and here it is not rate-preserving.
2. **CCQE.** Proposed: MaCCQE (spline) + NormCCQE (norm). MaCCQE's rate effect (±2.7% of CCQE) is small next to NormCCQE (+20/−15%), so the double counting is minor. Drop MaCCQEshape.
3. **MaCCQE prior.** GENIE's table gives ±3% (0.0298 GeV around 0.995; confirmed in the calculator, see the checks above). NEUT uses ±16.2% (0.196 GeV around 1.21). Proposed: GENIE's value.
4. **Branch/parameter names.** Proposed prefix `GENIE_` (e.g. `GENIE_MaCCQE`) to avoid clashing with NEUT's `MaCCQE`.
5. **BY and COH.** Proposed: include them only after their fixes give w(±0.01σ) = 1.
6. **Knot range.** Should the knots in physical units span NEUT's σ as well, so a widened prior needs no new splines?
   - Compare **relative** widths, each around its own generator's nominal (GENIE MaCCQE 0.995 GeV, NEUT 1.21 GeV), not NEUT's absolute range.
   - MaCCQE: NEUT ±2σ is ±32.4%, which around GENIE's 0.995 GeV means 0.673–1.317 GeV, against GENIE's own ±2σ of 0.935–1.055 GeV (±6%).
   - RES axial mass: GENIE's ±20% already exceeds NEUT's ±15.8%, so nothing changes. NEUT has no Mv dial to compare with.
   - The cost of the wider MaCCQE range is a coarser spline near nominal unless more knots are used, for example 9 uniform knots at ±0.5σ_NEUT steps.

## Open questions for the mentor

1. **FSI gap.** With hN there are no usable GENIE FSI dials, while the NEUT fit has 10. Should we accept the gap, switch to an hA sample or tune, or add an effective parameter?
2. **Broken BY/COH dials.** Both fixes are known and need no GENIE source change. Is it worth applying them, given DIS/COH are small in these samples, or should DIS/COH be covered by norms only?
3. **Priors without a covariance.** GENIE ships only a ±1σ table with no correlations. Is a diagonal prior from that table acceptable, or should priors come from elsewhere (for example NEUT-like widths or external fits)?
4. **GENIE pipeline.** GENIE has not yet been run through FakeRecon/FakeSelect (`genie_lbl_data/stock` holds only symlinks). That step is needed before any spline can be used.

## Assumptions / not verified

- Upstream GENIE Reweight R-1_04_02 was used, because the mastbaum fork is not on disk.
- σ values are GENIE's `GSystUncertaintyTable.xml` defaults.
- One-at-a-time scan, true Eν only (no E_rec shapes), and FakeSelect expectation values instead of the random draw.
- Five shape dials were subsampled.
- The cause of the negative MaCCQE weights (ν̄ Nieves tensor near the kinematic edge) was not traced in detail.

## Run commands

```bash
scripts/genie/build.sh
validation/genie/run_checks.sh scan [sample]   # one sample, or all 12 in turn
validation/genie/run_checks.sh summarize       # summary only
validation/genie/run_checks.sh aux [sample]    # optional checks (zero counts, ±0.01σ offsets, subsample errors)
```

Run from `$working_dir` with `data_dir` set; see [README.md](README.md) for the environment and for running the samples as SLURM jobs.
