# GENIE spline production plan (plan only; nothing run)

> Status note (2026-10-08): this is the plan as written on 2026-10-06, kept for the record. It has been carried out; the paths in §3 refer to the development area at the time. The implementation and the decisions taken are in [README.md](README.md) and in the decision log of `src/genie/genie_splines.cxx`.

This follows [report.md](report.md). The dial *choice* follows GENIE's models. Only the format, knots and GUNDAM interface copy the NEUT chain ([neut_knot_settings.md](neut_knot_settings.md)).

Every step that submits a job or runs a long computation needs your OK, and the open decisions in `report.md` come first.

## 1. Dials and knots

**Splines.**
- Knots at −2, −1, 0, +1, +2σ (as NEUT's model parameters), converted to physical values with GENIE's σ.
- GUNDAM `dialType: Spline`, `catmull-rom, uniformity(1e-6)`.
- Bounds are the outer knots.
- Nominals are from the tune config. The producer must read them from the calculators at run time and stop if they differ.

| parameter | nominal | σ | knots |
|---|---|---|---|
| MaCCQE | 0.995 GeV | 0.030 | 0.935, 0.965, 0.995, 1.025, 1.055 |
| CCQEPauliSupViaKF | k_F scale 1 | 0.3 | 0.4, 0.7, 1.0, 1.3, 1.6 |
| CoulombCCQE | potential scale 1 | 0.2 | 0.6, 0.8, 1.0, 1.2, 1.4 |
| MaCCRES, MaNCRES | 1.089 GeV | 0.218 | 0.653, 0.871, 1.089, 1.307, 1.525 |
| MvCCRES, MvNCRES | 0.840 GeV | 0.084 | 0.672, 0.756, 0.840, 0.924, 1.008 |
| MaNCEL | 0.995 GeV | 0.249 | 0.498, 0.746, 0.995, 1.244, 1.493 |
| EtaNCEL | 0.12 | 0.036 | 0.048, 0.084, 0.120, 0.156, 0.192 |

If decision 6 is taken, MaCCQE's knots would span NEUT's *relative* ±2σ (±32.4%) around GENIE's 0.995 GeV, i.e. 0.673–1.317 GeV, preferably with 9 uniform knots. The RES masses need no change.

**Normalisations.**
- Prior 1, with σ from GENIE's table:
  - NormCCQE +0.20/−0.15; a Gaussian prior needs one value, proposed 0.20;
  - MEC 0.5;
  - COH 1.0;
  - NonRESBG 0.5.
- Lower bound 0. This is where the σ = 50% norms sit at −2σ. The COH norms, whose −2σ value would be −1, are bounded at 0.

**Unphysical −2σ.** Only NormCCCOH/NormNCCOH go below 0. All spline parameters stay positive at −2σ. MaCCQE negative weights are clamped to 0 and counted.

## 2. GUNDAM interface (checked in `~/work/projects/gundam/main/src`)

- **Per-event dials** (`DialCollection.cpp`). `treeExpression` names a TGraph (or TSpline3) leaf on the event tree.
  - An empty graph, or one flat at 1, makes no dial; a graph flat at another value becomes a `Shift`.
  - `catmull-rom` assumes uniform knots, so use `akima` or `Graph` (linear) for uneven knots.
  - The minimum response defaults to 0.
- **Norm dials.** A `Norm` is one dial whose weight is the parameter value. It applies to the events passing `applyCondition`, a TFormula on tree branches.
  - **Implementation issue:** the FakeSelect output tree has no GENIE mode or NonRESBG-class branch. It has only `cc`, `PDGnu`, `RespCat`, flavour, sample and the reco/truth kinematics, and `FakeSelect.C` copies only that fixed list plus TGraph branches.
  - Norm parameters therefore need either:
    - a mode/class flag carried through a *new copy* of the select step (FakeSelect.C itself unchanged), used as `applyCondition`; or
    - as a fallback, a 2-point linear TGraph (`Graph` dial) on the affected events only.
  - The first is closer to "normalisation parameter"; to settle together with the decisions.
- **Priors** (`ParameterSet.cpp`). `parameterDefinitionFilePath` points to a ROOT file holding:
  - names (TObjArray);
  - priors, lower and upper bounds (TVectorD);
  - `covarianceMatrix` (TMatrixD).

  Write `xsecCovGENIE.root` in the `MakeCovariances.C` layout: diagonal, physical units.

## 3. Pipeline

New outputs go under `/storage/shared/xtao/genie_dial_sensitivity/lbl_data/`, in the `<flav>/<stem>` layout.

1. **FakeRecon for GENIE (not yet done).** Run `GENIE_DATA=<new dir> STEP=import`, then `STEP=recon bash /storage/shared/xtao/run_genie_phase1.sh`. This writes `<stem>_fakerecon.root`; check that entry counts equal the GHEP files.
2. **Spline producer.** A new `macros/genie_splines.cxx`, with the dial_scan setup.
   - Per sample: read GHEP and fakerecon **by entry**, asserting equal counts and matching Eν, PDGnu and CC/NC.
   - Compute all knot weights in one pass.
   - Write a clone of `FlatTree_VARS` plus one TGraph branch per dial (split level 0, named as the `treeExpression`, x = physical value, y = weight, nominal knot = 1, empty if unaffected), and the class flags if that route is chosen.
3. **FakeSelect** with input suffix `_spline.root`. It copies TGraph branches by class.
4. **`xsecCovGENIE.root`** and a new GENIE xsec parameter-set yaml.

## 4. Parallelisation and time

- One task per sample (12), with all dials in one pass.
- Fast dials cost ~10–30 s per 200k events for 4 points.
- Expect ~15–20 min per 200k sample, about 30 min wall for all 12.
- Only option B (shape-mode RES) would need event-chunk splitting, at ~490 CPU-h per CC dial.

## 5. First test (when approved)

1. Run on numu_x_numu_fhc only: FakeRecon, producer, FakeSelect on that one file.
2. Check:
   - nominal knot = 1;
   - graph values at ±1σ and ±2σ reproduce the scan's per-sample rates;
   - every branch reaches the FakeSelect tree as TGraph.
3. Run a GUNDAM dry run with a minimal config and check the per-parameter dial counts, then inject +1σ on one spline and one norm and compare with the scan.
4. Only then produce all 12.
