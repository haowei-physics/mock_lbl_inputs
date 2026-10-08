# NEUT knot settings and spline format in the mock_lbl_inputs → GUNDAM chain

Compiled 2026-10-05, read-only survey. Abbreviations used below:

- `MLI` = `/storage/shared/xtao/mock_lbl_inputs` (same content as `/storage/shared/haowei/neut_output/mock_lbl_inputs`; `diff -rq` excluding `.git` shows no differences)
- `T2KRW` = `/home/haowei/projects/T2KReWeight`
- `GLC` = `/storage/shared/haowei/neut_output/gundam_lbl_config`
- `DATA` = `/storage/shared/haowei/neut_output/Mock_FD_Input_Data_2` (the existing `$data_dir` outputs I inspected)
- `GUNDAM` = `/home/haowei/projects/gundam/src`

## 1. Dial/knot configuration

The knot config is **not** in MLI. `MLI/scripts/weight.sh` passes
`C=${T2KReWeight_ROOT}/var/T2KReWeight/config/CI/neut_dial_tweaks.toml` to `T2KReWeightCI`. That file is
`T2KRW/build/Linux/var/T2KReWeight/config/CI/neut_dial_tweaks.toml`, identical to `T2KRW/config/CI/neut_dial_tweaks.toml`.
Its tweak values are **in σ units**. T2KReWeight converts them to physical values and stores them in the
`TweakMetaData` tree (`<dial>` = physical tweak values, `<dial>_from` = nominal). I read the physical values below from
`DATA/numu/numu_x_numu_200K_FHC_Weight.flat.root`. σ is the ratio (tweak − from)/(σ multiple), and it matches the σ that
`MLI/src/MakeCovariances.C:169-246` writes to `DATA/parameters/xsec/xsecCovHK.root`.

"Spline knots" counts the TGraph points per event after `MakeSplines.C` adds the nominal (weight = 1) as a knot. It does
not add one when the nominal is already a tweak point (`MLI/src/MakeSplines.C:43-56,166-172`).

| Dial | Tweaks (σ) | Nominal | σ (phys.) | Physical tweak knots | Spline knots | GUNDAM range (bounds) | In GUNDAM? | Type |
|---|---|---|---|---|---|---|---|---|
| MaCCQE | −2,−1,1,2 | 1.21 | 0.196 | 0.818, 1.014, 1.406, 1.602 | 5 (uniform) | [0.818, 1.602] | yes | Spline, catmull-rom |
| MaRES | −2,−1,1,2 | 0.95 | 0.15 | 0.65, 0.80, 1.10, 1.25 | 5 (uniform) | [0.65, 1.25] | yes | Spline, catmull-rom |
| CA5RES | −2,−1,1,2 | 1.01 | 0.15 | 0.71, 0.86, 1.16, 1.31 | 5 (uniform) | [0.71, 1.31] | yes | Spline, catmull-rom |
| BgSclRES | −2,−1,1,2 | 1.30 | 0.15 | 1.00, 1.15, 1.45, 1.60 | 5 (uniform) | [1.00, 1.60] | yes | Spline, catmull-rom |
| MDLSPiEj | −3,−2 | 3 | 1 | 0, 1 | 3 (0,1,3; non-uniform) | — | **no** | TGraph written, unused |
| NucleonFSI_Double | −1,−0.5,0.5,1 | 1 | 0.5 | 0.5, 0.75, 1.25, 1.5 | 5 (uniform) | [0.5, 1.5] | yes | Spline, catmull-rom |
| NucleonFSI_Elastic | −1,−0.5,0.5,1 | 1 | 0.5 | 0.5, 0.75, 1.25, 1.5 | 5 | [0.5, 1.5] | yes | Spline, catmull-rom |
| NucleonFSI_Single | −1,−0.5,0.5,1 | 1 | 0.5 | 0.5, 0.75, 1.25, 1.5 | 5 | [0.5, 1.5] | yes | Spline, catmull-rom |
| NucleonFSI_Total | −1,−0.5,0.5,1 | 1 | 0.5 | 0.5, 0.75, 1.25, 1.5 | 5 | [0.5, 1.5] | yes | Spline, catmull-rom |
| NucleonFSI_PB | −1,−0.5,0.5,1 | 1 | 0.5 | 0.5, 0.75, 1.25, 1.5 | 5 | — | **no** | TGraph written, unused |
| PionFSI_AbsProb | −1,−0.5,0.5,1 | 1.404 | 0.432 | 0.972, 1.188, 1.620, 1.836 | 5 | [0.972, 1.836] | yes | Spline, catmull-rom |
| PionFSI_QELowMomProb | −1,−0.5,0.5,1 | 1.069 | 0.313 | 0.756, 0.9125, 1.2255, 1.382 | 5 | [0.756, 1.382] | yes | Spline, catmull-rom |
| PionFSI_QEHighMomProb | −1,−0.5,0.5,1 | 1.824 | 0.859 | 0.965, 1.3945, 2.2535, 2.683 | 5 | [0.965, 2.683] | yes | Spline, catmull-rom |
| PionFSI_CExLowMomProb | −1,−0.5,0.5,1 | 0.697 | 0.305 | 0.392, 0.5445, 0.8495, 1.002 | 5 | [0.392, 1.002] | yes | Spline, catmull-rom |
| PionFSI_CExHighMomProb | −1,−0.5,0.5,1 | 1.80 | 0.288 | 1.512, 1.656, 1.944, 2.088 | 5 | [1.512, 2.088] | yes | Spline, catmull-rom |
| PionFSI_InelProb | **−0.9**,−0.5,0.5,1 | 1.002 | 1.101 | 0.0111, 0.4515, 1.5525, 2.103 | 5 (**non-uniform**) | [0.0111, 2.103] | yes | Spline, **akima** |
| SRCFrac_C | 0.25, 0.5 (one-sided) | 0.34 | 0.33 | 0.4225, 0.505 | 3 (uniform) | — | **no** | TGraph written, unused |
| SRCFrac_O | 0.25, 0.5 (one-sided) | 0.34 | 0.33 | 0.4225, 0.505 | 3 (uniform) | **[0.34, 0.505]** (cannot go below nominal) | yes | Spline, catmull-rom |
| RadCorr_rate | 0, 0.25, 0.75 | 0 | 1 | 0, 0.25, 0.75 | 3 (nominal is already a knot) | — | **no** | TGraph written, unused |

Sources: `T2KRW/config/CI/neut_dial_tweaks.toml` (σ multiples), `TweakMetaData` in
`DATA/numu/numu_x_numu_200K_FHC_Weight.flat.root` (physical values), `DATA/parameters/xsec/xsecCovHK.root` (priors, σ and
bounds, dumped with ROOT), and `GLC/inputs/parameters/xsec/configXsecParSet.yaml` (which dials are used and their dial type).

**Normalisation vs spline.** All 15 cross-section parameters used in GUNDAM are per-event `Spline` dials. There are **no
cross-section normalisation dials**: no MEC norm, CCRES norm, NC norm or ν/ν̄ norm (`GLC/inputs/parameters/xsec/configXsecParSet.yaml`).
The only normalisation (Norm) dials in the fit are:
- flux bins: `GLC/inputs/parameters/flux/configFluxParSet.yaml` (`parametersBinningPath`, `parametersRange [0,2]`)
- per-sample detector norms: `GLC/inputs/parameters/det/configDetParSet.yaml`, with σ values from `MLI/src/MakeCovariances.C:88-94`

## 2. Spline-production chain

The stages run in this order (`MLI/scripts/run_all.sh`, `MLI/README.md`):

1. **`gen.sh`**: `neutroot2 gen/neut/<flav>/NEUT.*.toml` → `$data_dir/<flav>/<sample>.root` (12 samples, e.g. `numu_x_numu_200K_FHC`).
2. **`weight.sh`**: `T2KReWeightCI --no-plot -c neut_dial_tweaks.toml --neut-card <card> -i <sample>.root -o <sample>_Weight.root`.
   The files on disk are named `<sample>_Weight.flat.root` (see §5). They contain:
   - tree `T2KRWFlatTree`: one `vector<double>` branch per dial, holding the weight at each tweak in TweakMetaData order, plus `cv_weight`, `Mode`, `Enu`-type kinematics, etc.
   - tree `TweakMetaData`: `<dial>` (`vector<double>`, physical tweak values) and `<dial>_from` (double, nominal).
   - Pauli-blocked events are dropped (`MLI/src/MakeSplines.C:76-82`).
3. **`flat.sh`**: `nuisflat -i NEUT:<sample>.root` → `<sample>_flattree.root`.
4. **`recon.sh`**: `FakeRecon.C` → `<sample>_fakerecon.root`, tree `FlatTree_VARS`. It drops the same events (Q2 = −999).
5. **`splines.sh`**: `MLI/src/MakeSplines.C(fakerecon, Weight.flat, out)` → `<sample>_spline.root`, tree **`FlatTree_VARS`**.
   This is a clone of the fakerecon tree plus **one `TGraph` branch per dial, named exactly as the dial** (e.g. `MaCCQE`),
   written with split level 0 (`MLI/src/MakeSplines.C:204-219`). For each event and dial:
   - x = physical parameter value, y = weight.
   - Knots are sorted, and the nominal is added with y = 1 (`MakeSplines.C:150-183`).
   - All weights exactly 1 → empty TGraph, which GUNDAM treats as no dial (`MakeSplines.C:58-71,288-293`).
   - All weights 0, or any weight negative/non-finite, or weight ≠ 1 at a nominal knot → empty TGraph, counted as "rejected" (`MakeSplines.C:239-286`).
   - A response that is flat at a constant other than 1 is written in full, and GUNDAM turns it into a Shift dial (`MakeSplines.C:66-70`).
6. **`select.sh`**: `MLI/src/FakeSelect.C` chains `FlatTree_VARS` from every `*_spline.root` of one beam mode
   (`FakeSelect.C:28,103-112`). It writes **`$data_dir/fakeselect_{FHC,RHC}.root`, tree `FakeSelect`**. It copies every
   branch that inherits from `TGraph`, again with split level 0, found by class rather than by name (`FakeSelect.C:191-224`).
   The FakeSelect tree in `DATA/fakeselect_FHC.root` has 594769 entries and all 19 dial TGraph branches (checked with ROOT).
   Example, entry 4: `MaCCQE` = (0.818, 0.706) (1.014, 0.872) (1.21, 1) (1.406, 1.097) (1.602, 1.170).
7. **`make_covariances.sh`**: `MakeCovariances.C` reads one `_Weight.flat.root` and writes
   `parameters/xsec/xsecCovHK.root`, which holds `xsec_param_names/priors/lower/upper/cov`. The covariance is diagonal.
   σ is recovered from the knots, and the bounds are the outer knots, widened to include the nominal (`MakeCovariances.C:169-274`).

**Event matching.**
- Weight tree ↔ recon tree: matched by **row position (entry index)**, because neither tree has an event ID. A mismatch in
  entry count is fatal (`MakeSplines.C:73-82,110-118`).
- After that, each spline is a branch on its own event row, so FakeSelect and GUNDAM need no further matching. GUNDAM reads
  the TGraph from the same tree it reads events from, and the dataset loader has no friend trees (`MakeSplines.C:21-35`).

## 3. How GUNDAM reads them

- Top level: `GLC/config.yaml` → `GLC/inputs/propagator/configPropagator.yaml` → `parameterSetList: ./inputs/parameters/configParSet.yaml`. That file includes `./inputs/parameters/xsec/configXsecParSet.yaml`.
- Dataset: `GLC/inputs/datasets/configDatasets.yaml` sets `tree: "FakeSelect"` and reads `${data_dir}/fakeselect_FHC.root` and `${data_dir}/fakeselect_RHC.root`.
- Parameters: `parameterDefinitionFilePath: ${data_dir}/parameters/xsec/xsecCovHK.root`, with `covarianceMatrix: xsec_param_cov`, `parameterNameList: xsec_param_names`, `parameterPriorValueList: xsec_param_priors`, and `parameterLower/UpperBoundsList: xsec_param_lower/upper`. **The parameters are in physical units, not σ units.**
- Each parameter has `dialType: Spline`, `treeExpression: "<DialName>"` (the TGraph leaf), `applyOnDataSets: ["HK"]`, and either `options: catmull-rom, uniformity(1e-6)` or, for `PionFSI_InelProb` only, `options: akima`. The yaml comment explains why: catmull-rom becomes a CompactSpline, which assumes evenly spaced knots, so non-uniform knots would give a silently wrong response.
- Clamping:
  - Parameters are limited to `[xsec_param_lower, xsec_param_upper]` = outer knots, applied as minimiser bounds (`GUNDAM/ParametersManager/src/ParameterSet.cpp:1036,1150`).
  - Weights: the xsec yaml sets no `minDialResponse`/`maxDialResponse`. GUNDAM defaults to a minimum response of 0 and no maximum (`GUNDAM/DialDictionary/DialEngine/src/DialCollection.cpp:463-464,967-970`).
  - No extrapolation should happen, because the bounds equal the outer knots.
  - Flux/det norms use `parametersRange: [0.0, 2.0]`.
- GUNDAM accepts only `TGraph` or `TSpline3` leaves. Anything else, such as `vector<double>`, silently gives no dial (`MLI/src/MakeSplines.C:21-35`, quoting `DialUtils.cpp:14-19`).

## 4. Implications for GENIE (suggestions)

For format compatibility, which is a hard requirement if the same GUNDAM machinery is reused:
- Write one **`TGraph` branch per GENIE knob**, split level 0, named exactly as the GUNDAM `treeExpression`. Put it on the
  **same per-event row** of the tree GUNDAM reads: either `FlatTree_VARS` in a `<sample>_spline.root` that then goes through
  `FakeSelect.C`, which copies any TGraph branch automatically, or directly in the final selected tree.
- Use x = **physical parameter value**, not the GENIE twk-dial value in σ units, because `xsecCovHK.root` priors, σ and
  bounds are physical. You could instead build a GENIE covariance in σ units (prior 0, σ 1), but then do not mix the two
  conventions in one file.
- Always include the nominal as a knot with weight exactly 1.
- Write an empty TGraph when there is no response (all weights 1) and when the calculation failed (all 0, negative or NaN).
  Never write an all-zero graph, because GUNDAM would turn it into a Shift = 0 dial.
- Keep knots **evenly spaced, including the nominal**, if you use `catmull-rom, uniformity(1e-6)`. Use `akima` for any
  dial whose knots are uneven, as is done for PionFSI_InelProb.
- Join GENIE reweight output to recon by entry index only if both trees are guaranteed to drop the same events. An explicit
  event index branch would be safer.
- Make bounds equal the outer knots, so GUNDAM never extrapolates.

Knot settings that match NEUT:
- Model-parameter dials (MA, axial form factors, background scale): **±1, ±2 σ plus nominal = 5 uniform knots**.
- FSI dials: **±0.5, ±1 σ plus nominal = 5 uniform knots**.
- In GENIE Reweight, set the twk dial to {−2, −1, 0, 1, 2} or {−1, −0.5, 0, 0.5, 1}, then convert to physical units with GENIE's own σ.
- Keep the lowest knot away from unphysical values, such as a negative FSI fraction. That problem is why NEUT uses −0.9σ for PionFSI_InelProb.

Possible NEUT → GENIE counterparts (**suggestions only; the models differ and the σ definitions are not equivalent**):

| NEUT dial | Suggested GENIE knob(s) | Comment |
|---|---|---|
| MaCCQE | `MaCCQE` | direct analogue (dipole axial mass); GENIE nominal/σ differ |
| MaRES | `MaCCRES` (and `MaNCRES`) | Rein-Sehgal vs NEUT Graczyk-Sobczyk RES |
| CA5RES | none directly; `NormCCRES` as a rough proxy | C5A form factor has no GENIE RS counterpart |
| BgSclRES | `NonRESBGvpCC1pi`, `NonRESBGvnCC1pi`, … | loose: NEUT I½ background scale vs GENIE non-resonant DIS-region knobs |
| (no NEUT dial in chain) | `NormCCMEC` / MEC shape knobs | the NEUT GUNDAM config has **no MEC norm**; a GENIE MEC norm would be an addition with no NEUT twin |
| NucleonFSI_* (Total/Elastic/Single/Double) | `MFP_N`, `FrElas_N` (if available), `FrCEx_N`, `FrInel_N`, `FrAbs_N`, `FrPiProd_N` | NEUT cascade vs GENIE hA; the categories don't map one-to-one |
| PionFSI_AbsProb / QE* / CEx* / InelProb | `FrAbs_pi`, `MFP_pi`, `FrCEx_pi`, `FrInel_pi`, `FrPiProd_pi` | NEUT splits QE/CEx by momentum; GENIE hA does not |
| SRCFrac_O | none obvious | depends on GENIE's nuclear model (e.g. LFG/SuSAv2 tune) |

## 5. Not found or not verified

- No file in MLI holds the dial config. It is the T2KReWeight CI file named in `MLI/scripts/weight.sh`. I did not trace the T2KReWeight source that converts σ to physical units; the σ values above are derived from `TweakMetaData`.
- `weight.sh` names its output `*_Weight.root`, but `splines.sh` and `make_covariances.sh` read `*_Weight.flat.root`, and only `.flat.root` files exist in `DATA/`. I assume T2KReWeightCI adds the `.flat` suffix, but did not verify this.
- `MDLSPiEj`, `NucleonFSI_PB`, `SRCFrac_C` and `RadCorr_rate` are reweighted and written as TGraph branches but are not GUNDAM parameters (`configXsecParSet.yaml`). I found no stated reason.
- I inspected only one weight file (`numu_x_numu_200K_FHC`) and only `fakeselect_FHC.root`. I assume the other samples share the same TweakMetaData, since they use the same toml.
- I did not check which T2KReWeight engine produces the all-zero events (`MakeSplines.C:262-268` says this is undiagnosed).
- I did not verify GENIE knob names against the installed GENIE Reweight version.
- I did not check the GUNDAM line numbers quoted in MLI/GLC comments against the GUNDAM source checkout, except `DialCollection.cpp` (response limits) and `ParameterSet.cpp` (bounds).
