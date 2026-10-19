# mock_lbl_inputs

This repository contains the MC input generators needed to run a fake MC sensitivity study.

## Layout

    scripts/        scripts to run to generate the MC
    src/            ROOT macros for simulated selection/weights; NuFast_LBL.cpp is
                    the NuFast oscillation code (P. Denton, MIT licence), used by FakeRecon.C
    gen/            generator cards for generating neutrino events: gen/neut/ holds
                    the NEUT cards, one per sample; gen/genie/ holds the GENIE sample
                    manifest and message thresholds
    lbl_config      smearing/selection and beam flux files for different detectors

Stages that every generator uses sit at the top level of `scripts/` and `src/`
(`init_data_dir.sh`, `recon.sh`, `select.sh`; `FakeRecon.C`, `FakeSelect.C`,
`MakeCovariances.C`, `MakeOscGrid.C`, `NuFast_LBL.cpp`). Generator-specific
files live in a subfolder named after the generator: `scripts/neut/`,
`src/neut/` and `gen/neut/` for NEUT, `scripts/genie/`, `src/genie/` and
`gen/genie/` for GENIE. A new generator gets its own `scripts/<generator>/`,
`src/<generator>/` and `gen/<generator>/`.

## NEUT

Before running, export where this repository is and where the generated files go:

    export working_dir=/path/to/mock_lbl_inputs
    export data_dir=/path/to/output_data




### Generating the input

    scripts/neut/run_all.sh                        every stage below, in order, stopping at the first failure

or one at a time:

    scripts/init_data_dir.sh                      makes sub-directories under $data_dir

    scripts/neut/gen.sh                            NEUT events               numu/*.root ...
    scripts/neut/weight.sh                         T2KReWeight spline trees  *_Weight.flat.root
    scripts/neut/flat.sh                           flatten NEUT outputs      *_flattree.root
    scripts/recon.sh                               FakeRecon.C               *_fakerecon.root
    scripts/neut/splines.sh                        MakeSplines.C             *_spline.root
    scripts/select.sh                              FakeSelect.C              fakeselect_FHC.root, fakeselect_RHC.root
    scripts/neut/make_covariances.sh               MakeCovariances.C         parameters/{oscProb,det,flux,xsec}/*.root, *Binning.txt
                                                   MakeOscGrid.C             parameters/oscProb/grid_HK.root

MakeCovariances.C writes the oscillation, detector, cross-section and flux
covariances; it reads one *_Weight.flat.root for the cross-section sigmas and
the T2K public flux covariance for the flux block. The gundam configs read all
of these as ${data_dir}/parameters/..., which gundamFitter expands.




## GENIE

The same inputs can be made with GENIE (R-3_06_02, tune G18_10b_02_11b, GENIE
Reweight R-1_04_02) instead of NEUT and T2KReWeight. The GENIE stages are in
`scripts/genie/`, with the samples listed in `gen/genie/samples.txt`. They
write the same `<stem>_{flattree,fakerecon,spline}.root` files as NEUT, so
`scripts/init_data_dir.sh`, `scripts/recon.sh`, `scripts/select.sh` and the
macros `src/FakeRecon.C` (with `src/NuFast_LBL.cpp`), `src/FakeSelect.C`,
`src/MakeCovariances.C` (oscillation, detector and flux parts) and
`src/MakeOscGrid.C` are shared. Because the file names are the same, use a
`$data_dir` of its own for GENIE.

    scripts/genie/run_all.sh                       the whole GENIE chain, in order
    src/genie/          GENIE spline producer and xsec covariance macro

### Software

| package | version |
|---|---|
| GENIE Generator | R-3_06_02 |
| GENIE Reweight | R-1_04_02 |
| GENIE tune | G18_10b_02_11b (stock, no `--xml-path` / `GXMLPATH`) |
| GUNDAM | 2.1.1-110-g4446b740 (commit 4446b740, 2026-06-12) |
| ROOT | 6.30.04 |
| NUISANCE (built with GENIE) | 2.9.9 |
| LHAPDF | 5.9.1 |
| PYTHIA6 | 6.4.28 |
| log4cpp | 1.1.6 |

### Environment

```bash
export working_dir=/path/to/mock_lbl_inputs
export data_dir=/path/to/genie_output_data      # a data_dir of its own, not the NEUT one
```

The scripts do not set up software: source your GENIE and GENIE Reweight
environment first. `src/genie/genie_splines.cxx` is compiled against your local
ROOT, GENIE and GENIE Reweight (the libraries of GENIE Reweight's `grwght1p`);
put the executable at `$working_dir/build/genie_splines`, where
`scripts/genie/splines.sh` expects it.

### Pipeline

| step | script | tool | output in `$data_dir` |
|---|---|---|---|
| directories | `scripts/init_data_dir.sh` (shared) | mkdir | `numu/ ... parameters/` |
| cross sections | `scripts/genie/xsec.sh [flavour\|merge]` | `gmkspl`, `gspladd` | `genie/xsec/hk_v1_stock.xml` |
| events | `scripts/genie/gen.sh [sample]` | `gevgen` + `gen/genie/Messenger_production.xml` | `genie/<sample>/<sample>.ghep.root`, `input-flux.root` |
| flat trees | `scripts/genie/flat.sh [sample]` | `PrepareGENIE`, `nuisflat` | `<stem>_flattree.root` |
| fake recon | `scripts/recon.sh` (shared) | `src/FakeRecon.C` | `<stem>_fakerecon.root` |
| splines | `scripts/genie/splines.sh [sample]` | `src/genie/genie_splines.cxx` | `<stem>_spline.root` (+ `GenieDialMeta` tree) |
| selection | `scripts/select.sh` (shared) | `src/FakeSelect.C` | `fakeselect_{FHC,RHC}.root` |
| parameters | `scripts/genie/make_covariances.sh [all\|FHC\|RHC]` | `src/MakeCovariances.C` (osc/det/flux), `src/MakeOscGrid.C`, `src/genie/MakeGenieXsecCov.C` | `parameters/{oscProb,det,flux}/`, `parameters/xsec/xsecCovGENIE.root` |

Every GENIE script except `xsec.sh` reads the samples from the manifest
`gen/genie/samples.txt`.
`<stem>` is the NEUT name (e.g. `numu/numu_x_numu_200K_FHC`); `[sample]` is a
manifest name or 0-based row index, and without it all 12 samples run.

### Running on a cluster

The GENIE scripts have no SLURM settings; add your own. `sbatch` passes the environment.

```bash
cd $working_dir
for f in numu numubar nue nuebar; do       # cross sections: about a day per flavour
  sbatch --job-name=genie_xsec_$f --wrap="scripts/genie/xsec.sh $f"
done
scripts/genie/xsec.sh merge                # once the four jobs are done
sbatch --array=0-11 --wrap='scripts/genie/gen.sh $SLURM_ARRAY_TASK_ID && scripts/genie/flat.sh $SLURM_ARRAY_TASK_ID'
scripts/recon.sh
sbatch --array=0-11 --mem=8G --wrap='scripts/genie/splines.sh $SLURM_ARRAY_TASK_ID'
scripts/select.sh
scripts/genie/make_covariances.sh
```

### Model

- 9 spline dials (MaCCQE, CCQEPauliSupViaKF, CoulombCCQE, Ma/MvCCRES,
  Ma/MvNCRES, MaNCEL, EtaNCEL), each with 5 knots spanning −2 to +2σ.
- 21 normalisation dials (NormCCQE, NormCC/NCMEC, NormCC/NCCOH, 16 NonRESBG),
  each a linear graph y = x through −1, 0, +1σ.
- σ are the GENIE Reweight defaults, which are not any tune's systematic error.
- The prior is diagonal, in `parameters/xsec/xsecCovGENIE.root`.
- Parameter and branch names carry a `GENIE_` prefix, so a NEUT config cannot
  silently bind NEUT priors to GENIE graphs.

The decision log is in the header of `src/genie/genie_splines.cxx`. The GUNDAM
configuration is the GENIE option of haowei-physics/gundam_lbl_config
(`inputs/parameters/useGENIE.yaml`).

### Known limitations and open items

- The MaCCQE knot range is still under discussion.
- Interpolation was checked only on the manifest sample numu_x_numu_fhc; the antineutrino manifest samples (numubar_x_*, nuebar_x_*) in particular have not been checked.
- CCQEPauliSupViaKF has a step response, so a 5-knot interpolation cannot follow it exactly between knots.
