# GENIE path

Per-event GENIE cross-section response graphs ("splines") for the stock
G18_10b_02_11b HK samples, produced in the same layout as the NEUT chain so
that the shared `scripts/recon.sh`, `scripts/select.sh`, `src/FakeRecon.C` and
`src/FakeSelect.C` run on them unchanged and the output plugs into the same
GUNDAM machinery.

## Software

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

The splines in production were made with these versions on 2026-10-06.

## Environment

As for NEUT, every script checks two variables and works relative to them:

```bash
export working_dir=/path/to/mock_lbl_inputs
export data_dir=/path/to/genie_output_data      # a data_dir of its own, not the NEUT one
```

The scripts do not set up software. Before running, source the environment of
your GENIE and GENIE Reweight installation, so that `root`, `gmkspl`, `gspladd`,
`gevgen`, `PrepareGENIE` and `nuisflat` are on `PATH` and the GENIE and GENIE
Reweight libraries are found at run time.

## Pipeline

The samples are listed once, in `gen/genie/samples.txt`: GENIE name, NEUT
stem, probe PDG, flux file and histogram, number of events, gevgen run number
and seed, POT and flux flavour. Every script and program reads it from there.
Within each beam mode the row order is FakeSelect.C's file order.

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

`src/genie/genie_splines.cxx` is a compiled program. Build it against your
local ROOT, GENIE and GENIE Reweight installation (the same libraries as GENIE
Reweight's `grwght1p`) and place the executable at
`$working_dir/build/genie_splines`, where `scripts/genie/splines.sh` expects it.

`<stem>` is the NEUT name (e.g. `numu/numu_x_numu_200K_FHC`). `[sample]` is a
name or 0-based row index of the manifest; without it a script runs all 12
samples in turn. Events are matched by entry: the producer applies FakeRecon's
drop rules and checks every row against the fakerecon tree (Enu_true, PDGnu,
cc), and any mismatch is fatal. The GENIE programs run in `$data_dir/genie/run/`,
because the calculators may write debug files into their working directory.

`scripts/genie/run_all.sh` runs every step in order. `root -q` exits with the
macro's return value, and FakeRecon.C and FakeSelect.C return 1 on success, so
the shared `recon.sh` and `select.sh` exit with 1 when they succeed.
`run_all.sh` therefore removes their outputs first, ignores their exit status
and checks that every output exists afterwards. The GENIE scripts themselves
treat a return value of 1 as success.

The GUNDAM configuration for the GENIE parameters is not in this repository
yet (its location is still to be agreed).

## Running on a cluster

The stage scripts have no SLURM settings. Use one job per flavour or sample,
for example (add your partition, memory and time options):

```bash
cd $working_dir
# cross sections: one flavour per job (about a day each), then the merge
for f in numu numubar nue nuebar; do
  sbatch --job-name=genie_xsec_$f --output=xsec_$f.log --wrap="scripts/genie/xsec.sh $f"
done
scripts/genie/xsec.sh merge

# events, flat trees, splines: one array task per manifest row
sbatch --array=0-11 --job-name=genie_gen --output=gen_%a.log \
       --wrap='scripts/genie/gen.sh $SLURM_ARRAY_TASK_ID && scripts/genie/flat.sh $SLURM_ARRAY_TASK_ID'
scripts/recon.sh
sbatch --array=0-11 --mem=8G --job-name=genie_spl --output=spl_%a.log \
       --wrap='scripts/genie/splines.sh $SLURM_ARRAY_TASK_ID'
scripts/select.sh
scripts/genie/make_covariances.sh
```

`sbatch` passes the environment (`working_dir`, `data_dir`, `PATH`, ...) to
the job by default. A spline task takes about 1.5 min for a 50k sample and
2.5 min for a 200k sample (about 1 min of it is GENIE initialisation).

## Main decisions

The full decision log is in the header of `src/genie/genie_splines.cxx`; the
dial study behind it is in [report.md](report.md) and
[spline_plan.md](spline_plan.md).

- **GENIE's own settings.** Central values are the stock tune, read from the
  model registry at run time. Knots are at −2, −1, 0, +1, +2σ with the GENIE
  Reweight default σ, converted to physical values (`x = nominal·(1 + nσ)`).
  The prior is diagonal, with the same σ. The GENIE Reweight README says these
  defaults are not any tune's systematic error, so they are a starting point
  ([genie_default_dials.md](genie_default_dials.md)).
- **9 spline dials:** MaCCQE, CCQEPauliSupViaKF, CoulombCCQE, MaCCRES, MvCCRES,
  MaNCRES, MvNCRES, MaNCEL, EtaNCEL. GUNDAM reads them as
  `Spline, catmull-rom, uniformity(1e-6)` with bounds at the outer knots.
- **CCQEPauliSupViaKF also stays plain catmull-rom**, after comparing three
  interpolations (see "Known limitation").
  - Its per-event response is a step: `GReWeightFGM` returns exactly 0 or 1 for
    QELEventGenerator events. So no 5-knot interpolation is exact between
    knots, but all are exact at the knots.
  - **Flagged for mentor review:** the dial's physical relevance is uncertain.
    Nieves + LFG applies its own Pauli blocking in the cross section, while
    this dial takes k_F from GENIE's generic `PauliBlocker`.
- **21 normalisations:** NormCCQE, NormCCMEC, NormNCMEC, NormCCCOH, NormNCCOH
  and the 16 NonRESBG dials.
  - FakeSelect has no GENIE interaction-mode branch, so each affected event
    carries a graph y = x through −1, 0, +1σ, taken from GENIE's own calculator.
    GUNDAM reads it as a linear `Graph` dial.
  - Bounds are [0, 1+3σ]. The NormCCQE prior σ is 0.20 (GENIE gives +0.20/−0.15).
- **Dropped:**
  - Bodek-Yang and MaCOHpi/R0COHpi are broken for this tune.
  - AGKY, RDecBR and FracPN/FracDelta only change the final state, which
    FakeSelect cannot see.
  - MaCCQEshape and the RES shape dials are alternative parametrisations.
  - NormCCRES/NormNCRES would double count MaRES.
- **Per-event rules as in `src/neut/MakeSplines.C`:**
  - the nominal knot has weight exactly 1;
  - all weights equal to 1 gives an empty graph;
  - NaN, negative or all-zero weights give an empty graph, counted as rejected
    (only MaCCQE has any: 1–12 ν̄ CCQE events per ν̄-containing sample);
  - TGraph branches have split level 0.
- **Names** carry a `GENIE_` prefix, so a NEUT config cannot silently bind
  NEUT priors to GENIE graphs.
- **Knot settings** follow the NEUT chain's format
  ([neut_knot_settings.md](neut_knot_settings.md)); the MaCCQE knot range is
  still under discussion (see below).
- **Workarounds, with nothing fixed upstream:**
  - `root -q` exit codes, handled inside `scripts/genie` (see "Pipeline").
  - The "neutlike" config variant (not in this repository) sets EL-Axial-Eta
    = 0. This chain uses the stock config, and the producer stops unless
    EtaNCEL = 0.12.
  - `--one-sigma` stops on the fixed PMNS_SIGN_MASS_SQUARED_32, so the one-sigma
    check runs with `override_noOsc.yaml`, which disables the oscillation set.

## Checks done on the production (all passed unless noted)

| check | result |
|---|---|
| File level (all 12 samples): row match, TGraph split 0, knots, weight 1 at the nominal, no NaN or negative weights, norm y = x | passed |
| File level: ±1σ/±2σ per-sample rates vs the dial scan | agree to ≤ 4×10⁻¹¹ %, except MaCCQE in samples with rejected events (≤ 0.1%) |
| FakeSelect FHC and RHC | sample draw identical to the no-dial FakeSelect; graphs identical to their inputs |
| osc/det/flux parameter files | identical to the NEUT ones |
| GUNDAM dry run: dial counts per parameter and sample vs an independent recount from the FakeSelect output | 30 × 6 columns, 0 differences |
| GUNDAM +1σ per Eν bin vs the same variation computed from the graphs in the FakeSelect output | 1650 values, agree to ≤ 6×10⁻¹⁰ % |
| Interpolation at ±0.5σ and ±1.5σ: GUNDAM's interpolation vs GENIE computed directly (numu_x_numu_fhc) | within about 1% of each dial's effect for the catmull-rom dials; CCQEPauliSupViaKF see below |

These were run with the development scripts this code was restructured from.
The restructured spline producer was then checked against them: on
numu_x_nue_50K_FHC its output is bitwise identical to the production file
(every `GENIE_*` graph and every `GenieDialMeta` column).

**Known limitation: CCQEPauliSupViaKF interpolation.** A blocked event has
knots (1, 1, 1, 0, 0). Three GUNDAM interpolations were compared with GENIE
on numu_x_numu_fhc. The table gives the largest rate difference over the 3 FHC
samples, in % of the rate:

| twk | catmull-rom (in use) | linear (Graph) | catmull-rom, monotonic |
|---|---|---|---|
| −1.5σ | 2e-5 | 5e-5 | 5e-5 |
| −0.5σ | 0.177 | 1.2e-4 | 9e-6 |
| +0.5σ | 0.018 | 0.266 | 0.266 |
| +1.5σ | 0.068 | 0.211 | 0.346 |
| largest | **0.177** | 0.266 | 0.346 |

The rule agreed beforehand was to use monotonic only if it beat both on the
largest gap, and otherwise keep plain catmull-rom. Catmull-rom's error is an
overshoot below nominal (+0.18% in FHC 1Rμ at −0.5σ, where GENIE gives ~0).
Linear and monotonic remove that, but they underestimate the curvature above
nominal, where the dial acts. CoulombCCQE is mildly affected the same way:
about 2% of its small effect at +0.5σ.

The interpolation check has been done for numu_x_numu_fhc only; the
antineutrino samples have not been checked.

## Open items

- MaCCQE knot range: still under discussion (the GENIE prior is 3%).
- Interpolation check for the antineutrino samples: not done.
- GUNDAM config for the GENIE parameters: location to be agreed.

## Files

- `gen/genie/samples.txt`: the sample manifest; `gen/genie/Messenger_production.xml`: GENIE message thresholds.
- `scripts/genie/`: the stage scripts above.
- `src/genie/genie_splines.cxx`: producer, with the decision log in its header.
- `src/genie/MakeGenieXsecCov.C`: writes `xsecCovGENIE.root` (same keys as `xsecCovHK.root`).
- `docs/genie/`: this file; [report.md](report.md) (dial scan and dial choice), [spline_plan.md](spline_plan.md) (the plan before production), [neut_knot_settings.md](neut_knot_settings.md) (the NEUT chain's knots and spline format), [genie_default_dials.md](genie_default_dials.md) (GENIE Reweight's default dials and σ). The last two record the paths of the development area where they were examined.
