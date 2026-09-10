# mock_lbl_inputs

This repository contains the MC input generators needed to run a fake MC sensitivity study.

    scripts/        scripts to run to generate the MC
    src/            ROOT macros for simulated selection/weights
    gen/            generator cards, one per sample for generating neutrino events
    lbl_config      smearing/selection and beam flux files for different detectors


Before running, export where this repository is and where the generated files go:

    export working_dir=/path/to/mock_lbl_inputs
    export data_dir=/path/to/output_data




## Generating the input

    scripts/run_all.sh                             every stage below, in order, stopping at the first failure

or one at a time:

    scripts/init_data_dir.sh                      makes sub-directories under $data_dir

    scripts/gen.sh                                 NEUT events               numu/*.root ...
    scripts/weight.sh                              T2KReWeight spline trees  *_Weight.flat.root
    scripts/flat.sh                                flatten NEUT outputs      *_flattree.root
    scripts/recon.sh                               FakeRecon.C               *_fakerecon.root
    scripts/splines.sh                             MakeSplines.C             *_spline.root
    scripts/select.sh                              FakeSelect.C              fakeselect_FHC.root, fakeselect_RHC.root
    scripts/make_covariances.sh                    MakeCovariances.C         parameters/{oscProb,det,flux,xsec}/*.root, *Binning.txt
                                                   MakeOscGrid.C             parameters/oscProb/grid_HK.root

MakeCovariances.C writes the oscillation, detector, cross-section and flux
covariances; it reads one *_Weight.flat.root for the cross-section sigmas and
the T2K public flux covariance for the flux block. The gundam configs read all
of these as ${data_dir}/parameters/..., which gundamFitter expands.


