#!/bin/bash
# Writes the parameter covariances, their binning files and the oscillation
# energy grid into $data_dir/parameters/. The cross-section sigmas are read from
# one weight file.

[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }
[ -n "$working_dir" ] || { echo "working_dir is not set"; exit 1; }
cd "$working_dir"

root -l -b -q "src/MakeCovariances.C(\"$data_dir/numu/numu_x_numu_200K_FHC_Weight.flat.root\",\"$data_dir/parameters\")"
root -l -b -q "src/MakeOscGrid.C(\"$data_dir/parameters/oscProb/grid_HK.root\")"
