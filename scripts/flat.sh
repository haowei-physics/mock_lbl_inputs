#!/bin/bash
# Flatten to GenericFlux flattrees. Needs nuisflat on PATH.
[ -n "$working_dir" ] || { echo "working_dir is not set"; exit 1; }
cd "$working_dir"
[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }

nuisflat -i NEUT:$data_dir/numu/numu_x_numu_200K_FHC.root          -o $data_dir/numu/numu_x_numu_200K_FHC_flattree.root
nuisflat -i NEUT:$data_dir/numu/numu_x_numu_50K_RHC.root           -o $data_dir/numu/numu_x_numu_50K_RHC_flattree.root
nuisflat -i NEUT:$data_dir/numubar/numubar_x_numubar_200K_RHC.root -o $data_dir/numubar/numubar_x_numubar_200K_RHC_flattree.root
nuisflat -i NEUT:$data_dir/numubar/numubar_x_numubar_50K_FHC.root  -o $data_dir/numubar/numubar_x_numubar_50K_FHC_flattree.root
nuisflat -i NEUT:$data_dir/nue/nue_x_nue_200K_FHC.root             -o $data_dir/nue/nue_x_nue_200K_FHC_flattree.root
nuisflat -i NEUT:$data_dir/nue/nue_x_nue_50K_RHC.root              -o $data_dir/nue/nue_x_nue_50K_RHC_flattree.root
nuisflat -i NEUT:$data_dir/nuebar/nuebar_x_nuebar_200K_RHC.root    -o $data_dir/nuebar/nuebar_x_nuebar_200K_RHC_flattree.root
nuisflat -i NEUT:$data_dir/nuebar/nuebar_x_nuebar_50K_FHC.root     -o $data_dir/nuebar/nuebar_x_nuebar_50K_FHC_flattree.root
nuisflat -i NEUT:$data_dir/numu/numu_x_nue_50K_FHC.root            -o $data_dir/numu/numu_x_nue_50K_FHC_flattree.root
nuisflat -i NEUT:$data_dir/numu/numu_x_nue_50K_RHC.root            -o $data_dir/numu/numu_x_nue_50K_RHC_flattree.root
nuisflat -i NEUT:$data_dir/numubar/numubar_x_nuebar_50K_RHC.root   -o $data_dir/numubar/numubar_x_nuebar_50K_RHC_flattree.root
nuisflat -i NEUT:$data_dir/numubar/numubar_x_nuebar_50K_FHC.root   -o $data_dir/numubar/numubar_x_nuebar_50K_FHC_flattree.root
