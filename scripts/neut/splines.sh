#!/bin/bash
# Cross-section response splines, one file per sample. MakeSplines.C reads the
# fakerecon tree and the T2KReWeight flat tree and writes the graphs alongside.

[ -n "$working_dir" ] || { echo "working_dir is not set"; exit 1; }
cd "$working_dir"
[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }

root -l -b -q "src/neut/MakeSplines.C(\"$data_dir/numu/numu_x_numu_200K_FHC_fakerecon.root\",\"$data_dir/numu/numu_x_numu_200K_FHC_Weight.flat.root\",\"$data_dir/numu/numu_x_numu_200K_FHC_spline.root\")"
root -l -b -q "src/neut/MakeSplines.C(\"$data_dir/numu/numu_x_numu_50K_RHC_fakerecon.root\",\"$data_dir/numu/numu_x_numu_50K_RHC_Weight.flat.root\",\"$data_dir/numu/numu_x_numu_50K_RHC_spline.root\")"
root -l -b -q "src/neut/MakeSplines.C(\"$data_dir/numubar/numubar_x_numubar_200K_RHC_fakerecon.root\",\"$data_dir/numubar/numubar_x_numubar_200K_RHC_Weight.flat.root\",\"$data_dir/numubar/numubar_x_numubar_200K_RHC_spline.root\")"
root -l -b -q "src/neut/MakeSplines.C(\"$data_dir/numubar/numubar_x_numubar_50K_FHC_fakerecon.root\",\"$data_dir/numubar/numubar_x_numubar_50K_FHC_Weight.flat.root\",\"$data_dir/numubar/numubar_x_numubar_50K_FHC_spline.root\")"
root -l -b -q "src/neut/MakeSplines.C(\"$data_dir/nue/nue_x_nue_200K_FHC_fakerecon.root\",\"$data_dir/nue/nue_x_nue_200K_FHC_Weight.flat.root\",\"$data_dir/nue/nue_x_nue_200K_FHC_spline.root\")"
root -l -b -q "src/neut/MakeSplines.C(\"$data_dir/nue/nue_x_nue_50K_RHC_fakerecon.root\",\"$data_dir/nue/nue_x_nue_50K_RHC_Weight.flat.root\",\"$data_dir/nue/nue_x_nue_50K_RHC_spline.root\")"
root -l -b -q "src/neut/MakeSplines.C(\"$data_dir/nuebar/nuebar_x_nuebar_200K_RHC_fakerecon.root\",\"$data_dir/nuebar/nuebar_x_nuebar_200K_RHC_Weight.flat.root\",\"$data_dir/nuebar/nuebar_x_nuebar_200K_RHC_spline.root\")"
root -l -b -q "src/neut/MakeSplines.C(\"$data_dir/nuebar/nuebar_x_nuebar_50K_FHC_fakerecon.root\",\"$data_dir/nuebar/nuebar_x_nuebar_50K_FHC_Weight.flat.root\",\"$data_dir/nuebar/nuebar_x_nuebar_50K_FHC_spline.root\")"
root -l -b -q "src/neut/MakeSplines.C(\"$data_dir/numu/numu_x_nue_50K_FHC_fakerecon.root\",\"$data_dir/numu/numu_x_nue_50K_FHC_Weight.flat.root\",\"$data_dir/numu/numu_x_nue_50K_FHC_spline.root\")"
root -l -b -q "src/neut/MakeSplines.C(\"$data_dir/numu/numu_x_nue_50K_RHC_fakerecon.root\",\"$data_dir/numu/numu_x_nue_50K_RHC_Weight.flat.root\",\"$data_dir/numu/numu_x_nue_50K_RHC_spline.root\")"
root -l -b -q "src/neut/MakeSplines.C(\"$data_dir/numubar/numubar_x_nuebar_50K_RHC_fakerecon.root\",\"$data_dir/numubar/numubar_x_nuebar_50K_RHC_Weight.flat.root\",\"$data_dir/numubar/numubar_x_nuebar_50K_RHC_spline.root\")"
root -l -b -q "src/neut/MakeSplines.C(\"$data_dir/numubar/numubar_x_nuebar_50K_FHC_fakerecon.root\",\"$data_dir/numubar/numubar_x_nuebar_50K_FHC_Weight.flat.root\",\"$data_dir/numubar/numubar_x_nuebar_50K_FHC_spline.root\")"
