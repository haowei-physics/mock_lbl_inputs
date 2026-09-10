#!/bin/bash
# NEUT generation. Needs neutroot2 on PATH.
[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }
[ -n "$working_dir" ] || { echo "working_dir is not set"; exit 1; }
cd "$working_dir"

neutroot2 gen/neut/numu/NEUT.numu.200k.toml                    $data_dir/numu/numu_x_numu_200K_FHC.root
neutroot2 gen/neut/numu/NEUT.numu_x_numu_50K_RHC.toml          $data_dir/numu/numu_x_numu_50K_RHC.root
neutroot2 gen/neut/numubar/NEUT.numubar.200k.toml              $data_dir/numubar/numubar_x_numubar_200K_RHC.root
neutroot2 gen/neut/numubar/NEUT.numubar_x_numubar_50K_FHC.toml $data_dir/numubar/numubar_x_numubar_50K_FHC.root
neutroot2 gen/neut/nue/NEUT.nue.200k.toml                      $data_dir/nue/nue_x_nue_200K_FHC.root
neutroot2 gen/neut/nue/NEUT.nue_x_nue_50K_RHC.toml             $data_dir/nue/nue_x_nue_50K_RHC.root
neutroot2 gen/neut/nuebar/NEUT.nuebar.200k.toml                $data_dir/nuebar/nuebar_x_nuebar_200K_RHC.root
neutroot2 gen/neut/nuebar/NEUT.nuebar.50k.FHC.toml             $data_dir/nuebar/nuebar_x_nuebar_50K_FHC.root
neutroot2 gen/neut/numu/NEUT.numu_x_nue_50K_FHC.toml           $data_dir/numu/numu_x_nue_50K_FHC.root
neutroot2 gen/neut/numu/NEUT.numu_x_nue_50K_RHC.toml           $data_dir/numu/numu_x_nue_50K_RHC.root
neutroot2 gen/neut/numubar/NEUT.numubar_x_nuebar_50K_RHC.toml  $data_dir/numubar/numubar_x_nuebar_50K_RHC.root
neutroot2 gen/neut/numubar/NEUT.numubar_x_nuebar_50K_FHC.toml  $data_dir/numubar/numubar_x_nuebar_50K_FHC.root
