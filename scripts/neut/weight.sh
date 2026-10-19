#!/bin/bash

#SBATCH --job-name=HL_MC
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --array=0
#SBATCH --output=T2KRW_%A.out
#SBATCH --mem-per-cpu=4G
[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }
[ -n "$working_dir" ] || { echo "working_dir is not set"; exit 1; }
cd "$working_dir"

source "$HOME/setup_neut.sh"

C=${T2KReWeight_ROOT}/var/T2KReWeight/config/CI/neut_dial_tweaks.toml

T2KReWeightCI --no-plot -c $C --neut-card gen/neut/numu/NEUT.numu.200k.toml                    -i $data_dir/numu/numu_x_numu_200K_FHC.root          -o $data_dir/numu/numu_x_numu_200K_FHC_Weight.root
T2KReWeightCI --no-plot -c $C --neut-card gen/neut/numu/NEUT.numu_x_numu_50K_RHC.toml          -i $data_dir/numu/numu_x_numu_50K_RHC.root           -o $data_dir/numu/numu_x_numu_50K_RHC_Weight.root
T2KReWeightCI --no-plot -c $C --neut-card gen/neut/numubar/NEUT.numubar.200k.toml              -i $data_dir/numubar/numubar_x_numubar_200K_RHC.root -o $data_dir/numubar/numubar_x_numubar_200K_RHC_Weight.root
T2KReWeightCI --no-plot -c $C --neut-card gen/neut/numubar/NEUT.numubar_x_numubar_50K_FHC.toml -i $data_dir/numubar/numubar_x_numubar_50K_FHC.root  -o $data_dir/numubar/numubar_x_numubar_50K_FHC_Weight.root
T2KReWeightCI --no-plot -c $C --neut-card gen/neut/nue/NEUT.nue.200k.toml                      -i $data_dir/nue/nue_x_nue_200K_FHC.root             -o $data_dir/nue/nue_x_nue_200K_FHC_Weight.root
T2KReWeightCI --no-plot -c $C --neut-card gen/neut/nue/NEUT.nue_x_nue_50K_RHC.toml             -i $data_dir/nue/nue_x_nue_50K_RHC.root              -o $data_dir/nue/nue_x_nue_50K_RHC_Weight.root
T2KReWeightCI --no-plot -c $C --neut-card gen/neut/nuebar/NEUT.nuebar.200k.toml                -i $data_dir/nuebar/nuebar_x_nuebar_200K_RHC.root    -o $data_dir/nuebar/nuebar_x_nuebar_200K_RHC_Weight.root
T2KReWeightCI --no-plot -c $C --neut-card gen/neut/nuebar/NEUT.nuebar.50k.FHC.toml             -i $data_dir/nuebar/nuebar_x_nuebar_50K_FHC.root     -o $data_dir/nuebar/nuebar_x_nuebar_50K_FHC_Weight.root
T2KReWeightCI --no-plot -c $C --neut-card gen/neut/numu/NEUT.numu_x_nue_50K_FHC.toml           -i $data_dir/numu/numu_x_nue_50K_FHC.root            -o $data_dir/numu/numu_x_nue_50K_FHC_Weight.root
T2KReWeightCI --no-plot -c $C --neut-card gen/neut/numu/NEUT.numu_x_nue_50K_RHC.toml           -i $data_dir/numu/numu_x_nue_50K_RHC.root            -o $data_dir/numu/numu_x_nue_50K_RHC_Weight.root
T2KReWeightCI --no-plot -c $C --neut-card gen/neut/numubar/NEUT.numubar_x_nuebar_50K_RHC.toml  -i $data_dir/numubar/numubar_x_nuebar_50K_RHC.root   -o $data_dir/numubar/numubar_x_nuebar_50K_RHC_Weight.root
T2KReWeightCI --no-plot -c $C --neut-card gen/neut/numubar/NEUT.numubar_x_nuebar_50K_FHC.toml  -i $data_dir/numubar/numubar_x_nuebar_50K_FHC.root   -o $data_dir/numubar/numubar_x_nuebar_50K_FHC_Weight.root
