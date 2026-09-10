#!/bin/bash
# FakeRecon. Needs root on PATH. Run MakeResponse.C first.
#
# Arguments are: flattree, output, smearing, response, seed, flux file, flux hist,
# POT, nucleons, flux PDG.
#   POT       6.75e21 for FHC, 20.25e21 for RHC - the 1:3 split of 27e21 total.
#   nucleons  1.132e35, the 188 kt HK fiducial volume of water.
#   flux PDG  the flavour of the *flux*. For the appearance samples that is the
#             parent, not the interaction flavour, because an appeared nue and
#             an intrinsic beam nue both have PDGnu 12.
#   flux file the uniform-binned copies in flux/. NEUT samples a
#             histogram by bin content, so variable-width bins must not be used.
[ -n "$working_dir" ] || { echo "working_dir is not set"; exit 1; }
cd "$working_dir"
[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }
P=lbl_config/hk/t2kflux_2020_plus250kA_uniform.root
M=lbl_config/hk/t2kflux_2020_minus250kA_uniform.root

root -l -b -q "src/FakeRecon.C(\"$data_dir/numu/numu_x_numu_200K_FHC_flattree.root\",\"$data_dir/numu/numu_x_numu_200K_FHC_fakerecon.root\",\"lbl_config/hk/smearing.root\",\"lbl_config/hk/response.root\",1234,\"$P\",\"enu_sk_numu\",6.75e21,1.132e35,14)"
root -l -b -q "src/FakeRecon.C(\"$data_dir/numu/numu_x_numu_50K_RHC_flattree.root\",\"$data_dir/numu/numu_x_numu_50K_RHC_fakerecon.root\",\"lbl_config/hk/smearing.root\",\"lbl_config/hk/response.root\",1234,\"$M\",\"enu_sk_numu\",2.025e22,1.132e35,14)"
root -l -b -q "src/FakeRecon.C(\"$data_dir/numubar/numubar_x_numubar_200K_RHC_flattree.root\",\"$data_dir/numubar/numubar_x_numubar_200K_RHC_fakerecon.root\",\"lbl_config/hk/smearing.root\",\"lbl_config/hk/response.root\",1234,\"$M\",\"enu_sk_numub\",2.025e22,1.132e35,-14)"
root -l -b -q "src/FakeRecon.C(\"$data_dir/numubar/numubar_x_numubar_50K_FHC_flattree.root\",\"$data_dir/numubar/numubar_x_numubar_50K_FHC_fakerecon.root\",\"lbl_config/hk/smearing.root\",\"lbl_config/hk/response.root\",1234,\"$P\",\"enu_sk_numub\",6.75e21,1.132e35,-14)"
root -l -b -q "src/FakeRecon.C(\"$data_dir/nue/nue_x_nue_200K_FHC_flattree.root\",\"$data_dir/nue/nue_x_nue_200K_FHC_fakerecon.root\",\"lbl_config/hk/smearing.root\",\"lbl_config/hk/response.root\",1234,\"$P\",\"enu_sk_nue\",6.75e21,1.132e35,12)"
root -l -b -q "src/FakeRecon.C(\"$data_dir/nue/nue_x_nue_50K_RHC_flattree.root\",\"$data_dir/nue/nue_x_nue_50K_RHC_fakerecon.root\",\"lbl_config/hk/smearing.root\",\"lbl_config/hk/response.root\",1234,\"$M\",\"enu_sk_nue\",2.025e22,1.132e35,12)"
root -l -b -q "src/FakeRecon.C(\"$data_dir/nuebar/nuebar_x_nuebar_200K_RHC_flattree.root\",\"$data_dir/nuebar/nuebar_x_nuebar_200K_RHC_fakerecon.root\",\"lbl_config/hk/smearing.root\",\"lbl_config/hk/response.root\",1234,\"$M\",\"enu_sk_nueb\",2.025e22,1.132e35,-12)"
root -l -b -q "src/FakeRecon.C(\"$data_dir/nuebar/nuebar_x_nuebar_50K_FHC_flattree.root\",\"$data_dir/nuebar/nuebar_x_nuebar_50K_FHC_fakerecon.root\",\"lbl_config/hk/smearing.root\",\"lbl_config/hk/response.root\",1234,\"$P\",\"enu_sk_nueb\",6.75e21,1.132e35,-12)"
root -l -b -q "src/FakeRecon.C(\"$data_dir/numu/numu_x_nue_50K_FHC_flattree.root\",\"$data_dir/numu/numu_x_nue_50K_FHC_fakerecon.root\",\"lbl_config/hk/smearing.root\",\"lbl_config/hk/response.root\",1234,\"$P\",\"enu_sk_numu\",6.75e21,1.132e35,14)"
root -l -b -q "src/FakeRecon.C(\"$data_dir/numu/numu_x_nue_50K_RHC_flattree.root\",\"$data_dir/numu/numu_x_nue_50K_RHC_fakerecon.root\",\"lbl_config/hk/smearing.root\",\"lbl_config/hk/response.root\",1234,\"$M\",\"enu_sk_numu\",2.025e22,1.132e35,14)"
root -l -b -q "src/FakeRecon.C(\"$data_dir/numubar/numubar_x_nuebar_50K_RHC_flattree.root\",\"$data_dir/numubar/numubar_x_nuebar_50K_RHC_fakerecon.root\",\"lbl_config/hk/smearing.root\",\"lbl_config/hk/response.root\",1234,\"$M\",\"enu_sk_numub\",2.025e22,1.132e35,-14)"
root -l -b -q "src/FakeRecon.C(\"$data_dir/numubar/numubar_x_nuebar_50K_FHC_flattree.root\",\"$data_dir/numubar/numubar_x_nuebar_50K_FHC_fakerecon.root\",\"lbl_config/hk/smearing.root\",\"lbl_config/hk/response.root\",1234,\"$P\",\"enu_sk_numub\",6.75e21,1.132e35,-14)"
