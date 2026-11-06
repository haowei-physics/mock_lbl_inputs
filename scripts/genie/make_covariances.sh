#!/bin/bash
# GENIE counterpart of scripts/neut/make_covariances.sh. Writes into $data_dir/parameters/:
#   oscProb, det, flux   the generator-independent files of src/MakeCovariances.C
#                        (MakeOscCov, MakeDetCov, MakeFluxCov with the arguments
#                        MakeCovariances() uses; MakeXsecCov is not called, since it
#                        reads a T2KReWeight weight file)
#   oscProb/grid_HK.root src/MakeOscGrid.C
#   xsec/xsecCovGENIE.root src/genie/MakeGenieXsecCov.C, from the GenieDialMeta of
#                        the spline files of the given beam mode(s)
#
#   scripts/genie/make_covariances.sh [all|FHC|RHC]    (default all)
[ -n "$working_dir" ] || { echo "working_dir is not set"; exit 1; }
cd "$working_dir"
[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }
MODE="${1:-all}"
case "$MODE" in all|FHC|RHC) ;; *) echo "usage: scripts/genie/make_covariances.sh [all|FHC|RHC]"; exit 1 ;; esac

# root -q exits with the macro's (or the last -e expression's) return value;
# these macros return 1 on success.
root_ok() { local rc=0; root -l -b -q "$@" || rc=$?; [ "$rc" -eq 1 ] || { echo "ERROR: root exited with $rc"; return 1; }; }

P="$data_dir/parameters"
mkdir -p "$P"/{oscProb,det,flux,xsec}
OUTS=(oscProb/oscCovHK.root det/detCovHK.root det/detBinning.txt flux/fluxCovHK.root
      flux/fluxCovHK_ND.root flux/fluxBinning.txt oscProb/grid_HK.root xsec/xsecCovGENIE.root)
# No output of an earlier run can pass the checks below.
for f in "${OUTS[@]}"; do rm -f "${P:?}/${f:?}"; done

root_ok -e '.L src/MakeCovariances.C' \
  -e "MakeOscCov(\"$P/oscProb/oscCovHK.root\") & MakeDetCov(0.0,\"$P/det/detCovHK.root\",\"$P/det/detBinning.txt\") & MakeFluxCov(0.8,\"lbl_config/hk/t2kflux_2020_covariance.root\",\"$P/flux/fluxCovHK.root\",\"$P/flux/fluxCovHK_ND.root\",\"$P/flux/fluxBinning.txt\")" \
  || exit 1
root_ok "src/MakeOscGrid.C(\"$P/oscProb/grid_HK.root\")" || exit 1

FILES=""
while IFS='|' read -r name stem _; do
  case "$MODE:$name" in all:*|FHC:*_fhc|RHC:*_rhc) FILES="$FILES $data_dir/${stem}_spline.root" ;; esac
done < <(sed -e 's/#.*//' -e 's/[[:space:]]//g' -e '/^$/d' gen/genie/samples.txt)
root_ok "src/genie/MakeGenieXsecCov.C(\"$FILES\",\"$P/xsec/xsecCovGENIE.root\")" || exit 1

for f in "${OUTS[@]}"; do
  [ -s "$P/$f" ] || { echo "ERROR: missing $P/$f"; exit 1; }
done
echo "parameter files written to $P"
