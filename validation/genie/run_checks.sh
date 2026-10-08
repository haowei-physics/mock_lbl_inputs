#!/bin/bash
# Validation of the GENIE chain. Every output goes under $data_dir/validation/
# (logs in $data_dir/validation/logs/); nothing is written to the repository
# except ACLiC products in build/aclic/. Run scripts/genie/build.sh first.
#
# Dial scan (inputs of the "check" step, and of docs/genie/report.md):
#   scan [sample]           dial_scan, one-at-a-time scan at -2,-1,+1,+2 sigma
#                           -> dial_scan/part_<sample>.root (NEV, NSHAPE: event caps)
#   summarize               dial_scan summary -> dial_scan/genie_dial_sensitivity.root,
#                           ratio_plots.pdf, summary.log
#   aux [sample]            dial_aux -> dial_scan/aux_<sample>.txt (NEV: event cap)
#   modestats               mode_stats -> mode_stats/mode_stats.{root,pdf,md} (NEV)
# Spline files (after scripts/genie/splines.sh):
#   check [sample]          file-level checks of <stem>_spline.root vs the dial scan
#   interp [sample]         GENIE weights at -1.5,-0.5,+0.5,+1.5 sigma
#                           (interp/interp_<sample>.root) vs GUNDAM's interpolation
#   interpcheck [sample]    the comparison only, plus linear and monotonic variants
#                           for GENIE_CCQEPauliSupViaKF
# FakeSelect files (after scripts/select.sh):
#   checksel FHC|RHC        compare with a FakeSelect of the same events without dials
# GUNDAM (needs gundamFitter on PATH and GUNDAM_CONFIG = the GENIE GUNDAM config
# directory, with config.yaml, config_FHC.yaml and override_noOsc.yaml):
#   dryrun FHC|all          gundamFitter --dry-run
#   onesigma FHC|all        prefit one-sigma variations, oscillation set switched off
#   gundamcheck FHC|all     count_dials.C and bin_check.C against those outputs
#
# [sample]: a name or 0-based row index of gen/genie/samples.txt; all samples if omitted.
set -o pipefail

[ -n "$working_dir" ] || { echo "working_dir is not set"; exit 1; }
cd "$working_dir"
[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }
STEP="${1:-}"
ARG="${2:-}"
NEV="${NEV:--1}"
NSHAPE="${NSHAPE:-5000}"   # event cap per sample for the 5 slow shape-mode dials

BIN="$working_dir/build"
V="$data_dir/validation"
SCAN="$V/dial_scan"
LOGS="$V/logs"
RESP="$working_dir/lbl_config/hk/response.root"
mkdir -p "$SCAN" "$LOGS" "$V/interp" "$V/select" "$V/mode_stats"

# root -q exits with the macro's return value; FakeSelect.C returns 1 on success.
root_ok() { local rc=0; root -l -b -q "$@" || rc=$?; [ "$rc" -eq 1 ] || { echo "ERROR: root exited with $rc"; return 1; }; }

# Manifest rows without comments and blanks; samples() picks one by name or index.
ROWS=()
while read -r l; do ROWS+=("$l"); done < <(sed -e 's/#.*//' -e 's/[[:space:]]//g' -e '/^$/d' gen/genie/samples.txt)
samples() {
  SEL=()
  for i in "${!ROWS[@]}"; do
    IFS='|' read -r name _ <<< "${ROWS[$i]}"
    [ -z "$1" ] || [ "$1" = "$name" ] || [ "$1" = "$i" ] && SEL+=("${ROWS[$i]}")
  done
  [ ${#SEL[@]} -gt 0 ] || { echo "ERROR: no sample $1 in gen/genie/samples.txt"; exit 1; }
}
# Spline files of one beam mode, in manifest order, which within each mode is
# FakeSelect.C's file order (check_genie_splines select relies on it).
spline_files() {
  FILES=""
  for row in "${ROWS[@]}"; do
    IFS='|' read -r name stem _ <<< "$row"
    case "$1:$name" in all:*|FHC:*_fhc|RHC:*_rhc) FILES="$FILES $data_dir/${stem}_spline.root" ;; esac
  done
}
need() { [ -x "$BIN/$1" ] || { echo "ERROR: $BIN/$1 missing, run scripts/genie/build.sh"; exit 1; }; }
gundam_mode() {
  case "$1" in FHC) CFG=config_FHC ;; all) CFG=config ;; *) echo "usage: run_checks.sh $STEP FHC|all"; exit 1 ;; esac
  [ -n "$GUNDAM_CONFIG" ] || { echo "GUNDAM_CONFIG is not set"; exit 1; }
}
# The GENIE calculators may drop debug files in the working directory.
rundir() { mkdir -p "$data_dir/genie/run"; cd "$data_dir/genie/run"; }

case "$STEP" in
  scan)
    need dial_scan; samples "$ARG"; rundir
    for row in "${SEL[@]}"; do
      IFS='|' read -r S _ <<< "$row"
      "$BIN/dial_scan" scan "$working_dir" "$data_dir" "$S" "$SCAN" "$NEV" "$NSHAPE" 2>&1 \
        | tee "$LOGS/scan_$S.log" || exit 1
    done ;;
  summarize)
    need dial_scan
    "$BIN/dial_scan" summarize "$working_dir" "$data_dir" "$SCAN" | tee "$SCAN/summary.log" ;;
  aux)
    need dial_aux; samples "$ARG"; rundir
    for row in "${SEL[@]}"; do
      IFS='|' read -r S _ <<< "$row"
      SUFFIX=""; [ "$NEV" -gt 0 ] && SUFFIX="_test"
      "$BIN/dial_aux" "$working_dir" "$data_dir" "$S" "$SCAN/aux_${S}${SUFFIX}.txt" 5000 "$NEV" 2>&1 \
        | tee "$LOGS/aux_$S.log" || exit 1
    done ;;
  modestats)
    need mode_stats; rundir
    "$BIN/mode_stats" "$working_dir" "$data_dir" "$V/mode_stats" "$NEV" 2>&1 | tee "$LOGS/mode_stats.log" ;;
  check)
    need check_genie_splines; samples "$ARG"; fail=0
    for row in "${SEL[@]}"; do
      IFS='|' read -r S STEM _ <<< "$row"
      "$BIN/check_genie_splines" "$S" "$data_dir/${STEM}_spline.root" "$data_dir/${STEM}_fakerecon.root" \
          "$SCAN/part_$S.root" "$RESP" | tee "$LOGS/check_$S.log" || fail=1
    done
    exit $fail ;;
  interp)
    need genie_splines; need check_genie_splines; samples "$ARG"
    XSEC="$data_dir/genie/xsec/hk_v1_stock.xml"; MSG="$working_dir/gen/genie/Messenger_production.xml"
    for row in "${SEL[@]}"; do
      IFS='|' read -r S STEM _ <<< "$row"
      I="$V/interp/interp_$S.root"
      ( rundir
        "$BIN/genie_splines" --interp "$data_dir/genie/$S/$S.ghep.root" "$data_dir/${STEM}_flattree.root" \
            "$data_dir/${STEM}_fakerecon.root" "$XSEC" "$MSG" "$I" ) || exit 1
      "$BIN/check_genie_splines" interp "$S" "$I" "$RESP" | tee "$LOGS/interp_$S.log" || exit 1
    done ;;
  interpcheck)
    # CCQEPauliSupViaKF with plain catmull-rom (in use), linear Graph and
    # "catmull-rom, monotonic" (see the decision log in src/genie/genie_splines.cxx).
    need check_genie_splines; samples "$ARG"
    for row in "${SEL[@]}"; do
      IFS='|' read -r S _ <<< "$row"
      I="$V/interp/interp_$S.root"
      "$BIN/check_genie_splines" interp "$S" "$I" "$RESP" > "$LOGS/interp_$S.log"
      "$BIN/check_genie_splines" interp "$S" "$I" "$RESP" GENIE_CCQEPauliSupViaKF > "$LOGS/interp_${S}_pauliGraph.log"
      "$BIN/check_genie_splines" interp "$S" "$I" "$RESP" "" GENIE_CCQEPauliSupViaKF > "$LOGS/interp_${S}_pauliMonotonic.log"
      grep -h "^GENIE_CCQEPauliSupViaKF" "$LOGS/interp_$S.log" "$LOGS/interp_${S}_pauliGraph.log" "$LOGS/interp_${S}_pauliMonotonic.log"
    done ;;
  checksel)
    # The same events without dials (fakerecon inputs, same seed as scripts/select.sh),
    # only to check that the graphs do not change the sample draw.
    need check_genie_splines
    case "$ARG" in FHC) RHC=false ;; RHC) RHC=true ;; *) echo "usage: run_checks.sh checksel FHC|RHC"; exit 1 ;; esac
    NODIAL="$V/select/fakeselect_${ARG}_nodial.root"
    root_ok "src/FakeSelect.C($RHC,\"$data_dir\",\"$NODIAL\",\"lbl_config/hk/response.root\",20250731,\"_fakerecon.root\")" || exit 1
    spline_files "$ARG"
    "$BIN/check_genie_splines" select "$data_dir/fakeselect_$ARG.root" "$NODIAL" $FILES \
        | tee "$LOGS/checksel_$ARG.log" ;;
  dryrun)
    # Run from the config directory: its paths are relative. (--one-sigma is a
    # separate step: with the oscillation set on it stops at the fixed
    # PMNS_SIGN_MASS_SQUARED_32.)
    gundam_mode "$ARG"
    OUT="$V/gundam/dryrun_$ARG"; mkdir -p "$OUT"
    cd "$GUNDAM_CONFIG"
    gundamFitter -c $CFG.yaml --dry-run -t 8 --out-dir "$OUT/" 2>&1 | tee "$LOGS/dryrun_$ARG.log" ;;
  onesigma)
    # Oscillation set switched off (override_noOsc.yaml explains why); only for
    # checking the GENIE dials.
    gundam_mode "$ARG"
    OUT="$V/gundam/onesigma_$ARG"; mkdir -p "$OUT"
    cd "$GUNDAM_CONFIG"
    gundamFitter -c $CFG.yaml -of override_noOsc.yaml --dry-run --one-sigma -t 8 --out-dir "$OUT/" 2>&1 \
        | tee "$LOGS/onesigma_$ARG.log" ;;
  gundamcheck)
    gundam_mode "$ARG"
    SEL_FILES="$data_dir/fakeselect_FHC.root"; [ "$ARG" = all ] && SEL_FILES="$SEL_FILES $data_dir/fakeselect_RHC.root"
    B="$GUNDAM_CONFIG/inputs/samples/binning"
    ONE="$V/gundam/onesigma_$ARG/gundamFitter_${CFG}_With_override_noOsc_OneSigma_DryRun.root"
    ACLIC='gSystem->SetBuildDir("build/aclic", true);'
    root -l -b -q -e "$ACLIC" "validation/genie/count_dials.C+(\"$SEL_FILES\",\"$B\")" \
        | tee "$LOGS/count_dials_$ARG.log"
    root -l -b -q -e "$ACLIC" "validation/genie/bin_check.C+(\"$SEL_FILES\",\"$B\",\"$data_dir/parameters/xsec/xsecCovGENIE.root\",\"$ONE\")" \
        | tee "$LOGS/bin_check_$ARG.log" ;;
  *)
    echo "usage: validation/genie/run_checks.sh scan|aux [sample] | summarize | modestats |"
    echo "       check|interp|interpcheck [sample] | checksel FHC|RHC | dryrun|onesigma|gundamcheck FHC|all"
    exit 1 ;;
esac
