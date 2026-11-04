#!/bin/bash
# GENIE cross-section splines for the HK targets (O16 + H1), one gmkspl run per
# neutrino flavour, merged with gspladd into $data_dir/genie/xsec/hk_v1_stock.xml.
#  G18_10b_02_11b tune, no --xml-path. Needs gmkspl and gspladd on PATH.
#
#   scripts/genie/xsec.sh              all four flavours, then the merge
#   scripts/genie/xsec.sh <flavour>    one flavour (numu, numubar, nue, nuebar), no merge
#   scripts/genie/xsec.sh merge        the merge only
#
# One flavour takes about a day; README.md shows how to run the four
# as separate jobs.
[ -n "$working_dir" ] || { echo "working_dir is not set"; exit 1; }
cd "$working_dir"
[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }

TAG="stock"
OUTDIR="$data_dir/genie/xsec"
mkdir -p "$OUTDIR/parts"

FLAVOURS=(numu numubar nue nuebar)
declare -A PDG=([numu]=14 [numubar]=-14 [nue]=12 [nuebar]=-12)
declare -A SEED=([numu]=314159 [numubar]=314160 [nue]=314161 [nuebar]=314162)

one() {
  NAME=$1
  OUT="$OUTDIR/parts/hk_v1_${TAG}_${NAME}.xml"
  INFO="$OUTDIR/parts/RUN_INFO_${NAME}.txt"

  {
    echo "tag           : $TAG"
    echo "description   : unmodified GENIE tune, no --xml-path override"
    echo "genie         : 3.06.02"
    echo "tune          : G18_10b_02_11b"
    echo "xml-path      : (none)"
    echo "flavour       : $NAME (pdg ${PDG[$NAME]})"
    echo "targets       : 1000080160,1000010010"
    echo "knots / emax  : 300 / 30 GeV"
    echo "seed          : ${SEED[$NAME]}"
    echo "output        : $OUT"
    echo "host          : $(hostname)"
    echo "slurm         : ${SLURM_JOB_ID:-none}"
    echo "started       : $(date -Is)"
  } > "$INFO"
  cat "$INFO"

  start=$(date +%s)

  gmkspl \
    -p "${PDG[$NAME]}" \
    -t 1000080160,1000010010 \
    -n 300 \
    -e 30 \
    --tune G18_10b_02_11b \
    --event-generator-list Default \
    --seed "${SEED[$NAME]}" \
    -o "$OUT"

  status=$?
  end=$(date +%s)

  {
    echo "finished      : $(date -Is)"
    echo "elapsed       : $(( (end-start)/60 )) min"
    echo "exit status   : $status"
    [ -f "$OUT" ] && echo "sha256        : $(sha256sum "$OUT" | cut -d' ' -f1)"
  } >> "$INFO"
  cat "$INFO"

  return $status
}

merge() {
  INPUTS=""
  echo "Checking input spline files..."
  for f in "${FLAVOURS[@]}"; do
    p="$OUTDIR/parts/hk_v1_${TAG}_${f}.xml"
    if [ ! -s "$p" ]; then
      echo "ERROR: missing or empty: $p"
      return 1
    fi
    ls -lh "$p"
    INPUTS="$INPUTS${INPUTS:+,}$p"
  done

  echo "Merging..."
  rm -f "${OUTDIR:?}/hk_v1_${TAG:?}.xml"
  gspladd -f "$INPUTS" -o "$OUTDIR/hk_v1_${TAG}.xml" || return 1
  [ -s "$OUTDIR/hk_v1_${TAG}.xml" ] || { echo "ERROR: gspladd wrote nothing"; return 1; }

  echo "Done:"
  ls -lh "$OUTDIR/hk_v1_${TAG}.xml"
}

case "${1:-all}" in
  all)
    for f in "${FLAVOURS[@]}"; do one "$f" || exit 1; done
    merge || exit 1 ;;
  merge)
    merge || exit 1 ;;
  numu|numubar|nue|nuebar)
    one "$1" || exit 1 ;;
  *)
    echo "usage: scripts/genie/xsec.sh [all|merge|numu|numubar|nue|nuebar]"; exit 1 ;;
esac
