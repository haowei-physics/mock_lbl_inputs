#!/bin/bash
# GENIE cross-section response graphs, one file per sample. build/genie_splines
# reads the GHEP file (scripts/genie/gen.sh), the flat tree (scripts/genie/flat.sh)
# and the fakerecon tree (scripts/recon.sh) and writes <stem>_spline.root: the
# fakerecon tree plus one GENIE_* TGraph branch per dial, and the GenieDialMeta tree.
#
#   scripts/genie/splines.sh             all samples, in manifest order
#   scripts/genie/splines.sh <sample>    one sample, by name or 0-based row index
#   NEV=<n> scripts/genie/splines.sh ... first n events only (a test)
#
# Run scripts/genie/build.sh first.
[ -n "$working_dir" ] || { echo "working_dir is not set"; exit 1; }
cd "$working_dir"
[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }
NEV="${NEV:--1}"

BIN="$working_dir/build/genie_splines"
[ -x "$BIN" ] || { echo "ERROR: $BIN missing, run scripts/genie/build.sh"; exit 1; }
XSEC="$data_dir/genie/xsec/hk_v1_stock.xml"
MSG="$working_dir/gen/genie/Messenger_production.xml"

# Manifest rows without comments and blanks; $1 picks one by name or index.
ROWS=()
while read -r l; do ROWS+=("$l"); done < <(sed -e 's/#.*//' -e 's/[[:space:]]//g' -e '/^$/d' gen/genie/samples.txt)
SEL=()
for i in "${!ROWS[@]}"; do
  IFS='|' read -r name _ <<< "${ROWS[$i]}"
  [ -z "$1" ] || [ "$1" = "$name" ] || [ "$1" = "$i" ] && SEL+=("${ROWS[$i]}")
done
[ ${#SEL[@]} -gt 0 ] || { echo "ERROR: no sample $1 in gen/genie/samples.txt"; exit 1; }

# The calculators may drop debug files in the working directory.
mkdir -p "$data_dir/genie/run"
cd "$data_dir/genie/run"

for row in "${SEL[@]}"; do
  IFS='|' read -r S STEM _ <<< "$row"
  echo "--- genie_splines $S"
  "$BIN" "$data_dir/genie/$S/$S.ghep.root" "$data_dir/${STEM}_flattree.root" \
         "$data_dir/${STEM}_fakerecon.root" "$XSEC" "$MSG" "$data_dir/${STEM}_spline.root" "$NEV" \
    || { echo "ERROR: genie_splines failed for $S"; exit 1; }
done
