#!/bin/bash
# GENIE event generation with gevgen, one run per row of gen/genie/samples.txt.
# Needs gevgen on PATH and $data_dir/genie/xsec/hk_v1_stock.xml (scripts/genie/xsec.sh).
#
#   scripts/genie/gen.sh             all samples, in manifest order
#   scripts/genie/gen.sh <sample>    one sample, by name or 0-based row index
#
# Each sample runs in $data_dir/genie/<sample>/ and writes <sample>.ghep.root
# there, plus the input-flux.root that gevgen leaves in its working directory
# and scripts/genie/flat.sh reads.
set -o pipefail

[ -n "$working_dir" ] || { echo "working_dir is not set"; exit 1; }
cd "$working_dir"
[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }

XSEC="$data_dir/genie/xsec/hk_v1_stock.xml"
MSG="$working_dir/gen/genie/Messenger_production.xml"

TARGET="1000080160[0.8888888888888889],1000010010[0.1111111111111111]"

# Manifest rows without comments and blanks; $1 picks one by name or index.
ROWS=()
while read -r l; do ROWS+=("$l"); done < <(sed -e 's/#.*//' -e 's/[[:space:]]//g' -e '/^$/d' gen/genie/samples.txt)
SEL=()
for i in "${!ROWS[@]}"; do
  IFS='|' read -r name _ <<< "${ROWS[$i]}"
  [ -z "$1" ] || [ "$1" = "$name" ] || [ "$1" = "$i" ] && SEL+=("${ROWS[$i]}")
done
[ ${#SEL[@]} -gt 0 ] || { echo "ERROR: no sample $1 in gen/genie/samples.txt"; exit 1; }

for row in "${SEL[@]}"; do
  IFS='|' read -r NAME STEM PDG FLUX HIST NEVENTS RUN SEED _ _ <<< "$row"
  FLUX="$working_dir/$FLUX"

  RUNDIR="$data_dir/genie/$NAME"
  OUT="$RUNDIR/${NAME}.ghep.root"

  mkdir -p "$RUNDIR"
  cd "$RUNDIR" || exit 1

  echo "=================================================="
  echo "GENIE HK production"
  echo "Sample       : $NAME"
  echo "Host         : $(hostname)"
  echo "Start        : $(date)"
  echo "Probe PDG    : $PDG"
  echo "Events       : $NEVENTS"
  echo "Flux file    : $FLUX"
  echo "Flux hist    : $HIST"
  echo "Target       : $TARGET"
  echo "Run number   : $RUN"
  echo "Seed         : $SEED"
  echo "XSec         : $XSEC"
  echo "Output       : $OUT"
  echo "=================================================="

  if [ ! -s "$FLUX" ]; then
    echo "ERROR: flux file missing"
    exit 1
  fi

  if [ ! -s "$XSEC" ]; then
    echo "ERROR: xsec file missing"
    exit 1
  fi

  if [ ! -s "$MSG" ]; then
    echo "ERROR: messenger config missing"
    exit 1
  fi

  start_seconds=$(date +%s)

  gevgen \
    -r "$RUN" \
    -n "$NEVENTS" \
    -e 0,30 \
    -p "$PDG" \
    -t "$TARGET" \
    -f "$FLUX,$HIST,NOWIDTH" \
    -F 0,0,1,0,0,0,0 \
    -o "$OUT" \
    --seed "$SEED" \
    --cross-sections "$XSEC" \
    --tune G18_10b_02_11b \
    --event-generator-list Default \
    --message-thresholds "$MSG"

  status=$?

  end_seconds=$(date +%s)
  elapsed=$((end_seconds-start_seconds))

  echo "=================================================="
  echo "gevgen exit code: $status"
  echo "Elapsed seconds: $elapsed"
  echo "Elapsed minutes: $((elapsed/60))"
  echo "End: $(date)"

  if [ -f "$OUT" ]; then
    ls -lh "$OUT"
  fi

  echo "=================================================="

  [ "$status" -eq 0 ] || exit $status
done
