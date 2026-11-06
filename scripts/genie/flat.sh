#!/usr/bin/env bash
# Convert the GENIE GHEP samples to NUISANCE GenericVectors flat trees:
#   GHEP -> PrepareGENIE -> nuisflat -> $data_dir/<stem>_flattree.root
# so that the shared scripts/recon.sh reads them under their NEUT names.
# Needs PrepareGENIE and nuisflat (NUISANCE built with GENIE) on PATH.
#
#   scripts/genie/flat.sh             all samples, in manifest order
#   scripts/genie/flat.sh <sample>    one sample, by name or 0-based row index

set -eo pipefail

[ -n "$working_dir" ] || { echo "working_dir is not set"; exit 1; }
cd "$working_dir"
[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }

TARGET="1000080160[0.8888888888888889],1000010010[0.1111111111111111]"
FLUX_HIST="spectrum"

# Manifest rows without comments and blanks; $1 picks one by name or index.
ROWS=()
while read -r l; do ROWS+=("$l"); done < <(sed -e 's/#.*//' -e 's/[[:space:]]//g' -e '/^$/d' gen/genie/samples.txt)
SEL=()
for i in "${!ROWS[@]}"; do
  IFS='|' read -r name _ <<< "${ROWS[$i]}"
  [ -z "${1:-}" ] || [ "$1" = "$name" ] || [ "$1" = "$i" ] && SEL+=("${ROWS[$i]}")
done
[ ${#SEL[@]} -gt 0 ] || { echo "ERROR: no sample $1 in gen/genie/samples.txt"; exit 1; }


# ------------------------------------------------------------
# Environment checks
# ------------------------------------------------------------

echo
echo "============================================================"
echo " GENIE -> NUISANCE FlatTree Production"
echo "============================================================"
echo "ROOT       : $(root-config --version)"
echo "GENIE      : ${GENIE:-not set}"
echo "NUISANCE   : ${NUISANCE:-not set}"
echo "nuisflat   : $(command -v nuisflat)"
echo "PrepareGENIE: $(command -v PrepareGENIE)"
echo "Target     : $TARGET"
echo "============================================================"
echo


# ------------------------------------------------------------
# Conversion loop
# ------------------------------------------------------------

total=${#SEL[@]}
index=0

for row in "${SEL[@]}"; do

    IFS='|' read -r sample stem _ <<< "$row"
    index=$((index + 1))

    sample_dir="$data_dir/genie/$sample"

    ghep="$sample_dir/$sample.ghep.root"
    flux="$sample_dir/input-flux.root"

    prepared="$sample_dir/$sample.prepared.ghep.root"
    flat="$data_dir/${stem}_flattree.root"

    echo
    echo "============================================================"
    echo " [$index/$total] $sample"
    echo "============================================================"

    # --------------------------------------------------------
    # Check required input files
    # --------------------------------------------------------

    if [[ ! -f "$ghep" ]]; then
        echo "[ERROR] Missing GHEP file:"
        echo "        $ghep"
        exit 1
    fi

    if [[ ! -f "$flux" ]]; then
        echo "[ERROR] Missing flux file:"
        echo "        $flux"
        exit 1
    fi


    # --------------------------------------------------------
    # Step 1: PrepareGENIE
    # --------------------------------------------------------

    echo
    echo "[1/2] PrepareGENIE"

    rm -f "${prepared:?}"

    PrepareGENIE \
        -i "$ghep" \
        -f "$flux,$FLUX_HIST" \
        -t "$TARGET" \
        -o "$prepared"

    if [[ ! -s "$prepared" ]]; then
        echo "[ERROR] PrepareGENIE failed to create:"
        echo "        $prepared"
        exit 1
    fi


    # --------------------------------------------------------
    # Step 2: nuisflat
    # --------------------------------------------------------

    echo
    echo "[2/2] nuisflat"

    mkdir -p "$(dirname "$flat")"
    rm -f "${flat:?}"

    nuisflat \
        -i "GENIE:$prepared" \
        -f GenericVectors \
        -o "$flat"

    if [[ ! -s "$flat" ]]; then
        echo "[ERROR] nuisflat failed to create:"
        echo "        $flat"
        exit 1
    fi


    echo
    echo "[OK] $sample"
    echo "     FlatTree: $flat"
    ls -lh "$flat"

done


echo
echo "============================================================"
echo " $total SAMPLE(S) FINISHED"
echo "============================================================"
echo
