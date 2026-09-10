#!/bin/bash
# The whole chain, in order, into $data_dir. Stops at the first stage that fails.
# Source your setup first: neutroot2, nuisflat, T2KReWeightCI, root on PATH.

set -e

[ -n "$working_dir" ] || { echo "working_dir is not set"; exit 1; }
cd "$working_dir"
[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }

scripts/init_data_dir.sh
scripts/gen.sh
scripts/weight.sh
scripts/flat.sh
scripts/recon.sh
scripts/splines.sh
scripts/select.sh
scripts/make_covariances.sh
