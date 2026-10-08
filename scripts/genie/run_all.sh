#!/bin/bash
# The GENIE chain, in order, into $data_dir. Stops at the first stage that fails.
# Source your setup first: gmkspl, gspladd, gevgen, PrepareGENIE, nuisflat,
# root, genie-config on PATH and the variables scripts/genie/build.sh needs
# (docs/genie/README.md). FakeRecon and FakeSelect run through the shared
# scripts/recon.sh and scripts/select.sh, so use a data_dir of its own for GENIE.

set -e

[ -n "$working_dir" ] || { echo "working_dir is not set"; exit 1; }
cd "$working_dir"
[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }

# root -q exits with the macro's return value, and FakeRecon.C and FakeSelect.C
# return 1 on success, so the shared scripts exit with 1 when their last call
# succeeds, and their status says nothing about the earlier calls. It is not
# used here: their outputs are removed first and must all exist afterwards.
expect() {
  for f; do [ -s "$f" ] || { echo "ERROR: $f was not written"; exit 1; }; done
}
RECON=()
while IFS='|' read -r name stem _; do
  RECON+=("$data_dir/${stem}_fakerecon.root")
done < <(sed -e 's/#.*//' -e 's/[[:space:]]//g' -e '/^$/d' gen/genie/samples.txt)
SELECT=("$data_dir/fakeselect_FHC.root" "$data_dir/fakeselect_RHC.root")

scripts/init_data_dir.sh
scripts/genie/build.sh
scripts/genie/xsec.sh
scripts/genie/gen.sh
scripts/genie/flat.sh
for f in "${RECON[@]}"; do rm -f "${f:?}"; done
scripts/recon.sh || echo "scripts/recon.sh exited with $? (not used, outputs checked)"
expect "${RECON[@]}"
scripts/genie/splines.sh
for f in "${SELECT[@]}"; do rm -f "${f:?}"; done
scripts/select.sh || echo "scripts/select.sh exited with $? (not used, outputs checked)"
expect "${SELECT[@]}"
scripts/genie/make_covariances.sh
