#!/bin/bash
# FakeSelect, once every sample has a _fakerecon.root. Needs root on PATH.
# false = FHC, true = RHC.
[ -n "$working_dir" ] || { echo "working_dir is not set"; exit 1; }
cd "$working_dir"
[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }

root -l -b -q "src/FakeSelect.C(false,\"$data_dir\",\"$data_dir/fakeselect_FHC.root\")"
root -l -b -q "src/FakeSelect.C(true,\"$data_dir\",\"$data_dir/fakeselect_RHC.root\")"
