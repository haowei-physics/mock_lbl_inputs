#!/bin/bash
# Creates the empty data tree that every other script writes into.

[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }

mkdir -p $data_dir/numu
mkdir -p $data_dir/numubar
mkdir -p $data_dir/nue
mkdir -p $data_dir/nuebar
mkdir -p $data_dir/parameters/oscProb
mkdir -p $data_dir/parameters/det
mkdir -p $data_dir/parameters/flux
mkdir -p $data_dir/parameters/xsec
