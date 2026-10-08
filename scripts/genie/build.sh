#!/bin/bash
# Compiles the GENIE programs into build/ (ignored by git):
#   build/genie_splines, build/check_genie_splines      from src/genie/
#   build/dial_scan, build/dial_aux, build/mode_stats     from validation/genie/
# Needs root-config, genie-config and xml2-config on PATH, and
#   GENIE, GENIE_REWEIGHT   source trees (headers under src/; Reweight libraries under lib/)
#   LOG4CPP, LHAPDF5        install prefixes (include/ and lib/ below them)
set -e

[ -n "$working_dir" ] || { echo "working_dir is not set"; exit 1; }
cd "$working_dir"
for v in GENIE GENIE_REWEIGHT LOG4CPP LHAPDF5; do
  [ -n "${!v}" ] || { echo "$v is not set"; exit 1; }
done

mkdir -p build

# The link line of grwght1p in GENIE Reweight R-1_04_02. -Isrc is for
# NuFast_LBL.cpp, which the dial scan programs include.
CXX=(g++ -O2 -std=c++17 -pthread -w
     $(root-config --cflags) $(xml2-config --cflags)
     -I"$LOG4CPP/include" -I"$GENIE/src" -I"$GENIE_REWEIGHT/src" -Isrc)
LIBS=(-Wl,--no-as-needed
      $(root-config --libs) -lGeom -lEG -lEGPythia6 -lGenVector -lMathMore -lMinuit
      -L"$LHAPDF5/lib" -lLHAPDF -lgfortran -lxml2
      -L"$LOG4CPP/lib" -llog4cpp -lgsl -lgslcblas
      $(genie-config --libs) -L"$GENIE_REWEIGHT/lib" -lGRwFwk -lGRwIO -lGRwClc)

"${CXX[@]}" src/genie/genie_splines.cxx -o build/genie_splines "${LIBS[@]}"
for p in dial_scan dial_aux mode_stats; do
  "${CXX[@]}" validation/genie/$p.cxx -o build/$p "${LIBS[@]}"
done
# The checker needs ROOT only.
g++ -O2 -std=c++17 -w $(root-config --cflags) src/genie/check_genie_splines.cxx \
    -o build/check_genie_splines $(root-config --libs)

echo "built build/genie_splines build/check_genie_splines build/dial_scan build/dial_aux build/mode_stats"
