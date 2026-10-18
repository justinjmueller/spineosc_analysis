#!/bin/bash
# Process (cache) and run the AMR sensitivity surface for one binning variant of one configuration, then delete the caches.
#   scripts/amr_binning/run_variant.sh <config incl1D|excl1D> <variant> [threads]
# Same PROfit options as runs/amr_compare/run_stage01.sh (Asimov data, AMR, delta-chi2 = 1.6 refinement). Outputs go to
# runs/amr_binning/<config>_<variant>/ ; the *_prop.bin/_syst.bin caches (1-3 GB) are removed on success to save disk.
set -u
P=${PROFIT:-/Users/mueller/Projects/profit-asymmetric/build/bin/PROfit}
cfg=$1; var=$2; nt=${3:-12}
repo="$(cd "$(dirname "$0")/../.." && pwd)"
d="$repo/runs/amr_binning/${cfg}_${var}"
mkdir -p "$d/logs"; cd "$d" || exit 1
[ -e work ] || ln -s "$repo/runs/amr_compare/work" work
xml="$repo/runs/amr_binning/xml/${cfg}_${var}.xml"
tag="ab_${cfg}_${var}"
printf 'xml = %s\ntag = %s\nnthread = %s\nverbosity = 3\nfile-verbosity = 2\n' "$xml" "$tag" "$nt" > run.ini
if [ ! -f done_process ]; then
  /usr/bin/time -p $P --config run.ini -n "$nt" -o proc --use-fake-data -l logs/process.log process > logs/process.out 2> logs/process.time && echo ok > done_process
fi
[ -f done_process ] || { echo "process failed for $cfg $var"; exit 1; }
if [ ! -f done_surface ]; then
  /usr/bin/time -l $P --config run.ini -n "$nt" -o surf --use-fake-data -s 405 -l logs/surface.log \
    surface --surface-amr --amr-initial 12 --amr-levels 3 --amr-levels-chi2 1.6 \
            --xvar sinsq2thmm --yvar dmsq --xlo 1e-3 --xhi 1 --ylo 1e-2 --yhi 1e2 \
    > logs/surface.out 2> logs/surface.time && [ -s "${tag}_surf_surface_amr.txt" ] && echo ok > done_surface
fi
# (the fine_* variants have empty bins and abort in `process`; their population is made by run_population.sh)
if [ -f done_surface ]; then rm -f "${tag}"_prop.bin "${tag}"_syst.bin "${tag}"_data.bin; fi
