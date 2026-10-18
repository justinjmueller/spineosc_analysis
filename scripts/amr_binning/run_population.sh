#!/bin/bash
# Fine-binned population of a configuration, for the equal-population binnings (make_binnings.py --population).
#   scripts/amr_binning/run_population.sh <incl1D|excl1D>
# The fine_0.025 binning has empty bins (the reconstructed energy starts at ~0.33 GeV and the tail is empty), so PROfit's empty-bin
# check aborts it: process and plot are run with --force, only to read the populations. It is not a valid fit configuration.
set -u
P=${PROFIT:-/Users/mueller/Projects/profit-asymmetric/build/bin/PROfit}
cfg=$1
repo="$(cd "$(dirname "$0")/../.." && pwd)"
python3 "$repo/scripts/amr_binning/make_binnings.py" --base "$repo/runs/amr_compare/xml/contours_$cfg.xml" --name "$cfg" --out "$repo/runs/amr_binning/xml" --variants fine_0.025
d="$repo/runs/amr_binning/${cfg}_fine_0.025"; mkdir -p "$d/logs"; cd "$d" || exit 1
[ -e work ] || ln -s "$repo/runs/amr_compare/work" work
printf 'xml = %s\ntag = ab_%s_fine_0.025\nnthread = 12\nverbosity = 3\nfile-verbosity = 2\n' "$repo/runs/amr_binning/xml/${cfg}_fine_0.025.xml" "$cfg" > run.ini
$P --config run.ini -n 12 -o proc --use-fake-data --force -l logs/process.log process > logs/process.out 2> logs/process.time
$P --config run.ini -n 12 -o fine --use-fake-data --force -l logs/plot.log plot > logs/plot.out 2> logs/plot.time   # single-threaded covariance step, ~25 min
echo "equal-population edges: python3 scripts/amr_binning/make_binnings.py --base runs/amr_compare/xml/contours_$cfg.xml --name $cfg --out runs/amr_binning/xml --population $d/ab_${cfg}_fine_0.025_fine_PROplot_histograms.root --variants equalpop_9 equalpop_17"
