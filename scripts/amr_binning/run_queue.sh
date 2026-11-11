#!/bin/bash
# The binning matrix, sequentially (clean timings; every variant processes, scans, then deletes its caches). Variants that are
# already done are skipped. The equal-population variants need run_population.sh first (see the README of docs/amr_binning).
#   scripts/amr_binning/run_queue.sh
cd "$(dirname "$0")/../.." || exit 1
for cfg in incl1D excl1D; do
  for var in merge2 merge4 uniform_0.10 uniform_0.05 equalpop_5 equalpop_9 equalpop_17 equalpop_8 equalpop_15; do
    [ -f runs/amr_binning/xml/${cfg}_${var}.xml ] && scripts/amr_binning/run_variant.sh $cfg $var 12
  done
done
