#!/usr/bin/env python3
"""Q2: joint hashed-Poisson bootstrap (2000 replicates, salt 0) of N_b,V_b for every pair table; saved in the scratchpad (res/boot_<name>.npz)."""
from q2_lib import *
for n in NAMES:
    d = load(n); N, V = boot_counts(d, 2000, 0); np.savez_compressed(RES / f"boot_{n}.npz", N=N, V=V); print(n, flush=True)
