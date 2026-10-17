import pathlib as _pl, sys as _sys; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from spineosc_paths import REPO  # noqa: E402
import sys, json, numpy as np
sys.path.insert(0, str(REPO / "plotting/syst_study")); from q3_native import *
name, var, key = sys.argv[1:4]; r, nom, inj = response(name, var); e = json.load(open(Q3 / "emulation.json"))[name]
print(name, var, "vs", key, "max|native-emul|", np.abs(r - np.array(e[key])).max(), " vs prod", np.abs(r - np.array(e["prod"])).max())
