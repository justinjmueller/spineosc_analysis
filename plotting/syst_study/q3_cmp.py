import pathlib as _pl, sys as _sys; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from spineosc_paths import REPO  # noqa: E402
import sys, json, numpy as np
sys.path.insert(0, str(REPO / "plotting/syst_study")); from q3_native import *
name, var = sys.argv[1], sys.argv[2]; r, nom, inj = response(name, var); e = json.load(open(Q3 / "emulation.json"))[name]
d = r - np.array(e["pair"]); print(name, var, "max|native-pair|", np.abs(d).max(), "max|native-prod|", np.abs(r - np.array(e["prod"])).max()); print(np.round(r[:6], 5), np.round(e["pair"][:6], 5), nom[:3], e["Hn"][:3])
