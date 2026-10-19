"""Per-event quantities PROfit derives from a branch's XML definition, evaluated with ROOT (same expression syntax).

For each MC file: aw (the branch's additional_weight: selection x cvwgt; 0 when the event is not in the branch),
rule (the branch's model_rule; -1 when in no branch), branch (its index, -1 when none) and true_e (the truth-energy variable) and le (the model variable, float32 as PROfit stores it),
and x (the histogram-lookup variable `var3` of the hist1d detector systematics, float32).
"""
import numpy as np


def evaluate(filename, treename, branches, selection=None, xvar=3):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    df = ROOT.RDataFrame(treename, filename)
    cols = []
    for i, b in enumerate(branches):
        df = df.Define(f"_aw{i}", b["weight"])
        cols.append(f"_aw{i}")
    exprs = sorted({b["le_expr"] for b in branches})
    for j, e in enumerate(exprs):
        df = df.Define(f"_le{j}", e)
        cols.append(f"_le{j}")
    texprs = sorted({b["te_expr"] for b in branches})
    for j, e in enumerate(texprs):
        df = df.Define(f"_te{j}", e)
        cols.append(f"_te{j}")
    xexprs = sorted({b["variables"][xvar] for b in branches})       # the variable PROfit's hist1d systematics look up (xvar="var<k>")
    for j, e in enumerate(xexprs):
        df = df.Define(f"_x{j}", e)
        cols.append(f"_x{j}")
    if selection:
        df = df.Define("_sel", f"(bool)({selection})")
        cols.append("_sel")
    a = df.AsNumpy(cols)
    n = len(next(iter(a.values())))
    aw = np.zeros(n)
    rule = np.full(n, -1, dtype=np.int8)
    le = np.full(n, np.nan, dtype=np.float32)
    nbr = np.zeros(n, dtype=np.int8)
    branch = np.full(n, -1, dtype=np.int16)
    true_e = np.full(n, np.nan, dtype=np.float32)
    x = np.full(n, np.nan, dtype=np.float32)
    for i, b in enumerate(branches):
        w = np.asarray(a[f"_aw{i}"], dtype=np.float64)
        hit = w != 0
        nbr += hit
        aw[hit] = w[hit]
        rule[hit] = b["rule"]
        branch[hit] = i
        true_e[hit] = np.asarray(a[f"_te{texprs.index(b['te_expr'])}"], dtype=np.float32)[hit]
        le[hit] = np.asarray(a[f"_le{exprs.index(b['le_expr'])}"], dtype=np.float32)[hit]
        x[hit] = np.asarray(a[f"_x{xexprs.index(b['variables'][xvar])}"], dtype=np.float32)[hit]
    if (nbr > 1).any():
        raise ValueError(f"{treename}: {(nbr > 1).sum()} events fall in more than one branch")
    out = dict(aw=aw, rule=rule, le=le, branch=branch, true_e=true_e, x=x, n=n)
    if selection:
        out["sel"] = np.asarray(a["_sel"], dtype=bool)
    return out


def evaluate_formula(filename, treename, formula):
    """Boolean mask of a TFormula-style cut over every entry of a tree (same expression syntax as the XML weights)."""
    import ROOT
    ROOT.gROOT.SetBatch(True)
    df = ROOT.RDataFrame(treename, filename).Define("_cut", f"(bool)({formula})")
    return np.asarray(df.AsNumpy(["_cut"])["_cut"], dtype=bool)
