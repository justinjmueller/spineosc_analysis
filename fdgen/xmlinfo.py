"""What fdgen needs from the PROfit MC XML: per-file POT, per-branch model rule and weight, the model's L/E binning."""
import re

DETECTORS = {"SBND": "SBND", "ICARUSRun2": "ICARUSRun2", "ICARUSRun4": "ICARUSRun4"}


def attrs(tag):
    return dict(re.findall(r'(\w+)\s*=\s*"([^"]*)"', tag))


def load(xml_path):
    s = open(xml_path).read()
    det_pot = {m.group(1): float(m.group(2))
               for m in re.finditer(r'<detector name="([^"]+)"\s+pot="([^"]+)"\s*/>', s)}
    # the model variable is the variable with this index in every branch's <variable> list (variable_index)
    vi = int(re.search(r'<parameter name="L/E" variable_index="(\d+)"', s).group(1))
    le_bins = {(float(a), float(b), int(n)) for a, b, n in
               re.findall(r'<bins unit="True L/E"\s+min="([^"]+)" max="([^"]+)" nbins="([^"]+)"', s)}
    if len(le_bins) != 1:
        raise ValueError(f"expected one True L/E binning, found {le_bins}")
    te_bins = {(float(a), float(b), int(n)) for a, b, n in
               re.findall(r'<bins unit="True E \[GeV\]"\s+min="([^"]+)" max="([^"]+)" nbins="([^"]+)"', s)}
    if len(te_bins) != 1:
        raise ValueError(f"expected one True E binning, found {te_bins}")
    files = []
    for m in re.finditer(r'<MCFile .*?</MCFile>', s, re.S):
        block = m.group(0)
        a = attrs(re.search(r'<MCFile[^>]*>', block).group(0))
        branches = []
        for b in re.finditer(r'<branch\s(?:[^>"]|"[^"]*")*>.*?</branch>', block, re.S):
            bt = b.group(0)
            ba = attrs(re.match(r'<branch\s(?:[^>"]|"[^"]*")*>', bt).group(0))
            variables = re.findall(r'<variable(2D)?>(.*?)</variable(?:2D)?>', bt, re.S)
            # a <variable2D> occupies one index, like in PROfit
            le_expr = variables[vi][1].strip()
            te_expr = variables[2][1].strip()          # index 2 = the truth-energy variable the dial splines are binned in (binning="var2")
            det = re.match(r'nu_([A-Za-z0-9]+)_(?:ccqelike|nproton|inclusive)_', ba["associated_subchannel"]).group(1)
            branches.append(dict(subchannel=ba["associated_subchannel"], detector=det, rule=int(ba["model_rule"]),
                                 weight=ba["additional_weight"].strip(), le_expr=le_expr, te_expr=te_expr,
                                 variables=[v[1].strip() for v in variables]))
        files.append(dict(filename=a["filename"], treename=a["treename"], mc_pot=float(a["pot"]),
                          scale=float(a.get("scale", 1.0)), branches=branches))
    return dict(det_pot=det_pot, le_bins=next(iter(le_bins)), te_bins=next(iter(te_bins)), files=files, xml=xml_path)


def pot_scale(info, f):
    """PROfit's POT factor for a file: target POT of the (single) detector its branches belong to over the MC POT."""
    dets = {b["detector"] for b in f["branches"]}
    if len(dets) != 1:
        raise ValueError(f"{f['treename']}: branches belong to several detectors {dets}")
    return info["det_pot"][dets.pop()] / f["mc_pot"]


def load_detvars(xml_path):
    """The hist1d systematics of the XML (commented-out entries ignored): name -> dict(pattern, xvar, hist=(filename, histname) or
    sections=[(knobval, filename, histname), ...]). A name matching a <HistVarSection> carries that section's variations; otherwise the
    allowlist's own `filename` and `name` are the single (knob +1) ratio histogram."""
    s = re.sub(r"<!--.*?-->", "", open(xml_path).read(), flags=re.S)
    sections = {}
    for m in re.finditer(r'<HistVarSection\s+name="([^"]+)"\s*>(.*?)</HistVarSection>', s, re.S):
        sections[m.group(1)] = [(float(a["knobval"]), a["filename"], a["histname"])
                                for a in (attrs(v) for v in re.findall(r"<variation\s[^>]*/>", m.group(2)))]
    out = {}
    for m in re.finditer(r"<allowlist\s[^>]*>", s):
        a = attrs(m.group(0))
        if a.get("type") != "hist1d":
            continue
        entry = dict(pattern=a.get("apply_to_subchannel", ""), xvar=int(a["xvar"].replace("var", "")), plotname=a.get("plotname", ""))
        if a["name"] in sections:
            entry["sections"] = sections[a["name"]]
        else:
            entry["hist"] = (a["filename"], a["name"])
        out[a["name"]] = entry
    return out
