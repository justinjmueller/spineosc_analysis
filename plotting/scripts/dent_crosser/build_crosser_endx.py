import pathlib as _pl, sys as _sys; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))
from spineosc_paths import REPO  # noqa: E402
import re

ROOT = str(REPO)
SRC_MC = f'{ROOT}/xml/detector_data_access/sbnd_run1/muon_endx_zoom_inclusive.xml'
SRC_DATA = f'{ROOT}/xml/detector_data_access/sbnd_run1/muon_endx_zoom_inclusive_data.xml'
OUT_MC = f'{ROOT}/xml/detector_data_access/sbnd_run1/muon_endx_zoom_crosser.xml'
OUT_DATA = f'{ROOT}/xml/detector_data_access/sbnd_run1/muon_endx_zoom_crosser_data.xml'

CROSSER = {
    'stay': '(reco_vertex_x * reco_leading_muon_end_x) >= 0',
    'xing': '(reco_vertex_x * reco_leading_muon_end_x) < 0',
}

HEADER_NOTE = '''<!--
DETECTOR-ACCESS comparison: SBND Run 1, muon_endx, split by cathode-crosser
status, area normalized. Built from muon_endx_zoom_inclusive.xml (same
cathode-zoom binning, ±10cm merged dead-zone bin) by duplicating the single
"numu" channel into "stay" (reco_vertex_x and reco_leading_muon_end_x on the
same side: (reco_vertex_x * reco_leading_muon_end_x) >= 0) and "xing"
(opposite sides: < 0). Purpose: test whether the data/CV deficit seen near
the cathode in the inclusive muon_endx plot concentrates in the crossing
population, as expected if SPINE is failing to stitch some cathode-crossing
tracks (track breaking) -- raised by Gray/Tom re: the DENT force meeting.

Detector-response systematics (ChargeScale15, WireModxThetaXW, MattModYZ,
DENT) are NOT re-split by stay/xing population here (still apply_to_subchannel
="_SBND_", same single HistVarSection as the inclusive file) -- this is a
fast first look at the central-value shape, not a rescoped chi^2 test. If
this shows something, the systematics should be re-split per population
(as in dent_crosser/muon_score.xml) before trusting a chi^2 number.

POT is BOGUS (1.0e+19) on both sides: the comparison is area normalized.
-->

'''

HEADER_NOTE_DATA = '''<!--
DETECTOR-ACCESS comparison: SBND Run 1, muon_endx, split by cathode-crosser
status, area normalized. Data companion to muon_endx_zoom_crosser.xml -- see
that file for the stay/xing definition and purpose.
-->

'''

def load(path):
    return open(path).read()

def write(path, text):
    with open(path, 'w') as f:
        f.write(text)
    print('wrote', path)

# ---------- MC file ----------
text = load(SRC_MC)

# Strip the original leading header comments (everything up to and including
# the first <mode name="nu" /> line's preceding comment block) and replace
# with our new note, keeping everything from <mode onward.
mode_idx = text.index('<mode name="nu" />')
# find start of the comment block immediately preceding <mode ...>
pre = text[:mode_idx]
# Keep the XML declaration line, drop old comments.
decl_end = text.index('?>') + 2
new_text = text[:decl_end] + '\n\n' + HEADER_NOTE + text[mode_idx:]
text = new_text

# Duplicate the <channel name="numu" ...> ... </channel> block into stay/xing.
chan_pat = re.compile(r'<channel name="numu"(.*?)</channel>\n', re.S)
m = chan_pat.search(text)
assert m, 'channel block not found'
chan_body = m.group(0)

def make_channel(pop):
    block = chan_body.replace('name="numu"', f'name="{pop}"', 1)
    return block

stay_block = make_channel('stay')
xing_block = make_channel('xing')
text = text[:m.start()] + stay_block + xing_block + text[m.end():]

# Duplicate every <branch ... associated_subchannel="nu_SBND_numu_X" ...>...</branch>
# into a stay and xing copy, appending the crosser condition to additional_weight
# and renaming numu -> stay/xing in associated_subchannel.
branch_pat = re.compile(r'(    <branch\n(?:[^\n]*\n)*?    </branch>\n)', re.M)

def split_branch(m):
    block = m.group(1)
    if 'associated_subchannel' not in block or '_numu_' not in block:
        return block
    out = []
    for pop, cond in CROSSER.items():
        b = block.replace('_numu_', f'_{pop}_', 1)
        # Append crosser condition inside the additional_weight string.
        aw_pat = re.compile(r'(additional_weight\s*=\s*")([^"]*)(")')
        def add_cond(am, cond=cond):
            return f'{am.group(1)}{am.group(2)} && {cond}{am.group(3)}'
        b = aw_pat.sub(add_cond, b, count=1)
        out.append(b)
    return ''.join(out)

text = branch_pat.sub(split_branch, text)

write(OUT_MC, text)

# ---------- Data file ----------
dtext = load(SRC_DATA)
mode_idx = dtext.index('<mode name="nu" />')
decl_end = dtext.index('?>') + 2
dtext = dtext[:decl_end] + '\n\n' + HEADER_NOTE_DATA + dtext[mode_idx:]

m = chan_pat.search(dtext)
assert m, 'data channel block not found'
dchan_body = m.group(0)
dtext = dtext[:m.start()] + dchan_body.replace('name="numu"', 'name="stay"', 1) \
      + dchan_body.replace('name="numu"', 'name="xing"', 1) + dtext[m.end():]

dtext = branch_pat.sub(split_branch, dtext)

write(OUT_DATA, dtext)
