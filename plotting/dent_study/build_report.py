#!/usr/bin/env python3
"""Builds plotting/dent_study.ipynb and docs/dent_study/report.html from report_cells.CELLS.
Run in the mutools environment:  /Users/mueller/Projects/mutools/.venv/bin/python build_report.py [--no-html]"""
import base64, contextlib, html, io, pathlib, re, sys, traceback
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import nbformat
import report_cells, dentstudy as ds

REPO = ds.REPO
NB = REPO / "plotting" / "dent_study.ipynb"
HTML = REPO / "docs" / "dent_study" / "report.html"

def write_notebook():
    nb = nbformat.v4.new_notebook()
    nb.metadata["kernelspec"] = {"name": "mutools-dev", "display_name": "Python (mutools-dev)", "language": "python"}
    nb.cells = [nbformat.v4.new_markdown_cell(t) if k == "md" else nbformat.v4.new_code_cell(t) for k, t in report_cells.CELLS]
    nbformat.write(nb, NB); print("wrote", NB)

# ---- minimal markdown -> html (headings, paragraphs, lists, pipe tables, **bold**, *italic*, `code`)
def inline(s):
    s = html.escape(s, quote=False)
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s); s = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", s); s = re.sub(r"(?<!\*)\*([^*\n]+)\*(?!\*)", r"<i>\1</i>", s)
    return s

def md_to_html(text):
    out, lines, i = [], text.split("\n"), 0
    while i < len(lines):
        ln = lines[i]
        if not ln.strip(): i += 1; continue
        m = re.match(r"^(#{1,4})\s+(.*)", ln)
        if m: out.append(f"<h{len(m.group(1))}>{inline(m.group(2))}</h{len(m.group(1))}>"); i += 1; continue
        if ln.lstrip().startswith("|"):
            rows = []
            while i < len(lines) and lines[i].lstrip().startswith("|"):
                cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-{2,}:?", c) for c in cells): rows.append(cells)
                i += 1
            head, body = rows[0], rows[1:]
            out.append("<div class='scroll'><table class='tbl'><thead><tr>" + "".join(f"<th>{inline(c)}</th>" for c in head) + "</tr></thead><tbody>" +
                       "".join("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in r) + "</tr>" for r in body) + "</tbody></table></div>"); continue
        if re.match(r"^\s*[-*]\s+", ln) or re.match(r"^\s*\d+\.\s+", ln):
            ordered = bool(re.match(r"^\s*\d+\.\s+", ln)); items = []
            while i < len(lines) and (re.match(r"^\s*[-*]\s+", lines[i]) or re.match(r"^\s*\d+\.\s+", lines[i]) or (lines[i].startswith("  ") and items)):
                if re.match(r"^\s*([-*]|\d+\.)\s+", lines[i]): items.append(re.sub(r"^\s*([-*]|\d+\.)\s+", "", lines[i]))
                else: items[-1] += " " + lines[i].strip()
                i += 1
            tag = "ol" if ordered else "ul"; out.append(f"<{tag}>" + "".join(f"<li>{inline(x)}</li>" for x in items) + f"</{tag}>"); continue
        para = [ln]; i += 1
        while i < len(lines) and lines[i].strip() and not re.match(r"^(#{1,4}\s|\s*[-*]\s|\s*\d+\.\s|\s*\|)", lines[i]): para.append(lines[i]); i += 1
        out.append("<p>" + inline(" ".join(x.strip() for x in para)) + "</p>")
    return "\n".join(out)

def img_tag(path, caption, maxw=1100):
    from PIL import Image
    im = Image.open(path); im.thumbnail((maxw, maxw * 4)); buf = io.BytesIO(); im.save(buf, format="PNG", optimize=True)
    return f"<figure><img alt='' src='data:image/png;base64,{base64.b64encode(buf.getvalue()).decode()}'><figcaption>{inline(caption)}</figcaption></figure>"

CSS = """
/* layout: single reading column (74ch) for prose, wide scrolling tables and full-width figures; orange accent taken from the SBND plot style */
:root{--bg:#f6f7f8;--ink:#14202b;--ink2:#51606d;--rule:#d9dee3;--card:#eceff2;--accent:#b8520a;--accent-ink:#8f3f06;--fig:#ffffff;
 --serif:"Source Serif 4",Georgia,"Times New Roman",serif;--sans:"IBM Plex Sans",system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;--mono:"IBM Plex Mono",ui-monospace,Menlo,Consolas,monospace}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){--bg:#12171c;--ink:#e6eaee;--ink2:#a2afbb;--rule:#2b353e;--card:#1b232a;--accent:#f0923f;--accent-ink:#f5aa68;--fig:#ffffff;color-scheme:dark}}
:root[data-theme=dark]{--bg:#12171c;--ink:#e6eaee;--ink2:#a2afbb;--rule:#2b353e;--card:#1b232a;--accent:#f0923f;--accent-ink:#f5aa68;--fig:#ffffff;color-scheme:dark}
body{background:var(--bg);color:var(--ink);font:16px/1.6 var(--sans);margin:0;padding:0 18px;-webkit-font-smoothing:antialiased}
main{max-width:1060px;margin:0 auto;padding-block:40px 96px}
h1,h2,h3{font-family:var(--serif);text-wrap:balance;font-weight:600;letter-spacing:-.01em}
h1{font-size:clamp(30px,5vw,44px);line-height:1.1;margin:0 0 14px;max-width:20ch}
main>h1+p{font-family:var(--serif);font-size:18px;color:var(--ink2);font-style:italic;margin:0 0 28px;max-width:60ch}
h2{font-size:26px;line-height:1.2;margin:56px 0 12px;padding-top:18px;border-top:2px solid var(--ink)}
h2:first-of-type{margin-top:34px}
h3{font-size:19px;margin:32px 0 8px}
p,li{max-width:74ch}p{margin:0 0 14px}li{margin-bottom:8px}ul,ol{padding-left:22px}
strong{font-weight:600}
code{font-family:var(--mono);background:var(--card);padding:1px 5px;border-radius:3px;font-size:.86em}
.scroll{overflow-x:auto;margin:8px 0 22px;border:1px solid var(--rule);border-radius:4px}
table.tbl{border-collapse:collapse;font:12.5px/1.4 var(--mono);width:max-content;min-width:100%;font-variant-numeric:tabular-nums}
.tbl th,.tbl td{padding:5px 10px;border-bottom:1px solid var(--rule);text-align:right;white-space:nowrap}
.tbl th{color:var(--ink2);font-weight:500;background:var(--card);position:sticky;top:0}
.tbl td:first-child,.tbl th:first-child{text-align:left}
.tbl tr:last-child td{border-bottom:0}
figure{margin:20px 0 28px}figure img{max-width:100%;height:auto;background:var(--fig);border:1px solid var(--rule);border-radius:3px;display:block}
figcaption{font-size:13.5px;color:var(--ink2);margin-top:8px;max-width:80ch}
pre.out{font:12.5px/1.5 var(--mono);background:var(--card);padding:10px 14px;border-radius:3px;overflow-x:auto;white-space:pre-wrap;border-left:3px solid var(--accent);margin:12px 0 20px}
.cap{font:600 12px/1.4 var(--sans);letter-spacing:.06em;text-transform:uppercase;color:var(--accent-ink);margin:26px 0 4px;max-width:90ch}
a{color:var(--accent-ink)}:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
@media (max-width:600px){body{font-size:15px}h2{font-size:22px}}
"""

def build_html():
    ns = {"__name__": "report"}; parts = []; failures = 0
    for kind, text in report_cells.CELLS:
        if kind == "md": parts.append(md_to_html(text)); continue
        ds._SINK = []; buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf): exec(compile(text, "<cell>", "exec"), ns)
        except Exception:
            failures += 1; buf.write("\nCELL FAILED:\n" + traceback.format_exc())
        printed = buf.getvalue().strip()
        if printed: parts.append(f"<pre class='out'>{html.escape(printed)}</pre>")
        for item in ds._SINK:
            if item[0] == "table": parts.append((f"<div class='cap'>{inline(item[1])}</div>" if item[1] else "") + f"<div class='scroll'>{item[2]}</div>")
            else: parts.append(img_tag(item[2], item[1]))
        ds._SINK = None
    HTML.parent.mkdir(parents=True, exist_ok=True)
    HTML.write_text(f"<title>SBND DENT Study</title><link rel='stylesheet' href='https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&family=Source+Serif+4:ital,wght@0,500;0,600;1,400&display=swap'><style>{CSS}</style><main>\n" + "\n".join(parts) + "\n</main>")
    print("wrote", HTML, f"({HTML.stat().st_size/1e6:.1f} MB, {failures} failed cells)")

if __name__ == "__main__":
    write_notebook()
    if "--no-html" not in sys.argv: build_html()
