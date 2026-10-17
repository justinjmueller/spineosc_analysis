#!/usr/bin/env python3
"""Write the reco-energy closure study: the muon-score study XML with reco_visible_energy as var0.

The detector-access XMLs of the muon-score study already carry the production MC block, the full
systematic suite and the S2 selection. This swaps the observable for the analysis one (reconstructed
visible energy, the 17-bin axis of the detector-systematic ratio files) and the ratio file for the
production one (work/detsys_main.root). The pseudo-data template is derived the same way and is
only ever turned into a *simulation* sample by split_closure.py (its data file is replaced); no
real data XML is written, since reconstructed energy is not allowed for the detector data access
sample. The files go under runs/smoothing/studies (untracked).

    python scripts/detvar_smoothing/make_energy_study.py
"""
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SOURCE = REPO / "xml" / "dent_improvement" / "sbnd_run1"
OUT = REPO / "runs" / "smoothing" / "studies"
NAME = "reco_energy_inclusive_s2"
EDGES = "0.3 0.4 0.46 0.52 0.58 0.64 0.7 0.76 0.82 0.88 0.94 1.0 1.06 1.13 1.2 1.3 1.5 2.0"


def convert(text: str) -> str:
    """The muon-score XML with the energy observable, axis and production ratio file."""
    text = text.replace("<variable>reco_leading_muon_muon_softmax</variable>",
                        "<variable>reco_visible_energy</variable>")
    text = re.sub(r'<channel name="numu" unit="[^"]*">',
                  '<channel name="numu" unit="Reconstructed Visible Energy [GeV]">', text)
    text = re.sub(r'<bins unit="Reco Muon Score"[^>]*/>',
                  f'<bins unit="RecoE" edges="{EDGES}"/>', text, count=1)
    return text.replace("work/detsys_v2_muon_score.root", "work/detsys_main.root")


def main() -> None:
    """Write the three study files."""
    OUT.mkdir(parents=True, exist_ok=True)
    for suffix in (".xml", "_data.xml", "_plot.ini"):
        text = (SOURCE / f"muon_score_inclusive_s2{suffix}").read_text()
        text = text.replace("muon_score_inclusive_s2", NAME)
        if suffix != "_plot.ini":
            text = convert(text)
        else:
            text = re.sub(r"^(xml|data) = .*/", rf"\1 = {OUT}/", text, flags=re.M)
        (OUT / f"{NAME}{suffix}").write_text(text)
        print("wrote", OUT / f"{NAME}{suffix}")


if __name__ == "__main__":
    sys.exit(main())
