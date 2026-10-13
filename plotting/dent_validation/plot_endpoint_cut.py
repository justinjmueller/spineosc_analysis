import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
d=pd.read_csv("flip.csv",dtype={'cls':str})
BLUE,ORANGE,AQUA="#2a78d6","#eb6834","#1baf7a"; INK,INK2,GRID,SURF="#0b0b0b","#52514e","#e4e3df","#fcfcfb"
plt.rcParams.update({"font.size":10,"axes.edgecolor":INK2,"axes.labelcolor":INK,"xtick.color":INK2,"ytick.color":INK2,"axes.facecolor":SURF,"figure.facecolor":SURF,
  "axes.spines.top":False,"axes.spines.right":False,"axes.grid":True,"grid.color":GRID,"grid.linewidth":0.6,"axes.axisbelow":True})
fl=d[d.cls=="10"]; nc=d[d.cls=="00"]
fig,ax=plt.subplots(1,2,figsize=(11.5,4.6),constrained_layout=True)
a=ax[0]; bins=np.arange(0,62,2.0)
for x,col,c,l in ((fl,"excv",BLUE,"flipped events, CV endpoint"),(fl,"exdt",ORANGE,"flipped events, DENT endpoint"),(nc,"exdt",AQUA,"non-crossers (both samples), DENT endpoint")):
    v=x[col].abs(); a.hist(v,bins=bins,weights=np.full(len(v),1/len(v)/2.0),histtype="step",lw=1.7,color=c,label=l)
a.set_yscale("log"); a.set_xlabel("muon endpoint distance from the cathode  |x$_{end}$|  [cm]"); a.set_ylabel("fraction of events / cm")
a.set_title("(a) endpoints pile up at the cathode in the DENT sample",loc="left",color=INK,fontsize=10.5)
a.legend(loc="upper right",frameon=False,fontsize=8.8,labelcolor=INK)
a.axvline(5,color=INK2,lw=0.9,ls="--"); a.text(5.8,0.03,"5 cm",color=INK2,fontsize=8.5)
b=ax[1]; E=np.array([0.4,0.6,0.8,1.0,1.3,1.75]); xerr=np.array([.1,.1,.1,.1,.2,.25])
rows=[("no crosser cut (v0)",[0.9807,0.9832,0.9774,0.9780,0.9696,0.9700],AQUA,"o"),
      ("crosser == 0 (v1)",[0.9993,1.0132,1.0190,1.0279,1.0280,1.0330],BLUE,"s"),
      ("crosser == 0 and |x$_{end}$| ≥ 5 cm",[0.9839,0.9914,0.9928,1.0001,0.9993,0.9932],ORANGE,"D")]
for l,y,c,mk in rows: b.plot(E,y,marker=mk,ms=5.5,lw=1.6,color=c,mec=SURF,mew=1,label=l)
b.axhline(1,color=INK2,lw=0.8); b.set_ylim(0.955,1.045); b.set_xlabel("reconstructed visible energy  [GeV]"); b.set_ylabel("DENT / CV (selected events)")
b.set_title("(b) the endpoint cut restores a flat DENT effect",loc="left",color=INK,fontsize=10.5)
b.legend(loc="upper left",frameon=False,fontsize=8.8,labelcolor=INK)
fig.text(0.005,-0.02,"SBND, events passing cathode + high-y/z cuts, matched between CV and DENT samples. (a) normalised to each population; 2 cm bins. (b) statistical scatter is about 1% per point.",color=INK2,fontsize=8.3,ha="left")
fig.savefig("/Users/mueller/Projects/spineosc/dent_endpoint_cut.png",dpi=160,facecolor=SURF,bbox_inches="tight"); print("saved")
