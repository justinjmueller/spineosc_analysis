import pathlib as _pl, sys as _sys; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from spineosc_paths import REPO  # noqa: E402
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
d=pd.read_csv("flip.csv",dtype={'cls':str})
# validated categorical slots 1-3 (light); text stays in ink
BLUE,ORANGE,AQUA="#2a78d6","#eb6834","#1baf7a"
INK,INK2,GRID,SURF="#0b0b0b","#52514e","#e4e3df","#fcfcfb"
plt.rcParams.update({"font.size":10,"axes.edgecolor":INK2,"axes.labelcolor":INK,"xtick.color":INK2,"ytick.color":INK2,
    "axes.facecolor":SURF,"figure.facecolor":SURF,"axes.spines.top":False,"axes.spines.right":False,"axes.grid":True,"grid.color":GRID,"grid.linewidth":0.6,"axes.axisbelow":True})
fl=d[d.cls=="10"].copy(); st=d[d.cls=="11"].copy(); nc=d[d.cls=="00"].copy()
for x in (fl,st,nc):
    x["dL"]=x.Lcv-x.Ldt
    # straight-line extent of the CV track beyond the cathode (x=0): |end_x| for a track whose end is on the far side
    x["post_x"]=x.excv.abs()
    x["Lpost"]=x.Lcv*x.excv.abs()/(x.excv-x.sxcv).abs()      # length past the cathode if the track is straight
fig,ax=plt.subplots(2,2,figsize=(11,8.2),constrained_layout=True)

# (a) missing length
a=ax[0,0]; bins=np.linspace(-20,300,81)
over={}
for x,c,l in ((nc,AQUA,"non-crosser in both"),(st,BLUE,"crosser in both"),(fl,ORANGE,"flipped: crosser → non-crosser")):
    w=np.full(len(x),1.0/len(x)/np.diff(bins)[0])        # density normalised over ALL events
    a.hist(x.dL,bins=bins,weights=w,histtype="step",color=c,lw=1.6,label=l)
    over[l]=(x.dL>bins[-1]).sum()
a.set_yscale("log"); a.set_ylim(1e-6,0.4); a.set_xlabel("missing track length  L$_{CV}$ − L$_{DENT}$  [cm]"); a.set_ylabel("fraction of events / cm (density)")
a.set_title("(a) leading-muon length lost in the DENT sample",loc="left",color=INK,fontsize=10.5)
a.legend(loc="upper right",frameon=False,fontsize=9,labelcolor=INK,handlelength=1.6)

# (b) missing length vs length beyond cathode
b=ax[0,1]
b.scatter(fl.Lpost,fl.dL,s=7,color=ORANGE,alpha=0.35,linewidths=0)
m=max(fl.Lpost.quantile(.99),fl.dL.quantile(.99)); b.plot([0,m],[0,m],color=INK2,lw=1,ls="--"); b.set_xlim(0,m); b.set_ylim(-20,m)
b.set_xlabel("CV track length beyond the cathode (straight-line estimate)  [cm]"); b.set_ylabel("missing length  L$_{CV}$ − L$_{DENT}$  [cm]")
b.set_title("(b) flipped events: is it the post-cathode segment?",loc="left",color=INK,fontsize=10.5)
b.text(0.97,0.06,"dashed: missing length = length beyond cathode",transform=b.transAxes,ha="right",color=INK2,fontsize=8.5)
b.annotate("no length lost\n(muon unchanged)",xy=(m*0.55,0),xytext=(m*0.55,m*0.17),color=INK,fontsize=8.5,ha="center",arrowprops=dict(arrowstyle="-",color=INK2,lw=0.8))
b.annotate("whole post-cathode\nsegment lost",xy=(m*0.42,m*0.42),xytext=(m*0.08,m*0.62),color=INK,fontsize=8.5,ha="left",arrowprops=dict(arrowstyle="-",color=INK2,lw=0.8))
r=np.corrcoef(fl.Lpost,fl.dL)[0,1]

# (c) depth past the cathode in the CV track, flipped vs not
c=ax[1,0]; bc=np.linspace(0,200,41)
c.hist(st.post_x,bins=bc,weights=np.full(len(st),1/len(st)/5.0),histtype="step",color=BLUE,lw=1.6)
c.hist(fl.post_x,bins=bc,weights=np.full(len(fl),1/len(fl)/5.0),histtype="step",color=ORANGE,lw=1.6)
c.set_xlabel("CV muon endpoint distance beyond the cathode  |x$_{end}$|  [cm]"); c.set_ylabel("fraction of crossers / cm (density)")
c.set_title("(c) where the CV track ends",loc="left",color=INK,fontsize=10.5)
c.legend([plt.Line2D([],[],color=ORANGE,lw=1.8),plt.Line2D([],[],color=BLUE,lw=1.8)],["flipped","crosser in both"],loc="upper right",frameon=False,fontsize=9,labelcolor=INK)

# (d) flip probability vs depth
e=ax[1,1]; allx=pd.concat([fl,st]); allx["flip"]=allx.index.isin(fl.index)
edges=np.array([0,5,10,20,30,45,60,80,110,150,250]); allx["bin"]=pd.cut(allx.post_x,edges,include_lowest=True)
g=allx.groupby("bin",observed=True).flip.agg(["sum","count"]); p=g["sum"]/g["count"]; err=np.sqrt(p*(1-p)/g["count"])
ctr=0.5*(edges[:-1]+edges[1:])
e.errorbar(ctr,100*p,yerr=100*err,fmt="o",ms=5,color=ORANGE,ecolor=ORANGE,elinewidth=1.2,capsize=0,mec=SURF,mew=1)
e.set_xlabel("CV muon endpoint distance beyond the cathode  |x$_{end}$|  [cm]"); e.set_ylabel("crossers reconstructed as non-crossers [%]")
e.set_title("(d) probability of flipping",loc="left",color=INK,fontsize=10.5); e.set_ylim(0,None)
fig.suptitle("SBND DENT vs CV: cathode crossers that lose track length  (events passing cathode + high-y/z cuts, matched by event)",color=INK,fontsize=11.5,x=0.01,ha="left")
fig.text(0.01,-0.012,f"Panel (b) correlation r = {r:.2f}; points beyond the 99th percentile not shown.",color=INK2,fontsize=8.5,ha="left")
out=str(REPO / "dent_missing_length.png"); fig.savefig(out,dpi=160,facecolor=SURF,bbox_inches="tight")
print("saved",out)
# numbers for the write-up
print("flipped: median dL %.1f cm, mean %.1f, 25/75%% %.1f/%.1f, frac with dL>20cm %.2f, >50cm %.2f"%(fl.dL.median(),fl.dL.mean(),fl.dL.quantile(.25),fl.dL.quantile(.75),(fl.dL>20).mean(),(fl.dL>50).mean()))
print("flipped: median Lpost %.1f cm; median (dL/Lpost)=%.2f; r(Lpost,dL)=%.2f"%(fl.Lpost.median(),(fl.dL/fl.Lpost).median(),r))
print("crossers(both): |end_x| median %.1f, flipped median %.1f cm"%(st.post_x.median(),fl.post_x.median()))
print(g.assign(p=(100*p).round(1)).to_string())
