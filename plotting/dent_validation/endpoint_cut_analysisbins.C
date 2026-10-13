// DENT/CV ratio of selected events in the analysis reco-energy binning, per proton-multiplicity channel, for three selections:
//   0: no crosser cut (v0)   1: reco_cathode_crosser == 0 (v1)   2: crosser == 0 && |reco_leading_muon_end_x| >= 4 cm
// Events are matched by (Run, Subrun, Evt, event_max_neutrino_energy) against the generated-event list of the other sample, as in
// the detsys calculator. Each sample's selection uses its own reco quantities.
#include <map>
#include <set>
#include <tuple>
#include <fstream>
using K=std::tuple<int,int,int,double>;
struct R{double E,ca,hy,cx,ex,pm;};
std::set<K> evkeys(TFile*f,const char*k){std::set<K> s; ROOT::RDataFrame d(Form("events/%s/events",k),f->GetName()); d.Foreach([&s](int r,int sr,int e,double m){s.insert(K(r,sr,e,m));},{"Run","Subrun","Evt","event_max_neutrino_energy"}); return s;}
std::multimap<K,R> load(TFile*f,const char*k){std::multimap<K,R> m; ROOT::RDataFrame d(Form("events/%s/selected",k),f->GetName());
  d.Foreach([&](int r,int sr,int e,double mm,double E,double ca,double hy,double cx,double ex,double pm){m.insert({K(r,sr,e,mm),{E,ca,hy,cx,ex,pm}});},
   {"Run","Subrun","Evt","event_max_neutrino_energy","reco_visible_energy","reco_fiducialize_cathode","reco_veto_sbnd_highy_highz","reco_cathode_crosser","reco_leading_muon_end_x","reco_proton_multiplicity"}); return m;}
void endpoint_cut_analysisbins(double XCUT=4.0){
  TFile f("/Users/mueller/data/spineosc/v1.1.0/spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root");
  auto kc=evkeys(&f,"sbnd_detvar_newcv"), kd=evkeys(&f,"sbnd_detvar_newdent");
  auto cv=load(&f,"sbnd_detvar_newcv"), dt=load(&f,"sbnd_detvar_newdent");
  std::vector<double> e1={0.3,0.4,0.46,0.52,0.58,0.64,0.7,0.76,0.82,0.88,0.94,1.0,1.06,1.13,1.2,1.3,1.5};
  std::vector<double> e2={0.4,0.5,0.59,0.67,0.75,0.83,0.91,0.99,1.07,1.16,1.26,1.37,1.5,1.7,2.0};
  std::vector<double>* ed[2]={&e1,&e2};
  std::map<std::tuple<int,int,int,int>,double> N; // (chan, mode, sample, bin)
  auto fill=[&](std::multimap<K,R>&m,std::set<K>&other,int s){ for(auto&[k,r]:m){ if(!other.count(k)) continue; if(!(r.ca==1&&r.hy==1)) continue;
      int ch = r.pm==1 ? 0 : (r.pm>1 ? 1 : -1); if(ch<0) continue; auto&E=*ed[ch]; int b=-1; for(size_t i=0;i+1<E.size();i++) if(r.E>=E[i]&&r.E<E[i+1]) b=i; if(b<0) continue;
      N[{ch,0,s,b}]++; if(r.cx==0){ N[{ch,1,s,b}]++; if(fabs(r.ex)>=XCUT) N[{ch,2,s,b}]++; } } };
  fill(cv,kd,0); fill(dt,kc,1);
  std::ofstream o("endcut_bins.csv"); o<<"channel,mode,lo,hi,ncv,ndent\n";
  for(int ch=0;ch<2;ch++) for(int m=0;m<3;m++) for(size_t i=0;i+1<ed[ch]->size();i++) o<<(ch?"nproton":"ccqelike")<<","<<m<<","<<(*ed[ch])[i]<<","<<(*ed[ch])[i+1]<<","<<N[{ch,m,0,(int)i}]<<","<<N[{ch,m,1,(int)i}]<<"\n";
  printf("done\n");
}
