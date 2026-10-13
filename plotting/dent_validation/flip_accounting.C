#include <map>
#include <set>
#include <tuple>
using K=std::tuple<int,int,int,double>;
struct R{double E,ca,hy,cx;};
std::set<K> evkeys(TFile*f,const char*k){std::set<K> s; ROOT::RDataFrame d(Form("events/%s/events",k),f->GetName()); d.Foreach([&s](int r,int sr,int e,double m){s.insert(K(r,sr,e,m));},{"Run","Subrun","Evt","event_max_neutrino_energy"}); return s;}
std::multimap<K,R> load(TFile*f,const char*k){std::multimap<K,R> m; ROOT::RDataFrame d(Form("events/%s/selected",k),f->GetName());
  d.Foreach([&](int r,int sr,int e,double mm,double E,double ca,double hy,double cx){m.insert({K(r,sr,e,mm),{E,ca,hy,cx}});},
   {"Run","Subrun","Evt","event_max_neutrino_energy","reco_visible_energy","reco_fiducialize_cathode","reco_veto_sbnd_highy_highz","reco_cathode_crosser"}); return m;}
void mig2(){
  TFile f("/Users/mueller/data/spineosc/v1.1.0/spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root");
  auto kc=evkeys(&f,"sbnd_detvar_newcv"), kd=evkeys(&f,"sbnd_detvar_newdent");
  auto cv=load(&f,"sbnd_detvar_newcv"), dt=load(&f,"sbnd_detvar_newdent");
  double ed[]={0.3,0.5,0.7,0.9,1.1,1.5,2.0}; const int nb=6;
  auto bin=[&](double E){for(int i=0;i<nb;i++) if(E>=ed[i]&&E<ed[i+1]) return i; return -1;};
  double Ncv[nb]={},Ndt[nb]={},flip[nb]={},Ncv_all[nb]={},Ndt_all[nb]={};
  // CV-side crosser lookup (any selected CV row with same key that passes v0 selection and is a crosser)
  std::set<K> cvCross; for(auto&[k,r]:cv) if(r.ca==1&&r.hy==1&&r.cx!=0) cvCross.insert(k);
  for(auto&[k,r]:cv){ if(!kd.count(k)) continue; int b=bin(r.E); if(b<0||!(r.ca==1&&r.hy==1)) continue; Ncv_all[b]++; if(r.cx==0) Ncv[b]++; }
  for(auto&[k,r]:dt){ if(!kc.count(k)) continue; int b=bin(r.E); if(b<0||!(r.ca==1&&r.hy==1)) continue; Ndt_all[b]++; if(r.cx==0){ Ndt[b]++; if(cvCross.count(k)) flip[b]++; } }
  printf("%-8s | %8s %8s %7s | %8s %8s %7s | %6s %20s\n","E [GeV]","CV all","DENT all","ratio","CV x=0","DENT x=0","ratio","flips","ratio w/o flip-ins");
  for(int i=0;i<nb;i++) printf("%.1f-%-4.1f | %8.0f %8.0f %7.4f | %8.0f %8.0f %7.4f | %6.0f %20.4f\n",ed[i],ed[i+1],Ncv_all[i],Ndt_all[i],Ndt_all[i]/Ncv_all[i],Ncv[i],Ndt[i],Ndt[i]/Ncv[i],flip[i],(Ndt[i]-flip[i])/Ncv[i]);
}
