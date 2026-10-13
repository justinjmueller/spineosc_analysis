#include <map>
#include <tuple>
#include <fstream>
using K=std::tuple<int,int,int,double>;
struct R{double E,ca,hy,cx,L,sx,ex,vx,tcx,tsx,tex; int n;};
std::map<K,R> load(TFile*f,const char*k){std::map<K,R> m; ROOT::RDataFrame d(Form("events/%s/selected",k),f->GetName());
  d.Foreach([&](int r,int sr,int e,double mm,double E,double ca,double hy,double cx,double L,double sx,double ex,double vx,double tcx,double tsx,double tex){
    K key(r,sr,e,mm); auto it=m.find(key); if(it!=m.end()){it->second.n++;return;} m[key]={E,ca,hy,cx,L,sx,ex,vx,tcx,tsx,tex,1};},
   {"Run","Subrun","Evt","event_max_neutrino_energy","reco_visible_energy","reco_fiducialize_cathode","reco_veto_sbnd_highy_highz","reco_cathode_crosser","reco_leading_muon_length","reco_leading_muon_start_x","reco_leading_muon_end_x","reco_vertex_x","true_cathode_crosser","true_leading_muon_start_x","true_leading_muon_end_x"}); return m;}
void dump(){
  TFile f("/Users/mueller/data/spineosc/v1.1.0/spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root");
  auto cv=load(&f,"sbnd_detvar_newcv"), dt=load(&f,"sbnd_detvar_newdent");
  std::ofstream o("/private/tmp/claude-550/-Users-mueller-Projects-spineosc/32857cee-ac1d-4b8e-b98a-74f527609260/scratchpad/flip.csv"); o<<"cls,Ecv,Edt,Lcv,Ldt,sxcv,excv,sxdt,exdt,vxcv,vxdt,tcx,tsx,tex\n"; long n=0,dup=0;
  for(auto&[k,a]:cv){ auto it=dt.find(k); if(it==dt.end()) continue; auto&b=it->second; if(a.n>1||b.n>1){dup++;continue;}
    if(!(a.ca==1&&a.hy==1&&b.ca==1&&b.hy==1)) continue;
    o<<(a.cx!=0)<<(b.cx!=0)<<","<<a.E<<","<<b.E<<","<<a.L<<","<<b.L<<","<<a.sx<<","<<a.ex<<","<<b.sx<<","<<b.ex<<","<<a.vx<<","<<b.vx<<","<<a.tcx<<","<<a.tsx<<","<<a.tex<<"\n"; n++; }
  printf("wrote %ld matched events, skipped %ld duplicate-key events\n",n,dup);
}
