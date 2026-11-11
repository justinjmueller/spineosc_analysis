#include <map>
#include <set>
#include <tuple>
using K=std::tuple<int,int,int,double>;
struct R{double E,ca,hy,cx,ex;};
std::set<K> evkeys(TFile*f,const char*k){std::set<K> s; ROOT::RDataFrame d(Form("events/%s/events",k),f->GetName()); d.Foreach([&s](int r,int sr,int e,double m){s.insert(K(r,sr,e,m));},{"Run","Subrun","Evt","event_max_neutrino_energy"}); return s;}
std::multimap<K,R> load(TFile*f,const char*k){std::multimap<K,R> m; ROOT::RDataFrame d(Form("events/%s/selected",k),f->GetName());
  d.Foreach([&](int r,int sr,int e,double mm,double E,double ca,double hy,double cx,double ex){m.insert({K(r,sr,e,mm),{E,ca,hy,cx,ex}});},
   {"Run","Subrun","Evt","event_max_neutrino_energy","reco_visible_energy","reco_fiducialize_cathode","reco_veto_sbnd_highy_highz","reco_cathode_crosser","reco_leading_muon_end_x"}); return m;}
#include <TString.h>
#include <TSystem.h>

// Release ntuple directory: $SPINEOSC_DATA/v1.1.0 (default $HOME/data/spineosc/v1.1.0), with a trailing slash.
static TString spineosc_dir()
{
  const char *root = gSystem->Getenv("SPINEOSC_DATA");
  TString base = (root && *root) ? TString(root) : TString(gSystem->Getenv("HOME")) + "/data/spineosc";
  return base + "/v1.1.0/";
}

void cut(){
  TFile f(spineosc_dir() + "spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root");
  auto kc=evkeys(&f,"sbnd_detvar_newcv"), kd=evkeys(&f,"sbnd_detvar_newdent");
  auto cv=load(&f,"sbnd_detvar_newcv"), dt=load(&f,"sbnd_detvar_newdent");
  double ed[]={0.3,0.5,0.7,0.9,1.1,1.5,2.0}; const int nb=6; double X[]={0,3,5,10,20}; const int nx=5;
  auto bin=[&](double E){for(int i=0;i<nb;i++) if(E>=ed[i]&&E<ed[i+1]) return i; return -1;};
  // mode 0: no crosser cut (v0). modes 1..nx: crosser==0 && |end_x| >= X[j]
  double N[nx+1][2][nb]={};
  auto fill=[&](std::multimap<K,R>&m,std::set<K>&other,int s){ for(auto&[k,r]:m){ if(!other.count(k)) continue; int b=bin(r.E); if(b<0||!(r.ca==1&&r.hy==1)) continue;
      N[0][s][b]++; if(r.cx==0) for(int j=0;j<nx;j++) if(fabs(r.ex)>=X[j]) N[j+1][s][b]++; } };
  fill(cv,kd,0); fill(dt,kc,1);
  printf("DENT/CV ratio by reco-energy bin\n%-30s","selection");
  for(int i=0;i<nb;i++) printf(" %.1f-%.1f",ed[i],ed[i+1]); printf(" | kept(CV)  slope[%%/GeV]\n");
  for(int m=0;m<=nx;m++){ char lab[64]; if(m==0) snprintf(lab,64,"no crosser cut (v0)"); else if(X[m-1]==0) snprintf(lab,64,"crosser==0 (v1)"); else snprintf(lab,64,"crosser==0 && |end_x|>=%.0f cm",X[m-1]);
    printf("%-30s",lab); double xs[nb],ys[nb],tot=0,tot0=0; for(int i=0;i<nb;i++){ double r=N[m][1][i]/N[m][0][i]; printf("   %.4f",r); xs[i]=0.5*(ed[i]+ed[i+1]); ys[i]=r; tot+=N[m][0][i]; tot0+=N[0][0][i]; }
    double sx=0,sy=0,sxx=0,sxy=0; for(int i=0;i<nb;i++){sx+=xs[i];sy+=ys[i];sxx+=xs[i]*xs[i];sxy+=xs[i]*ys[i];} double slope=(nb*sxy-sx*sy)/(nb*sxx-sx*sx);
    printf(" | %6.1f%%   %+6.2f\n",100*tot/tot0,100*slope); }
  // stat error on the high-energy bin for reference
  printf("\nevents in last two bins (CV) for each selection: "); for(int m=0;m<=nx;m++) printf("%.0f ",N[m][0][4]+N[m][0][5]); printf("\n");
}
