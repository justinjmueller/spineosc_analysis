// Copies the BindingEnergy_* histograms from the file written by make_be_v1.py into a detsys file, replacing the existing ones
// (kOverwrite). Usage: root -l -b -q 'inject_be.C("new.root", "detvar/detsys-v1.1.0v1.root")'
void inject_be(const char* src, const char* dst){
  TFile s(src); TFile d(dst, "UPDATE"); int n = 0;
  for(auto k : *s.GetListOfKeys()){
    TString name = ((TKey*)k)->GetName(); if(!name.BeginsWith("BindingEnergy_")) continue;
    TH1* h = (TH1*)s.Get(name); TH1* old = (TH1*)d.Get(name);
    if(!old){ printf("MISSING in destination: %s\n", name.Data()); continue; }
    if(old->GetNbinsX() != h->GetNbinsX() || fabs(old->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmax()) > 1e-9){ printf("BINNING MISMATCH: %s\n", name.Data()); continue; }
    d.cd(); TH1* c = (TH1*)h->Clone(name); c->SetDirectory(&d); c->Write(name, TObject::kOverwrite); n++;
  }
  printf("replaced %d histograms\n", n);
}
