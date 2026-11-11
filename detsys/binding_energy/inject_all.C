// Copies every histogram of a source file into a destination detsys file, replacing the same-named object (kOverwrite).
// Refuses on a missing name or a binning mismatch. Usage: root -l -b -q 'inject_all.C("src.root","dst.root")'
void inject_all(const char* src, const char* dst){
  TFile s(src); TFile d(dst, "UPDATE"); int n = 0;
  for(auto k : *s.GetListOfKeys()){
    TString name = ((TKey*)k)->GetName(); TH1* h = dynamic_cast<TH1*>(s.Get(name)); if(!h) continue;
    TH1* old = dynamic_cast<TH1*>(d.Get(name));
    if(!old){ printf("MISSING in destination: %s\n", name.Data()); continue; }
    if(old->GetNbinsX() != h->GetNbinsX() || fabs(old->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmax()) > 1e-9){ printf("BINNING MISMATCH: %s\n", name.Data()); continue; }
    d.cd(); TH1* c = (TH1*)h->Clone(name); c->SetDirectory(&d); c->Write(name, TObject::kOverwrite); n++;
  }
  printf("replaced %d histograms\n", n);
}
