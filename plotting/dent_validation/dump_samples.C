// Dumps the energy-free, scalar branches shared by the SBND detector-access data and the DENT/CV detvar MC
// into CSV files, after one common selection. Run:  root -l -b -q 'dump_samples.C("<outdir>")'
//
// Selection (applied identically to all three samples, matching xml/detector_data_access/sbnd_run1):
//   reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1
//   && 0.647 * event_largest_flash_pe > 2000   (MC ONLY: trigger emulation; NOT applied to data)
// No energy or momentum branch is read: the data file does not carry them and none is requested here.
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <cstdio>
#include <string>
#include <vector>


#include <TString.h>
#include <TSystem.h>

// Release ntuple directory: $SPINEOSC_DATA/v1.1.0 (default $HOME/data/spineosc/v1.1.0), with a trailing slash.
static TString spineosc_dir()
{
  const char *root = gSystem->Getenv("SPINEOSC_DATA");
  TString base = (root && *root) ? TString(root) : TString(gSystem->Getenv("HOME")) + "/data/spineosc";
  return base + "/v1.1.0/";
}

void dump_one(const char *file, const char *tree_path, const std::string &out,
              const std::vector<std::string> &cols, bool flash_cut = true)
{
    TFile f(Form("%s%s", spineosc_dir().Data(), file));
    TTree *t = static_cast<TTree *>(f.Get(tree_path));
    if (!t) { printf("missing %s\n", tree_path); return; }
    std::vector<TLeaf *> lv;
    for (auto &c : cols) { TLeaf *l = t->GetLeaf(c.c_str()); if (!l) { printf("missing leaf %s\n", c.c_str()); return; } lv.push_back(l); }
    TLeaf *mult = t->GetLeaf("reco_proton_multiplicity"), *fid = t->GetLeaf("reco_fiducialize_cathode"),
          *hy = t->GetLeaf("reco_veto_sbnd_highy_highz"), *pe = t->GetLeaf("event_largest_flash_pe");
    FILE *o = fopen(out.c_str(), "w");
    for (size_t i = 0; i < cols.size(); ++i) fprintf(o, "%s%s", cols[i].c_str(), i + 1 < cols.size() ? "," : "\n");
    long n = 0, tot = t->GetEntries();
    for (long e = 0; e < tot; ++e)
    {
        t->GetEntry(e);
        if (!(mult->GetValue() >= 1 && fid->GetValue() == 1 && hy->GetValue() == 1 && (!flash_cut || 0.647 * pe->GetValue() > 2000))) continue;
        for (size_t i = 0; i < lv.size(); ++i) fprintf(o, "%.7g%s", lv[i]->GetValue(), i + 1 < lv.size() ? "," : "\n");
        ++n;
    }
    fclose(o);
    printf("%s: %ld of %ld entries selected -> %s\n", tree_path, n, tot, out.c_str());
}

void dump_samples(const char *outdir = ".")
{
    std::vector<std::string> cols = {
        "Run", "Subrun", "Evt", "event_largest_flash_pe", "event_largest_flash_time",
        "reco_cathode_crosser", "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz", "reco_neutrino2026_veto",
        "reco_electron_multiplicity", "reco_muon_multiplicity", "reco_photon_multiplicity", "reco_pion_multiplicity", "reco_proton_multiplicity",
        "reco_vertex_x", "reco_vertex_y", "reco_vertex_z",
        "reco_leading_muon_start_x", "reco_leading_muon_start_y", "reco_leading_muon_start_z",
        "reco_leading_muon_end_x", "reco_leading_muon_end_y", "reco_leading_muon_end_z",
        "reco_leading_muon_start_dir_x", "reco_leading_muon_start_dir_y", "reco_leading_muon_start_dir_z",
        "reco_leading_muon_muon_softmax", "reco_leading_muon_pion_softmax", "reco_leading_muon_primary_softmax", "reco_leading_muon_proton_softmax",
        "reco_leading_muon_calo_csda_diff",
        "reco_leading_proton_start_x", "reco_leading_proton_start_y", "reco_leading_proton_start_z",
        "reco_leading_proton_end_x", "reco_leading_proton_end_y", "reco_leading_proton_end_z",
        "reco_leading_proton_proton_softmax", "reco_leading_proton_muon_softmax", "reco_leading_proton_pion_softmax",
        "reco_leading_proton_hadron_softmax", "reco_leading_proton_mip_softmax", "reco_leading_proton_primary_softmax",
        "reco_leading_proton_secondary_softmax", "reco_leading_proton_calo_csda_diff"};
    std::string o = outdir;
    dump_one("spineosc_numu_sbnd_release-v1.1.0_detector_access_DONOTLOOK.root", "events/sbnd_bnblight_full/selected", o + "/data.csv", cols, /*flash_cut=*/false);   // the flash threshold is MC trigger emulation; data is already triggered
    dump_one("spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root", "events/sbnd_detvar_newcv/selected", o + "/cv.csv", cols);
    dump_one("spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root", "events/sbnd_detvar_newdent/selected", o + "/dent.csv", cols);
}
