/*
 * @file main.cc
 * @brief Main entry point for the detector systematics assessment code.
 * @details This code is intended to be used to assess the impact of detector
 * systematics and enable the implementation of the detector systematic in the
 * PROfit framework. This is accomplished by calculating the ratio of selected
 * events between the nominal and variation sample using the list of common
 * events between the two samples. The ratio is then used to create a histogram
 * that can be used to reweight the nominal sample to account for the detector
 * systematic.
 * @author Justin Mueller (FNAL)
 */
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TTree.h>

#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <algorithm>

#include "classes.h"
#include "configuration.h"

namespace detsys
{
    /**
     * @brief Fetches an open TFile for the given path, opening it if needed.
     * @details Caches open TFiles by path so a file shared by multiple
     * systematics (e.g. the SBND detvar file, or any nosyst file) is only
     * opened once. Exits the process if the file cannot be opened.
     * @param file_cache The path -> TFile* cache to look up/insert into.
     * @param path The path of the ROOT file to open.
     * @return A pointer to the (possibly newly-opened) TFile.
     */
    TFile *get_file(std::map<std::string, TFile *> &file_cache, const std::string &path)
    {
        auto it = file_cache.find(path);
        if (it != file_cache.end()) return it->second;

        TFile *file = TFile::Open(path.c_str());
        if (!file || file->IsZombie())
        {
            std::cerr << "ERROR: could not open " << path << std::endl;
            std::exit(1);
        }
        file_cache[path] = file;
        return file;
    }

    /**
     * @brief The tuple used to uniquely match an event across samples.
     * @details Run/Subrun/Evt alone is not guaranteed unique (e.g. multiple
     * neutrino interactions can share a triggered event), so
     * event_max_neutrino_energy is included as a tiebreaker. It is read from
     * both the events and selected trees.
     */
    using EventKey = std::tuple<int, int, int, double>;

    /**
     * @brief Loads the unique event key set for a sample's raw events.
     * @details Reads the events/{key}/events tree, which enumerates every
     * event that was processed (independent of whether it passed any
     * selection), and returns its identifying keys as a set.
     * @param file The ROOT file to read from.
     * @param key The tree key to read, e.g. "sbnd" or "icarus_run4".
     * @return The set of EventKeys present in the tree.
     */
    std::set<EventKey> load_event_keys(TFile *file, const std::string &key)
    {
        const std::string tree_path = "events/" + key + "/events";
        TTree *tree = static_cast<TTree *>(file->Get(tree_path.c_str()));
        if (!tree)
        {
            std::cerr << "ERROR: missing tree " << tree_path << std::endl;
            std::exit(1);
        }

        ROOT::RDataFrame df(*tree);
        std::set<EventKey> keys;
        df.Foreach(
            [&keys](int run, int subrun, int evt, double event_max_neutrino_energy) {
                keys.insert({run, subrun, evt, event_max_neutrino_energy});
            },
            {"Run", "Subrun", "Evt", "event_max_neutrino_energy"});
        return keys;
    }

    /**
     * @brief A central-value/variation pair of RDataFrame nodes restricted
     * to their common processed events.
     */
    struct CommonCollection
    {
        ROOT::RDF::RNode variation;  ///< The variation sample, restricted to events also processed by the CV sample.
        ROOT::RDF::RNode cv;         ///< The CV sample, restricted to events also processed by the variation sample.
    };

    /**
     * @brief Restricts a CV/variation sample pair to their common events.
     * @details Restricts each side's `selected` RDataFrame to rows whose
     * EventKey is present in the *other* side's raw `events` key set. This
     * does not discard events whose selection outcome differs between the
     * CV and variation samples -- only events that were never processed at
     * all on one side.
     * @param file_cache The path -> TFile* cache to fetch samples through.
     * @param variation_file Path to the ROOT file holding the varied sample.
     * @param variation_key Tree key for the varied sample.
     * @param cv_file Path to the ROOT file holding the central-value sample.
     * @param cv_key Tree key for the central-value sample.
     * @return The variation/CV RDataFrame nodes, each restricted to events
     * common to both samples.
     */
    CommonCollection build_common_collection(std::map<std::string, TFile *> &file_cache,
                                              const std::string &variation_file,
                                              const std::string &variation_key,
                                              const std::string &cv_file,
                                              const std::string &cv_key)
    {
        TFile *variation_tfile = get_file(file_cache, variation_file);
        TFile *cv_tfile = get_file(file_cache, cv_file);

        auto variation_keys = load_event_keys(variation_tfile, variation_key);
        auto cv_keys = load_event_keys(cv_tfile, cv_key);

        TTree *variation_sel_tree =
            static_cast<TTree *>(variation_tfile->Get(("events/" + variation_key + "/selected").c_str()));
        TTree *cv_sel_tree =
            static_cast<TTree *>(cv_tfile->Get(("events/" + cv_key + "/selected").c_str()));

        ROOT::RDataFrame variation_df(*variation_sel_tree);
        ROOT::RDataFrame cv_df(*cv_sel_tree);

        // Captured by value (shared_ptr-backed set copies are cheap relative to
        // the tree sizes here, and this keeps each RNode self-contained).
        auto variation_restricted = variation_df.Filter(
            [cv_keys](int run, int subrun, int evt, double event_max_neutrino_energy) {
                return cv_keys.count({run, subrun, evt, event_max_neutrino_energy}) > 0;
            },
            {"Run", "Subrun", "Evt", "event_max_neutrino_energy"});
        auto cv_restricted = cv_df.Filter(
            [variation_keys](int run, int subrun, int evt, double event_max_neutrino_energy) {
                return variation_keys.count({run, subrun, evt, event_max_neutrino_energy}) > 0;
            },
            {"Run", "Subrun", "Evt", "event_max_neutrino_energy"});

        return {ROOT::RDF::RNode(variation_restricted), ROOT::RDF::RNode(cv_restricted)};
    }

    /**
     * @brief The ratio and nominal histograms produced for one systematic/subchannel.
     */
    struct RatioHistograms
    {
        TH1D *ratio;    ///< Per-bin ratio (variation / CV), fallback 1.0 where CV is zero.
        TH1D *nominal;  ///< The CV (nominal) histogram itself.
    };

    /**
     * @brief Computes the per-bin CV/variation ratio and nominal histograms for one subchannel.
     * @details Applies the subchannel selection to both sides of a common
     * collection and histograms `branch` with the given edges. The ratio
     * histogram holds the per-bin ratio (variation / CV), with a
     * zero-CV-denominator fallback of 1.0; the nominal histogram is the CV
     * histogram itself.
     * @param common The CV/variation RDataFrame node pair to histogram.
     * @param subchannel The subchannel selection, branch, and binning to use.
     * @param hist_name The name to give the returned ratio histogram.
     * @return The newly-allocated ratio and nominal histograms (caller takes ownership of both).
     */
    RatioHistograms compute_ratio_hist(CommonCollection common, const SubchannelConfig &subchannel,
                                        const std::string &hist_name)
    {
        ROOT::RDF::RNode cv_filtered = common.cv.Filter(subchannel.selection_expr);
        ROOT::RDF::RNode variation_filtered = common.variation.Filter(subchannel.selection_expr);

        // A plain column name is histogrammed directly. Anything else, e.g.
        // "std::abs(reco_vertex_x)", is JIT-compiled into a temporary column
        // first; the expression must evaluate to double.
        std::string column = subchannel.branch;
        if (!cv_filtered.HasColumn(column))
        {
            column = "detsys_subchannel_value";
            cv_filtered = cv_filtered.Define(column, subchannel.branch);
            variation_filtered = variation_filtered.Define(column, subchannel.branch);
        }

        const int nbins = static_cast<int>(subchannel.edges.size()) - 1;
        auto cv_hist = cv_filtered.Histo1D<double>(
            ROOT::RDF::TH1DModel((hist_name + "_cv").c_str(), "", nbins, subchannel.edges.data()),
            column);
        auto variation_hist = variation_filtered.Histo1D<double>(
            ROOT::RDF::TH1DModel((hist_name + "_var").c_str(), "", nbins, subchannel.edges.data()),
            column);

        TH1D *ratio_hist = static_cast<TH1D *>(cv_hist->Clone(hist_name.c_str()));
        ratio_hist->SetDirectory(nullptr);
        TH1D *nominal_hist = static_cast<TH1D *>(cv_hist->Clone((hist_name + "_nominal").c_str()));
        nominal_hist->SetDirectory(nullptr);

        for (int bin = 1; bin <= nbins; ++bin)
        {
            const double cv_content = cv_hist->GetBinContent(bin);
            const double var_content = variation_hist->GetBinContent(bin);
            const double ratio = cv_content > 0 ? var_content / cv_content : 1.0;
            ratio_hist->SetBinContent(bin, ratio);
            ratio_hist->SetBinError(bin, 0.0);
            nominal_hist->SetBinError(bin, 0.0);
        }
        return {ratio_hist, nominal_hist};
    }

    /**
     * @brief Builds the nominal histogram scaled by a per-bin ratio.
     * @details Multiplies each bin of `nominal` by the corresponding bin of
     * `ratio`, giving the nominal shape reweighted by the systematic's
     * per-bin effect. `nominal` and `ratio` must share the same binning.
     * @param nominal The nominal (CV) histogram.
     * @param ratio The per-bin ratio histogram.
     * @param hist_name The name to give the returned histogram.
     * @return The newly-allocated scaled histogram (caller takes ownership).
     */
    TH1D *scale_hist_by_ratio(const TH1D *nominal, const TH1D *ratio, const std::string &hist_name)
    {
        TH1D *scaled_hist = static_cast<TH1D *>(nominal->Clone(hist_name.c_str()));
        scaled_hist->SetDirectory(nullptr);
        for (int bin = 1; bin <= nominal->GetNbinsX(); ++bin)
        {
            scaled_hist->SetBinContent(bin, nominal->GetBinContent(bin) * ratio->GetBinContent(bin));
            scaled_hist->SetBinError(bin, 0.0);
        }
        return scaled_hist;
    }

    /**
     * @brief Builds a ratio histogram for a flat (hand-set) systematic.
     * @details Every bin of every subchannel gets the same hand-set
     * ratio_at_one_sigma value -- unlike compute_ratio_hist(), there is no
     * CV or variation sample to read, so no ROOT file/TTree access occurs.
     * @param subchannel The subchannel binning to use.
     * @param ratio_at_one_sigma The flat ratio value to fill every bin with.
     * @param hist_name The name to give the returned histogram.
     * @return The newly-allocated ratio histogram (caller takes ownership).
     */
    TH1D *build_flat_ratio_hist(const SubchannelConfig &subchannel, double ratio_at_one_sigma,
                                 const std::string &hist_name)
    {
        const int nbins = static_cast<int>(subchannel.edges.size()) - 1;
        TH1D *ratio_hist = new TH1D(hist_name.c_str(), "", nbins, subchannel.edges.data());
        ratio_hist->SetDirectory(nullptr);
        for (int bin = 1; bin <= nbins; ++bin)
        {
            ratio_hist->SetBinContent(bin, ratio_at_one_sigma);
            ratio_hist->SetBinError(bin, 0.0);
        }
        return ratio_hist;
    }

    /**
     * @brief Writes one or more histograms to a single (freshly-created) ROOT file.
     * @details Each histogram is renamed to its paired in-file name before
     * being written. All histograms destined for the same out_path must be
     * passed in one call, since the file is opened in RECREATE mode.
     * @param hists The (histogram, in-file name) pairs to write.
     * @param out_path The path of the ROOT file to create.
     */
    void write_systematics(const std::vector<std::pair<TH1D *, std::string>> &hists, const std::string &out_path)
    {
        TFile out(out_path.c_str(), "RECREATE");
        for (const auto &[hist, name] : hists)
        {
            hist->SetName(name.c_str());
            hist->Write();
        }
        out.Close();
    }

    /**
     * @brief Writes a single histogram to its own ROOT file.
     * @param hist The histogram to write.
     * @param hist_name The in-file name to give the histogram.
     * @param out_path The path of the ROOT file to create.
     */
    void write_systematic(TH1D *hist, const std::string &hist_name, const std::string &out_path)
    {
        write_systematics({{hist, hist_name}}, out_path);
    }

    /**
     * @brief Builds the in-file/output histogram name for a systematic.
     * @param syst_name The systematic's name, e.g. "SCE_0x".
     * @param detector The detector configuration the systematic belongs to.
     * @return "<systname>_<detector>_<run>", e.g. "SCE_0x_sbnd_Run1".
     */
    std::string build_hist_name(const std::string &syst_name, const DetectorConfig &detector)
    {
        return syst_name + "_" + detector.detector_label + "_" + detector.run_label;
    }

    /**
     * @brief Builds the output ROOT file path for one (systematic, subchannel) pair.
     * @param output_dir The directory to write output ROOT files into.
     * @param detector_name The detector's output-filename name, e.g. "icarus_run2".
     * @param syst_name The systematic's name, e.g. "SCE_0x".
     * @param subchannel_name The subchannel's name, e.g. "by_vis_energy".
     * @return "<output_dir>/<detector_name>_<syst_name>_<subchannel_name>.root".
     */
    std::string build_output_path(const std::string &output_dir, const std::string &detector_name,
                                   const std::string &syst_name, const std::string &subchannel_name)
    {
        return output_dir + "/" + detector_name + "_" + syst_name + "_" + subchannel_name + ".root";
    }

    /**
     * @brief Formats a histogram's bin contents as a "[a, b, c]" string.
     * @param hist The histogram to read bin contents from.
     * @return The formatted bin-contents string.
     */
    std::string format_bin_contents(const TH1D *hist)
    {
        std::ostringstream out;
        out << "[";
        for (int bin = 1; bin <= hist->GetNbinsX(); ++bin)
        {
            out << hist->GetBinContent(bin);
            if (bin != hist->GetNbinsX()) out << ", ";
        }
        out << "]";
        return out.str();
    }

    /**
     * @brief Writes one or more related histograms, prints a summary line, and frees them.
     * @param hists The (histogram, in-file name) pairs to emit (each deleted before returning).
     * @param out_path The path of the ROOT file to create.
     * @param summary The per-subchannel summary text to print after the first histogram's name.
     */
    void emit_systematics(const std::vector<std::pair<TH1D *, std::string>> &hists, const std::string &out_path,
                           const std::string &summary)
    {
        write_systematics(hists, out_path);
        std::cout << hists.front().second << " " << summary << std::endl;
        for (const auto &[hist, name] : hists) delete hist;
    }

    /**
     * @brief Writes a single histogram, prints its summary line, and frees it.
     * @param hist The histogram to emit (deleted before returning).
     * @param hist_name The in-file name to give the histogram.
     * @param out_path The path of the ROOT file to create.
     * @param summary The per-subchannel summary text to print after hist_name.
     */
    void emit_systematic(TH1D *hist, const std::string &hist_name, const std::string &out_path,
                          const std::string &summary)
    {
        emit_systematics({{hist, hist_name}}, out_path, summary);
    }

    /**
     * @brief Builds the in-file histogram name for a cut systematic universe.
     * @details Cut systematics follow "<name><sign>_<detector name>" (e.g.
     * "TriggerEmulationM1_sbnd_run1") rather than build_hist_name's
     * "<name>_<label>_<run>", preserving the names the XML configs'
     * <HistVarSection> entries already reference.
     * @param syst_name The systematic's base name, e.g. "TriggerEmulation".
     * @param sign The universe tag, "M1" or "P1".
     * @param detector The detector configuration the systematic belongs to.
     * @return "<systname><sign>_<detector name>".
     */
    std::string build_cut_hist_name(const std::string &syst_name, const std::string &sign,
                                     const DetectorConfig &detector)
    {
        return syst_name + sign + "_" + detector.name;
    }

    /**
     * @brief Formats a cut expression "(scale * branch) > threshold".
     * @details Written at full double precision so the emitted RDataFrame
     * expression reproduces the configured calibration exactly.
     * @param branch The branch the cut is applied to.
     * @param scale The calibration factor to multiply the branch by.
     * @param threshold The threshold the scaled branch must exceed.
     * @return The filter expression string.
     */
    std::string build_cut_expression(const std::string &branch, double scale, double threshold)
    {
        std::ostringstream out;
        out.precision(17);
        out << "(" << scale << " * " << branch << ") > " << threshold;
        return out.str();
    }

    /**
     * @brief Computes the ratio and nominal histograms for one cut-systematic universe.
     * @details Both numerator and denominator come from the same `selected`
     * tree: the denominator is the nominal cut (scale), the numerator the
     * shifted cut (scale +/- sigma). Because it is one sample, no common-event
     * restriction is needed -- the two sides are subsets of the same rows. The
     * pair is handed to compute_ratio_hist, which applies the subchannel
     * selection and performs the division, so the zero-CV fallback and the
     * nominal/scaled companions behave identically to the two-sample path.
     * @param file_cache The path -> TFile* cache to fetch the sample through.
     * @param syst The cut systematic being processed.
     * @param shifted_scale The calibration factor for this universe.
     * @param subchannel The subchannel selection, branch, and binning to use.
     * @param hist_name The name to give the returned ratio histogram.
     * @return The newly-allocated ratio and nominal histograms (caller owns both).
     */
    RatioHistograms compute_cut_ratio_hist(std::map<std::string, TFile *> &file_cache,
                                            const CutSystematicConfig &syst, double shifted_scale,
                                            const SubchannelConfig &subchannel, const std::string &hist_name)
    {
        TFile *file = get_file(file_cache, syst.file);
        const std::string tree_path = "events/" + syst.key + "/selected";
        TTree *tree = static_cast<TTree *>(file->Get(tree_path.c_str()));
        if (!tree)
        {
            std::cerr << "ERROR: missing tree " << tree_path << std::endl;
            std::exit(1);
        }

        ROOT::RDataFrame df(*tree);
        auto cv_node = df.Filter(build_cut_expression(syst.branch, syst.scale, syst.threshold));
        auto variation_node = df.Filter(build_cut_expression(syst.branch, shifted_scale, syst.threshold));

        CommonCollection common{ROOT::RDF::RNode(variation_node), ROOT::RDF::RNode(cv_node)};
        return compute_ratio_hist(common, subchannel, hist_name);
    }

    /**
     * @brief Processes every systematic configured for a single detector.
     * @details For each data-driven systematic, computes the ratio, nominal,
     * and ratio-scaled-nominal histograms for every subchannel and writes
     * all three to one output file; for each flat systematic, computes just
     * the ratio histogram; for each cut systematic, emits the same three
     * histograms per one-sided universe (M1 and P1). Every output file gets
     * a summary line of its ratio histogram's bin contents.
     * @param detector The detector configuration to process.
     * @param output_dir The directory to write output ROOT files into.
     * @param file_cache The path -> TFile* cache to fetch samples through.
     */
    /**
     * @brief Selects the subchannels that should be computed for one detector.
     * @details A subchannel with an empty SubchannelConfig::detectors list is
     * detector-independent and is always kept. A non-empty list restricts the
     * subchannel to the named detector configurations, so that a variable
     * whose binning follows one detector's geometry -- a vertex or endpoint
     * coordinate, say -- is not computed against another detector's data,
     * where the resulting histogram could never be used.
     *
     * The "generic" configuration is exempt: it writes the flat systematics
     * once as *_generic_RunX, and every detector's ToroidCalibration_* and
     * ProtonMiss_* HistVarSection in the XMLs resolves to those same
     * histograms. Skipping generic for a restricted subchannel would leave
     * those entries dangling for every detector in that subchannel's XML.
     *
     * @param all Every configured subchannel.
     * @param detector The detector configuration being processed.
     * @return The subset to process for this detector.
     */
    std::vector<SubchannelConfig> select_subchannels(const std::vector<SubchannelConfig> &all,
                                                     const DetectorConfig &detector)
    {
        std::vector<SubchannelConfig> kept;
        for (const auto &subchannel : all)
        {
            const bool unrestricted = subchannel.detectors.empty();
            const bool is_generic = detector.name == "generic";
            const bool named = std::find(subchannel.detectors.begin(), subchannel.detectors.end(),
                                         detector.name) != subchannel.detectors.end();
            if (unrestricted || is_generic || named)
                kept.push_back(subchannel);
        }
        return kept;
    }

    void process_detector(const DetectorConfig &detector, const std::string &output_dir,
                           std::map<std::string, TFile *> &file_cache)
    {
        // Filtered once here rather than guarded inside each of the three
        // systematic loops below, so the data-driven, flat and cut paths
        // cannot drift apart in which subchannels they honour.
        const auto subchannels = select_subchannels(build_subchannels(), detector);
        if (subchannels.size() != build_subchannels().size())
        {
            std::cout << "[" << detector.name << "] processing " << subchannels.size() << " of "
                      << build_subchannels().size() << " subchannels (detector-restricted ones skipped)"
                      << std::endl;
        }

        for (const auto &syst : detector.systematics)
        {
            CommonCollection common = build_common_collection(
                file_cache, syst.variation_file, syst.variation_key, syst.cv_file, syst.cv_key);
            const std::string hist_name = build_hist_name(syst.name, detector);

            for (const auto &subchannel : subchannels)
            {
                RatioHistograms hists = compute_ratio_hist(common, subchannel, hist_name);
                TH1D *scaled_hist = scale_hist_by_ratio(hists.nominal, hists.ratio, hist_name + "_scaled");

                const std::string out_path = build_output_path(output_dir, detector.name, syst.name, subchannel.name);
                const std::string summary = subchannel.name + " " + format_bin_contents(hists.ratio);

                emit_systematics(
                    {
                        {hists.ratio, hist_name},
                        {hists.nominal, hist_name + "_nominal"},
                        {scaled_hist, hist_name + "_scaled"},
                    },
                    out_path, summary);
            }
        }

        for (const auto &flat_syst : detector.flat_systematics)
        {
            const std::string hist_name = build_hist_name(flat_syst.name, detector);

            for (const auto &subchannel : subchannels)
            {
                TH1D *ratio_hist = build_flat_ratio_hist(subchannel, flat_syst.ratio_at_one_sigma, hist_name);
                const std::string out_path =
                    build_output_path(output_dir, detector.name, flat_syst.name, subchannel.name);

                std::ostringstream summary;
                summary << subchannel.name << " (flat ratio = " << flat_syst.ratio_at_one_sigma << ")";
                emit_systematic(ratio_hist, hist_name, out_path, summary.str());
            }
        }

        for (const auto &cut_syst : detector.cut_systematics)
        {
            // Each cut systematic yields two one-sided universes: the calibration
            // shifted down (M1) and up (P1) by one sigma, each divided by the
            // nominal-cut counts.
            const std::vector<std::pair<std::string, double>> universes = {
                {"M1", cut_syst.scale - cut_syst.sigma},
                {"P1", cut_syst.scale + cut_syst.sigma},
            };

            for (const auto &[sign, shifted_scale] : universes)
            {
                const std::string syst_name = cut_syst.name + sign;
                const std::string hist_name = build_cut_hist_name(cut_syst.name, sign, detector);

                for (const auto &subchannel : subchannels)
                {
                    RatioHistograms hists =
                        compute_cut_ratio_hist(file_cache, cut_syst, shifted_scale, subchannel, hist_name);
                    TH1D *scaled_hist = scale_hist_by_ratio(hists.nominal, hists.ratio, hist_name + "_scaled");

                    const std::string out_path =
                        build_output_path(output_dir, detector.name, syst_name, subchannel.name);
                    const std::string summary = subchannel.name + " " + format_bin_contents(hists.ratio);

                    emit_systematics(
                        {
                            {hists.ratio, hist_name},
                            {hists.nominal, hist_name + "_nominal"},
                            {scaled_hist, hist_name + "_scaled"},
                        },
                        out_path, summary);
                }
            }
        }
    }
}  // namespace detsys

/**
 * @brief Program entry point.
 * @details Builds the requested detector configuration(s) and writes a
 * ratio histogram ROOT file per (detector, systematic, subchannel) triple
 * into the output directory.
 * @param argc Argument count.
 * @param argv Argument values: [detector selection] [output directory].
 * Detector selection is one of "sbnd", "icarus", "generic", or "all"
 * (default "all"); output directory defaults to "./profit_systematics".
 * @return 0 on success, 1 if the detector selection is not recognized.
 */
int main(int argc, char **argv)
{
    const std::string which = argc > 1 ? argv[1] : "all";
    const std::string output_dir = argc > 2 ? argv[2] : "../profit_systematics";

    gSystem->mkdir(output_dir.c_str(), /*recursive=*/true);

    std::map<std::string, TFile *> file_cache;

    if (which == "sbnd" || which == "all")
    {
        detsys::process_detector(detsys::build_sbnd_config(), output_dir, file_cache);
    }
    if (which == "icarus" || which == "all")
    {
        detsys::process_detector(detsys::build_icarus_run2_config(), output_dir, file_cache);
        detsys::process_detector(detsys::build_icarus_run4_config(), output_dir, file_cache);
    }
    if (which == "generic" || which == "all")
    {
        detsys::process_detector(detsys::build_generic_config(), output_dir, file_cache);
    }
    if (which != "sbnd" && which != "icarus" && which != "generic" && which != "all")
    {
        std::cerr << "Unknown detector selection: " << which << " (expected sbnd, icarus, generic, or all)"
                  << std::endl;
        return 1;
    }

    for (auto &[path, file] : file_cache)
    {
        file->Close();
    }

    return 0;
}
