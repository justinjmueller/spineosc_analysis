/**
 * @file classes.h
 * @brief Configuration types for the detector systematics ratio writer.
 * @details Defines the plain-data structs used to describe systematics,
 * subchannel binning, and per-detector configuration for the PROfit
 * hist1d systematic ratio pipeline implemented in main.cc.
 * @author Justin Mueller (FNAL)
 */
#ifndef DETSYS_CLASSES_H
#define DETSYS_CLASSES_H

#include <string>
#include <vector>

namespace detsys
{
    /**
     * @brief Configuration for a single data-driven detector systematic.
     * @details Describes the pair of ROOT ntuples (central-value and
     * variation) that a systematic's per-bin event-count ratio is computed
     * from.
     */
    struct SystematicConfig
    {
        std::string name;            ///< Output/histogram name, e.g. "SCE_0x".
        std::string variation_file;  ///< Path to the ROOT file holding the varied sample.
        std::string variation_key;   ///< Tree key for the varied sample -> events/{key}/{selected,events}.
        std::string cv_file;         ///< Path to the ROOT file holding the central-value sample.
        std::string cv_key;          ///< Tree key for the central-value sample.
    };

    /**
     * @brief Configuration for a single reconstruction subchannel.
     * @details Defines the event selection, histogrammed branch, and bin
     * edges used to build a ratio histogram for one subchannel of a
     * systematic.
     */
    struct SubchannelConfig
    {
        std::string name;            ///< Subchannel name, e.g. "by_vis_energy".
        std::string selection_expr;  ///< RDataFrame filter expression selecting this subchannel.
        std::string branch;          ///< Branch/expression to histogram, e.g. "reco_visible_energy".
        std::vector<double> edges;   ///< Histogram bin edges (size = nbins + 1).

        /**
         * @brief Detector configurations this subchannel is computed for.
         * @details Empty (the default) means every detector, which is what
         * all the detector-independent variables want. Naming detectors
         * restricts the subchannel to those configs only -- for variables
         * that are meaningful in one detector alone, such as a vertex or
         * track-endpoint coordinate, whose binning must follow that
         * detector's geometry. Computing an SBND-ranged coordinate against
         * ICARUS data produces a histogram nothing can use.
         *
         * Names match DetectorConfig::name, i.e. "sbnd_run1",
         * "icarus_run2", "icarus_run4".
         *
         * The "generic" config is ALWAYS processed regardless of this list.
         * Its flat systematics are written once as *_generic_RunX and every
         * detector's ToroidCalibration_* / ProtonMiss_* HistVarSection in
         * the XMLs points at those same histograms, so skipping generic for
         * a restricted subchannel would break those entries for every
         * detector in that subchannel's XML.
         */
        std::vector<std::string> detectors;
    };

    /**
     * @brief Configuration for a flat (hand-set) systematic.
     * @details Represents a systematic with no CV/variation ntuples at
     * all -- just a fixed multiplicative ratio (e.g. a flat normalization
     * uncertainty) applied uniformly to every bin of every subchannel.
     */
    struct FlatSystematicConfig
    {
        std::string name;           ///< Output/histogram name, e.g. "ToroidCalibration".
        double ratio_at_one_sigma;  ///< Ratio value representing the +1 sigma variation.
    };

    /**
     * @brief Configuration for a cut-threshold (single-sample) systematic.
     * @details Unlike SystematicConfig, which divides one production by a
     * separate varied production, this varies a selection cut *within* one
     * sample: the ratio is the count passing a shifted cut over the count
     * passing the nominal cut, on the same events. Used for trigger
     * emulation, where an event fires when scale * branch > threshold and
     * the PE-to-threshold calibration `scale` carries an uncertainty
     * `sigma`. Each entry yields two one-sided universes, M1 (scale -
     * sigma) and P1 (scale + sigma).
     */
    struct CutSystematicConfig
    {
        std::string name;       ///< Base output/histogram name, e.g. "TriggerEmulation".
        std::string file;       ///< Path to the ROOT file holding the sample.
        std::string key;        ///< Tree key for the sample -> events/{key}/selected.
        std::string branch;     ///< Branch the cut is applied to, e.g. "event_largest_flash_pe".
        double scale;           ///< Nominal calibration factor applied to `branch`.
        double sigma;           ///< One-sigma uncertainty on `scale`.
        double threshold;       ///< An event passes when scale * branch > threshold.
    };

    /**
     * @brief Configuration for a single detector (or detector/run pair).
     * @details Bundles the data-driven, flat, and cut-threshold systematics
     * to be processed for one detector configuration, along with the labels
     * used to build output filenames and in-file histogram names.
     */
    struct DetectorConfig
    {
        std::string name;            ///< Used for output filenames, e.g. "icarus_run2".
        std::string detector_label;  ///< Used in the in-file histogram name, e.g. "icarus".
        std::string run_label;       ///< Used in the in-file histogram name, e.g. "Run2".
        std::vector<SystematicConfig> systematics;           ///< Data-driven systematics to process.
        std::vector<FlatSystematicConfig> flat_systematics;  ///< Hand-set systematics (no ntuples read).
        std::vector<CutSystematicConfig> cut_systematics;    ///< Cut-threshold systematics (one sample, shifted cut).
    };
}  // namespace detsys

#endif  // DETSYS_CLASSES_H
