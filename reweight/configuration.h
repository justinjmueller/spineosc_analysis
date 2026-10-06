/**
 * @file configuration.h
 * @brief Configurable inputs for the event reweighting code.
 * @details Defines the per-sample post-selection cuts whose combined value is
 * added to the copied `selected` tree by main.cc. This is the file to edit
 * when workshopping what the selection means: the PE calibration and the
 * thresholds live in one table, and each cut is a named, separately
 * toggleable expression.
 * @author Justin Mueller (FNAL)
 */
#ifndef REWEIGHT_CONFIGURATION_H
#define REWEIGHT_CONFIGURATION_H

#include <sstream>
#include <string>
#include <vector>

namespace reweight
{
    /**
     * @brief The flash-PE trigger-emulation cut for one sample.
     * @details The cut is `scale * event_largest_flash_pe > threshold`.
     *
     * The threshold belongs to the detector and never changes between
     * samples of it. The scale is the calibration that puts simulated PE on
     * the same footing as data, so it applies to the MC samples and is 1.0
     * for the off-beam data ones, where the PE is already on that footing.
     */
    struct PeCut
    {
        double scale;      ///< PE calibration, 1.0 for off-beam data.
        double threshold;  ///< Detector's PE threshold, the same for all its samples.
    };

    /**
     * @brief Everything about one sample the selection depends on.
     * @details Keyed by the tree key the tool is run against, i.e. the {key}
     * in events/{key}/selected.
     *
     * WHY THIS IS PER-SAMPLE RATHER THAN PER-EVENT. The obvious alternative
     * is one expression that decides the detector from true_detector and Run,
     * as the XML's `additional_weight` does. That cannot work here: in every
     * off-beam tree true_detector is NaN for all events, so a
     * `true_detector == 1` style predicate is false everywhere and would
     * reject the entire sample, while an EE filter written as an implication
     * over those predicates would silently not apply. Run numbers cannot
     * stand in either -- they identify the sample, not the detector (SBND
     * beam runs 3-383, SBND off-beam 18250-18503).
     *
     * The tool is pointed at one sample at a time and is told which, so the
     * detector and the sample kind are known before an event is read. This
     * table records that, and the cut expressions stay free of any per-event
     * detector inference.
     */
    struct SampleConfig
    {
        std::string key;             ///< Tree key under events/, e.g. "sbnd_offbeam".
        PeCut pe;                    ///< The sample's flash-PE cut.
        bool apply_filter_ee;        ///< Whether the east-east containment filter applies.
        bool expects_weight_trees;   ///< Whether multisim/multisigma trees should be present.
    };

    /**
     * @brief One named post-selection cut.
     * @details The expression is handed to a TTreeFormula built against the
     * input `selected` tree, so it may use any branch of that tree and any
     * function TFormula understands. An event survives the cut when the
     * expression evaluates nonzero.
     */
    struct SelectionCut
    {
        std::string name;        ///< Short label, used in the per-cut summary.
        std::string expression;  ///< TFormula expression; nonzero means the event survives.
        bool enabled;            ///< Whether the cut participates at all.
    };

    /**
     * @brief The selection branch written into the copied tree.
     * @details The branch holds 1.0 for an event passing every enabled cut and
     * 0.0 otherwise, so it can be applied downstream as a weight that zeroes
     * the events the cuts reject. It is stored as a double rather than a flag
     * so a cut list that grows a genuinely weight-valued term later does not
     * need the branch type to change underneath it.
     */
    struct SelectionConfig
    {
        std::string name;                ///< Name of the branch added to the copied tree.
        std::vector<SelectionCut> cuts;  ///< Cuts, combined with a logical AND.
    };

    /**
     * @brief The samples the tool knows how to build a selection for.
     * @details Thresholds are per detector: SBND 2000, ICARUS Run 2 6000,
     * ICARUS Run 4 3000. Within a detector, the beam MC, dirt MC and off-beam
     * data samples share that threshold and differ only in the calibration
     * scale, which is 1.0 for the off-beam data.
     *
     * The EE containment filter is an ICARUS Run 2 matter and applies to its
     * data and MC alike, off-beam included, which is what the XML does.
     *
     * The off-beam samples carry no multisim/multisigma trees -- they are
     * data-driven cosmics, so there are no generator weights to vary -- and
     * are marked accordingly so their absence is expected rather than an
     * error.
     *
     * NOT COVERED YET: the `selected_nonmatched` cosmic trees. They live
     * inside a beam sample's key rather than having one of their own, and the
     * tool copies `selected` alone, so they need the copier to grow a notion
     * of which tree carries the selection before they can be configured here.
     *
     * @return The per-sample configurations.
     */
    inline std::vector<SampleConfig> build_samples()
    {
        return
        {
            //  key                    PE scale  threshold   EE filter  weight trees
            // v1.0.10 names the SBND beam sample "sbnd"; v1.1.0 splits it into
            // "sbnd_large" and "sbnd_small". All three are the same detector
            // and take the same calibration and threshold.
            {"sbnd",                 {0.647,     2000.0},    false,     true},
            {"sbnd_large",           {0.647,     2000.0},    false,     true},
            {"sbnd_small",           {0.647,     2000.0},    false,     true},
            {"sbnd_dirt",            {0.647,     2000.0},    false,     true},
            {"sbnd_offbeam",         {1.0,       2000.0},    false,     false},

            {"icarus_run2",          {0.615,     6000.0},    true,      true},
            {"icarus_dirt_run2",     {0.615,     6000.0},    true,      true},
            {"icarus_offbeam_run2",  {1.0,       6000.0},    true,      false},

            {"icarus_run4",          {0.361,     3000.0},    false,     true},
            {"icarus_dirt_run4",     {0.361,     3000.0},    false,     true},
            {"icarus_offbeam_run4",  {1.0,       3000.0},    false,     false},
        };
    }

    /**
     * @brief Formats a sample's PE cut as a TFormula expression.
     * @details Written at full double precision so the emitted expression
     * reproduces the configured calibration exactly rather than a rounded
     * rendering of it. A unit scale is left out of the expression entirely,
     * both because multiplying by one earns nothing and because seeing the
     * bare branch in the log is how one confirms at a glance that a data
     * sample is being cut on uncalibrated PE.
     * @param pe The sample's PE cut.
     * @return The filter expression string.
     */
    inline std::string format_pe_cut(const PeCut &pe)
    {
        std::ostringstream out;
        out.precision(17);
        out << "(";
        if (pe.scale != 1.0) out << pe.scale << " * ";
        out << "event_largest_flash_pe > " << pe.threshold << ")";
        return out.str();
    }

    /**
     * @brief Builds the post-selection cuts for one sample.
     * @details These reproduce the post-selection cuts that
     * xml/contours/contours_base.xml applies inside every `additional_weight`
     * expression, so that they can eventually be replaced there by this one
     * branch. The truth-level and proton-multiplicity terms in those
     * expressions are NOT included: they are what separates one subchannel
     * from another, not a cut every event must pass.
     *
     * The `enabled` flags are the workshopping knob -- turning one off leaves
     * the others untouched and reports the difference in the per-cut summary.
     *
     * @param sample The sample being processed.
     * @return The selection configuration to evaluate.
     */
    inline SelectionConfig build_selection(const SampleConfig &sample)
    {
        SelectionConfig selection;
        selection.name = "selected";

        // Flash-PE trigger emulation, against this sample's own calibration.
        selection.cuts.push_back({"pe_cut", format_pe_cut(sample.pe), true});

        // Containment: the interaction is inside the fiducial volume,
        // cathode-crossers included. Applied to every sample.
        selection.cuts.push_back({"fiducialize_cathode", "(reco_fiducialize_cathode == 1)", true});

        // East-east containment filter, ICARUS Run 2 only. Because the sample
        // is already known to be Run 2, this is the bare cut -- no predicate
        // guarding it, and so nothing for a NaN true_detector to defeat.
        if (sample.apply_filter_ee)
            selection.cuts.push_back({"filter_ee", "(reco_user_containment_cut_filterEE == 1)", true});

        return selection;
    }

    /**
     * @brief A systematic whose central value is folded into the event weight.
     * @details The multisigma trees store one weight per knob per event, and
     * the knob values themselves in a parallel `<branch>_sigma` vector. A
     * systematic that carries a genuine central-value correction has a knob at
     * zero sigma whose weight is NOT 1: that entry is the correction, and the
     * +-sigma entries are variations around it.
     *
     * Folding that entry into the event weight is what `reweight_aFF` does in
     * the GUMP loader (analysis_village/gump/loaddf.py), where cvwgt is the
     * product of one `cv` column per systematic in `xsec_cv_rwgt`. Doing it
     * here, once, in the file itself means the fits inherit the corrected CV
     * without any of them having to know about it.
     */
    struct CvReweight
    {
        std::string branch;  ///< Branch in the multisigma tree holding the per-knob weights.
        double knob;         ///< Knob value whose weight is the central value, normally 0.
        bool required;       ///< Whether an absent branch or knob is fatal rather than a warning.
    };

    /**
     * @brief A systematic whose stored universes are divided by their CV knob.
     * @details Applying a systematic's CV correction through the weight (above)
     * and ALSO leaving it in the universes would count it twice. Dividing every
     * knob of the branch by its own zero-sigma entry leaves the spline at
     * exactly 1 there, so the component still varies -- it stays fittable --
     * but no longer moves the central value.
     *
     * This is what the four ZExpPCA MvA components need: they share one
     * bit-identical CV column, so PROfit currently applies the same axial
     * form-factor correction once per component, i.e. to the fourth power.
     */
    struct CvNormalise
    {
        std::string branch;  ///< Branch in the multisigma tree to renormalise.
        double knob;         ///< Knob whose weight becomes the divisor, normally 0.
    };

    /**
     * @brief How the central value is rebuilt and which systematics survive it.
     */
    struct ReweightConfig
    {
        std::string weight_branch;             ///< Name of the weight column added to `selected`.
        std::vector<CvReweight> reweights;     ///< Folded into the weight; drop these from the fit.
        std::vector<CvNormalise> normalise;    ///< Renormalised in place; keep these in the fit.
    };

    /**
     * @brief The central-value reweighting applied to a copied sample.
     * @details The three entries under `reweights` are the GUMP loader's
     * `xsec_cv_rwgt` list. Folding them in here and dropping them from the XML
     * allowlist reproduces that analysis's central value while leaving the
     * ZExpPCA components free to be fitted.
     *
     * MEASURED AGAINST spineosc_numu_sbnd_v1.1.0_wsyst.root, the three do not
     * store their central value in the same place, which is why the knob
     * differs per entry:
     *   ZExpPCAWeighter_SBN_v3_MvA_b1   knobs [1,-1,2,-2,3,-3,0]. The zero
     *                                   entry is the correction: mean 1.085,
     *                                   63% of events not 1.
     *   CCQEXSecCorr_SBN_v3_CCQEXSecCorr  a SINGLE knob, labelled +1 sigma.
     *                                   With nothing to interpolate between it
     *                                   is not a variation at all but a flat
     *                                   correction, so that one entry is the
     *                                   central value. Every weight is exactly
     *                                   1 in this input, making it a no-op
     *                                   here; it is configured so it starts
     *                                   contributing as soon as an input
     *                                   carries real values.
     *   GENIEReWeight_SBN_v3_FrKin_PiProFix_N  knobs [1,-1,2,-2,3,-3], no zero
     *                                   entry, so the +1 sigma point is taken
     *                                   as the centre. Mean 1.0003 with 0.5%
     *                                   of events not 1.
     *
     * Taking a +1 sigma entry as a central value is a deliberate choice, not
     * an approximation of a missing zero knob: these two are reweighted and
     * NOT fitted, so there is no spline through them whose centre could
     * disagree. Only the ZExpPCA components stay fittable, and only they are
     * renormalised below.
     *
     * All four ZExpPCA components are renormalised, b1 included: its CV is
     * taken through the weight instead, so leaving it in the universes as well
     * would reapply it.
     * @return The reweighting configuration.
     */
    inline ReweightConfig build_reweight()
    {
        return
        {
            "cvwgt",
            {
                //  branch                                   knob   required
                {"ZExpPCAWeighter_SBN_v3_MvA_b1",            0.0,   true},
                {"CCQEXSecCorr_SBN_v3_CCQEXSecCorr",         1.0,   true},
                {"GENIEReWeight_SBN_v3_FrKin_PiProFix_N",    1.0,   true},
            },
            {
                {"ZExpPCAWeighter_SBN_v3_MvA_b1", 0.0},
                {"ZExpPCAWeighter_SBN_v3_MvA_b2", 0.0},
                {"ZExpPCAWeighter_SBN_v3_MvA_b3", 0.0},
                {"ZExpPCAWeighter_SBN_v3_MvA_b4", 0.0},
            },
        };
    }

}  // namespace reweight

#endif  // REWEIGHT_CONFIGURATION_H
