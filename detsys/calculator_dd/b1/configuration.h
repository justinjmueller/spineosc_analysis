/**
 * @file configuration.h
 * @brief Configurable inputs for the detector systematics ratio writer.
 * @details Defines the ntuple base paths, subchannel binning, and
 * per-detector systematic lists consumed by the processing pipeline in
 * main.cc. This is the file to edit when adding/removing a systematic,
 * changing binning, or pointing at a different production's ntuples.
 * @author Justin Mueller (FNAL)
 */
#ifndef DETSYS_CONFIGURATION_H
#define DETSYS_CONFIGURATION_H

#include <string>
#include <vector>

#include "classes.h"

namespace detsys
{
    constexpr const char *k_sbnd_path = "/Users/mueller/data/spineosc/v1.1.0";
    constexpr const char *k_icarus_path = "/Users/mueller/data/spineosc/v1.1.0";

    // v1.1.0 regained the SBND WireMod samples, so the v1.0.10 fallback that
    // stood in for them is retired.

    /**
     * @brief Joins a base directory and a filename into a single path.
     * @param base_dir The base directory.
     * @param file_name The filename to append.
     * @return "<base_dir>/<file_name>".
     */
    inline std::string build_file_path(const std::string &base_dir, const std::string &file_name)
    {
        return base_dir + "/" + file_name;
    }

    /**
     * @brief Builds the list of reconstruction subchannels to process.
     * @details Subchannel definitions are identical across every detector
     * configuration.
     * @return The list of subchannel configurations.
     */
    inline std::vector<SubchannelConfig> build_subchannels()
    {
        return
        {
            {"dd_muon_calo_csda", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_leading_muon_calo_csda_diff", {-0.321, -0.182, -0.164, -0.153, -0.146, -0.139, -0.134, -0.129, -0.125, -0.121, -0.117, -0.113, -0.109, -0.105, -0.101, -0.0957, -0.0903, -0.0838, -0.0749, -0.0596, 0.00676}, {"sbnd_run1"}},
            {"dd_muon_abs_endx", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "std::abs(reco_leading_muon_end_x)", {0.75, 7.35, 14.9, 23, 30.8, 38.9, 47, 55, 63.2, 71.9, 80.5, 89.2, 98.6, 108, 118, 128, 139, 150, 163, 176, 188}, {"sbnd_run1"}},
            {"dd_crosser_flag", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_cathode_crosser", {-0.5, 0.5, 1.5}, {"sbnd_run1"}},
            {"dd_muon_endx_full", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_leading_muon_end_x", {-187, -167, -147, -128, -111, -95, -79.7, -64.7, -50.5, -36.8, -23, -9.45, 2.55, 15.1, 27.5, 41, 53.9, 67.4, 81.5, 95.9, 112, 128, 146, 167, 187}, {"sbnd_run1"}},
            {"dd_muon_proton_score", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_leading_muon_proton_softmax", {0, 1.85e-05, 0.000129, 0.000372, 0.000809, 0.00155, 0.00281, 0.00494, 0.0086, 0.0146, 0.0242, 0.04, 0.0671, 0.118, 0.219, 0.438, 1}, {"sbnd_run1"}}
        };
    }

    /**
     * @brief Builds the detector configuration for SBND.
     * @return The SBND detector configuration.
     */
    inline DetectorConfig build_sbnd_config()
    {
        const std::string detvar_file =
            build_file_path(k_sbnd_path, "spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root");
        const std::string nosyst_file =
            build_file_path(k_sbnd_path, "spineosc_numu_sbnd_v1.1.0_nosyst.root");
        // The v1.1.0 SBND CV is split into sbnd_small/sbnd_large trees;
        // sbnd_large alone is used as the go-forward CV/nosyst reference.
        const std::string nosyst_key = "sbnd_large";

        // The v1.0.10 WireMod fallback is retired: v1.1.0 has both samples.

        // Every SBND systematic but DENT and the legacy WireMod pair compares
        // its detvar_file variation against the same nosyst_file/nosyst_key
        // central value.
        auto detvar_syst = [&](const std::string &name, const std::string &detvar_key) {
            return SystematicConfig{name, detvar_file, detvar_key, nosyst_file, nosyst_key};
        };

        return {
            "sbnd_run1",
            "sbnd",
            "Run1",
            {
                detvar_syst("SCE0x", "sbnd_detvar_sce0x"),
                detvar_syst("SCE2x", "sbnd_detvar_sce2x"),
                detvar_syst("ChargeScale15", "sbnd_detvar_smear15"),
                // v1.1.0 now reproduces both WireMod samples, so these no longer fall
                // back to v1.0.10. WireModYZ is still COMPUTED (keeping the detvar
                // files complete and allowing a comparison against MattModYZ) but the
                // XMLs allowlist MattModYZ in its place; allowlisting both would
                // double-count one physical effect.
                detvar_syst("WireModxThetaXW", "sbnd_detvar_wiremod_xthetaxw"),
                detvar_syst("MattModYZ", "sbnd_detvar_mattmod_yz"),
                {"DENT", detvar_file, "sbnd_detvar_newdent", detvar_file, "sbnd_detvar_newcv"},
            },
            { },
            {
                // Trigger emulation: PE-to-threshold calibration, its one-sigma
                // uncertainty, and the PE-equivalent threshold. Detector/run
                // specific and supplied directly -- not derivable from the ntuples.
                {"TriggerEmulation", nosyst_file, nosyst_key, "event_largest_flash_pe", 0.647, 0.004, 2000},
            },
        };
    }

    /**
     * @brief Builds the detector configuration for one ICARUS run.
     * @details ICARUS Run 2 and Run 4 read from separate, run-specific
     * detvar/nosyst files and tree keys, so both configs are generated from
     * this single template parameterized on run_suffix.
     * @param run_name Used for output filenames, e.g. "icarus_run2".
     * @param run_label Used in the in-file histogram name, e.g. "Run2".
     * @param run_suffix Used to build file names and tree keys, e.g. "run2".
     * @param trigger_scale Nominal trigger PE-to-threshold calibration factor.
     * @param trigger_sigma One-sigma uncertainty on trigger_scale.
     * @param trigger_threshold PE-equivalent trigger threshold.
     * @return The requested ICARUS run's detector configuration.
     */
    inline DetectorConfig build_icarus_run_config(const std::string &run_name, const std::string &run_label,
                                                   const std::string &run_suffix, double trigger_scale,
                                                   double trigger_sigma, double trigger_threshold)
    {
        const std::string detvar_file =
            build_file_path(k_icarus_path, "spineosc_numu_icarus_" + run_suffix + "_v1.1.0_detvar_nosyst.root");
        const std::string nosyst_file =
            build_file_path(k_icarus_path, "spineosc_numu_icarus_" + run_suffix + "_v1.1.0_nosyst.root");
        const std::string nosyst_key = "icarus_" + run_suffix;

        // Every ICARUS systematic compares its detvar_file variation against
        // the same nosyst_file/nosyst_key central value.
        auto detvar_syst = [&](const std::string &name, const std::string &detvar_key_prefix) {
            return SystematicConfig{name, detvar_file, detvar_key_prefix + "_" + run_suffix, nosyst_file, nosyst_key};
        };

        return {
            run_name,
            "icarus",
            run_label,
            {
                detvar_syst("SCE0x", "icarus_detvar_sce0x"),
                detvar_syst("SCE2x", "icarus_detvar_sce2x"),
                detvar_syst("WireModxThetaXW", "icarus_wiremod_xthetaxw"),
            },
            { },
            {
                // Trigger emulation: PE-to-threshold calibration, its one-sigma
                // uncertainty, and the PE-equivalent threshold. Detector/run
                // specific and supplied directly -- not derivable from the ntuples.
                {"TriggerEmulation", nosyst_file, nosyst_key, "event_largest_flash_pe",
                 trigger_scale, trigger_sigma, trigger_threshold},
            },
        };
    }

    /**
     * @brief Builds the detector configuration for ICARUS Run 2.
     * @return The ICARUS Run 2 detector configuration.
     */
    inline DetectorConfig build_icarus_run2_config()
    {
        return build_icarus_run_config("icarus_run2", "Run2", "run2", 0.615, 0.018, 6000);
    }

    /**
     * @brief Builds the detector configuration for ICARUS Run 4.
     * @return The ICARUS Run 4 detector configuration.
     */
    inline DetectorConfig build_icarus_run4_config()
    {
        return build_icarus_run_config("icarus_run4", "Run4", "run4", 0.361, 0.013, 3000);
    }

    /**
     * @brief Builds the generic (detector-independent) configuration.
     * @details Holds systematics that are hand-set flat ratios rather than
     * derived from any detector's ntuples (e.g. POT normalization
     * uncertainties), so a single configuration covers every detector/run.
     * @return The generic detector configuration.
     */
    inline DetectorConfig build_generic_config()
    {
        return {
            "generic",
            "generic",
            "RunX",
            { },
            {
                {"ToroidCalibration", 1.005},
                {"ProtonMissLow", 0.9975},
                {"ProtonMissHigh", 1.000}
            },
            { },
        };
    }
}  // namespace detsys

#endif  // DETSYS_CONFIGURATION_H
