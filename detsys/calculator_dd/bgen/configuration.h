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

#include <cstdlib>
#include <string>
#include <vector>

#include "classes.h"

namespace detsys
{
    /**
     * @brief Directory holding the v1.1.0 release ntuples.
     * @details Read at run time from the SPINEOSC_DATA environment variable
     * (the directory that contains the v1.1.0 folder); defaults to
     * $HOME/data/spineosc. No rebuild is needed to point at another machine's data.
     * @return "<SPINEOSC_DATA>/v1.1.0".
     */
    inline std::string release_dir()
    {
        const char *root = std::getenv("SPINEOSC_DATA");
        if (root != nullptr && *root != '\0')
        {
            return std::string(root) + "/v1.1.0";
        }
        const char *home = std::getenv("HOME");
        return std::string(home != nullptr ? home : ".") + "/data/spineosc/v1.1.0";
    }

    inline const std::string k_sbnd_path = release_dir();
    inline const std::string k_icarus_path = release_dir();

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
            {"dd_muon_len3d", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "sqrt(pow(reco_leading_muon_end_x-reco_leading_muon_start_x,2)+pow(reco_leading_muon_end_y-reco_leading_muon_start_y,2)+pow(reco_leading_muon_end_z-reco_leading_muon_start_z,2))", {49.9, 57.1, 65, 72.8, 80.5, 88.3, 96.3, 105, 113, 122, 132, 143, 154, 167, 181, 198, 217, 240, 271, 319, 430}, {"sbnd_run1"}},
            {"dd_muon_abs_endx", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "std::abs(reco_leading_muon_end_x)", {0.75, 7.35, 14.9, 23, 30.8, 38.9, 47, 55, 63.2, 71.9, 80.5, 89.2, 98.6, 108, 118, 128, 139, 150, 163, 176, 188}, {"sbnd_run1"}},
            {"dd_crosser_flag", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_cathode_crosser", {-0.5, 0.5, 1.5}, {"sbnd_run1"}},
            {"dd_muon_endx_full", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_leading_muon_end_x", {-187, -167, -147, -128, -111, -95, -79.7, -64.7, -50.5, -36.8, -23, -9.45, 2.55, 15.1, 27.5, 41, 53.9, 67.4, 81.5, 95.9, 112, 128, 146, 167, 187}, {"sbnd_run1"}},
            {"dd_vertexx_full", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_vertex_x", {-182, -159, -141, -124, -109, -94.5, -81.2, -68.2, -56, -44, -31.9, -20.2, -8.7, 12.6, 24, 35.6, 47.7, 60.1, 73, 86.7, 101, 117, 135, 156, 182}, {"sbnd_run1"}},
            {"dd_proton_endx_full", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_leading_proton_end_x", {-186, -166, -147, -129, -113, -98, -83.9, -70.7, -57.8, -45.1, -32.5, -20, -5.55, 11.5, 24.1, 36.5, 49, 62, 75.5, 89.9, 105, 121, 140, 162, 186}, {"sbnd_run1"}},
            {"dd_muon_startx_full", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_leading_muon_start_x", {-182, -159, -141, -124, -108, -94.4, -81.2, -68.2, -56, -44, -32, -20.2, -8.55, 12.8, 24.1, 35.5, 47.5, 60.1, 73, 86.6, 101, 117, 135, 156, 181}, {"sbnd_run1"}},
            {"dd_proton_calo_csda", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_leading_proton_calo_csda_diff", {-0.765, -0.644, -0.598, -0.568, -0.543, -0.522, -0.501, -0.48, -0.461, -0.443, -0.425, -0.408, -0.391, -0.374, -0.357, -0.34, -0.323, -0.304, -0.282, -0.255, -0.201}, {"sbnd_run1"}},
            {"dd_muon_proton_score", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_leading_muon_proton_softmax", {0, 1.85e-05, 0.000129, 0.000372, 0.000809, 0.00155, 0.00281, 0.00494, 0.0086, 0.0146, 0.0242, 0.04, 0.0671, 0.118, 0.219, 0.438, 1}, {"sbnd_run1"}},
            {"dd_muon_pion_score", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_leading_muon_pion_softmax", {0, 0.00377, 0.0132, 0.0208, 0.0284, 0.0365, 0.0456, 0.0555, 0.0667, 0.0801, 0.0969, 0.118, 0.146, 0.187, 0.255, 0.421, 1}, {"sbnd_run1"}},
            {"dd_proton_mip_score", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_leading_proton_mip_softmax", {0, 7.86e-05, 0.000219, 0.000396, 0.000642, 0.000998, 0.00152, 0.00232, 0.00359, 0.00571, 0.0095, 0.0169, 0.0333, 0.0718, 0.172, 0.409, 1}, {"sbnd_run1"}},
            {"dd_proton_pion_score", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_leading_proton_pion_softmax", {0, 6.4e-05, 0.000177, 0.000321, 0.000515, 0.000788, 0.00117, 0.00174, 0.0026, 0.00396, 0.00623, 0.0104, 0.0191, 0.0398, 0.104, 0.346, 1}, {"sbnd_run1"}},
            {"dd_proton_muon_score", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_leading_proton_muon_softmax", {0, 3.23e-06, 2.12e-05, 4.39e-05, 7.53e-05, 0.000121, 0.000192, 0.000305, 0.000496, 0.000825, 0.00145, 0.00271, 0.00557, 0.0129, 0.0367, 0.159, 1}, {"sbnd_run1"}},
            {"dd_muon_endy", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_leading_muon_end_y", {-187, -165, -142, -122, -104, -86.9, -70.4, -54.5, -38.9, -23.5, -8.25, 6.75, 22, 37.6, 53.5, 70.4, 87.8, 107, 127, 150, 186}, {"sbnd_run1"}},
            {"dd_vertexy", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_vertex_y", {-182, -154, -133, -114, -96.5, -80.2, -64.8, -50.2, -35.9, -21.8, -7.72, 6.19, 20.4, 35, 50, 65.7, 82.3, 100, 120, 144, 180}, {"sbnd_run1"}},
            {"dd_proton_endy", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_leading_proton_end_y", {-186, -161, -138, -118, -99.8, -83, -67, -51.8, -36.8, -22, -7.65, 6.75, 21.5, 36.5, 52, 68.2, 85.6, 104, 125, 150, 184}, {"sbnd_run1"}},
            {"dd_proton_endz", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_leading_proton_end_z", {16.9, 40.5, 60.4, 79, 96.7, 114, 131, 149, 167, 185, 203, 222, 242, 263, 285, 308, 333, 360, 390, 425, 474}, {"sbnd_run1"}},
            {"dd_muon_dir_x", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_leading_muon_start_dir_x", {-0.986, -0.88, -0.769, -0.663, -0.56, -0.462, -0.364, -0.267, -0.167, -0.0652, 0.0363, 0.135, 0.231, 0.324, 0.414, 0.503, 0.593, 0.686, 0.783, 0.887, 0.987}, {"sbnd_run1"}},
            {"dd_muon_dir_y", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_leading_muon_start_dir_y", {-0.987, -0.882, -0.771, -0.665, -0.566, -0.472, -0.38, -0.288, -0.196, -0.104, -0.00862, 0.0883, 0.184, 0.279, 0.373, 0.466, 0.563, 0.662, 0.769, 0.881, 0.987}, {"sbnd_run1"}},
            {"dd_muon_dir_z", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_leading_muon_start_dir_z", {-0.956, -0.644, -0.4, -0.215, -0.0662, 0.061, 0.168, 0.26, 0.343, 0.418, 0.488, 0.552, 0.611, 0.666, 0.717, 0.765, 0.81, 0.853, 0.896, 0.94, 0.991}, {"sbnd_run1"}},
            {"dd_pion_mult", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_pion_multiplicity", {-0.5, 0.5, 1.5, 2.5, 3.5}, {"sbnd_run1"}},
            {"dd_photon_mult", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_photon_multiplicity", {-0.5, 0.5, 1.5, 2.5, 3.5}, {"sbnd_run1"}},
            {"dd_electron_mult", "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000", "reco_electron_multiplicity", {-0.5, 0.5, 1.5, 2.5}, {"sbnd_run1"}}
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
