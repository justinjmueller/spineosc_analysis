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
            {
                "by_vis_energy",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_visible_energy",
                {0.3, 0.4, 0.46, 0.52, 0.58, 0.64, 0.70, 0.76, 0.82, 0.88, 0.94, 1.0, 1.06, 1.13, 1.2, 1.3, 1.5, 2.0}
            },
            {
                "by_muon_score",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_leading_muon_muon_softmax",
                {0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00}
            },
            // The remaining score variables share the muon-score axis: a uniform
            // 0.00-1.00 covers a softmax's full domain, so no event can fall
            // outside the ratio histogram and be silently forced to weight 1.
            // This is deliberately decoupled from the XML channel binning, which
            // may start well above zero without affecting the lookup.
            {
                "by_proton_score",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_leading_proton_proton_softmax",
                {0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00}
            },
            {
                "by_muon_primary_score",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_leading_muon_primary_softmax",
                {0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00}
            },
            {
                "by_proton_primary_score",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_leading_proton_primary_softmax",
                {0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00}
            },
            // The range/momentum agreement variables are signed residuals, not
            // softmaxes, so the "0-1 covers the whole domain" argument above does
            // NOT transfer: each needs an axis wide enough to contain its own
            // tails, or events falling outside are silently forced to weight 1.
            // The calorimetric-vs-CSDA residuals are well behaved and sit inside
            // [-1, 1] for >99.99% of events in every detector/run.
            {
                "by_muon_calo_csda",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_leading_muon_calo_csda_diff",
                {-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}
            },
            {
                "by_proton_calo_csda",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_leading_proton_calo_csda_diff",
                {-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}
            },
            // MCS-vs-CSDA needs a wider, asymmetric axis than the calorimetric
            // pair: its positive tail is long (SBND p99 ~ 2.8, max ~ 700), and a
            // [-1, 1] axis would hold only 97.2% of SBND events. [-1, 3] holds
            // 99.47% (SBND) / 99.87% (ICARUS Run 2) / 99.86% (ICARUS Run 4).
            {
                "by_muon_mcs_csda",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_leading_muon_mcs_csda_diff",
                {-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0}
            },
            // Transverse-kinematic-imbalance and single-particle kinematics.
            // These are kept in the ntuples' NATIVE units (MeV, MeV/c, rad)
            // rather than the GeV convention of xml/varset1_devsample.xml:
            // the branch below is handed straight to RDataFrame's
            // Histo1D<double> as a column name, so an expression such as
            // "reco_pn/1000" is not available here. The XML channel binning
            // is kept in the same native units so that the fit variable and
            // this ratio axis cannot disagree.
            {
                "by_proton_ke",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_leading_proton_ke",
                {0.0, 60.0, 120.0, 180.0, 240.0, 300.0, 360.0, 420.0, 480.0, 540.0, 600.0, 660.0, 720.0, 780.0, 840.0, 900.0, 960.0, 1020.0, 1080.0, 1140.0, 1200.0}
            },
            // Leading-proton momentum, in the same native MeV/c as the KE axis
            // above. The analysis binning (xml/sbnd_data_access/proton_p.xml)
            // is 300-1200 MeV/c; this ratio axis runs to 1800 so the sparse
            // upper tail still lands on a defined bin rather than in overflow.
            {
                "by_proton_p",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_leading_proton_p",
                {300.0, 375.0, 450.0, 525.0, 600.0, 675.0, 750.0, 825.0, 900.0, 975.0, 1050.0, 1125.0, 1200.0, 1275.0, 1350.0, 1425.0, 1500.0, 1575.0, 1650.0, 1725.0, 1800.0}
            },
            // ---- phase1-tier2: single-particle muon kinematics -------------
            // Ratio axes are deliberately wider than the analysis binning so no
            // event falls outside and is silently given weight 1. Verified on
            // the nominal sample: no empty bins on any of these.
            {
                "by_muon_polar_angle",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_leading_muon_polar_angle",
                {0.0, 0.16, 0.32, 0.48, 0.64, 0.80, 0.96, 1.12, 1.28, 1.44, 1.60, 1.76, 1.92, 2.08, 2.24, 2.40, 2.56, 2.72, 2.88, 3.04, 3.20}
            },
            {
                "by_muon_p",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_leading_muon_p",
                {225.0, 300.0, 375.0, 450.0, 525.0, 600.0, 675.0, 750.0, 825.0, 900.0, 975.0, 1050.0, 1125.0, 1200.0, 1275.0, 1350.0, 1425.0, 1500.0, 1575.0, 1650.0, 1725.0}
            },
            {
                "by_opening_angle",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_bivar_muon_proton_opening_angle",
                {0.0, 0.16, 0.32, 0.48, 0.64, 0.80, 0.96, 1.12, 1.28, 1.44, 1.60, 1.76, 1.92, 2.08, 2.24, 2.40, 2.56, 2.72, 2.88, 3.04, 3.20}
            },
            // ---- phase1-tier3: energy estimators ----------------------------
            // The visible-energy subchannel already exists above. energy_qel has
            // pathological tails (min -10142, max +2433 GeV), so its axis starts
            // at 0.15 (a 0.0 start leaves the first bin empty); ~1.5% of events
            // fall outside it, which is the price of a usable axis.
            {
                "by_energy_qel",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_energy_qel",
                {0.15, 0.30, 0.45, 0.60, 0.75, 0.90, 1.05, 1.20, 1.35, 1.50, 1.65, 1.80, 1.95, 2.10, 2.25, 2.40, 2.55, 2.70, 2.85, 3.00}
            },
            // reco_pn is the full momentum imbalance (delta p); it agrees
            // bin-for-bin with sqrt(dpT^2 + dpL^2) computed from the ntuples.
            {
                "by_dp",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_pn",
                {0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0}
            },
            {
                "by_dpt",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_dpT",
                {0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0}
            },
            // The angular variables span [0, pi] (or [-pi, pi]) exactly. The
            // axis is carried out to 3.2 instead of pi so the bin edges stay
            // round while still containing pi, and a handful of events with a
            // non-finite dalphaT/dphiT (< 0.1%) simply fail to fill.
            {
                "by_dalphat",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_dalphaT",
                {0.0, 0.16, 0.32, 0.48, 0.64, 0.80, 0.96, 1.12, 1.28, 1.44, 1.60, 1.76, 1.92, 2.08, 2.24, 2.40, 2.56, 2.72, 2.88, 3.04, 3.20}
            },
            {
                "by_dphit",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_dphiT",
                {0.0, 0.16, 0.32, 0.48, 0.64, 0.80, 0.96, 1.12, 1.28, 1.44, 1.60, 1.76, 1.92, 2.08, 2.24, 2.40, 2.56, 2.72, 2.88, 3.04, 3.20}
            },
            {
                "by_muon_azimuth",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_leading_muon_azimuthal_angle",
                {-3.20, -2.88, -2.56, -2.24, -1.92, -1.60, -1.28, -0.96, -0.64, -0.32, 0.0, 0.32, 0.64, 0.96, 1.28, 1.60, 1.92, 2.24, 2.56, 2.88, 3.20}
            },
            // ----------------------------------------------------------------
            // Position variables: interaction vertex and muon track endpoint.
            //
            // These are the first DETECTOR-RESTRICTED subchannels. A coordinate
            // only means something in the geometry it was measured in, so each
            // is computed for one detector's configs alone via
            // SubchannelConfig::detectors -- an SBND-ranged coordinate run
            // against ICARUS data produces a histogram nothing can use. ICARUS
            // Run 2 and Run 4 share a geometry (their distributions agree to a
            // few tenths of a cm), so they share one subchannel each.
            //
            // Ranges verified against the v1.0.10 ntuples with the
            // reco_fiducialize_cathode == 1 selection applied; every range
            // below holds 100% of selected events in its detector.
            //
            // Note ICARUS z spans BOTH cryostats, -885 to +845, straddling the
            // beam origin -- a [0, 850] axis would silently discard the entire
            // upstream half (only 44.6% retained). Likewise ICARUS y runs
            // -171.5 to +120.4, so [-150, 100] would clip both tails (93.1%).
            // ----------------------------------------------------------------

            // --- SBND: x,y span the full -200..200 cathode-to-cathode extent;
            //     z runs 0..450 (vertex) and 0..500 (muon endpoint, which
            //     exits downstream of its own vertex and reaches 490).
            {
                "by_vertexx_sbnd",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_vertex_x",
                {-200.0, -180.0, -160.0, -140.0, -120.0, -100.0, -80.0, -60.0, -40.0, -20.0, 0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0},
                {"sbnd_run1"}
            },
            {
                "by_vertexy_sbnd",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_vertex_y",
                {-200.0, -180.0, -160.0, -140.0, -120.0, -100.0, -80.0, -60.0, -40.0, -20.0, 0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0},
                {"sbnd_run1"}
            },
            {
                "by_vertexz_sbnd",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_vertex_z",
                {0.0, 25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0},
                {"sbnd_run1"}
            },
            {
                "by_muon_endx_sbnd",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_leading_muon_end_x",
                {-200.0, -180.0, -160.0, -140.0, -120.0, -100.0, -80.0, -60.0, -40.0, -20.0, 0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0},
                {"sbnd_run1"}
            },
            {
                "by_muon_endy_sbnd",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_leading_muon_end_y",
                {-200.0, -180.0, -160.0, -140.0, -120.0, -100.0, -80.0, -60.0, -40.0, -20.0, 0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0},
                {"sbnd_run1"}
            },
            {
                "by_muon_endz_sbnd",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_leading_muon_end_z",
                {0.0, 25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 500.0},
                {"sbnd_run1"}
            },

            // --- ICARUS (Run 2 and Run 4 share a geometry, so one subchannel
            //     each): |x| 70..350, y -185..125, z -900..900 across both
            //     cryostats.
            //
            //     x is FOLDED: the two cryostats sit either side of x = 0 with
            //     no TPC in |x| < 70, so a signed axis has empty central bins
            //     that PROfit's empty-bin check rejects. The branch is an
            //     expression, JIT-compiled by compute_ratio_hist; the XMLs use
            //     the TTreeFormula spelling abs(...) for the same quantity.
            {
                "by_vertexx_icarus",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "std::abs(reco_vertex_x)",
                {70.0, 105.0, 140.0, 175.0, 210.0, 245.0, 280.0, 315.0, 350.0},
                {"icarus_run2", "icarus_run4"}
            },
            {
                "by_vertexy_icarus",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_vertex_y",
                {-185.0, -169.5, -154.0, -138.5, -123.0, -107.5, -92.0, -76.5, -61.0, -45.5, -30.0, -14.5, 1.0, 16.5, 32.0, 47.5, 63.0, 78.5, 94.0, 109.5, 125.0},
                {"icarus_run2", "icarus_run4"}
            },
            {
                "by_vertexz_icarus",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_vertex_z",
                {-900.0, -810.0, -720.0, -630.0, -540.0, -450.0, -360.0, -270.0, -180.0, -90.0, 0.0, 90.0, 180.0, 270.0, 360.0, 450.0, 540.0, 630.0, 720.0, 810.0, 900.0},
                {"icarus_run2", "icarus_run4"}
            },
            {
                "by_muon_endx_icarus",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "std::abs(reco_leading_muon_end_x)",
                {70.0, 105.0, 140.0, 175.0, 210.0, 245.0, 280.0, 315.0, 350.0},
                {"icarus_run2", "icarus_run4"}
            },
            {
                "by_muon_endy_icarus",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_leading_muon_end_y",
                {-185.0, -169.5, -154.0, -138.5, -123.0, -107.5, -92.0, -76.5, -61.0, -45.5, -30.0, -14.5, 1.0, 16.5, 32.0, 47.5, 63.0, 78.5, 94.0, 109.5, 125.0},
                {"icarus_run2", "icarus_run4"}
            },
            {
                "by_muon_endz_icarus",
                "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1",
                "reco_leading_muon_end_z",
                {-900.0, -810.0, -720.0, -630.0, -540.0, -450.0, -360.0, -270.0, -180.0, -90.0, 0.0, 90.0, 180.0, 270.0, 360.0, 450.0, 540.0, 630.0, 720.0, 810.0, 900.0},
                {"icarus_run2", "icarus_run4"}
            }
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
#if 0  // TOGGLE: already computed, not re-run
                detvar_syst("SCE0x", "sbnd_detvar_sce0x"),
#endif  // TOGGLE
#if 0  // TOGGLE: already computed, not re-run
                detvar_syst("SCE2x", "sbnd_detvar_sce2x"),
#endif  // TOGGLE
#if 0  // TOGGLE: already computed, not re-run
                detvar_syst("ChargeScale", "sbnd_detvar_charge_scale"),
#endif  // TOGGLE
#if 0  // TOGGLE: already computed, not re-run
                detvar_syst("ChargeScale15", "sbnd_detvar_smear15"),
#endif  // TOGGLE
                // v1.1.0 now reproduces both WireMod samples, so these no longer fall
                // back to v1.0.10. WireModYZ is still COMPUTED (keeping the detvar
                // files complete and allowing a comparison against MattModYZ) but the
                // XMLs allowlist MattModYZ in its place; allowlisting both would
                // double-count one physical effect.
                detvar_syst("WireModxThetaXW", "sbnd_detvar_wiremod_xthetaxw"),
                detvar_syst("WireModYZ", "sbnd_detvar_wiremod_yz"),
#if 0  // TOGGLE: already computed, not re-run
                detvar_syst("MattModYZ", "sbnd_detvar_mattmod_yz"),
#endif  // TOGGLE
#if 0  // TOGGLE: already computed, not re-run
                detvar_syst("GainHigh", "sbnd_detvar_p5_gain"),
#endif  // TOGGLE
#if 0  // TOGGLE: already computed, not re-run
                detvar_syst("GainLow", "sbnd_detvar_m5_gain"),
#endif  // TOGGLE
#if 0  // TOGGLE: already computed, not re-run
                {"DENT", detvar_file, "sbnd_detvar_newdent", detvar_file, "sbnd_detvar_newcv"},
#endif  // TOGGLE
            },
            { },
            {
                // Trigger emulation: PE-to-threshold calibration, its one-sigma
                // uncertainty, and the PE-equivalent threshold. Detector/run
                // specific and supplied directly -- not derivable from the ntuples.
#if 0  // TOGGLE: already computed, not re-run
                {"TriggerEmulation", nosyst_file, nosyst_key, "event_largest_flash_pe", 0.647, 0.004, 2000},
#endif  // TOGGLE
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
#if 0  // TOGGLE: already computed, not re-run
                detvar_syst("SCE0x", "icarus_detvar_sce0x"),
#endif  // TOGGLE
#if 0  // TOGGLE: already computed, not re-run
                detvar_syst("SCE2x", "icarus_detvar_sce2x"),
#endif  // TOGGLE
#if 0  // TOGGLE: already computed, not re-run
                detvar_syst("ChargeScale", "icarus_detvar_charge_scale"),
#endif  // TOGGLE
                detvar_syst("WireModxThetaXW", "icarus_wiremod_xthetaxw"),
                detvar_syst("WireModYZ", "icarus_wiremod_yz"),
            },
            { },
            {
                // Trigger emulation: PE-to-threshold calibration, its one-sigma
                // uncertainty, and the PE-equivalent threshold. Detector/run
                // specific and supplied directly -- not derivable from the ntuples.
#if 0  // TOGGLE: already computed, not re-run
                {"TriggerEmulation", nosyst_file, nosyst_key, "event_largest_flash_pe",
                 trigger_scale, trigger_sigma, trigger_threshold},
#endif  // TOGGLE
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
