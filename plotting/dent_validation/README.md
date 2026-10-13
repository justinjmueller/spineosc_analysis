# DENT data/MC validation (SBND)

Compares the SBND detector-access data (`spineosc_numu_sbnd_release-v1.1.0_detector_access_DONOTLOOK.root`,
which has no energy or momentum column) with the detvar CV and DENT MC (`sbnd_detvar_newcv`, `sbnd_detvar_newdent`),
after one common selection: >=1 reco proton, `reco_fiducialize_cathode == 1`, `reco_veto_sbnd_highy_highz == 1`,
`0.647 * event_largest_flash_pe > 2000`. Shapes are compared at unit area; the cathode-fraction plots use fractions of
all selected events. `event_largest_flash_pe` is used only in the selection and is never plotted.

Regenerate (the CSV dumps are large and live outside the repo):

    root -l -b -q 'dump_samples.C("<outdir>")'      # writes data.csv, cv.csv, dent.csv
    # point SCRATCH in dentval.py at <outdir>, then
    python3 fig2_cathode.py   # cathode-region features + cathode_feature_numbers.csv
    python3 fig134.py         # ranked scan (needs scan_all_variables.csv), calorimetry/PID, other variables

Outputs: `fig1_scan_summary.png`, `fig2_cathode_features.png`, `fig3_calorimetry_pid.png`, `fig4_other_variables.png`,
`scan_all_variables.csv`, `cathode_feature_numbers.csv`.

Caveats: the detvar samples are neutrino-only (the selection is >99% pure in data); the detvar CV is not the production
nominal MC, so some chi2 values reflect data/MC differences unrelated to DENT. Chi2 values are shape-only.
