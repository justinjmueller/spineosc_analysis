# Which cuts cause selection churn in a detector variation (ChargeScale)

Script: `plotting/studies/detvar_churn.py`. Simulation only. ChargeScale (`sbnd_detvar_charge_scale`) against the production nominal
(`sbnd_large`), matched event by event, analysis selection S2 (crosser veto + 4 cm end-point cut), all multiplicities.

## Method and its limit
Events selected in one sample but not the other ("lost" and "gained") are attributed to the reconstruction criteria that fail in the sample
where they are not selected. The per-criterion flags exist only in the `signal` tree, i.e. for **true-signal events**; for those the flags
reproduce the `selected` tree for 99.9% of events. Selected events that are not true signal have no per-criterion flags, so they are broken
down by true category instead. That says what kind of event leaves the selection, not which cut removed it.

## Results
All energies (0.3-3 GeV): 22,916 events lost, 19,331 gained; about half of each are not true signal and cannot be attributed.
Attributed signal events that are lost fail (an event can fail several): at least one proton 5,916 (3,797 as the only failing criterion),
no charged pions 4,918 (2,021), single muon 2,174 (263), end point >= 4 cm 1,438 (17), high-y/z veto 478, containment 437.

The fraction of selected nominal events that leave the selection depends strongly on true category:

| True category | Meaning | Lost fraction |
|---|---|---|
| 0, 2 | signal (1p, Np) | 4.6%, 5.4% |
| 4 | contained, 1 muon, any protons, >0 charged pions | 23.7% |
| 5 | same, non-contained | 24.8% |
| 6, 7 | CC, non-signal (fiducial / non-fiducial) | 20.1%, 13.6% |
| 8 | NC | 12.5% |
| 9 | cosmic | 23.6% |

**The 1.5-2.0 GeV bin** (the one behind the large multi-proton ChargeScale uncertainty): 88 events lost, 27 gained, net -61 of the -65
change in the bin; energy migration across the bin edges cancels (32 out, 28 in). Of the 88 lost events, 20 are true signal (13 fail the
no-charged-pion cut) and 68 are not; of those, 47 are category 4 (a true charged pion that the nominal reconstruction missed). At >= 1.5 GeV
27% of the selected category-4 events leave the selection, against 0.6% of the category-2 signal.

## Interpretation (inference, not directly observed)
The variation moves events across the pion and proton criteria. The high-energy loss is dominated by background events with true charged
pions that pass the reco pion veto in the nominal sample and fail it in the variation. The per-cut flags for these events are not in the ntuples.
Confirming it needs the per-criterion flags for non-signal events, or the pion's reconstructed kinetic energy in both samples.
