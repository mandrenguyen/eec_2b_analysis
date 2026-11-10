# Local jet calibration inputs

`Prompt24HIpp_Placeholder1Percent_JESUncertainty_AK4PF.txt` defines the symmetric 1% JES uncertainty agreed upon for the 2024 pp-reference AK4PF analysis. It is an analysis placeholder, rather than an official JEC uncertainty payload.

The file uses the text format read by `Run3HeavyIonJetCalibrations/include/JetUncertainty.h`. It covers `-5.191 <= eta <= 5.191` and is constant from 0 to 10000 GeV; the reader clamps values outside the listed pT range to its endpoint uncertainty.

`Prompt24HIpp_Derived_MC_PtResolution_AK4PF.txt` is the nominal fractional MC
jet-pT resolution used by the hybrid JER smearing. It was derived with
`../derive_jer_resolution.cpp` from the first two million events of QCD block
0000 after applying the analysis L1FastJet and L2Relative corrections. Reco/gen
pairs require `refdrjt < 0.2`. In each `(pT_gen, |eta|)` bin, the response core
is fit iteratively with a Gaussian over `mean +/- 1.5 sigma`; the resulting
widths are fit with the standard JME four-parameter resolution form. A separate
two-million-event derivation from block 0001 was used as a stability check.

`Prompt24HIpp_Estimated_MC_SF_AK4PF.txt` contains the analysis JER scale-factor
estimate. In the table format expected by `JetSmearer.h`, the three values are
nominal, down, and up: `1.10`, `1.05`, and `1.15`.
