#ifndef RESULT_PATHS_H
#define RESULT_PATHS_H



// apply_unfolding_2d.C WRITES these paths and apply_weights_and_systematics.C READS them.
//
// Both macros include this header, so there is ONE definition of every mapping and the
// two cannot drift. 
// If you add a variation, add it in ONE place here and both macros follow.

#include "TString.h"
#include <cstdlib>
#include <iostream>

// A tag no path will ever contain, returned for a variation name that is not recognised.
// The point is that it must NOT be "": an empty tag is the nominal, so an unknown value
// silently becoming nominal is the bug this header exists to prevent.
inline const TString &kUnknownTag()
{
    static const TString t = "_UNKNOWN_VARIATION";
    return t;
}

inline void warnUnknown(const char *what, const TString &value, const char *allowed)
{
    std::cerr << "ERROR: unknown " << what << " '" << value << "' (use " << allowed << "). "
              << "Paths built from it carry " << kUnknownTag() << " and will not open -- "
              << "this is deliberate, so the value cannot alias onto the nominal result."
              << std::endl;
}

// ---- Known values -------------------------------------------------------------------
// The canonical lists. Validate against these before doing any work.
inline bool isKnownGenerator(const TString &g)
{
    return g == "pythia" || g == "herwig";
}

inline bool isKnownSample(const TString &s)
{
    return s == "qcd" || s == "bjet" || s == "both";
}

// Which of template_fit.cpp's fits supplies the signal fraction. These are the varNames[]
// strings of Help_Functions.h, which are also the prefixes on its per-variation output
// files, so a fit file can be found from the variation name alone.
inline bool isKnownTfVariation(const TString &v)
{
    return v == "nominal" || v == "var0B_2" || v == "var0B_0";
}

// Which UParT b-tag efficiency scale factor is divided into the reco-level data.
// "off" applies no SF at all -- a comparison, not a systematic.
// Two DIFFERENT ways to use the JP HF calibration, deliberately both available:
//   jpcalib_hf            applies the alternative SF curve itself, signs and all. The shift
//                         flips sign across dr (up in bins 1-3, down from 4 on), so it
//                         RESHAPES the distribution and survives the unit-area normalisation.
//   jpsyst_up / jpsyst_dn shift the CENTRAL SF by +/- |jphf - central| coherently in every
//                         dr bin. Same magnitude per bin, sign discarded -- the usual
//                         "propagate the quoted uncertainty" treatment. A coherent shift is
//                         much closer to a normalisation, so most of it cancels in the
//                         unit-area renormalisation.
// They are the same numbers used two ways, so they must never both enter the band.
inline bool isKnownSfupartVariation(const TString &v)
{
    return v == "nominal" || v == "off"
        || v == "jpcalib_hf" || v == "qqrate_up" || v == "qqrate_down"
        || v == "statup"     || v == "statdn"
        || v == "jpsyst_up"  || v == "jpsyst_dn";
}

// ---- Path tags ----------------------------------------------------------------------
// "" means nominal in every one of these, so that nominal paths carry no tag at all and
// are unchanged from before any variation existed.

// TF_GENERATOR: which generator's template fit supplied the signal fraction.
inline TString tfTag(const TString &tfGenerator)
{
    if (!isKnownGenerator(tfGenerator)) {
        warnUnknown("TF_GENERATOR", tfGenerator, "pythia | herwig");
        return kUnknownTag();
    }
    return (tfGenerator != "pythia") ? ("_tf" + tfGenerator) : TString("");
}

// TRACK_EFF_UNC: the 3% track-drop MC production. The literal must stay in step with
// TrkEffSyst::kDropFraction in tracking_efficiency_syst.h.
inline TString trkTag(bool trackEffUnc)
{
    return trackEffUnc ? TString("_trkdrop030") : TString("");
}

// EEC_WEIGHT_OFF: the yield run. NOT a systematic variation -- the same chain measuring a
// different observable, with the (pt1*pt2)^n weight dropped everywhere, so the unfolded
// result is dN/dr instead of the EEC. Its own tag so a yield file can never be mistaken for,
// or overwrite, an EEC one. Must stay in step with EecWeight::tag() in
// create_files_for_template_fit.cpp, which is what the WRITER uses.
inline TString eecWeightTag(bool eecWeightOff)
{
    return eecWeightOff ? TString("_noeecw") : TString("");
}

// SV_SF (added 2026-10-12): the SV b scale-factor MC production (SV_SF=nominal in
// run_agg_ntuple_chunks.sh, see sv_sf.h). Read from the SAME environment variable the run
// script uses, so a chain run as `SV_SF=<variation> <driver>` reads that variation's MC and
// writes its own fits and results. Must stay in step with SvSF::tag(), which is what the
// WRITER uses. Data is never weighted, so dataTemplateFile() ignores it.
//
// svsfNominal(): what an UNSET SV_SF means. ON since 2026-10-12 (Zoe: the result WITH the SV
// b SF is the nominal, its uncertainty a source in the band) -- the default chain reads the
// "_svsf" MC and results. SV_SF=off reads the no-SF chain, whose results are untouched; set
// this false to make that the default again. Same pattern as jecInputs() below.
inline bool svsfNominal() { return true; }

inline TString svsfTag()
{
    const char *env = std::getenv("SV_SF");
    const TString v = (env && *env) ? TString(env) : TString(svsfNominal() ? "nominal" : "off");
    if (v == "off")             return "";
    if (v == "nominal")         return "_svsf";
    if (v == "up")              return "_svsfup";
    if (v == "down")            return "_svsfdn";
    // The per-contribution variations (2026-10-13), tag "_svsf<name>" -- SvSF::modeNames().
    for (const char *c : {"statup", "statdn", "nomuonup", "nomuondn", "bbup", "bbdn",
                          "jpup", "jpdn", "qqup", "qqdn"})
        if (v == c) return "_svsf" + v;
    warnUnknown("SV_SF", v, "off | nominal | up | down | stat|nomuon|bb|jp|qq + up|dn");
    return kUnknownTag();
}

// TRACK_BDT_SF: effective MC working points that reproduce the measured track-level SFb
// and SFl. The two sources are independent, while their in/out-SV values move together.
// Data never carries this tag. An unset variable means the corrected nominal production.
inline TString trackBdtSfTag()
{
    const char *env = std::getenv("TRACK_BDT_SF");
    const TString v = (env && *env) ? TString(env) : TString("nominal");
    if (v == "off")     return "";
    if (v == "nominal") return "_trkbdtsf";
    if (v == "boff")    return "_trkbdtsf_boff";
    if (v == "loff")    return "_trkbdtsf_loff";
    warnUnknown("TRACK_BDT_SF", v, "off | nominal | boff | loff");
    return kUnknownTag();
}

// TF_VARIATION: the template-fit variation supplying the signal fraction.
inline TString tfVarTag(const TString &tfVariation)
{
    if (!isKnownTfVariation(tfVariation)) {
        warnUnknown("TF_VARIATION", tfVariation, "nominal | var0B_2 | var0B_0");
        return kUnknownTag();
    }
    if (tfVariation == "var0B_2") return "_tf0Bx2";
    if (tfVariation == "var0B_0") return "_tf0Bx0";
    return "";
}

// SFUPART_VARIATION: the UParT SF variation. Namespaced to UParT because more scale
// factors are coming and each should get its own mapping rather than share a generic one.
inline TString sfupartVarTag(const TString &sfupartVariation)
{
    if (!isKnownSfupartVariation(sfupartVariation)) {
        warnUnknown("SFUPART_VARIATION", sfupartVariation,
                    "nominal | off | jpcalib_hf | qqrate_up | qqrate_down | statup | statdn "
                    "| jpsyst_up | jpsyst_dn");
        return kUnknownTag();
    }
    if (sfupartVariation == "off")         return "_sfupartoff";
    if (sfupartVariation == "jpcalib_hf")  return "_sfupartjphf";
    if (sfupartVariation == "qqrate_up")   return "_sfupartqqup";
    if (sfupartVariation == "qqrate_down") return "_sfupartqqdn";
    if (sfupartVariation == "statup")      return "_sfupartstatup";
    if (sfupartVariation == "statdn")      return "_sfupartstatdn";
    if (sfupartVariation == "jpsyst_up")   return "_sfupartjpsystup";
    if (sfupartVariation == "jpsyst_dn")   return "_sfupartjpsystdn";
    return "";
}

// ---- The UParT scale factor itself ---------------------------------------------------
inline TString sfupartFile()
{
    // "finalbins": source bins 8-9 merged into one output bin, so this file has 8 dr bins
    // against the analysis's 9. Its last bin spans 0.35-0.45 and therefore applies to BOTH
    // of the analysis's last two bins -- the application maps by bin centre so that happens
    // on its own. The merge also tamed the last bins: 1.149 +/- 0.020 and 1.086 +/- 0.029
    // became a single 1.010 +/- 0.017.
    // The 9-bin predecessor (..._eec_systematics.root, no "finalbins") is superseded.
    return "/home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/workflow/"
           "lifetime_jp_sfb_agg_dr_wp712_eec_finalbins_systematics.root";
}

// NOTE: these are the KEY names in that file. The histograms' own internal names differ
// (key "h_SFb_dr_central" holds a histogram named "h_SFb_dr_wp712_eec_central"), and
// TFile::Get() takes the key -- the internal name returns null.
// "off" returns the marker "(none)": no SF is applied, the caller checks for it.
inline TString sfupartHist(const TString &sfupartVariation)
{
    if (!isKnownSfupartVariation(sfupartVariation)) {
        warnUnknown("SFUPART_VARIATION", sfupartVariation,
                    "nominal | off | jpcalib_hf | qqrate_up | qqrate_down | statup | statdn "
                    "| jpsyst_up | jpsyst_dn");
        return "";
    }
    if (sfupartVariation == "off")         return "(none)";
    // statup/statdn shift the CENTRAL SF by +/- its own bin error, so they read the same
    // histogram; the shift is applied where the SF is used, not here. Coherent across all
    // dr bins, because the analysis's last two bins share ONE calibration bin and their
    // errors are therefore 100% correlated -- a coherent shift gets that right, and is the
    // conservative choice elsewhere.
    // jpsyst_up/jpsyst_dn likewise start from the CENTRAL SF and are shifted by the JP HF
    // systematic where the SF is used -- see sfupartSystHist() just below.
    if (sfupartVariation == "nominal"
     || sfupartVariation == "statup"
     || sfupartVariation == "statdn"
     || sfupartVariation == "jpsyst_up"
     || sfupartVariation == "jpsyst_dn")   return "h_SFb_dr_central";
    if (sfupartVariation == "jpcalib_hf")  return "h_SFb_dr_jpcalib_hf";
    if (sfupartVariation == "qqrate_up")   return "h_SFb_dr_qqrate_up";
    if (sfupartVariation == "qqrate_down") return "h_SFb_dr_qqrate_down";
    return "";
}

// The per-bin UNCERTAINTY histogram added to / subtracted from the central SF, for the
// variations that propagate a quoted error rather than swapping in another curve. "" means
// no such shift, which is every other variation.
//
// h_SFb_dr_syst_jpcalib_hf is |h_SFb_dr_jpcalib_hf - h_SFb_dr_central| bin by bin (checked
// 2026-09-14, all 8 bins). That is why this is NOT a duplicate of jpcalib_hf: the same
// magnitudes are applied with one coherent sign instead of the measured, sign-flipping one.
inline TString sfupartSystHist(const TString &sfupartVariation)
{
    if (sfupartVariation == "jpsyst_up" || sfupartVariation == "jpsyst_dn")
        return "h_SFb_dr_syst_jpcalib_hf";
    return "";
}

// +1 / -1 / 0: which way sfupartSystHist() is applied.
inline double sfupartSystShift(const TString &sfupartVariation)
{
    if (sfupartVariation == "jpsyst_up") return  1.;
    if (sfupartVariation == "jpsyst_dn") return -1.;
    return 0.;
}

// ---- Result paths --------------------------------------------------------------------
// The combined tag every path carries. One function, so a folder and the label inside it
// can never disagree about which variation they are.
// Which OBSERVABLE the result measures. "" for dr, so every existing dr path is unchanged;
// "_B" for the momentum balance. It is a different measurement, not a variation of the dr
// one, so it gets its own results folder -- a B result sitting in the dr folder would be
// picked up by apply_weights_and_systematics.C by name and silently compared against dr.
//
// Same rule as the rest of this header: ONE definition, both the writer
// (apply_unfolding_2d.C) and the reader (apply_weights_and_systematics.C) go through it.
//
// ⚠️ RENAMED 2026-09-22: the momentum balance is "B", tagged "_B". It was "z"/"_z", and
// nothing translates between the two -- results and MC written under the old spelling are
// not found by this code, and are not overwritten by it either. They sit alongside the new
// ones in the results area under their "_z" names.
//
// ⚠️ fnb, fb and lnfb are DISABLED (2026-09-22) together with their ObsDefs in
// observables.h and their binnings in binning_histos_small.h. Their branches are left
// commented rather than deleted, so re-enabling one is an uncomment here too.
inline bool isKnownObservable(const TString &observable)
{
    return observable == "dr" || observable == "B";
    //  || observable == "fnb" || observable == "fb" || observable == "lnfb";
}

inline TString observableTag(const TString &observable)
{
    if (observable == "dr")  return "";
    if (observable == "B")   return "_B";
    // DISABLED 2026-09-22:
    //   if (observable == "fnb")  return "_fnb";
    //   if (observable == "fb")   return "_fb";
    //   if (observable == "lnfb") return "_lnfb";
    return kUnknownTag();   // never "" -- an unknown value must not alias onto dr
}

// ---- Which MC PRODUCTION carries an observable ----------------------------------------
// This is the run scripts' OUT_TAG, i.e. an INPUT path tag, not a results-path tag -- it is
// here for the same reason trkTag() is: the literal is shared by the production filenames
// and by three macros that have to agree on it (template_fit.cpp, apply_unfolding_2d.C and
// the data-file choice in both). An observable only exists in blocks produced AFTER the code
// that fills it, so a run pointed at the wrong production fails on a missing histogram --
// the right way round, but it has to be pointed correctly from one place.
//
//   dr   "_upartv2"          the original production, dR only
//   B    "_upartv2_B"        NOT PRODUCED YET. The momentum balance, under its new name.
//
// ⚠️ WHY B NEEDS A NEW PRODUCTION, even though the observable itself did not change.
// The balance histograms used to be named "..._z" and are now "..._B" (ObsDef::suffix in
// observables.h). The earlier productions -- "_upartv2_zfirst", "_upartv2_3obs",
// "_upartv2_lnfblin", "_upartv2_fb" -- all carry the "_z" spelling, so this code opens them
// and finds nothing on the balance axis. Re-run step 1 with OUT_TAG=_upartv2_B
// (run_agg_ntuple_chunks.sh) before unfolding B.
//
// The old productions are NOT deleted and NOT renamed: the "_z" results already on disk stay
// reproducible from the files they were made from, by the version of the code that made them.
//
//   fnb  "_upartv2_3obs"     ) DISABLED 2026-09-22 with the observables themselves. All three
//   lnfb "_upartv2_lnfblin"  ) predate the rename, so re-enabling any of them means a fresh
//   fb   "_upartv2_fb"       ) production too, not just an uncomment.
// ---- Which INPUT ERA the whole chain runs on (added 2026-09-28) ----------------------
// Marc re-made every MC and data input chunk on 2026-09-22..24 as merged_block_*_rho.root;
// the inputs the "_upartv2" productions were made from are gone. Zoe chose (2026-09-28) to
// re-produce BOTH observables on the _rho inputs, so dR and B come from the same MC.
//
// The rho era gets its own tag at EVERY level -- MC/data production names, template-fit
// directories and result folders -- because the old productions cannot be overwritten in
// place: _rho has 8 Pythia blocks per sample where the old had 10 (QCD) / 9 (bJet), so
// writing the new blocks over the old names would leave old blocks 8-9 behind, and the
// block-summing in apply_unfolding_2d.C would silently add them in.
//
// Flip rhoInputs() to false to point the chain back at the pre-_rho results on disk.
inline bool rhoInputs() { return true; }

// ---- Jet energy corrections (added 2026-10-01) -----------------------------------------
// true = read the productions made with APPLY_JEC=true (2024 pp-reference JEC on jtpt, see
// the Jec namespace in create_files_for_template_fit.cpp). They are tagged "jec" at every
// level, so the JEC and no-JEC chains sit side by side. Only meaningful with rhoInputs():
// the JEC needs the event rho the _rho forests carry.
// ON since 2026-10-01 (the upartv2rhojec productions exist). Set false to read the
// uncorrected _upartv2rho results, which are untouched.
inline bool jecInputs() { return true; }

// "_upartv2rho" in the rho era, "_upartv2" before it, "_upartv2rhojec" with the JEC. Ends
// every production name, every template-fit directory and every result folder.
inline TString eraTag()
{
    if (!rhoInputs()) return TString("_upartv2");
    return jecInputs() ? TString("_upartv2rhojec") : TString("_upartv2rho");
}

inline TString obsProdTag(const TString &observable)
{
    if (observable == "dr")  return eraTag();
    if (observable == "B")   return eraTag() + "_B";
    // DISABLED 2026-09-22:
    //   if (observable == "fnb")  return "_upartv2_3obs";
    //   if (observable == "lnfb") return "_upartv2_lnfblin";
    //   if (observable == "fb")   return "_upartv2_fb";
    warnUnknown("OBSERVABLE", observable, "dr | B");
    return kUnknownTag();   // never "_upartv2": must not alias onto the dr production
}

// ---- The DATA file an observable is measured from -------------------------------------
// The reco-level data templates, as hadd'ed to the top of agg_template_chunks. Three macros
// need it -- template_fit.cpp fits it, apply_unfolding_2d.C unfolds it, plot_raw_yields.C
// draws it -- and they must all open the SAME file, or the raw yields shown do not belong
// to the fit that was run on them.
//
// Two things pick the file, and both are properties of the measurement, not of the run:
//   the EEC weight   h3D_data is filled with eec * weight_tree, so EEC-weighted data fitted
//                    against unweighted templates would mix two observables. The yield run
//                    takes its own data production.
//   the observable   only a production that filled this observable's axis has its
//                    histograms at all -- obsProdTag() above.
inline TString dataTemplateFile(bool eecWeightOff, const TString &observable = "dr")
{
    const TString hp = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/HardProbes/"
                       "agg_template_chunks/";
    // The rho era: every observable, EEC-weighted or not, from Zoe's own data production on
    // the _rho inputs, hadd'ed to the top of agg_template_chunks. No fallback onto Afnan's.
    if (rhoInputs())
        return hp + "Run3_btagWP0712_template_for_fit_histos_3D_data_fMCGEN"
                  + (eecWeightOff ? TString("_noeecw") : TString(""))
                  + obsProdTag(observable) + ".root";
    if (observable != "dr")
        return hp + "Run3_btagWP0712_template_for_fit_histos_3D_data_fMCGEN"
                  + (eecWeightOff ? TString("_noeecw") : TString(""))
                  + obsProdTag(observable) + ".root";
    if (eecWeightOff)
        return hp + "Run3_btagWP0712_template_for_fit_histos_3D_data_fMCGEN_noeecw_upartv2.root";
    // The dR EEC measurement still reads Afnan's data production, which predates all of this.
    return "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/HardProbes/agg_template_chunks/"
           "Run3_btagWP712_template_for_fit_histos_3D_data_f_80_9999_2MCGEN.root";
}

// ---- What the NOMINAL run of an observable is ----------------------------------------
// Two things differ between dr and B for reasons that are properties of the OBSERVABLE, not
// choices of a particular run. Both were being applied by hand in one macro each, which is
// the drift this header exists to stop: the writer forced them and the reader had to
// reproduce the same forcing from memory, or it would look for a path that is never written.
// One definition, both macros go through it.

// The UParT b-tag efficiency SF is measured in reco dR bins and has no B equivalent, so the
// nominal B run applies NO SF and says so in its path ("_sfupartoff"). That is a KNOWN GAP,
// not a silent skip: the dr result carries a correction (0.933 -> 1.149 across its range)
// that the B result does not.
//
// Written as "everything that is not dr", so the disabled fraction observables are still
// covered by it if one of them comes back: the SF exists in reco dR and nowhere else.
inline TString nominalSfupart(const TString &observable)
{
    return (observable == "dr") ? TString("nominal") : TString("off");
}

// The momentum balance is a YIELD measurement: B is measured with the (pt1*pt2)^n EEC weight
// OFF, so the unfolded result is dN/dB. Unlike the SF above this is a convention rather than
// an impossibility -- the EEC-weighted B production exists -- but it is THE convention, so
// both macros default to it and a B run that quietly used the weighted templates would be
// compared against an unweighted band without anything saying so.
// Flip this one function if the weighted B ever becomes the measurement.
inline bool nominalEecWeightOff(const TString &observable)
{
    return observable == "B";
    // DISABLED 2026-09-22. The fraction observables followed B here: they are
    // momentum-sharing fractions, not the angular structure of the pair, so weighting them
    // by (pt1*pt2)^n would measure something nobody asked for.
    //  || observable == "fnb" || observable == "fb" || observable == "lnfb";
}

inline TString variationTag(const TString &tfGenerator, bool trackEffUnc,
                            const TString &tfVariation, const TString &sfupartVariation,
                            bool eecWeightOff = false,
                            const TString &observable = "dr")
{
    return observableTag(observable) + tfTag(tfGenerator) + trkTag(trackEffUnc) + svsfTag()
         + trackBdtSfTag()
         + tfVarTag(tfVariation) + sfupartVarTag(sfupartVariation)
         + eecWeightTag(eecWeightOff);
}

// ---- Which unfolding produced it ------------------------------------------------------
// Matrix inversion results live in their own subdirectory of the results area, so a whole
// MI production -- every variation, every observable, plots and logs included -- is one
// directory to sync, compare or throw away.
//
// This is a DIRECTORY, unlike every other tag in this header, because it is not a variation
// of the measurement: the same variation unfolded two ways is the same systematic, and
// putting the two side by side under one parent would interleave two complete productions.
//
// Bayesian keeps the flat layout it always had -- "" here, the same rule as every other
// nominal in this file -- so every Bayesian path that exists on disk today is unchanged,
// including the ones already used in talks.
inline TString methodDir(bool unfoldBayes)
{
    return unfoldBayes ? TString("") : TString("matrix_inversion/");
}

// unfoldBayes sits third, before the defaulted arguments and in the same position it holds
// in resultLabel/resultFile, so it CANNOT be forgotten: a call written for the old signature
// passes a TString where a bool is expected and fails to compile, rather than quietly
// building a Bayesian path for a matrix-inversion run.
inline TString resultFolder(const TString &sample, const TString &generator, bool unfoldBayes,
                            const TString &tfGenerator = "pythia",
                            bool trackEffUnc = false,
                            const TString &tfVariation = "nominal",
                            const TString &sfupartVariation = "nominal",
                            bool eecWeightOff = false,
                            const TString &observable = "dr")
{
    return "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/"
         + methodDir(unfoldBayes)
         + "unfolding_"
         + sample + "_" + generator
         + variationTag(tfGenerator, trackEffUnc, tfVariation, sfupartVariation,
                        eecWeightOff, observable)
         + eraTag() + "/";
}

// <sample>_<generator><variation tags>_<bayesian|MI>_<mode>.
// This is the "dataset" string apply_unfolding_2d.C names its output histograms file with,
// so the reader reconstructs exactly what the writer wrote.
inline TString resultLabel(const TString &sample, const TString &generator,
                           int test_mode, bool unfoldBayes,
                           const TString &tfGenerator = "pythia",
                           bool trackEffUnc = false,
                           const TString &tfVariation = "nominal",
                           const TString &sfupartVariation = "nominal",
                           bool eecWeightOff = false,
                           const TString &observable = "dr")
{
    TString label = sample + "_" + generator
                  + variationTag(tfGenerator, trackEffUnc, tfVariation, sfupartVariation,
                                 eecWeightOff, observable);
    label += unfoldBayes ? "_bayesian" : "_MI";
    if (test_mode == 0) label += "_full_closure";
    if (test_mode == 1) label += "_split_test";
    if (test_mode == 2) label += "_data";
    return label;
}

inline TString resultFile(const TString &sample, const TString &generator,
                          int test_mode, bool unfoldBayes,
                          const TString &tfGenerator = "pythia",
                          bool trackEffUnc = false,
                          const TString &tfVariation = "nominal",
                          const TString &sfupartVariation = "nominal",
                          bool eecWeightOff = false,
                          const TString &observable = "dr")
{
    return resultFolder(sample, generator, unfoldBayes, tfGenerator, trackEffUnc, tfVariation,
                        sfupartVariation, eecWeightOff, observable)
         + "histos_" + resultLabel(sample, generator, test_mode, unfoldBayes, tfGenerator,
                                   trackEffUnc, tfVariation, sfupartVariation, eecWeightOff,
                                   observable)
         + "_after_unfolding_2D.root";
}

#endif // RESULT_PATHS_H
