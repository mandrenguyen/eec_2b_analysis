// Post-processing of the unfolded result.
//
// Takes the output of apply_unfolding_2d.C, applies the final weights, builds a
// systematic band from the variations (Pythia vs Herwig to start with) and writes the
// result carrying stat and syst errors separately plus their quadrature sum.
//
// Prerequisite: run apply_unfolding_2d.C once per (SAMPLE, GENERATOR) you want in here.
// A variation whose file is missing is reported and skipped, so the nominal alone still
// works -- you just get stat errors until the Herwig run exists.

#include "result_paths.h"   // the ONE definition of every result path and variation tag
#include "observables.h"    // ObsDef: the axis title and binning of the observable measured
#include "TBox.h"
#include "TSystem.h"         // Setenv/Unsetenv for the per-variation SV_SF
#include <vector>
#include <algorithm>
#include <cmath>
#include <map>

// ---- Paths ---------------------------------------------------------------------------
// resultFolder/resultLabel/resultFile and every variation tag now come from
// result_paths.h, which apply_unfolding_2d.C includes too. They used to be duplicated here
// and drifted -- see the note at the top of that header for what that cost.

// A colour blended towards white, i.e. what "alpha over a white canvas" looks like.
// SetFillColorAlpha() is honoured in the PDF but ignored by ROOT's PNG backend, which
// paints the fill solid; blending gives the same light tone in both.
Color_t blendWithWhite(const char *hex, double alpha)
{
    TColor *c = gROOT->GetColor(TColor::GetColor(hex));
    if (!c) return kWhite;
    Float_t r, g, b;
    c->GetRGB(r, g, b);
    return (Color_t) TColor::GetColor(Float_t(1. - alpha + alpha * r),
                                      Float_t(1. - alpha + alpha * g),
                                      Float_t(1. - alpha + alpha * b));
}

// Generator name as it should appear on a plot.
TString prettyGen(const TString &g)
{
    if (g == "herwig") return "Herwig";
    if (g == "pythia") return "Pythia";
    return g;
}

// A correction-level variation: no second unfolding, the alternative generator's correction
// is applied to the final result and the result re-normalised. inverse = the correction is
// DIVIDED into the data (an efficiency); false = it MULTIPLIES the data (the EEC weight).
struct CorrVariation {
    TString name;     // printout column and histogram name (no spaces)
    TString label;    // legend text
    TString effType;  // plotNice_UParT_efficiency effType, for the "run this first" message
    TString file;     // that macro's output .root
    bool    inverse;
    TString sysLabel; // name in the "systematics applied" tag on the data-vs-gen plot
};

// One systematic variation: the same unfolding redone with something changed.
struct Variation {
    TString name;        // printout column and output histogram name (no spaces)
    TString sample;
    TString generator;   // UNFOLDING_GENERATOR of that run
    TString tfGenerator; // TF_GENERATOR of that run ("pythia" unless the fit is varied)
    TString label;       // legend text
    bool    trackEffUnc; // TRACK_EFF_UNC of that run (false unless the MC itself is varied)
    TString tfVariation = "nominal"; // TF_VARIATION of that run (the template-fit variation)
    TString sfupartVariation = "nominal"; // SFUPART_VARIATION of that run (the UParT SF variation)
    // Entries sharing a non-empty group are the two SIDES OF ONE uncertainty (an up/down
    // pair). Their envelope -- the larger of the two |shifts| in each bin -- enters the
    // quadrature total once. An empty group means an independent source that enters on its
    // own, which is how every entry behaved before groups existed.
    TString group = "";
    // Name in the "systematics applied" tag on the data-vs-gen plot. The two sides of a
    // pair carry the SAME sysLabel and are listed once. Leave it empty to keep an entry
    // out of the tag while it still counts in the band.
    TString sysLabel = "";
    // false = DISPLAY ONLY. The curve is loaded, drawn and printed, but contributes nothing
    // to sigma+/sigma- and nothing to the tag. For showing the effect of a correction --
    // "here is the result without the UParT SF" -- which is a comparison, not a systematic.
    bool    inBand = true;
    // true = |delta| goes into BOTH sigma+ and sigma-. false (the default) = SIGNED: the
    // shift feeds only the side it actually pushes the result towards. Only the tracking
    // efficiency is symmetrised (Zoe, 2026-09-29); every other source is signed.
    bool    symmetrise = false;
    // SV_SF of that run (svsfTag() in result_paths.h, added 2026-10-12). "" = the same as the
    // nominal, i.e. whatever SV_SF this band is built with; "up" / "down" for the two sides of
    // the SV b scale-factor uncertainty, whose results live under _svsfup / _svsfdn.
    TString svsfVariation = "";
    // TRACK_BDT_SF of that run. The nominal applies both SFb and SFl; "boff" and "loff"
    // remove one source at a time while retaining the other.
    TString trackBdtSfVariation = "";
};

// ---- Weights = scale factors ---------------------------------------------------------
// Applied to the nominal AND to every variation, before any difference is taken: a shape
// systematic has to be compared after the same weighting, otherwise a pure normalisation
// difference between the samples leaks into the shape band.
//
// This is the GENERIC per-bin scale-factor hook, deliberately not named for any one SF --
// the UParT SF has its own sfupart* plumbing in apply_unfolding_2d.C (it is applied at reco
// level, before unfolding), and this stays free for a future SF that genuinely belongs on
// the final result.
//
// sf_file/sf_hist is the per-bin scale factor histogram, multiplied bin by bin (it must
// have the same binning as the result). global_scale is a single number for a flat SF.
struct Weights {
    TString sf_file             = "";    // per-bin scale factors, "" = off
    TString sf_hist             = "";
    bool    sf_errors           = true;  // false = treat the SF as exact (ignore its errors)
    // An EFFICIENCY scale factor is DIVIDED into the result, it does not multiply it. The
    // chain corrected the data with the efficiency measured in MC, and the true efficiency
    // is eps_data = SF * eps_MC, so
    //     corrected = raw / eps_MC   ->   true = raw / (SF * eps_MC) = corrected / SF.
    // Same convention as CorrVariation::inverse elsewhere in this file. Set false for a
    // factor that genuinely multiplies (a weight, a luminosity ratio).
    bool    sf_divide           = true;
    double  global_scale        = 1.0;   // flat scale factor / luminosity / prescale; 1 = off
    bool    normalise_unit_area = true;  // Scale(1/integral, "width") -- per bin width, unit area
};

bool applyWeights(TH1D *h, const Weights &w)
{
    if (!h) return false;

    if (w.sf_file.Length() > 0) {
        TFile *fw = TFile::Open(w.sf_file);
        if (!fw || fw->IsZombie()) {
            std::cerr << "ERROR: cannot open scale-factor file " << w.sf_file << std::endl;
            return false;
        }
        TH1D *hsf = dynamic_cast<TH1D *>(fw->Get(w.sf_hist));
        if (!hsf) {
            std::cerr << "ERROR: '" << w.sf_hist << "' not in " << w.sf_file << std::endl;
            fw->Close();
            return false;
        }
        if (hsf->GetNbinsX() != h->GetNbinsX()) {
            std::cerr << "ERROR: scale factors have " << hsf->GetNbinsX() << " bins, result has "
                      << h->GetNbinsX() << std::endl;
            fw->Close();
            return false;
        }
        // Applied MANUALLY, bin by bin, with GetBinContent/SetBinContent -- deliberately
        // NOT TH1::Divide/Multiply. Those require the two histograms to share a binning and
        // return silently when they do not, which is exactly how a scale factor binned
        // differently from the result would fail: quietly, leaving the result unchanged.
        // A loop cannot fail that way. (The UParT SF is applied the same way, at reco level
        // in apply_unfolding_2d.C, where the SF has 8 dr bins against the analysis's 9 --
        // a binning Divide() would refuse outright.)
        //
        // The SF's statistical error is folded into the result's bin errors, i.e. into its
        // STATISTICAL error; sf_errors = false treats the SF as exact.
        for (int i = 1; i <= h->GetNbinsX(); ++i) {
            const double sf = hsf->GetBinContent(i);
            if (sf == 0.) {
                std::cerr << "ERROR: scale factor is zero in bin " << i << " of " << w.sf_hist
                          << " -- refusing to " << (w.sf_divide ? "divide" : "multiply")
                          << std::endl;
                fw->Close();
                return false;
            }
            const double factor = w.sf_divide ? (1. / sf) : sf;
            const double rel_sf = w.sf_errors ? (hsf->GetBinError(i) / sf) : 0.;

            const double c = h->GetBinContent(i);
            const double e = h->GetBinError(i);
            const double c_new = c * factor;
            const double e_new = (c != 0.)
                ? std::fabs(c_new) * std::sqrt((e / c) * (e / c) + rel_sf * rel_sf)
                : e * factor;
            h->SetBinContent(i, c_new);
            h->SetBinError(i, e_new);
        }
        fw->Close();
    }

    if (w.global_scale != 1.0) h->Scale(w.global_scale);

    if (w.normalise_unit_area) {
        const double integral = h->Integral();
        if (integral <= 0.) {
            std::cerr << "WARNING: " << h->GetName() << " has integral " << integral
                      << ", not normalising" << std::endl;
            return false;
        }
        h->Scale(1. / integral, "width");
    }
    return true;
}

// Read one result histogram, detached from its file, with the weights already applied.
TH1D *loadWeighted(const TString &filename, const TString &histname,
                   const TString &newname, const Weights &w)
{
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cerr << "MISSING: " << filename << std::endl;
        return nullptr;
    }
    TH1D *h_in = dynamic_cast<TH1D *>(f->Get(histname));
    if (!h_in) {
        std::cerr << "ERROR: '" << histname << "' not found in " << filename << std::endl;
        f->Close();
        return nullptr;
    }
    TH1D *h = (TH1D *) h_in->Clone(newname);
    h->SetDirectory(nullptr);
    f->Close();
    if (!applyWeights(h, w)) return nullptr;
    return h;
}

// =====================================================================================
// GENERATOR / TF_GENERATOR are the UNFOLDING_GENERATOR and TF_GENERATOR of the NOMINAL
// unfolding run; the variations below say which of them to flip.
//
// OBSERVABLE: dr | B -- which measurement to build the band for. It is threaded into every
//             result path (result_paths.h), so a B run reads B results and only B results;
//             it also picks the variation list, because the two observables do not have the
//             same systematics available. unfoldBayes defaults to false (matrix inversion)
//             to match apply_unfolding_2d.C -- the two must agree or every file is "missing".
//
//   root -l -b -q 'apply_weights_and_systematics.C("both","pythia")'                      # EEC(dr)
//   root -l -b -q 'apply_weights_and_systematics.C("both","pythia",2,false,"h_data_fully_corrected_2D","pythia",2,false,"B")'   # dN/dB
void apply_weights_and_systematics(TString SAMPLE = "both", TString GENERATOR = "pythia",
                                   int test_mode = 2, bool unfoldBayes = false,
                                   TString histname = "h_data_fully_corrected_2D",
                                   TString TF_GENERATOR = "pythia",
                                   int ibin_pt = 2,
                                   bool EEC_WEIGHT_OFF = false,
                                   TString OBSERVABLE = "dr")
{
    if (SAMPLE != "qcd" && SAMPLE != "bjet" && SAMPLE != "both") {
        std::cerr << "ERROR: unknown SAMPLE '" << SAMPLE << "' (use qcd | bjet | both)" << std::endl;
        return;
    }
    if (GENERATOR != "pythia" && GENERATOR != "herwig") {
        std::cerr << "ERROR: unknown GENERATOR '" << GENERATOR << "' (use pythia | herwig)" << std::endl;
        return;
    }
    if (!isKnownObservable(OBSERVABLE)) {
        std::cerr << "ERROR: unknown OBSERVABLE '" << OBSERVABLE << "' (use dr | B)" << std::endl;
        return;
    }
    // The two per-observable conventions, read from result_paths.h rather than repeated here.
    // apply_unfolding_2d.C forces the same two when it WRITES, so these are what is actually
    // on disk; reproducing them by hand in this file is exactly the drift that header exists
    // to stop.
    //   B carries no UParT SF (it is measured in reco dR bins and has no B equivalent), so
    //   the nominal B result and every B variation live under "_sfupartoff".
    //   B is a yield measurement: no EEC weight, so the result is dN/dB.
    const TString SFUPART_NOM = nominalSfupart(OBSERVABLE);
    if (!EEC_WEIGHT_OFF && nominalEecWeightOff(OBSERVABLE)) {
        std::cout << "NOTE: OBSERVABLE '" << OBSERVABLE << "' is measured WITHOUT the EEC "
                  << "weight (dN/d" << OBSERVABLE << "). Reading the EEC_WEIGHT_OFF results -- "
                  << "see nominalEecWeightOff() in result_paths.h." << std::endl;
        EEC_WEIGHT_OFF = true;
    }
    // A generator systematic only means something when both unfoldings ran on the SAME
    // input, which is only true for data. In the closure modes each run unfolds its own MC,
    // so the difference is not a systematic on a measurement.
    if (test_mode != 2)
        std::cerr << "WARNING: test_mode " << test_mode << " unfolds MC, not data -- the "
                  << "variation differences below are NOT a systematic on a measurement."
                  << std::endl;

    // =========================== SET SCALE FACTORS HERE ===============================
    // The UParT b-tag efficiency SF is NOT applied here any more. It is measured in reco
    // R_BB bins, so it is divided into the RECO-level data inside apply_unfolding_2d.C,
    // before unfolding -- see SFUPART_VARIATION there. Applying it here would put a reco-binned
    // correction on particle-level bins and ignore the migration between them.
    // sf_file/sf_hist stay available for a genuinely final-state, per-bin factor.
    Weights w;
    w.sf_file             = "";    // per-bin scale factors: file and histogram name
    w.sf_hist             = "";
    w.sf_errors           = true;  // propagate the SF errors into the stat error
    w.sf_divide           = true;  // an efficiency SF divides into the data (see Weights)
    w.global_scale        = 1.0;   // flat scale factor; 1 = off
    w.normalise_unit_area = true;  // the EEC is a shape: unit area, divided by bin width
    // ==================================================================================

    // ============================= SET VARIATIONS HERE ================================
    // Each entry is re-read, re-weighted and compared bin by bin against the nominal.
    //
    // Nominal is (UNFOLDING_GENERATOR, TF_GENERATOR) = (pythia, pythia). Each variation
    // flips ONE of the two and is subtracted from it:
    //   template_fit    (pythia, herwig) -- the signal fraction from the Herwig fit
    //   unfolding_model (herwig, pythia) -- the whole unfolding MC: migration matrix AND
    //                                       purity, efficiency, truth, hgenjet_2b_*
    //   herwig_full     (herwig, herwig) -- both at once. Since 2026-09-29 this REPLACES the
    //                                       two above in the dr-rho and B lists: one source,
    //                                       so their partial cancellation is kept.
    // They are different sources, so their signed shifts are added in quadrature below,
    // upward and downward shifts separately. (If you ever want them treated as two sides of ONE uncertainty,
    // give them a shared Variation::group -- the total loop then takes their envelope once,
    // which is how the light-jet mistag pair below is combined.)
    const TString other_gen = (GENERATOR == "pythia") ? "herwig" : "pythia";
    // What the observable is called on every axis of every plot this macro draws. With the
    // EEC weight off the same chain measures a yield, so labelling it "EEC" would be wrong,
    // and with OBSERVABLE = B it is not dr either. One definition, used by every axis title
    // below.
    //
    // The x-axis title for dr stays the literal "#Delta r" it always was rather than
    // obs.axis ("#DeltaR"), so every existing dr plot is pixel-for-pixel unchanged; B takes
    // its full label from the ObsDef, which is the one place it is written down.
    const ObsDef  obs       = obsByName(OBSERVABLE);
    const TString obs_sym   = obsSymbol(OBSERVABLE);
    const TString obs_axis  = (OBSERVABLE == "dr") ? "#Delta r" : obs.axis;
    const TString obs_title = EEC_WEIGHT_OFF ? ("dN/d" + obs_sym) : ("EEC(" + obs_sym + ")");

    // An axis with no explicit range is padded by ROOT when it draws it: B lives on
    // [0.5, 1.0] in 5 bins and comes out drawn to 1.1, which puts empty frame where the
    // observable cannot go at all (B = 1 needs one B hadron with zero pT), and makes the last bin
    // look like it is missing. Asking for the bin range explicitly is what stops it --
    // SetNdivisions(..., kFALSE) does NOT, the padding is not the tick optimisation.
    //
    // dr already draws as exactly [0, 0.45], so it is left untouched and every existing dr
    // plot is unchanged.
    auto fixXRange = [&](TH1 *h) {
        if (OBSERVABLE != "dr") h->GetXaxis()->SetRange(1, h->GetNbinsX());
    };
    // Each observable gets the variations that EXIST for it, and only those. Booking an entry
    // whose production was never made just prints a wall of "skipped" lines and invites
    // someone to read a band built from whatever happened to load -- so a variation that is
    // not producible today is commented out here, next to what it would take to produce it.
    std::vector<Variation> variations;

    if (OBSERVABLE != "dr") {
      // ---- B (dN/dB) --------------------------------------------------------------------
      // The full dr set, minus the two sources that cannot exist for B (2026-09-24). Same
      // entries, same groups, same sysLabels as dr, so the two bands are directly comparable
      // source by source. All of them run off the "_upartv2_B" productions
      // (obsProdTag() in result_paths.h) with the EEC weight off:
      //   template_fit     Herwig B template fit           needs the Herwig B production
      //   unfolding_model  Herwig B unfolding MC           the same Herwig B production
      //   tracking_eff     3% track drop, MC AND refit     TRACK_EFF_UNC=true B production
      //   mistag_0B_*      0B x2 / x0 refits               no new MC, one unfolding each
      //
      // Deliberately NOT in this band, each for a different reason:
      //   UParT SF (all of it)  no equivalent exists -- the SF is measured in reco dR bins.
      //                         Every such path carries "_sfupartoff" and says so. KNOWN GAP:
      //                         the dr result carries a correction this result does not.
      //   EEC weight modelling  the weight is off here, so the correction is identically 1.
      variations = {
        // The two Herwig swaps as SEPARATE sources again (Zoe, 2026-10-01: "separate as we
        // were doing before"), replacing herwig_full below.
        { "template_fit",    SAMPLE, GENERATOR, other_gen,    "TF " + prettyGen(other_gen),        false,
          "nominal", SFUPART_NOM, "", "MC template modeling" },
        { "unfolding_model", SAMPLE, other_gen, TF_GENERATOR, "Unfolding " + prettyGen(other_gen), false,
          "nominal", SFUPART_NOM, "", "Detector response" },
        /* ---- disabled (kept for reference): herwig_full, 2026-09-29 -> 2026-10-01 ----
        // All the Herwig dependence in ONE source (Zoe, 2026-09-29): the whole chain redone
        // with Herwig -- template fit AND unfolding MC -- against the all-Pythia nominal.
        // The two separate swaps above partly cancel each other (opposite signs in most
        // bins), which adding them in quadrature ignores. Signed, like everything but tracking.
        { "herwig_full",     SAMPLE, other_gen, other_gen,    "Full " + prettyGen(other_gen) + " chain", false,
          "nominal", SFUPART_NOM, "", "MC generator modeling" },
        ---- end disabled: herwig_full ---- */
        { "tracking_eff",    SAMPLE, GENERATOR, TF_GENERATOR, "Tracking eff. (-3% tracks)", true,
          "nominal", SFUPART_NOM, "", "Tracking efficiency", true, /*symmetrise*/ true },
        { "mistag_0B_up",   SAMPLE, GENERATOR, TF_GENERATOR, "Light-jet mistag (0B #times 2)", false,
          "var0B_2", SFUPART_NOM, "mistag_0B", "Light jet mistagging" },
        { "mistag_0B_down", SAMPLE, GENERATOR, TF_GENERATOR, "Light-jet mistag (0B #times 0)", false,
          "var0B_0", SFUPART_NOM, "mistag_0B", "Light jet mistagging" },
      };
    } else if (rhoInputs() && !EEC_WEIGHT_OFF) {
      // ---- dr (EEC) on the _rho inputs: THE SAME FOUR SOURCES AS B (Zoe, 2026-09-28) -----
      // "the same systematics for both, even if they are not all of them": the dr band books
      // exactly the B list above -- same names, groups, sysLabels -- so the two measurements
      // carry the same uncertainty model and compare source by source.
      //
      // Left OUT of the dr band on purpose, although producible for dr:
      //   UParT SF stat / JP HF swap / qq rate   no B equivalent exists
      //   MC modeling of EEC weight              identically 1 for B (weight off)
      // The UParT SF CORRECTION itself is still applied to the dr nominal (SFUPART_NOM =
      // "nominal"): it is a correction, not a systematic -- only its uncertainty is dropped.
      // The full pre-_rho dr list is the next branch, untouched.
      variations = {
        // The two Herwig swaps as SEPARATE sources again (Zoe, 2026-10-01: "separate as we
        // were doing before"), replacing herwig_full below.
        { "template_fit",    SAMPLE, GENERATOR, other_gen,    "TF " + prettyGen(other_gen),        false,
          "nominal", SFUPART_NOM, "", "MC template modeling" },
        { "unfolding_model", SAMPLE, other_gen, TF_GENERATOR, "Unfolding " + prettyGen(other_gen), false,
          "nominal", SFUPART_NOM, "", "Detector response" },
        /* ---- disabled (kept for reference): herwig_full, 2026-09-29 -> 2026-10-01 ----
        // All the Herwig dependence in one source -- see the B branch above.
        { "herwig_full",     SAMPLE, other_gen, other_gen,    "Full " + prettyGen(other_gen) + " chain", false,
          "nominal", SFUPART_NOM, "", "MC generator modeling" },
        ---- end disabled: herwig_full ---- */
        { "tracking_eff",   SAMPLE, GENERATOR, TF_GENERATOR, "Tracking eff. (-3% tracks)", true,
          "nominal", SFUPART_NOM, "", "Tracking efficiency", true, /*symmetrise*/ true },
        { "mistag_0B_up",   SAMPLE, GENERATOR, TF_GENERATOR, "Light-jet mistag (0B #times 2)", false,
          "var0B_2", SFUPART_NOM, "mistag_0B", "Light jet mistagging" },
        { "mistag_0B_down", SAMPLE, GENERATOR, TF_GENERATOR, "Light-jet mistag (0B #times 0)", false,
          "var0B_0", SFUPART_NOM, "mistag_0B", "Light jet mistagging" },
      };
    } else if (EEC_WEIGHT_OFF) {
      // The dr YIELD run has only ONE variation produced so far: the 3% track drop, re-run
      // with the weight off. Every other entry below needs a production that does not exist
      // unweighted (Herwig, the 0B template refits, the UParT SF variations). So the yield
      // band is exactly the tracking uncertainty, and says so.
      variations = {
        { "tracking_eff", SAMPLE, GENERATOR, TF_GENERATOR, "Tracking eff. (-3% tracks)", true,
          "nominal", "nominal", "", "Tracking efficiency", true, /*symmetrise*/ true },
      };
    } else {
      variations = {
        { "template_fit",    SAMPLE, GENERATOR, other_gen,    "TF " + prettyGen(other_gen),        false,
          "nominal", "nominal", "", "MC template modeling" },
        { "unfolding_model", SAMPLE, other_gen, TF_GENERATOR, "Unfolding " + prettyGen(other_gen), false,
          "nominal", "nominal", "", "Detector response" },
        // Tracking efficiency: the same unfolding redone on MC in which 3% of the
        // reconstructed tracks were thrown away during the substructure extraction. That
        // value covers the residual data/MC difference for tracking in the jet core seen
        // in the jet energy scale determination. It is an independent source from the two
        // above -- a detector effect, not a generator choice -- and the ONLY symmetrised
        // one: its |nominal - variation| enters both sides of the quadrature sum.
        // Produce it with TRACK_EFF_UNC=true in run_agg_ntuple_chunks.sh, then
        //   apply_unfolding_2d.C(SAMPLE, GENERATOR, 2, true, false, TF_GENERATOR, true)
        { "tracking_eff",    SAMPLE, GENERATOR, TF_GENERATOR, "Tracking eff. (-3% tracks)", true,
          "nominal", "nominal", "", "Tracking efficiency", true, /*symmetrise*/ true },
        // Light-jet mistagging: the two sides of ONE uncertainty, so they share a group and
        // their envelope counts once in the total (see Variation::group).
        //
        // The 0B template is the jets that passed the b tag with NO gen b hadron in them --
        // mistagged light and charm. It is not a free component of the fit: it is folded
        // into the background PDF at the MC-predicted ratio c' = 0B/(0B+1B), so the mistag
        // rate reaches the measurement only through the SHAPE of the effective background,
        // and the fitted signal fraction is what moves. The two variations rebuild that
        // shape with the 0B admixture doubled and removed -- a conservative +/-100% on the
        // mistag rate, which brackets the measured light-flavour mistag scale factors
        // (20-50% at a tight WP) by a wide margin and needs no external SF input.
        //
        // These need no new MC and no refit: template_fit.cpp already writes one fit per
        // variation into every template-fit directory. Produce each with
        //   apply_unfolding_2d.C(SAMPLE, GENERATOR, 2, true, false, TF_GENERATOR, false, "var0B_2")
        // and the same with "var0B_0".
        { "mistag_0B_up",   SAMPLE, GENERATOR, TF_GENERATOR, "Light-jet mistag (0B #times 2)", false,
          "var0B_2", "nominal", "mistag_0B", "Light jet mistagging" },
        { "mistag_0B_down", SAMPLE, GENERATOR, TF_GENERATOR, "Light-jet mistag (0B #times 0)", false,
          "var0B_0", "nominal", "mistag_0B", "Light jet mistagging" },
        // UParT SF statistical precision. The SF's own error used to be folded into the
        // data's statistical error; it is now a systematic in its own right, so the data
        // error bar stays a data error bar and this uncertainty is separately attributable.
        // The two sides shift the SF coherently by +/- its bin error, so they are a genuine
        // pair: grouped, envelope counted once, and free to come out asymmetric.
        // Produce with SFUPART_VARIATION = "statup" / "statdn".
        // sysLabel deliberately EMPTY: this pair counts in the band and shows on the
        // per-source breakdown, but is kept off the tag on the data-vs-gen plot.
        { "sfupart_stat_up", SAMPLE, GENERATOR, TF_GENERATOR, "UParT SF (stat +1#sigma)", false,
          "nominal", "statup", "sfupart_stat", "" },
        { "sfupart_stat_dn", SAMPLE, GENERATOR, TF_GENERATOR, "UParT SF (stat #minus1#sigma)", false,
          "nominal", "statdn", "sfupart_stat", "" },

        // DISPLAY ONLY (inBand = false): the same unfolding with no UParT SF applied, so
        // the plot shows the result before and after the correction. It is a correction,
        // not a systematic, so it must not enter the band -- hence the trailing false.
        // Produce it with
        //   apply_unfolding_2d.C(SAMPLE, GENERATOR, 2, true, false, TF_GENERATOR, false,
        //                        "nominal", "off")
        { "sfupart_off", SAMPLE, GENERATOR, TF_GENERATOR, "Before UParT SF", false,
          "nominal", "off", "", "", false },

        // The UParT SF's two CALIBRATION systematics, alongside its statistical precision
        // above. Each is a separate unfolding with a different SF divided into the reco
        // data, so the response matrix carries the varied SF through the migration.
        //   JP HF swap : one alternative calibration -> ungrouped and signed: each bin's
        //                shift feeds only the side it pushes towards.
        //   qq rate    : a genuine +/-25% pair straddling the central -> grouped, so the
        //                total takes its envelope once and it can come out asymmetric.
        // Produce each with
        //   apply_unfolding_2d.C(SAMPLE, GENERATOR, 2, true, false, TF_GENERATOR, false,
        //                        "nominal", "<jpcalib_hf|qqrate_up|qqrate_down>")
        { "sfupart_jphf",   SAMPLE, GENERATOR, TF_GENERATOR, "UParT SF (JP HF swap)", false,
          "nominal", "jpcalib_hf",  "", "UParT SF: JP HF swap" },
        { "sfupart_qq_up",  SAMPLE, GENERATOR, TF_GENERATOR, "UParT SF (qq rate +25%)", false,
          "nominal", "qqrate_up",   "sfupart_qqrate", "UParT SF: qq rate #pm25%" },
        { "sfupart_qq_down",SAMPLE, GENERATOR, TF_GENERATOR, "UParT SF (qq rate #minus25%)", false,
          "nominal", "qqrate_down", "sfupart_qqrate", "UParT SF: qq rate #pm25%" },

        // THE SECOND WAY TO USE THE SAME JP HF MEASUREMENT -- for comparison against
        // "sfupart_jphf" above, which is the first.
        //
        //   sfupart_jphf        : the alternative SF CURVE, applied as measured. Its shift
        //                         flips sign across dr (+ in bins 1-3, - from 4 on), so it
        //                         reshapes the EEC and survives the unit-area normalisation.
        //   sfupart_jpsyst_*    : the quoted JP HF UNCERTAINTY, +/- |jphf - central| applied
        //                         coherently in every dr bin. Identical magnitudes, one sign.
        //                         A coherent shift is nearly a normalisation, so most of it
        //                         cancels when the result is renormalised to unit area.
        //
        // Both read the same two histograms in the calibration file, so booking both in the
        // band would count the JP HF calibration twice. The pair is therefore DISPLAY ONLY
        // (trailing false) while the two treatments are being compared -- flip that to true
        // and drop "sfupart_jphf" if this is the treatment you want in the band.
        // Produce with SFUPART_VARIATION = "jpsyst_up" / "jpsyst_dn".
        { "sfupart_jpsyst_up", SAMPLE, GENERATOR, TF_GENERATOR, "UParT SF (JP HF syst +1#sigma)",
          false, "nominal", "jpsyst_up", "sfupart_jpsyst", "", false },
        { "sfupart_jpsyst_dn", SAMPLE, GENERATOR, TF_GENERATOR, "UParT SF (JP HF syst #minus1#sigma)",
          false, "nominal", "jpsyst_dn", "sfupart_jpsyst", "", false },
        // { "sample", (SAMPLE == "both") ? "qcd" : "both", GENERATOR, TF_GENERATOR, "qcd only", false }, // MC composition
      };
    }

    // ---- SV b scale factor uncertainty, per contribution (Zoe, 2026-10-13) --------------
    // Booked whenever the nominal itself carries the SV b SF (svsfTag() == "_svsf", the
    // default since svsfNominal()), for dr and B alike. Each contribution in the SF file is
    // its own source -- the WHOLE chain redone with only that one shifted, coherent in SV pT
    // (sv_sf.h) -- and the five enter the quadrature total separately, each pair as its
    // larger shift per side (signed), like the mistag pair.
    //
    // JP HF and qq rate are calibrations SHARED with the UParT b-tag SF, so for dr each of
    // those chains moves BOTH SFs the same way: SV SF jp+/- with UParT jpsyst_up/dn (the same
    // +/- |JP HF swap - central| treatment), SV SF qq up/down with UParT qqrate_up/down. B has
    // no UParT SF, so there they move the SV SF alone.
    //
    // The 2026-10-12 version (one pair: central +/- stat (+) syst_total) is kept below, off.
    // Produce with run_svsf_components.sh.
    if (svsfTag() == "_svsf" && (OBSERVABLE != "dr" || (rhoInputs() && !EEC_WEIGHT_OFF))) {
        const bool    dr  = (OBSERVABLE == "dr");
        const TString jpU = dr ? "jpsyst_up"   : SFUPART_NOM, jpD = dr ? "jpsyst_dn"   : SFUPART_NOM;
        const TString qqU = dr ? "qqrate_up"   : SFUPART_NOM, qqD = dr ? "qqrate_down" : SFUPART_NOM;
        const TString jpLab = dr ? "JP HF swap (SV + UParT SF)" : "JP HF swap (SV SF)";
        const TString qqLab = dr ? "qq rate (SV + UParT SF)"    : "qq rate (SV SF)";
        struct P { TString name, label, sysLabel, svUp, svDn, sfU, sfD; };
        const std::vector<P> pairs = {
            { "sv_sf_stat",   "SV b SF stat",            "SV b SF: stat",     "statup",   "statdn",   SFUPART_NOM, SFUPART_NOM },
            { "sv_sf_nomuon", "SV b SF no-muon",         "SV b SF: no muon",  "nomuonup", "nomuondn", SFUPART_NOM, SFUPART_NOM },
            { "sv_sf_bb",     "SV b SF bb model",        "SV b SF: bb model", "bbup",     "bbdn",     SFUPART_NOM, SFUPART_NOM },
            { "jp_hf",        jpLab,                     "JP HF swap",        "jpup",     "jpdn",     jpU, jpD },
            { "qq_rate",      qqLab,                     "qq rate",           "qqup",     "qqdn",     qqU, qqD },
        };
        for (const P &q : pairs) {
            variations.push_back({ q.name + "_up",   SAMPLE, GENERATOR, TF_GENERATOR, q.label + " (+)", false,
                                   "nominal", q.sfU, q.name, q.sysLabel, true, false, q.svUp });
            variations.push_back({ q.name + "_down", SAMPLE, GENERATOR, TF_GENERATOR, q.label + " (#minus)", false,
                                   "nominal", q.sfD, q.name, q.sysLabel, true, false, q.svDn });
        }
        /* ---- disabled 2026-10-13 (kept for reference): the total SV SF as ONE pair ----
        variations.push_back({ "sv_sf_up",   SAMPLE, GENERATOR, TF_GENERATOR, "SV b SF (+1#sigma)", false,
                               "nominal", SFUPART_NOM, "sv_sf", "SV b scale factor", true, false, "up" });
        variations.push_back({ "sv_sf_down", SAMPLE, GENERATOR, TF_GENERATOR, "SV b SF (#minus1#sigma)", false,
                               "nominal", SFUPART_NOM, "sv_sf", "SV b scale factor", true, false, "down" });
        ---- end disabled ---- */
    }

    // ---- Track-BDT SF on/off systematics ---------------------------------------------
    // With no dedicated uncertainties, use the full correction as the uncertainty. SFb
    // and SFl are independent, while the in-SV and out-of-SV values move together.
    if (trackBdtSfTag() == "_trkbdtsf") {
        variations.push_back({ "track_bdt_sfb_off", SAMPLE, GENERATOR, TF_GENERATOR,
                               "Track-BDT SF_{b} off", false, "nominal", SFUPART_NOM,
                               "", "Track-BDT SF_{b}", true, true, "", "boff" });
        variations.push_back({ "track_bdt_sfl_off", SAMPLE, GENERATOR, TF_GENERATOR,
                               "Track-BDT SF_{l} off", false, "nominal", SFUPART_NOM,
                               "", "Track-BDT SF_{l}", true, true, "", "loff" });
    }

    // ---- Correction-level variations (no second unfolding; see the loop below) --------
    // NOTE ON DOUBLE COUNTING: "unfolding_model" above already swaps the whole unfolding MC
    // to the other generator, which includes these very corrections. Keeping both in the
    // quadrature total counts the correction difference twice. Use these when you want the
    // corrections booked on their own -- and then drop "unfolding_model", or read it as the
    // deliberately conservative option.
    const TString corr_dir = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/";
    const TString corr_tag = Form("_%s_pt%d.root", SAMPLE.Data(), ibin_pt);
    // Correction-level variations are all EEC-weight modelling, which is meaningless once
    // the weight is off -- the correction itself is identically 1 there.
    // Not booked in the rho era either: the band is the four sources B can have (see the dr
    // rho branch above), and the eec_weight input file was made from the pre-_rho MC.
    std::vector<CorrVariation> corrVariations = (EEC_WEIGHT_OFF || rhoInputs())
      ? std::vector<CorrVariation>{}
      : std::vector<CorrVariation>{
        /* ---- disabled (kept for reference): 2SV+UParT efficiency variation ----
        // Dropped on request 2026-09-09: not used in the result. Commenting the entry
        // out removes it from the quadrature total, the printout, both plots and the
        // written histograms in one go -- everything downstream iterates this vector.
        // Re-enable by uncommenting; nothing else has to change.
        { "eff_2sv_btag", "2SV+UParT eff. " + prettyGen(other_gen), "2sv_btag",
          corr_dir + "UParT_efficiency_2sv_btag" + corr_tag, true  },
        ---- end disabled: eff_2sv_btag ---- */
        { "eec_weight",   "EEC weight "     + prettyGen(other_gen), "eec_weight",
          corr_dir + "UParT_efficiency_eec_weight" + corr_tag,  false, "MC modeling of EEC weight" },
    };
    // ==================================================================================

    // ---- Validate every booked variation name BEFORE touching a file ------------------
    // An unrecognised name poisons its tag (result_paths.h) instead of aliasing onto the
    // nominal, so it would show up as a missing file. Checking up front turns that into one
    // clear message naming the offending entry, rather than a variation quietly dropping
    // out of the band.
    {
        bool ok = true;
        for (const Variation &v : variations) {
            if (!isKnownSample(v.sample)) {
                std::cerr << "ERROR: variation '" << v.name << "' has unknown sample '"
                          << v.sample << "'" << std::endl; ok = false;
            }
            if (!isKnownGenerator(v.generator) || !isKnownGenerator(v.tfGenerator)) {
                std::cerr << "ERROR: variation '" << v.name << "' has an unknown generator"
                          << std::endl; ok = false;
            }
            if (!isKnownTfVariation(v.tfVariation)) {
                std::cerr << "ERROR: variation '" << v.name << "' has unknown tfVariation '"
                          << v.tfVariation << "'" << std::endl; ok = false;
            }
            if (!isKnownSfupartVariation(v.sfupartVariation)) {
                std::cerr << "ERROR: variation '" << v.name << "' has unknown sfupartVariation '"
                          << v.sfupartVariation << "'" << std::endl; ok = false;
            }
            // An observable with no UParT SF has none to vary either: such an entry would
            // build a path apply_unfolding_2d.C refuses to write, so it could only ever be
            // reported as missing. Say why here instead.
            if (SFUPART_NOM == "off" && v.sfupartVariation != "off") {
                std::cerr << "ERROR: variation '" << v.name << "' asks for sfupartVariation '"
                          << v.sfupartVariation << "', but OBSERVABLE '" << OBSERVABLE
                          << "' has no UParT SF -- use \"off\"." << std::endl; ok = false;
            }
        }
        // Two entries writing to the same name would silently overwrite each other's
        // histograms in the output file.
        for (size_t a = 0; a < variations.size(); ++a)
            for (size_t b = a + 1; b < variations.size(); ++b)
                if (variations[a].name == variations[b].name) {
                    std::cerr << "ERROR: two variations are both named '"
                              << variations[a].name << "'" << std::endl; ok = false;
                }
        if (!ok) {
            std::cerr << "Refusing to run: fix the variation list above." << std::endl;
            return;
        }
    }

    // The nominal always has trackEffUnc = false: every variation above is measured
    // against the unvaried MC. SFUPART_NOM is "nominal" for dr and "off" for B -- the nominal
    // B run applies no UParT SF and its path says so.
    const TString label       = resultLabel(SAMPLE, GENERATOR, test_mode, unfoldBayes, TF_GENERATOR,
                                            false, "nominal", SFUPART_NOM, EEC_WEIGHT_OFF,
                                            OBSERVABLE);
    const TString nominalFile = resultFile(SAMPLE, GENERATOR, test_mode, unfoldBayes, TF_GENERATOR,
                                           false, "nominal", SFUPART_NOM, EEC_WEIGHT_OFF,
                                           OBSERVABLE);

    // The folder every plot and the final file are written to. One expression, used
    // everywhere below -- it used to be spelled out at each of the six call sites, and the
    // observable would have had to be added to all six.
    const TString outFolder   = resultFolder(SAMPLE, GENERATOR, unfoldBayes, TF_GENERATOR,
                                             false, "nominal", SFUPART_NOM, EEC_WEIGHT_OFF,
                                             OBSERVABLE);

    std::cout << "Observable: " << OBSERVABLE << "  (" << obs_title << ")" << std::endl;
    std::cout << "Unfolding : " << (unfoldBayes ? "Bayesian" : "matrix inversion") << std::endl;
    std::cout << "Nominal : " << nominalFile << std::endl;
    std::cout << "Histogram: " << histname << std::endl;

    TH1D *h_nom = loadWeighted(nominalFile, histname, "h_nominal", w);
    if (!h_nom) {
        // The full call, out to the observable, as for the variations below: a B nominal
        // rebuilt from a two-argument call would be the dr one.
        std::cerr << "ERROR: no nominal result -- run apply_unfolding_2d.C(\"" << SAMPLE
                  << "\",\"" << GENERATOR << "\"," << test_mode << ","
                  << (unfoldBayes ? "true" : "false") << ",false,\"" << TF_GENERATOR
                  << "\",false,\"nominal\",\"" << SFUPART_NOM << "\","
                  << (EEC_WEIGHT_OFF ? "true" : "false") << ",\"" << OBSERVABLE
                  << "\") first" << std::endl;
        return;
    }

    const int nbins = h_nom->GetNbinsX();

    // Per-variation shifts, kept so each one can be plotted and debugged: h_shifts holds
    // |delta| (written out only), h_deltas the signed one the band is built from, h_vars the
    // variation curve itself (for the comparison plot).
    std::vector<TH1D *> h_shifts;
    std::vector<TH1D *> h_deltas;
    std::vector<TH1D *> h_vars;
    std::vector<TString> shift_names;
    std::vector<TString> shift_labels;
    // Parallel to h_shifts: the source name shown on the breakdown plot (Variation::sysLabel,
    // shared by the two sides of a pair).
    std::vector<TString> shift_syslabels;
    // Parallel to h_shifts: the group id each entry belongs to ("" = independent source).
    std::vector<TString> shift_groups;
    // Parallel to h_shifts: false = drawn and printed but kept out of the band.
    std::vector<bool> shift_in_band;
    // Parallel to h_shifts: true = |delta| into both sides (Variation::symmetrise).
    std::vector<bool> shift_symm;
    // Particle-level MC curves, filled by the data-vs-gen plot below and written out with
    // the rest so the comparison can be redrawn from the final file alone.
    std::vector<TH1D *> h_gen_curves;
    std::vector<TString> gen_curve_names;
    // Display names for the "systematics applied" tag on the data-vs-gen plot, in booking
    // order and deduplicated -- an up/down pair shares one sysLabel and is listed once.
    // Built from the entries that actually made it into the band, so the tag cannot claim a
    // source whose result file was missing and got skipped.
    std::vector<TString> syst_tag_labels;
    auto addTagLabel = [&](const TString &lbl) {
        if (!lbl.Length()) return;
        for (const TString &t : syst_tag_labels) if (t == lbl) return;
        syst_tag_labels.push_back(lbl);
    };

    for (const Variation &v : variations) {
        // A variation with its own SV_SF builds its path under that setting, then the
        // band's own setting is put back for everything else.
        const char *svsf_env  = std::getenv("SV_SF");
        const bool  svsf_had  = (svsf_env != nullptr);
        const TString svsf_prev = svsf_had ? TString(svsf_env) : TString("");
        const char *trkbdt_env = std::getenv("TRACK_BDT_SF");
        const bool trkbdt_had = (trkbdt_env != nullptr);
        const TString trkbdt_prev = trkbdt_had ? TString(trkbdt_env) : TString("");
        if (v.svsfVariation.Length()) gSystem->Setenv("SV_SF", v.svsfVariation);
        if (v.trackBdtSfVariation.Length())
            gSystem->Setenv("TRACK_BDT_SF", v.trackBdtSfVariation);
        const TString varFile =
            resultFile(v.sample, v.generator, test_mode, unfoldBayes, v.tfGenerator,
                       v.trackEffUnc, v.tfVariation, v.sfupartVariation, EEC_WEIGHT_OFF,
                       OBSERVABLE);
        if (v.svsfVariation.Length()) {
            if (svsf_had) gSystem->Setenv("SV_SF", svsf_prev);
            else          gSystem->Unsetenv("SV_SF");
        }
        if (v.trackBdtSfVariation.Length()) {
            if (trkbdt_had) gSystem->Setenv("TRACK_BDT_SF", trkbdt_prev);
            else            gSystem->Unsetenv("TRACK_BDT_SF");
        }
        std::cout << "Variation '" << v.name << "' : " << varFile << std::endl;

        TH1D *h_var = loadWeighted(varFile, histname, "h_var_" + v.name, w);
        if (!h_var) {
            // The exact call that produces it, all the way out to the observable: the
            // trailing arguments are what a B or yield variation differs by, and a command
            // that stopped at the SF variation would rebuild the dr result instead.
            std::cerr << "   -> skipped (run "
                      << (v.svsfVariation.Length() ? "SV_SF=" + v.svsfVariation + " " : TString(""))
                      << (v.trackBdtSfVariation.Length()
                              ? "TRACK_BDT_SF=" + v.trackBdtSfVariation + " " : TString(""))
                      << "apply_unfolding_2d.C(\"" << v.sample << "\",\""
                      << v.generator << "\"," << test_mode << "," << (unfoldBayes ? "true" : "false")
                      << ",false,\"" << v.tfGenerator << "\","
                      << (v.trackEffUnc ? "true" : "false")
                      << ",\"" << v.tfVariation << "\",\"" << v.sfupartVariation << "\","
                      << (EEC_WEIGHT_OFF ? "true" : "false") << ",\"" << OBSERVABLE << "\""
                      << ") to include it)" << std::endl;
            continue;
        }
        if (h_var->GetNbinsX() != nbins) {
            std::cerr << "   -> skipped: " << h_var->GetNbinsX() << " bins vs " << nbins
                      << " in the nominal" << std::endl;
            continue;
        }

        // The signed shift is what the band is built from (see the total below); |shift| is
        // kept alongside it for the output file.
        TH1D *h_shift  = (TH1D *) h_nom->Clone("h_syst_" + v.name);
        TH1D *h_signed = (TH1D *) h_nom->Clone("h_delta_" + v.name);
        h_shift->SetDirectory(nullptr);  h_shift->Reset();
        h_signed->SetDirectory(nullptr); h_signed->Reset();
        for (int i = 1; i <= nbins; ++i) {
            const double d = h_nom->GetBinContent(i) - h_var->GetBinContent(i);
            h_signed->SetBinContent(i, d);
            h_shift->SetBinContent(i, std::fabs(d));
        }

        h_shifts.push_back(h_shift);
        h_deltas.push_back(h_signed);
        h_vars.push_back(h_var);
        shift_names.push_back(v.name);
        shift_labels.push_back(v.label.Length() ? v.label : v.name);
        shift_syslabels.push_back(v.sysLabel.Length() ? v.sysLabel : (v.label.Length() ? v.label : v.name));
        shift_groups.push_back(v.group);
        shift_in_band.push_back(v.inBand);
        shift_symm.push_back(v.symmetrise);
        if (v.inBand) addTagLabel(v.sysLabel);
        else std::cout << "   -> display only, not in the band" << std::endl;
    }

    // ---- Correction-level variations, evaluated on the final result -------------------
    // These do not need another unfolding run: the correction enters the result as a
    // per-dr factor, so swapping the generator's correction for the other one is a bin-by-bin
    // rescale of the nominal, followed by RE-NORMALISING to unit area.
    //
    // The re-normalisation is the point. The EEC is a shape, so a correction ratio that is
    // flat in dr cancels completely and only its SHAPE contributes. Comparing the raw
    // correction curves (as plotNice_UParT_efficiency draws them) overstates the effect --
    // e.g. the 2SV+UParT efficiency differs by 3-4%, but almost all of that is a constant
    // offset that self-normalisation removes.
    //
    //   efficiency (divided into the data) : varied = nominal * eff_nominalGen / eff_otherGen
    //   EEC weight (multiplied)            : varied = nominal * w_otherGen  / w_nominalGen
    //
    // Inputs are the histograms plotNice_UParT_efficiency() writes; run it first for each
    // correction, at the SAME pT bin as the unfolding.
    for (const CorrVariation &cv : corrVariations) {
        TFile *fc = TFile::Open(cv.file);
        if (!fc || fc->IsZombie()) {
            std::cerr << "Correction variation '" << cv.name << "': MISSING " << cv.file
                      << "\n   -> skipped (run plotNice_UParT_efficiency(\"" << SAMPLE << "\","
                      << ibin_pt << ",\"" << cv.effType << "\") to include it)" << std::endl;
            continue;
        }
        TH1D *h_gen_nom = dynamic_cast<TH1D *>(fc->Get(GENERATOR == "pythia" ? "h_eff_pythia"
                                                                             : "h_eff_herwig"));
        TH1D *h_gen_alt = dynamic_cast<TH1D *>(fc->Get(GENERATOR == "pythia" ? "h_eff_herwig"
                                                                             : "h_eff_pythia"));
        if (!h_gen_nom || !h_gen_alt) {
            std::cerr << "ERROR: efficiency histograms not found in " << cv.file << std::endl;
            fc->Close();
            continue;
        }
        if (h_gen_nom->GetNbinsX() != nbins) {
            std::cerr << "Correction variation '" << cv.name << "': " << h_gen_nom->GetNbinsX()
                      << " bins vs " << nbins << " in the result -- skipped" << std::endl;
            fc->Close();
            continue;
        }

        TH1D *h_var = (TH1D *) h_nom->Clone("h_var_" + cv.name);
        h_var->SetDirectory(nullptr);
        for (int i = 1; i <= nbins; ++i) {
            const double c_nom = h_gen_nom->GetBinContent(i);
            const double c_alt = h_gen_alt->GetBinContent(i);
            if (c_nom == 0. || c_alt == 0.) continue;
            const double scale = cv.inverse ? (c_nom / c_alt)   // correction divided into data
                                            : (c_alt / c_nom);  // correction multiplied in
            h_var->SetBinContent(i, h_var->GetBinContent(i) * scale);
            h_var->SetBinError(i, h_var->GetBinError(i) * scale);
        }
        fc->Close();

        // Self-normalise, exactly as the nominal is normalised.
        if (w.normalise_unit_area) {
            const double integral = h_var->Integral();
            if (integral > 0.) h_var->Scale(1. / integral, "width");
        }

        TH1D *h_shift  = (TH1D *) h_nom->Clone("h_syst_" + cv.name);
        TH1D *h_signed = (TH1D *) h_nom->Clone("h_delta_" + cv.name);
        h_shift->SetDirectory(nullptr);  h_shift->Reset();
        h_signed->SetDirectory(nullptr); h_signed->Reset();
        for (int i = 1; i <= nbins; ++i) {
            const double d = h_nom->GetBinContent(i) - h_var->GetBinContent(i);
            h_signed->SetBinContent(i, d);
            h_shift->SetBinContent(i, std::fabs(d));
        }

        std::cout << "Correction variation '" << cv.name << "' : " << cv.file << std::endl;
        h_shifts.push_back(h_shift);
        h_deltas.push_back(h_signed);
        h_vars.push_back(h_var);
        shift_names.push_back(cv.name);
        shift_labels.push_back(cv.label);
        shift_syslabels.push_back(cv.sysLabel.Length() ? cv.sysLabel : cv.label);
        shift_groups.push_back("");   // correction variations are independent sources
        shift_in_band.push_back(true);
        shift_symm.push_back(false);  // signed, like every source except tracking_eff
        addTagLabel(cv.sysLabel);
    }

    if (h_shifts.empty())
        std::cerr << "WARNING: no variation available -- the output carries stat errors only"
                  << std::endl;

    // ---- Total systematic: an ASYMMETRIC quadrature sum -----------------------------
    // Two bands, sigma+ and sigma-, because a genuine up/down variation does not have to
    // move the result by the same amount in each direction. Which side a source feeds
    // depends on whether it is really two-sided:
    //
    //   GROUPED entries (a shared, non-empty Variation::group) are the two sides of ONE
    //     uncertainty. Each member is signed: the ones that push the result UP set sigma+,
    //     the ones that push it DOWN set sigma-, and the source contributes its largest
    //     shift on each side -- once, never both members to the same side.
    //
    //   UNGROUPED entries have a single alternative (a Pythia<->Herwig swap, the 3% track
    //     drop, the EEC-weight swap). They are SIGNED too (Zoe, 2026-09-29): in each bin the
    //     shift goes into sigma+ if it pushes the result up and into sigma- if it pushes it
    //     down, never both. Upward shifts are summed in quadrature into sigma+, downward
    //     ones into sigma-, separately.
    //
    //   EXCEPTION: entries with Variation::symmetrise = true -- only tracking_eff -- put
    //     |delta| into BOTH sides.
    //
    // A pair whose two members push the SAME way in a bin contributes only to that side
    // there (its larger shift, once) and nothing to the other. This used to fall back to
    // the envelope on both sides; that symmetrisation is disabled below.
    TH1D *h_syst_up = (TH1D *) h_nom->Clone("h_syst_total_up");
    TH1D *h_syst_dn = (TH1D *) h_nom->Clone("h_syst_total_down");
    for (TH1D *h : {h_syst_up, h_syst_dn}) { h->SetDirectory(nullptr); h->Reset(); }
    std::vector<TString> one_sided_groups;
    for (int i = 1; i <= nbins; ++i) {
        double s2p = 0., s2m = 0.;
        std::map<TString, double> grp_up, grp_dn;
        for (size_t k = 0; k < h_shifts.size(); ++k) {
            if (!shift_in_band[k]) continue;   // display-only comparison, not a systematic
            // h_deltas is (nominal - variation); the shift OF THE RESULT is the negative.
            const double d = -h_deltas[k]->GetBinContent(i);
            if (shift_groups[k].Length()) {
                double &u = grp_up[shift_groups[k]];   // value-initialised to 0 on insert
                double &m = grp_dn[shift_groups[k]];
                if (d > 0. && d  > u) u =  d;
                if (d < 0. && -d > m) m = -d;
            } else if (shift_symm[k]) {
                s2p += d * d;
                s2m += d * d;
            } else if (d > 0.) {
                s2p += d * d;
            } else {
                s2m += d * d;
            }
        }
        for (auto &g : grp_up) {
            double up = g.second, dn = grp_dn[g.first];
            if (up == 0. || dn == 0.) {
                /* ---- disabled (kept for reference): same-side pair symmetrisation ----
                // Both sides of the pair went the same way in this bin: symmetrise it.
                const double env = std::max(up, dn);
                up = dn = env;
                ---- end disabled: same-side pair symmetrisation ---- */
                bool seen = false;
                for (const TString &t : one_sided_groups) if (t == g.first) seen = true;
                if (!seen) one_sided_groups.push_back(g.first);
            }
            s2p += up * up;
            s2m += dn * dn;
        }
        h_syst_up->SetBinContent(i, std::sqrt(s2p));
        h_syst_dn->SetBinContent(i, std::sqrt(s2m));
    }
    for (const TString &g : one_sided_groups)
        std::cerr << "NOTE: both sides of pair '" << g << "' shift the result the same way "
                  << "in at least one bin -- it contributes to that side only there."
                  << std::endl;

    // The symmetrised total, for anything that wants one number per bin (the TH1 bin errors
    // below, and any downstream script that predates the split). The larger side, so it is
    // never the optimistic one.
    TH1D *h_syst = (TH1D *) h_nom->Clone("h_syst_total");
    h_syst->SetDirectory(nullptr);
    h_syst->Reset();
    for (int i = 1; i <= nbins; ++i)
        h_syst->SetBinContent(i, std::max(h_syst_up->GetBinContent(i),
                                          h_syst_dn->GetBinContent(i)));

    // ---- Result with stat + syst ----
    // h_total keeps the central values with the quadrature error in the bin errors; the
    // graph is what you draw as the systematic box around the points.
    TH1D *h_total = (TH1D *) h_nom->Clone("h_result_stat_syst");
    h_total->SetDirectory(nullptr);
    // Asymmetric: this is the band you draw, and the two sides are independent numbers.
    TGraphAsymmErrors *g_syst = new TGraphAsymmErrors(nbins);
    g_syst->SetName("g_syst_band");

    // One signed-shift column per variation -- the shift OF THE RESULT (variation - nominal),
    // so the sign in the table is the direction that source actually pushes the point, and
    // it matches which of syst+/syst- the source feeds.
    // Column header, 15 chars wide. A longer name is elided in the MIDDLE, not chopped at
    // the end: the tail is what tells the two sides of a pair apart (sfupart_jpsyst_up vs
    // ..._dn), so a plain "%.15s" would print two identical headers over two different
    // columns -- which is how a comparison gets misread.
    auto colHeader = [](const TString &name, bool in_band) -> TString {
        const TString s = in_band ? name : ("*" + name);
        if (s.Length() <= 15) return s;
        return TString(s(0, 7)) + "~" + TString(s(s.Length() - 7, 7));
    };
    printf("\n bin |        value |     stat");
    for (size_t k = 0; k < shift_names.size(); ++k)
        printf(" | %15.15s", colHeader(shift_names[k], shift_in_band[k]).Data());
    printf(" |    syst+ |    syst- |    total |  syst/value\n");
    printf("-----+--------------+---------");
    for (size_t k = 0; k < shift_names.size(); ++k) printf("-+----------------");
    printf("-+----------+----------+----------+------------\n");

    for (int i = 1; i <= nbins; ++i) {
        const double v     = h_nom->GetBinContent(i);
        const double stat  = h_nom->GetBinError(i);
        const double sup   = h_syst_up->GetBinContent(i);
        const double sdn   = h_syst_dn->GetBinContent(i);
        const double syst  = h_syst->GetBinContent(i);           // the larger side
        const double tot   = std::sqrt(stat * stat + syst * syst);

        h_total->SetBinError(i, tot);
        g_syst->SetPoint(i - 1, h_nom->GetBinCenter(i), v);
        g_syst->SetPointError(i - 1, 0.5 * h_nom->GetBinWidth(i), 0.5 * h_nom->GetBinWidth(i),
                              sdn, sup);   // low, high

        printf(" %3d | %12.5g | %8.4g", i, v, stat);
        for (TH1D *h : h_deltas) printf(" | %+15.4g", -h->GetBinContent(i));
        printf(" | %8.4g | %8.4g | %8.4g | %10.3f\n",
               sup, sdn, tot, (v != 0. ? syst / std::fabs(v) : 0.));
    }
    {
        bool any_display_only = false;
        for (bool b : shift_in_band) if (!b) any_display_only = true;
        if (any_display_only)
            printf("\n * = shown for comparison only; NOT included in syst+/syst-.\n");
    }

    // ---- Plot: every variation's EEC curve against the nominal ----------------------
    // Top pad: the EEC itself, nominal with its stat errors plus one curve per variation.
    // Bottom pad: variation/nominal - 1, i.e. the fractional pull of each systematic, with
    // the total band drawn around zero so each source can be read against the total.
    if (!h_vars.empty()) {
        // Palette convention: red = the measurement, the rest for the variations.
        const Color_t col_nom  = (Color_t) TColor::GetColor("#C44E52");
        // Light enough to sit behind the curves, dark enough to actually read. 0.10 is
        // effectively white on paper; raise or lower this one number to taste.
        const Color_t col_band = blendWithWhite("#C44E52", 0.2); // systematic band
        // Red stays the measurement; the variations take blue, green, purple, then the
        // spares. k % size() only reads as distinct while size() >= the number of variations
        // booked, so this list has to grow with the band -- at six, the three UParT SF
        // curves wrapped onto the same colours as the generator swaps.
        const std::vector<Color_t> col_var = { (Color_t) TColor::GetColor("#4C72B0"),   // blue
                                               (Color_t) TColor::GetColor("#4F8F52"),   // green
                                               (Color_t) TColor::GetColor("#8C6BB1"),   // purple
                                               (Color_t) TColor::GetColor("#C97430"),   // orange
                                               (Color_t) TColor::GetColor("#3D9E9E"),   // teal
                                               (Color_t) TColor::GetColor("#937860"),   // brown
                                               (Color_t) TColor::GetColor("#CC79A7"),   // pink
                                               (Color_t) TColor::GetColor("#7F7F2E"),   // olive
                                               (Color_t) TColor::GetColor("#56638A"),   // slate
                                               (Color_t) TColor::GetColor("#B04A3A"),   // brick
                                               (Color_t) TColor::GetColor("#2F6F6F"),   // deep teal
                                               (Color_t) TColor::GetColor("#6A5ACD"),   // slate blue
                                               // 2026-10-13: five SV SF pairs -> 15 curves
                                               (Color_t) TColor::GetColor("#D4A017"),   // mustard
                                               (Color_t) TColor::GetColor("#8B3A62"),   // plum
                                               (Color_t) TColor::GetColor("#5F9EA0"),   // cadet
                                               (Color_t) TColor::GetColor("#A0522D") }; // sienna
        if (h_vars.size() > col_var.size())
            std::cerr << "NOTE: " << h_vars.size() << " variations but only " << col_var.size()
                      << " colours -- some curves repeat a colour and are told apart only by "
                      << "line style. Add entries to col_var." << std::endl;
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
        gStyle->SetLegendBorderSize(0);
        gStyle->SetLegendFillColor(0);

        // Same layout and fonts as the bottomline plot in apply_unfolding_2d.C.
        const Float_t font_scale   = 1200. / 800.;
        const Style_t font_code    = 43;
        const Float_t label_size   = 15. * font_scale;
        const Float_t title_size   = 15. * font_scale;
        const Float_t legend_size  = 15. * font_scale;
        const Float_t title_offset = 3.0;

        TCanvas *c_sys = new TCanvas("c_syst_curves", "", 800, 800);
        // 45/55 rather than the usual 30/70: on THIS plot the ratio pad is the content --
        // the per-source pulls are a few percent and all six curves sit on top of each
        // other in the top pad, so squeezing them into 0.3 of the canvas hides the answer.
        TPad *p_top = new TPad("p_top", "", 0., 0.45, 1., 1.);
        TPad *p_bot = new TPad("p_bot", "", 0., 0.,   1., 0.45);
        for (TPad *p : {p_top, p_bot}) {
            p->SetTicks(1, 0);
            p->SetFillColor(0);
        }
        p_top->SetMargin(0.1, 0.1, 0.0,  0.1);   // left, right, bottom, top
        p_bot->SetMargin(0.1, 0.1, 0.16, 0.0);
        c_sys->cd();
        p_top->Draw(); p_bot->Draw();

        // ---- top ----
        p_top->cd();
        TH1D *h_nom_draw = (TH1D *) h_nom->Clone("h_nom_draw");
        h_nom_draw->SetDirectory(nullptr);
        double ymax = 0.;
        for (int i = 1; i <= nbins; ++i)
            ymax = std::max(ymax, h_nom_draw->GetBinContent(i)
                                + h_nom_draw->GetBinError(i) + h_syst_up->GetBinContent(i));
        for (TH1D *h : h_vars)
            for (int i = 1; i <= nbins; ++i) ymax = std::max(ymax, h->GetBinContent(i));
        h_nom_draw->SetTitle("");
        h_nom_draw->GetYaxis()->SetRangeUser(0., ymax * 1.75);
        // TLatex (#Delta), not TMathText: TMathText axis titles do not render here.
        h_nom_draw->GetYaxis()->SetTitle(obs_title);
        h_nom_draw->GetYaxis()->CenterTitle(true);
        h_nom_draw->GetYaxis()->SetTitleFont(font_code); h_nom_draw->GetYaxis()->SetTitleSize(title_size);
        h_nom_draw->GetYaxis()->SetTitleOffset(1.5); // as the bottomline plot's main pad
        h_nom_draw->GetYaxis()->SetLabelFont(font_code); h_nom_draw->GetYaxis()->SetLabelSize(label_size);
        h_nom_draw->GetXaxis()->SetTitleSize(0);
        h_nom_draw->GetXaxis()->SetLabelSize(0);
        h_nom_draw->SetLineColor(col_nom);   h_nom_draw->SetMarkerColor(col_nom);
        h_nom_draw->SetMarkerStyle(kFullCircle); h_nom_draw->SetMarkerSize(1);
        h_nom_draw->SetLineWidth(2);
        fixXRange(h_nom_draw);
        h_nom_draw->Draw("AXIS");

        // Total systematic as a light red band on the nominal, drawn first so the points
        // and the variation curves stay readable on top of it. Same colour as the band in
        // the ratio pad below, so the two read as one quantity.
        g_syst->SetFillColor(col_band);
        g_syst->SetLineWidth(0);
        g_syst->Draw("2 same");

        // Two columns across the top: eight rows in one column reach down into the peak
        // whichever side they are parked on, and halving the row count clears it.
        // Three columns: at two, fourteen entries reach down into the peak whatever the
        // headroom. Rows = ceil(entries / columns).
        const int n_cols = 3;
        const int n_rows = (2 + (int) h_vars.size() + n_cols - 1) / n_cols;
        TLegend *leg = new TLegend(0.10, 0.87 - 0.040 * n_rows, 0.95, 0.87);
        leg->SetNColumns(n_cols);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetMargin(0.15);
        leg->SetTextFont(font_code);
        leg->SetTextSize(legend_size * 0.62);  // three columns: the long labels have a third of the width
        leg->AddEntry(h_nom_draw, "nominal (stat)", "pe1");
        leg->AddEntry(g_syst, "quadrature systematic", "f");

        for (size_t k = 0; k < h_vars.size(); ++k) {
            TH1D *h = h_vars[k];
            h->SetLineColor(col_var[k % col_var.size()]);
            h->SetMarkerColor(col_var[k % col_var.size()]);
            h->SetLineWidth(2);
            // ROOT only defines line styles 1..10; the old 2 + k ran past that once the
            // band grew, and an invalid style draws nothing at all -- including in the
            // legend, where the entry silently lost its swatch. Cycling a short set of
            // clearly distinct patterns keeps every curve drawable, and pairs with the
            // colour so a repeated colour still reads as a different curve.
            h->SetLineStyle(1 + (int) (k % 4));
            h->Draw("HIST same");
            leg->AddEntry(h, shift_labels[k], "l");
        }
        h_nom_draw->Draw("PE X0 same");
        leg->Draw();

        // Same banner as the bottomline plot.
        TLatex cms_bl;
        cms_bl.SetNDC();
        cms_bl.SetTextFont(62); cms_bl.SetTextSize(0.045); cms_bl.DrawLatex(0.10, 0.945, "CMS");
        cms_bl.SetTextFont(52); cms_bl.SetTextSize(0.037); cms_bl.DrawLatex(0.205, 0.945, "Work In Progress");
        cms_bl.SetTextFont(42); cms_bl.SetTextSize(0.037); cms_bl.DrawLatex(0.66, 0.945, "pp #sqrt{s} = 5.36 TeV");

        // ---- bottom: fractional pull of each variation ----
        p_bot->cd();
        TH1D *h_frame = (TH1D *) h_nom->Clone("h_ratio_frame");
        h_frame->SetDirectory(nullptr);
        h_frame->Reset();
        double rmax = 0.;
        for (int i = 1; i <= nbins; ++i) {
            const double v = h_nom->GetBinContent(i);
            if (v == 0.) continue;
            for (TH1D *h : h_vars)
                rmax = std::max(rmax, std::fabs(h->GetBinContent(i) / v - 1.));
            rmax = std::max(rmax, std::max(h_syst_up->GetBinContent(i),
                                           h_syst_dn->GetBinContent(i)) / std::fabs(v));
        }
        if (rmax <= 0.) rmax = 0.05;
        h_frame->SetTitle("");
        h_frame->GetYaxis()->SetRangeUser(-1.6 * rmax, 1.6 * rmax);
        h_frame->GetXaxis()->SetTitle(obs_axis);
        // Short enough to fit the pad height -- "variation / nominal - 1" ran off the end
        // and lost its "1", which made the panel look like it plotted a ratio, not a shift.
        h_frame->GetYaxis()->SetTitle("Variation / nominal #minus 1");
        h_frame->GetXaxis()->CenterTitle(true);
        h_frame->GetYaxis()->CenterTitle(true);
        h_frame->GetYaxis()->SetNdivisions(505);
        for (TAxis *ax : { h_frame->GetXaxis(), h_frame->GetYaxis() }) {
            ax->SetTitleFont(font_code);
            ax->SetTitleSize(title_size);
            ax->SetLabelFont(font_code);
            ax->SetLabelSize(label_size);
        }
        // title_offset (3.0) is what the bottomline ratio pad uses, but with these margins
        // it puts both titles outside the pad, so they are set to what actually renders.
        h_frame->GetXaxis()->SetTitleOffset(1.1);
        h_frame->GetYaxis()->SetTitleOffset(1.5);
        fixXRange(h_frame);
        h_frame->Draw("AXIS");

        // Total systematic as a light red band around zero, for scale.
        TGraphAsymmErrors *g_band = new TGraphAsymmErrors(nbins);
        g_band->SetName("g_syst_band_rel");
        for (int i = 1; i <= nbins; ++i) {
            const double v  = h_nom->GetBinContent(i);
            const double hw = 0.5 * h_nom->GetBinWidth(i);
            g_band->SetPoint(i - 1, h_nom->GetBinCenter(i), 0.);
            g_band->SetPointError(i - 1, hw, hw,
                                  (v != 0.) ? h_syst_dn->GetBinContent(i) / std::fabs(v) : 0.,
                                  (v != 0.) ? h_syst_up->GetBinContent(i) / std::fabs(v) : 0.);
        }
        g_band->SetFillColor(col_band);
        g_band->SetLineWidth(0);
        g_band->Draw("2 same");

        for (size_t k = 0; k < h_vars.size(); ++k) {
            TH1D *h_r = (TH1D *) h_vars[k]->Clone("h_relshift_" + shift_names[k]);
            h_r->SetDirectory(nullptr);
            for (int i = 1; i <= nbins; ++i) {
                const double v = h_nom->GetBinContent(i);
                h_r->SetBinContent(i, (v != 0.) ? h_vars[k]->GetBinContent(i) / v - 1. : 0.);
                h_r->SetBinError(i, 0.);
            }
            h_r->SetLineColor(col_var[k % col_var.size()]);
            h_r->SetLineWidth(2);
            h_r->SetLineStyle(1 + (int) (k % 4));   // same cycle as the top pad
            h_r->Draw("HIST same");
        }

        TLine *l0 = new TLine(h_frame->GetXaxis()->GetXmin(), 0.,
                              h_frame->GetXaxis()->GetXmax(), 0.);
        l0->SetLineStyle(2);
        l0->Draw();
        p_bot->RedrawAxis();

        const TString plot_stem = outFolder
                                + "systematics_curves_" + label;
        c_sys->Print(plot_stem + ".pdf");
        c_sys->Print(plot_stem + ".png");
        std::cout << "\nWrote " << plot_stem << ".{pdf,png}" << std::endl;
    }

    // ---- Plot: systematic breakdown per bin, CMS bar style (Zoe, 2026-10-13) -----------
    // In each bin, one bar per SOURCE side by side (an up/down pair is one source), going
    // from -sigma_down to +sigma_up of that source exactly as it enters the band: a pair
    // gives its larger shift on each side, a single variation only the side it pushes
    // towards, tracking efficiency both sides. The open black box is the total (syst+ up,
    // syst- down). Absolute, in the units of the result. Display-only entries are left out.
    {
        struct BarSrc { TString key, title; std::vector<double> up, dn; };
        std::vector<BarSrc> bsrc;
        for (size_t k = 0; k < h_deltas.size(); ++k) {
            if (!shift_in_band[k]) continue;
            const TString key = shift_groups[k].Length() ? shift_groups[k] : shift_names[k];
            BarSrc *b = nullptr;
            for (BarSrc &x : bsrc) if (x.key == key) b = &x;
            if (!b) { bsrc.push_back({ key, shift_syslabels[k], std::vector<double>(nbins + 1, 0.),
                                       std::vector<double>(nbins + 1, 0.) }); b = &bsrc.back(); }
            for (int i = 1; i <= nbins; ++i) {
                const double d = -h_deltas[k]->GetBinContent(i);   // shift of the result
                if (shift_symm[k] && !shift_groups[k].Length()) {
                    b->up[i] = std::max(b->up[i], std::fabs(d)); b->dn[i] = std::max(b->dn[i], std::fabs(d));
                } else if (d > 0.) b->up[i] = std::max(b->up[i], d);
                else               b->dn[i] = std::max(b->dn[i], -d);
            }
        }
        // The GROUPED view (Zoe, 2026-10-13): one bar per category, its sources added in
        // quadrature side by side. b tagging = every SV SF contribution (+ the UParT SF they
        // move with, for dr); anything not listed keeps its own bar.
        auto category = [&](const TString &key) -> TString {
            if (key.BeginsWith("sv_sf") || key == "jp_hf" || key == "qq_rate")
                return (OBSERVABLE == "dr") ? "b tagging (UParT + SV SF)" : "b tagging (SV SF)";
            return "";
        };
        std::vector<BarSrc> bcat;
        for (const BarSrc &b : bsrc) {
            const TString cat = category(b.key);
            const TString ck  = cat.Length() ? cat : b.key;
            BarSrc *c = nullptr;
            for (BarSrc &x : bcat) if (x.key == ck) c = &x;
            if (!c) { bcat.push_back({ ck, cat.Length() ? cat : b.title, std::vector<double>(nbins + 1, 0.),
                                       std::vector<double>(nbins + 1, 0.) }); c = &bcat.back(); }
            for (int i = 1; i <= nbins; ++i) {
                c->up[i] = std::hypot(c->up[i], b.up[i]);
                c->dn[i] = std::hypot(c->dn[i], b.dn[i]);
            }
        }
        // Display order (Zoe, 2026-10-13): b tagging, tracking efficiency, light-jet mistag,
        // MC template modeling, detector response, then anything else. Stable, so the b-tag
        // contributions keep their booking order (stat, no muon, bb, JP HF, qq) in "detailed".
        auto rank = [&](const BarSrc &b) -> int {
            if (b.key.BeginsWith("sv_sf") || b.key == "jp_hf" || b.key == "qq_rate" || b.key.BeginsWith("b tagging")) return 0;
            if (b.key == "tracking_eff")    return 1;
            if (b.key == "mistag_0B")       return 2;
            if (b.key == "template_fit")    return 3;
            if (b.key == "unfolding_model") return 4;
            return 5;
        };
        auto byRank = [&](const BarSrc &a, const BarSrc &b) { return rank(a) < rank(b); };
        std::stable_sort(bsrc.begin(), bsrc.end(), byRank);
        std::stable_sort(bcat.begin(), bcat.end(), byRank);
        for (const std::pair<const char *, std::vector<BarSrc> *> view :
             { std::pair<const char *, std::vector<BarSrc> *>{"detailed", &bsrc},
               std::pair<const char *, std::vector<BarSrc> *>{"grouped",  &bcat} }) {
        const std::vector<BarSrc> &bsrc_v = *view.second;
        const TString vtag = view.first;
        if (!bsrc_v.empty()) {
            double ymax = 0.;
            for (int i = 1; i <= nbins; ++i)
                ymax = std::max({ ymax, h_syst_up->GetBinContent(i), h_syst_dn->GetBinContent(i) });
            // Pastel palette (Zoe, 2026-10-13), in display order.
            const std::vector<Color_t> col_bar = {
                (Color_t) TColor::GetColor("#7FB3D5"), (Color_t) TColor::GetColor("#82C99A"),
                (Color_t) TColor::GetColor("#F5B26B"), (Color_t) TColor::GetColor("#B39DDB"),
                (Color_t) TColor::GetColor("#F28B8B"), (Color_t) TColor::GetColor("#F7DC6F"),
                (Color_t) TColor::GetColor("#76D7C4"), (Color_t) TColor::GetColor("#F1A7C6"),
                (Color_t) TColor::GetColor("#C5A880"), (Color_t) TColor::GetColor("#A9CCE3"),
                (Color_t) TColor::GetColor("#D7BDE2"), (Color_t) TColor::GetColor("#ABEBC6") };
            TCanvas *c_br = new TCanvas("c_breakdown_" + vtag + "_" + label, "", 800, 800);
            c_br->SetMargin(0.16, 0.04, 0.12, 0.07);
            c_br->SetTicks(1, 1);
            const double x0 = h_nom->GetXaxis()->GetBinLowEdge(1);
            const double x1 = h_nom->GetXaxis()->GetBinUpEdge(nbins);
            TH1F *fr = c_br->DrawFrame(x0, -1.4 * ymax, x1, 2.3 * ymax);
            fr->GetXaxis()->SetTitle(OBSERVABLE == "dr" ? TString("#Delta r") : obs_sym);
            fr->GetYaxis()->SetTitle("Systematic uncertainty (" + obs_title + ", normalised)");
            fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.045);
            fr->GetXaxis()->SetLabelSize(0.042); fr->GetYaxis()->SetLabelSize(0.040);
            fr->GetYaxis()->SetTitleOffset(1.6); fr->GetXaxis()->SetTitleOffset(1.1);
            fr->GetXaxis()->CenterTitle(false);
            TLine *l0 = new TLine(x0, 0., x1, 0.); l0->SetLineColor(kBlack); l0->Draw();
            TLegend *lg = new TLegend(0.55, 0.89 - 0.040 * (bsrc_v.size() + 1), 0.93, 0.89);
            lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.030);
            const int ns = (int) bsrc_v.size();
            for (int i = 1; i <= nbins; ++i) {
                const double lo = h_nom->GetXaxis()->GetBinLowEdge(i), w_ = h_nom->GetXaxis()->GetBinWidth(i);
                const double inner = 0.80 * w_, bw = inner / ns, start = lo + 0.10 * w_;
                for (int s = 0; s < ns; ++s) {
                    TBox *bx = new TBox(start + s * bw, -bsrc_v[s].dn[i], start + (s + 1) * bw, bsrc_v[s].up[i]);
                    bx->SetFillColor(col_bar[s % col_bar.size()]); bx->SetLineWidth(0);
                    bx->Draw();
                }
                TBox *tot = new TBox(lo + 0.05 * w_, -h_syst_dn->GetBinContent(i), lo + 0.95 * w_, h_syst_up->GetBinContent(i));
                tot->SetFillStyle(0); tot->SetLineColor(kBlack); tot->SetLineWidth(2);
                tot->Draw("l");
            }
            for (int s = 0; s < ns; ++s) {
                TH1F *dummy = new TH1F(Form("lgbar_%s_%s_%d", vtag.Data(), label.Data(), s), "", 1, 0, 1);
                dummy->SetFillColor(col_bar[s % col_bar.size()]); dummy->SetLineWidth(0);
                lg->AddEntry(dummy, bsrc_v[s].title, "f");
            }
            TH1F *dtot = new TH1F("lgbar_tot_" + vtag + "_" + label, "", 1, 0, 1);
            dtot->SetLineColor(kBlack); dtot->SetLineWidth(2);
            lg->AddEntry(dtot, "Total", "l");
            lg->Draw();
            TLatex t; t.SetNDC();
            t.SetTextFont(62); t.SetTextSize(0.060); t.DrawLatex(0.20, 0.84, "CMS");
            t.SetTextFont(52); t.SetTextSize(0.040); t.DrawLatex(0.20, 0.795, "Work In Progress");
            t.SetTextFont(42); t.SetTextSize(0.034);
            t.DrawLatex(0.20, 0.745, obs_title + ", 2B jets");
            t.DrawLatex(0.20, 0.705, Form("p_{T}^{jet} > %g GeV, |#eta_{jet}| < 2", jtpt_binsVector[ibin_pt - 1]));
            t.DrawLatex(0.20, 0.665, unfoldBayes ? "Bayesian unfolding" : "Matrix inversion");
            t.SetTextSize(0.034); t.DrawLatex(0.66, 0.945, "pp #sqrt{s} = 5.36 TeV");
            c_br->RedrawAxis();
            const TString stem_br = outFolder + "systematics_breakdown_" + vtag + "_" + label;
            c_br->Print(stem_br + ".pdf");
            c_br->Print(stem_br + ".png");
            std::cout << "Wrote " << stem_br << ".{pdf,png}" << std::endl;
        }
        }   // detailed / grouped
    }

    // ---- Plot: the final measurement against particle-level MC ----------------------
    // The money plot. Red = the measurement (unfolded data, stat bars + systematic box),
    // blue = particle-level Pythia, purple = particle-level Herwig when that unfolding
    // exists. h_mc_true_2D is hgenjet_2b_all in data mode: the gen reference the corrected
    // data is supposed to land on, written by apply_unfolding_2d.C already restricted to
    // the drawn dr range and unit-normalised, exactly like the result next to it.
    //
    // The gen curves are only RE-normalised here, never weighted: the scale factors in w
    // correct data for detector effects and have no meaning at particle level.
    {
        auto loadGen = [&](const TString &file, const TString &newname) -> TH1D * {
            TFile *f = TFile::Open(file);
            if (!f || f->IsZombie()) { if (f) delete f; return nullptr; }
            TH1D *h = dynamic_cast<TH1D *>(f->Get("h_mc_true_2D"));
            if (!h) { f->Close(); return nullptr; }
            TH1D *out = (TH1D *) h->Clone(newname);
            out->SetDirectory(nullptr);
            f->Close();
            if (out->GetNbinsX() != nbins) {
                std::cerr << "Gen MC '" << newname << "': " << out->GetNbinsX() << " bins vs "
                          << nbins << " in the result -- skipped" << std::endl;
                delete out;
                return nullptr;
            }
            if (w.normalise_unit_area) {
                const double integral = out->Integral();
                if (integral > 0.) out->Scale(1. / integral, "width");
            }
            return out;
        };

        // Draw the OTHER generator's particle level alongside the nominal one? Off: the
        // plot is the measurement against one theory curve. Flip to true to put Herwig
        // back next to Pythia -- nothing else has to change, every step below is guarded
        // on h_gen_alt being non-null.
        const bool show_alt_gen = false;

        TH1D *h_gen_nom = loadGen(nominalFile, "h_gen_" + GENERATOR);
        TH1D *h_gen_alt = show_alt_gen
                        ? loadGen(resultFile(SAMPLE, other_gen, test_mode, unfoldBayes,
                                             TF_GENERATOR, false, "nominal", SFUPART_NOM,
                                             EEC_WEIGHT_OFF, OBSERVABLE),
                                  "h_gen_" + other_gen)
                        : nullptr;

        if (!h_gen_nom) {
            std::cerr << "\nNo particle-level MC in " << nominalFile
                      << " (h_mc_true_2D) -- skipping the data-vs-gen plot.\n"
                      << "   That histogram's Write() was disabled until 2026-09-10, so a "
                      << "result unfolded before then does not carry it.\n"
                      << "   Re-run apply_unfolding_2d.C(\"" << SAMPLE << "\",\"" << GENERATOR
                      << "\"," << test_mode << "," << (unfoldBayes ? "true" : "false")
                      << ",false,\"" << TF_GENERATOR << "\",false,\"nominal\",\"" << SFUPART_NOM
                      << "\"," << (EEC_WEIGHT_OFF ? "true" : "false") << ",\"" << OBSERVABLE
                      << "\") to get it." << std::endl;
        } else {
            const Color_t col_data = (Color_t) TColor::GetColor("#C44E52");  // measurement
            const Color_t col_pl   = (Color_t) TColor::GetColor("#4C72B0");  // particle-level MC
            const Color_t col_pl2  = (Color_t) TColor::GetColor("#8C6BB1");  // second PL MC
            const Color_t col_box  = blendWithWhite("#C44E52", 0.2);
            gStyle->SetOptStat(0);
            gStyle->SetOptTitle(0);
            gStyle->SetLegendBorderSize(0);
            gStyle->SetLegendFillColor(0);

            const Float_t font_scale  = 1200. / 800.;
            const Style_t font_code   = 43;
            const Float_t label_size  = 15. * font_scale;
            const Float_t title_size  = 15. * font_scale;
            const Float_t legend_size = 15. * font_scale;

            TCanvas *c_gen = new TCanvas("c_data_vs_gen", "", 800, 800);
            TPad *q_top = new TPad("q_top", "", 0., 0.3, 1., 1.);
            TPad *q_bot = new TPad("q_bot", "", 0., 0.,  1., 0.3);
            for (TPad *q : {q_top, q_bot}) { q->SetTicks(1, 0); q->SetFillColor(0); }
            q_top->SetMargin(0.1, 0.1, 0.0,  0.1);
            q_bot->SetMargin(0.1, 0.1, 0.23, 0.0);
            c_gen->cd();
            q_top->Draw(); q_bot->Draw();

            // ---- top: EEC ----
            q_top->cd();
            TH1D *h_meas = (TH1D *) h_nom->Clone("h_meas_draw");
            h_meas->SetDirectory(nullptr);
            double ymax = 0.;
            for (int i = 1; i <= nbins; ++i)
                ymax = std::max(ymax, h_nom->GetBinContent(i) + h_nom->GetBinError(i)
                                    + h_syst_up->GetBinContent(i));
            for (TH1D *h : {h_gen_nom, h_gen_alt})
                if (h) for (int i = 1; i <= nbins; ++i)
                    ymax = std::max(ymax, h->GetBinContent(i) + h->GetBinError(i));
            h_meas->SetTitle("");
            // Headroom: the legend sits upper right, the systematics tag upper left, and
            // both have to clear the peak. 1.75 is what fits the tag's rows above it.
            h_meas->GetYaxis()->SetRangeUser(0., ymax * 1.75);
            h_meas->GetYaxis()->SetTitle(obs_title);
            h_meas->GetYaxis()->CenterTitle(true);
            h_meas->GetYaxis()->SetTitleFont(font_code); h_meas->GetYaxis()->SetTitleSize(title_size);
            h_meas->GetYaxis()->SetTitleOffset(1.5);
            h_meas->GetYaxis()->SetLabelFont(font_code); h_meas->GetYaxis()->SetLabelSize(label_size);
            h_meas->GetXaxis()->SetTitleSize(0);
            h_meas->GetXaxis()->SetLabelSize(0);
            h_meas->SetLineColor(col_data);   h_meas->SetMarkerColor(col_data);
            h_meas->SetMarkerStyle(kFullCircle); h_meas->SetMarkerSize(1);
            h_meas->SetLineWidth(2);
            fixXRange(h_meas);
            h_meas->Draw("AXIS");

            // Systematic as a box on each point, drawn under the MC so the curves stay read
            g_syst->SetFillColor(col_box);
            g_syst->SetLineWidth(0);
            g_syst->Draw("2 same");

            // Upper RIGHT, unlike the systematics plot: the EEC peaks at small dr and the
            // left half of this pad is where the curves and the biggest systematic box are.
            // Sized by entry count -- a TLegend spreads its rows over the whole box, so a
            // fixed height leaves three entries floating with too much air between them.
            const int n_leg = 3 + (h_gen_alt ? 1 : 0);
            TLegend *leg_g = new TLegend(0.55, 0.87 - 0.0625 * n_leg, 0.93, 0.87);
            leg_g->SetFillStyle(0);
            leg_g->SetBorderSize(0);
            leg_g->SetMargin(0.15);
            leg_g->SetTextFont(font_code);
            leg_g->SetTextSize(legend_size);
            leg_g->AddEntry(h_meas, "Data, unfolded (stat)", "pe1");
            leg_g->AddEntry(g_syst, "Systematic uncertainty", "f");

            h_gen_nom->SetLineColor(col_pl); h_gen_nom->SetMarkerColor(col_pl);
            h_gen_nom->SetLineWidth(2); h_gen_nom->SetLineStyle(1);
            h_gen_nom->Draw("HIST same");
            leg_g->AddEntry(h_gen_nom, prettyGen(GENERATOR) + " (particle level)", "l");
            if (h_gen_alt) {
                h_gen_alt->SetLineColor(col_pl2); h_gen_alt->SetMarkerColor(col_pl2);
                h_gen_alt->SetLineWidth(2); h_gen_alt->SetLineStyle(2);
                h_gen_alt->Draw("HIST same");
                leg_g->AddEntry(h_gen_alt, prettyGen(other_gen) + " (particle level)", "l");
            }
            h_meas->Draw("PE X0 same");
            leg_g->Draw();

            TLatex cms_g;
            cms_g.SetNDC();
            cms_g.SetTextFont(62); cms_g.SetTextSize(0.045); cms_g.DrawLatex(0.10, 0.945, "CMS");
            cms_g.SetTextFont(52); cms_g.SetTextSize(0.037); cms_g.DrawLatex(0.205, 0.945, "Work In Progress");
            cms_g.SetTextFont(42); cms_g.SetTextSize(0.037); cms_g.DrawLatex(0.66, 0.945, "pp #sqrt{s} = 5.36 TeV");

            // ---- "systematics applied" tag, upper left ----
            // Listed from syst_tag_labels, which is built from the entries that actually
            // entered the band -- so a source whose result file was missing and got skipped
            // cannot be claimed here.
            if (!syst_tag_labels.empty()) {
                TLatex tag;
                tag.SetNDC();
                tag.SetTextFont(42);
                tag.SetTextSize(0.030);
                tag.SetTextAlign(13);           // left, top
                double y_tag = 0.855;
                tag.SetTextFont(62);
                tag.DrawLatex(0.145, y_tag, "Systematic uncertainties");
                tag.SetTextFont(42);
                for (const TString &t : syst_tag_labels) {
                    y_tag -= 0.042;
                    tag.DrawLatex(0.160, y_tag, "#bullet " + t);
                }
            }

            // ---- bottom: MC / data ----
            // Ratio to the MEASUREMENT, so the data uncertainty shows as a band around 1 and
            // each MC reads directly against it. Each ratio keeps its numerator's colour.
            q_bot->cd();
            TH1D *q_frame = (TH1D *) h_nom->Clone("h_gen_ratio_frame");
            q_frame->SetDirectory(nullptr);
            q_frame->Reset();
            double qmax = 0.;
            for (int i = 1; i <= nbins; ++i) {
                const double v = h_nom->GetBinContent(i);
                if (v == 0.) continue;
                for (TH1D *h : {h_gen_nom, h_gen_alt})
                    if (h) qmax = std::max(qmax, std::fabs(h->GetBinContent(i) / v - 1.));
                qmax = std::max(qmax,
                                (std::max(h_syst_up->GetBinContent(i), h_syst_dn->GetBinContent(i))
                                 + h_nom->GetBinError(i)) / std::fabs(v));
            }
            if (qmax <= 0.) qmax = 0.05;
            q_frame->SetTitle("");
            // 1.7, not 1.5: with one MC curve qmax is set by the bin-1 systematic alone and
            // 1.5 puts the frame edge just under a tick, so ROOT draws the top label half
            // outside the pad. The extra headroom keeps every label inside.
            q_frame->GetYaxis()->SetRangeUser(1. - 1.7 * qmax, 1. + 1.7 * qmax);
            q_frame->GetXaxis()->SetTitle(obs_axis);
            q_frame->GetYaxis()->SetTitle("MC / Data");
            q_frame->GetXaxis()->CenterTitle(true);
            q_frame->GetYaxis()->CenterTitle(true);
            q_frame->GetYaxis()->SetNdivisions(505);
            for (TAxis *ax : { q_frame->GetXaxis(), q_frame->GetYaxis() }) {
                ax->SetTitleFont(font_code); ax->SetTitleSize(title_size);
                ax->SetLabelFont(font_code); ax->SetLabelSize(label_size);
            }
            q_frame->GetXaxis()->SetTitleOffset(1.1);
            q_frame->GetYaxis()->SetTitleOffset(1.5);
            fixXRange(q_frame);
            q_frame->Draw("AXIS");

            // Data uncertainty around 1: the systematic as the wide light box, the
            // statistical error as a narrower, more saturated box drawn on top of it. Both
            // solid -- a hatched fill at this width renders as a few stray ticks and reads
            // as nothing at all.
            TGraphAsymmErrors *g_dsyst = new TGraphAsymmErrors(nbins);
            TGraphAsymmErrors *g_dstat = new TGraphAsymmErrors(nbins);
            g_dsyst->SetName("g_data_syst_rel");
            g_dstat->SetName("g_data_stat_rel");
            for (int i = 1; i <= nbins; ++i) {
                const double v = h_nom->GetBinContent(i);
                const double x = h_nom->GetBinCenter(i);
                const double hw = 0.5 * h_nom->GetBinWidth(i);
                const double st = (v != 0.) ? h_nom->GetBinError(i) / std::fabs(v) : 0.;
                g_dsyst->SetPoint(i - 1, x, 1.);
                g_dstat->SetPoint(i - 1, x, 1.);
                // The systematic box is asymmetric about 1; the statistical one is not.
                g_dsyst->SetPointError(i - 1, hw, hw,
                                       (v != 0.) ? h_syst_dn->GetBinContent(i) / std::fabs(v) : 0.,
                                       (v != 0.) ? h_syst_up->GetBinContent(i) / std::fabs(v) : 0.);
                g_dstat->SetPointError(i - 1, hw * 0.30, hw * 0.30, st, st);
            }
            g_dsyst->SetFillColor(col_box);                        g_dsyst->SetLineWidth(0);
            g_dsyst->Draw("2 same");
            g_dstat->SetFillColor(blendWithWhite("#C44E52", 0.5)); g_dstat->SetLineWidth(0);
            g_dstat->Draw("2 same");

            for (TH1D *h : {h_gen_nom, h_gen_alt}) {
                if (!h) continue;
                TH1D *h_q = (TH1D *) h->Clone(TString(h->GetName()) + "_ratio");
                h_q->SetDirectory(nullptr);
                for (int i = 1; i <= nbins; ++i) {
                    const double v = h_nom->GetBinContent(i);
                    h_q->SetBinContent(i, (v != 0.) ? h->GetBinContent(i) / v : 0.);
                    h_q->SetBinError(i, 0.);
                }
                h_q->SetLineColor(h->GetLineColor());
                h_q->SetLineWidth(2);
                h_q->SetLineStyle(h->GetLineStyle());
                h_q->Draw("HIST same");
            }

            TLine *l1 = new TLine(q_frame->GetXaxis()->GetXmin(), 1.,
                                  q_frame->GetXaxis()->GetXmax(), 1.);
            l1->SetLineStyle(2);
            l1->Draw();
            q_bot->RedrawAxis();

            const TString gen_stem = outFolder
                                   + "data_vs_gen_" + label;
            c_gen->Print(gen_stem + ".pdf");
            c_gen->Print(gen_stem + ".png");
            std::cout << "\nWrote " << gen_stem << ".{pdf,png}" << std::endl;

            // Kept for the write block below.
            h_gen_curves.push_back(h_gen_nom);
            gen_curve_names.push_back("h_gen_" + GENERATOR);
            if (h_gen_alt) {
                h_gen_curves.push_back(h_gen_alt);
                gen_curve_names.push_back("h_gen_" + other_gen);
            }
        }
    }

    // ---- Plot: the UParT SF, before and after, and nothing else ----------------------
    // The systematics plot carries the before-SF curve too, but buried among every other
    // variation. This one shows the correction on its own: the two curves and their ratio,
    // which is what you look at to judge the SF itself.
    //
    // Driven by the display-only "sfupart_off" entry, so it appears only when that entry is
    // booked and its unfolding exists. Comment the entry out and this plot stops being
    // produced, along with the curve on the systematics plot.
    {
        int k_off = -1;
        for (size_t k = 0; k < shift_names.size(); ++k)
            if (shift_names[k] == "sfupart_off") k_off = (int) k;

        if (k_off >= 0) {
            TH1D *h_before = h_vars[k_off];   // no SF applied
            TH1D *h_after  = h_nom;           // the nominal: SF applied

            // Particle-level MC for reference, so the plot answers not just "how big is the
            // SF" but "does it move the data toward the theory". Loaded by the data-vs-gen
            // block above and already unit-area normalised, the same as the two data
            // curves. Null when that block did not run, and everything below is guarded.
            TH1D *h_gen_ref = nullptr;
            for (size_t g = 0; g < h_gen_curves.size(); ++g)
                if (gen_curve_names[g] == "h_gen_" + GENERATOR) h_gen_ref = h_gen_curves[g];

            // Palette: red is the measurement, and the result WITH the SF is the
            // measurement. Before-SF is an intermediate correction stage -> green.
            // Particle-level MC is blue, as everywhere else.
            const Color_t col_after  = (Color_t) TColor::GetColor("#C44E52");
            const Color_t col_before = (Color_t) TColor::GetColor("#4F8F52");
            const Color_t col_gen    = (Color_t) TColor::GetColor("#4C72B0");
            gStyle->SetOptStat(0);
            gStyle->SetOptTitle(0);
            gStyle->SetLegendBorderSize(0);
            gStyle->SetLegendFillColor(0);

            const Float_t font_scale  = 1200. / 800.;
            const Style_t font_code   = 43;
            const Float_t label_size  = 15. * font_scale;
            const Float_t title_size  = 15. * font_scale;
            const Float_t legend_size = 15. * font_scale;

            TCanvas *c_sf = new TCanvas("c_sfupart_before_after", "", 800, 800);
            TPad *s_top = new TPad("s_top", "", 0., 0.35, 1., 1.);
            TPad *s_bot = new TPad("s_bot", "", 0., 0.,   1., 0.35);
            for (TPad *q : {s_top, s_bot}) { q->SetTicks(1, 0); q->SetFillColor(0); }
            s_top->SetMargin(0.1, 0.1, 0.0,  0.1);
            s_bot->SetMargin(0.1, 0.1, 0.20, 0.0);
            c_sf->cd();
            s_top->Draw(); s_bot->Draw();

            // ---- top: the two EEC curves ----
            s_top->cd();
            TH1D *h_a = (TH1D *) h_after->Clone("h_sfupart_after_draw");
            TH1D *h_b = (TH1D *) h_before->Clone("h_sfupart_before_draw");
            for (TH1D *h : {h_a, h_b}) h->SetDirectory(nullptr);

            double ymax = 0.;
            for (int i = 1; i <= nbins; ++i)
                ymax = std::max(ymax, std::max(h_a->GetBinContent(i) + h_a->GetBinError(i),
                                               h_b->GetBinContent(i) + h_b->GetBinError(i)));
            if (h_gen_ref)
                for (int i = 1; i <= nbins; ++i)
                    ymax = std::max(ymax, h_gen_ref->GetBinContent(i));
            h_a->SetTitle("");
            h_a->GetYaxis()->SetRangeUser(0., ymax * 1.45);
            h_a->GetYaxis()->SetTitle(obs_title);
            h_a->GetYaxis()->CenterTitle(true);
            h_a->GetYaxis()->SetTitleFont(font_code); h_a->GetYaxis()->SetTitleSize(title_size);
            h_a->GetYaxis()->SetTitleOffset(1.5);
            h_a->GetYaxis()->SetLabelFont(font_code); h_a->GetYaxis()->SetLabelSize(label_size);
            h_a->GetXaxis()->SetTitleSize(0);
            h_a->GetXaxis()->SetLabelSize(0);

            h_b->SetLineColor(col_before); h_b->SetMarkerColor(col_before);
            h_b->SetLineWidth(2); h_b->SetLineStyle(2);
            h_a->SetLineColor(col_after);  h_a->SetMarkerColor(col_after);
            h_a->SetMarkerStyle(kFullCircle); h_a->SetMarkerSize(1); h_a->SetLineWidth(2);

            TH1D *h_g = nullptr;
            if (h_gen_ref) {
                h_g = (TH1D *) h_gen_ref->Clone("h_sfupart_gen_draw");
                h_g->SetDirectory(nullptr);
                h_g->SetLineColor(col_gen); h_g->SetMarkerColor(col_gen);
                h_g->SetLineWidth(2); h_g->SetLineStyle(1);
            }

            h_a->Draw("AXIS");
            if (h_g) h_g->Draw("HIST same");
            h_b->Draw("HIST same");
            h_a->Draw("PE X0 same");

            // Sized by row count so two and three entries both sit tight.
            const int n_sf_rows = 2 + (h_g ? 1 : 0);
            TLegend *leg_sf = new TLegend(0.52, 0.85 - 0.055 * n_sf_rows, 0.93, 0.85);
            leg_sf->SetFillStyle(0);
            leg_sf->SetBorderSize(0);
            leg_sf->SetMargin(0.15);
            leg_sf->SetTextFont(font_code);
            leg_sf->SetTextSize(legend_size);
            leg_sf->AddEntry(h_b, "Before UParT SF", "l");
            leg_sf->AddEntry(h_a, "After UParT SF",  "pe1");
            if (h_g) leg_sf->AddEntry(h_g, prettyGen(GENERATOR) + " (particle level)", "l");
            leg_sf->Draw();

            TLatex cms_sf;
            cms_sf.SetNDC();
            cms_sf.SetTextFont(62); cms_sf.SetTextSize(0.045); cms_sf.DrawLatex(0.10, 0.945, "CMS");
            cms_sf.SetTextFont(52); cms_sf.SetTextSize(0.037); cms_sf.DrawLatex(0.205, 0.945, "Work In Progress");
            cms_sf.SetTextFont(42); cms_sf.SetTextSize(0.037); cms_sf.DrawLatex(0.66, 0.945, "pp #sqrt{s} = 5.36 TeV");

            // ---- bottom: after / before ----
            // The effect of the SF, on its own. Not 1/SF: the unit-area renormalisation and
            // the migration through the unfolding both change it, which is the whole reason
            // the SF is applied at reco level rather than here.
            s_bot->cd();
            TH1D *h_r = (TH1D *) h_a->Clone("h_sfupart_ratio");
            h_r->SetDirectory(nullptr);
            double rlo = 1., rhi = 1.;
            for (int i = 1; i <= nbins; ++i) {
                const double b = h_b->GetBinContent(i);
                const double r = (b != 0.) ? h_a->GetBinContent(i) / b : 1.;
                h_r->SetBinContent(i, r);
                h_r->SetBinError(i, 0.);
                rlo = std::min(rlo, r);
                rhi = std::max(rhi, r);
            }
            const double pad = std::max(0.02, 0.25 * (rhi - rlo));
            h_r->SetTitle("");
            h_r->GetYaxis()->SetRangeUser(rlo - pad, rhi + pad);
            h_r->GetXaxis()->SetTitle(obs_axis);
            h_r->GetYaxis()->SetTitle("After / Before");
            h_r->GetXaxis()->CenterTitle(true);
            h_r->GetYaxis()->CenterTitle(true);
            h_r->GetYaxis()->SetNdivisions(505);
            for (TAxis *ax : { h_r->GetXaxis(), h_r->GetYaxis() }) {
                ax->SetTitleFont(font_code); ax->SetTitleSize(title_size);
                ax->SetLabelFont(font_code); ax->SetLabelSize(label_size);
            }
            h_r->GetXaxis()->SetTitleOffset(1.0);
            h_r->GetYaxis()->SetTitleOffset(1.5);
            h_r->SetLineColor(col_after);
            h_r->SetLineWidth(2);
            h_r->SetMarkerSize(0);
            h_r->Draw("HIST");

            TLine *l_one = new TLine(h_r->GetXaxis()->GetXmin(), 1.,
                                     h_r->GetXaxis()->GetXmax(), 1.);
            l_one->SetLineStyle(2);
            l_one->Draw();
            s_bot->RedrawAxis();

            const TString sfupart_stem = outFolder
                                  + "sfupart_before_after_" + label;
            c_sf->Print(sfupart_stem + ".pdf");
            c_sf->Print(sfupart_stem + ".png");
            std::cout << "\nWrote " << sfupart_stem << ".{pdf,png}" << std::endl;
        }
    }

    // ---- Plot: the two JP HF treatments, as two bands on the EEC -----------------------
    // The JP HF calibration can be used in two ways and the analysis must pick ONE. This
    // plot puts the resulting uncertainty band from each on the measured EEC, so the choice
    // is made by looking at the bands rather than at a column of numbers:
    //
    //   "SF variation"              (sfupart_jphf) the alternative SF CURVE is divided into
    //                               the data instead of the central one, and the difference
    //                               from the nominal EEC is SYMMETRISED -- one alternative
    //                               gives the size of the shift, not its direction.
    //   "SF uncertainty propagated" (sfupart_jpsyst_up/dn) the quoted per-bin SF UNCERTAINTY
    //                               is added to / subtracted from the central SF, coherently
    //                               in every dr bin, and the pair is combined signed: the
    //                               member that pushes the result up sets the upper side,
    //                               the one that pushes down sets the lower.
    //
    // Those two phrases are the legend text, so the plot and this comment cannot drift.
    //
    // Each band is therefore drawn exactly as that treatment would enter the total, so the
    // plot answers "which band do I get if I choose this one".
    //
    // The two carry the SAME per-bin magnitudes: h_SFb_dr_syst_jpcalib_hf is exactly
    // |h_SFb_dr_jpcalib_hf - h_SFb_dr_central| (checked 2026-09-14, all 8 bins). They differ
    // only in the SIGN PATTERN -- the SF variation keeps the measured one, which flips across
    // dr and so reshapes the EEC and survives the unit-area normalisation; the propagated
    // uncertainty applies one coherent sign, which is close to a pure normalisation and
    // largely cancels in it. That is why the bands come out different sizes, and why booking
    // both would count one measurement twice.
    //
    // Driven by the three entries being booked, like the before/after plot above: comment
    // them out and this plot stops being produced.
    //
    // Colours follow the house palette: red is the measurement, and the two bands are
    // secondary -- green for the SF variation, purple for the propagated uncertainty.
    {
        int k_curve = -1, k_up = -1, k_dn = -1;
        for (size_t k = 0; k < shift_names.size(); ++k) {
            if (shift_names[k] == "sfupart_jphf")      k_curve = (int) k;
            if (shift_names[k] == "sfupart_jpsyst_up") k_up    = (int) k;
            if (shift_names[k] == "sfupart_jpsyst_dn") k_dn    = (int) k;
        }

        if (k_curve >= 0 && k_up >= 0 && k_dn >= 0) {
            const Color_t col_meas = (Color_t) TColor::GetColor("#C44E52");  // measurement
            const Color_t col_1    = (Color_t) TColor::GetColor("#4F8F52");  // SF variation
            const Color_t col_2    = (Color_t) TColor::GetColor("#8C6BB1");  // SF unc. propagated
            gStyle->SetOptStat(0);
            gStyle->SetOptTitle(0);
            gStyle->SetLegendBorderSize(0);
            gStyle->SetLegendFillColor(0);

            const Float_t font_scale  = 1200. / 800.;
            const Style_t font_code   = 43;
            const Float_t label_size  = 15. * font_scale;
            const Float_t title_size  = 15. * font_scale;
            const Float_t legend_size = 14. * font_scale;

            TCanvas *c_jp = new TCanvas("c_jphf_treatments", "", 800, 850);
            TPad *j_top = new TPad("j_top", "", 0., 0.40, 1., 1.);
            TPad *j_bot = new TPad("j_bot", "", 0., 0.,   1., 0.40);
            for (TPad *q : {j_top, j_bot}) { q->SetTicks(1, 0); q->SetFillColor(0); }
            j_top->SetMargin(0.12, 0.05, 0.0,  0.08);
            j_bot->SetMargin(0.12, 0.05, 0.18, 0.0);
            c_jp->cd();
            j_top->Draw(); j_bot->Draw();

            // ---- the two bands, built the way each treatment is actually combined -------
            TGraphAsymmErrors *g1 = new TGraphAsymmErrors(nbins);
            TGraphAsymmErrors *g2 = new TGraphAsymmErrors(nbins);
            g1->SetName("g_jphf_band_curve");
            g2->SetName("g_jphf_band_syst");

            // The same two bands again, but centred on zero and in percent of the nominal,
            // for the lower pad.
            TGraphAsymmErrors *r1 = new TGraphAsymmErrors(nbins);
            TGraphAsymmErrors *r2 = new TGraphAsymmErrors(nbins);
            r1->SetName("g_jphf_rel_variation");
            r2->SetName("g_jphf_rel_uncprop");

            double ymax = 0., rmax = 0.;
            for (int i = 1; i <= nbins; ++i) {
                const double v  = h_nom->GetBinContent(i);
                const double hw = 0.5 * h_nom->GetBinWidth(i);
                const double x  = h_nom->GetBinCenter(i);

                // SF variation: one alternative -> symmetrised.
                const double a1 = std::fabs(h_deltas[k_curve]->GetBinContent(i));
                // SF uncertainty propagated: a pair -> signed, each side from the member
                // that pushes that way.
                // h_deltas is (nominal - variation), so the shift of the RESULT is minus it.
                const double u = -h_deltas[k_up]->GetBinContent(i);
                const double d = -h_deltas[k_dn]->GetBinContent(i);
                const double up2 = std::max(0., std::max(u, d));
                const double dn2 = std::max(0., std::max(-u, -d));

                g1->SetPoint(i - 1, x, v);
                g1->SetPointError(i - 1, hw, hw, a1, a1);
                g2->SetPoint(i - 1, x, v);
                g2->SetPointError(i - 1, hw, hw, dn2, up2);

                ymax = std::max(ymax, v + std::max(a1, up2) + h_nom->GetBinError(i));

                // As a RATIO to the nominal, so the lower pad reads 1.00 +/- the band and
                // needs no unit. Same two bands, just divided by the central value.
                const double p1  = (v != 0.) ? a1  / std::fabs(v) : 0.;
                const double p2u = (v != 0.) ? up2 / std::fabs(v) : 0.;
                const double p2d = (v != 0.) ? dn2 / std::fabs(v) : 0.;
                r1->SetPoint(i - 1, x, 1.);
                r1->SetPointError(i - 1, hw, hw, p1, p1);
                r2->SetPoint(i - 1, x, 1.);
                r2->SetPointError(i - 1, hw, hw, p2d, p2u);
                rmax = std::max(rmax, std::max(p1, std::max(p2u, p2d)));
            }

            // ---- Band styles, shared by both pads --------------------------------------
            // BOTH are drawn as bands, and they overlap in every bin, so they cannot both be
            // plain solid fills -- whichever is smaller would vanish underneath the other.
            // The variation is a light solid fill; the propagated uncertainty is HATCHED, so
            // the fill underneath shows through the gaps and the two stay readable wherever
            // one contains the other. Hatching rather than a second alpha because it survives
            // conversion to PDF and printing, which stacked transparencies do not reliably do.
            auto styleBand1 = [&](TGraphAsymmErrors *g) {
                g->SetFillColorAlpha(col_1, 0.45);
                g->SetFillStyle(1001);
                g->SetLineColor(col_1);
                g->SetLineWidth(1);
            };
            auto styleBand2 = [&](TGraphAsymmErrors *g) {
                g->SetFillColor(col_2);
                g->SetFillStyle(3354);
                g->SetLineColor(col_2);
                g->SetLineWidth(2);
            };
            styleBand1(g1); styleBand1(r1);
            styleBand2(g2); styleBand2(r2);

            // ---- top: the EEC with both bands ------------------------------------------
            j_top->cd();
            TH1D *h_frame = (TH1D *) h_nom->Clone("h_jphf_frame");
            h_frame->SetDirectory(nullptr);
            h_frame->SetTitle("");
            h_frame->GetYaxis()->SetRangeUser(0., ymax * 1.45);
            h_frame->GetYaxis()->SetTitle(obs_title);
            h_frame->GetYaxis()->CenterTitle(true);
            h_frame->GetYaxis()->SetTitleFont(font_code);
            h_frame->GetYaxis()->SetTitleSize(title_size);
            h_frame->GetYaxis()->SetTitleOffset(1.5);
            h_frame->GetYaxis()->SetLabelFont(font_code);
            h_frame->GetYaxis()->SetLabelSize(label_size);
            h_frame->GetXaxis()->SetTitleSize(0);
            h_frame->GetXaxis()->SetLabelSize(0);

            TH1D *h_pts = (TH1D *) h_nom->Clone("h_jphf_points");
            h_pts->SetDirectory(nullptr);
            h_pts->SetLineColor(col_meas); h_pts->SetMarkerColor(col_meas);
            h_pts->SetMarkerStyle(kFullCircle); h_pts->SetMarkerSize(1.0);
            h_pts->SetLineWidth(2);

            h_frame->Draw("AXIS");
            g1->Draw("2 same");
            g2->Draw("2 same");
            h_pts->Draw("PE X0 same");
            gPad->RedrawAxis();

            TLegend *leg_jp = new TLegend(0.42, 0.60, 0.93, 0.90);
            leg_jp->SetTextFont(font_code); leg_jp->SetTextSize(legend_size);
            leg_jp->AddEntry(h_pts, "Unfolded EEC (#pm stat)", "pe");
            leg_jp->AddEntry(g1, "SF variation", "f");
            leg_jp->AddEntry(g2, "SF uncertainty propagated", "f");
            leg_jp->Draw();

            // ---- bottom: the same two bands as a ratio to the nominal ------------------
            // The top pad shows the bands on the measurement; this one is where their SIZES
            // can actually be read off and compared bin by bin. Same two bands, same styles,
            // centred on 1.
            j_bot->cd();
            TH1D *h_rframe = (TH1D *) h_nom->Clone("h_jphf_rel_frame");
            h_rframe->SetDirectory(nullptr);
            h_rframe->Reset();
            h_rframe->SetTitle("");
            // 1.9, not 1.6: the largest band would otherwise run into the legend.
            h_rframe->GetYaxis()->SetRangeUser(1. - 1.9 * rmax, 1. + 1.9 * rmax);
            h_rframe->GetYaxis()->SetTitle("Ratio to nominal");
            h_rframe->GetXaxis()->SetTitle(obs_axis);
            h_rframe->GetXaxis()->CenterTitle(true);
            h_rframe->GetYaxis()->CenterTitle(true);
            h_rframe->GetYaxis()->SetNdivisions(505);
            for (TAxis *ax : { h_rframe->GetXaxis(), h_rframe->GetYaxis() }) {
                ax->SetTitleFont(font_code); ax->SetTitleSize(title_size);
                ax->SetLabelFont(font_code); ax->SetLabelSize(label_size);
            }
            h_rframe->GetXaxis()->SetTitleOffset(1.1);
            h_rframe->GetYaxis()->SetTitleOffset(1.5);

            h_rframe->Draw("AXIS");
            r1->Draw("2 same");
            r2->Draw("2 same");

            TLine *l_one = new TLine(h_rframe->GetXaxis()->GetXmin(), 1.,
                                     h_rframe->GetXaxis()->GetXmax(), 1.);
            l_one->SetLineColor(kGray + 1); l_one->SetLineStyle(2);
            l_one->Draw();
            j_bot->RedrawAxis();

            TLegend *leg_jp2 = new TLegend(0.16, 0.82, 0.88, 0.98);
            leg_jp2->SetNColumns(2);
            leg_jp2->SetTextFont(font_code); leg_jp2->SetTextSize(legend_size);
            leg_jp2->AddEntry(r1, "SF variation", "f");
            leg_jp2->AddEntry(r2, "SF uncertainty propagated", "f");
            leg_jp2->Draw();

            const TString jp_stem = outFolder
                                  + "jphf_treatments_" + label;
            c_jp->Print(jp_stem + ".pdf");
            c_jp->Print(jp_stem + ".png");
            std::cout << "\nWrote " << jp_stem << ".{pdf,png}" << std::endl;

            // The same numbers, flagging the bins where the choice actually matters.
            printf("\n  JP HF: SF variation vs SF uncertainty propagated, %% of the nominal\n");
            printf("  %3s %8s | %9s %8s | %9s %9s %8s\n",
                   "bin", "value", "var +/-", "rel%", "unc up", "unc dn", "rel%");
            printf("  ----------------+--------------------+------------------------------\n");
            for (int i = 1; i <= nbins; ++i) {
                const double v  = h_nom->GetBinContent(i);
                const double a1 = std::fabs(h_deltas[k_curve]->GetBinContent(i));
                const double u  = -h_deltas[k_up]->GetBinContent(i);
                const double d  = -h_deltas[k_dn]->GetBinContent(i);
                const double env = std::max(std::fabs(u), std::fabs(d));
                printf("  %3d %8.4f | %9.5f %7.2f%% | %+9.5f %+9.5f %7.2f%%%s\n",
                       i, v, a1, v ? 100. * a1 / v : 0., u, d, v ? 100. * env / v : 0.,
                       (env > 0. && a1 > 0. && (a1 / env > 2. || env / a1 > 2.))
                           ? "   <-- differ >2x" : "");
            }
            printf("\n  Same per-bin magnitudes, different sign pattern: the SF variation keeps\n"
                   "  the measured one and reshapes, the propagated uncertainty is coherent\n"
                   "  across dr and cancels in the unit-area normalisation. Pick one -- never\n"
                   "  book both.\n");
        }
    }

    // ---- Write ----
    const TString fout_name = outFolder
                            + "final_" + label + "_with_systematics.root";
    TFile fout(fout_name, "RECREATE");
    h_nom->Write("h_nominal");          // central values, stat errors
    h_syst_up->Write("h_syst_total_up");   // absolute systematic per bin, upward side
    h_syst_dn->Write("h_syst_total_down"); // ... and downward side
    h_syst->Write("h_syst_total");         // the larger of the two, for one-number consumers
    h_total->Write("h_result_stat_syst");  // central values, stat (+) the larger syst side
    g_syst->Write("g_syst_band");          // TGraphAsymmErrors: the band to draw
    for (size_t k = 0; k < h_shifts.size(); ++k) {
        h_shifts[k]->Write("h_syst_"  + shift_names[k]);   // symmetrised |delta|
        // NOTE the sign: the STORED convention is (nominal - variation), unchanged, because
        // plot_systematics_summary() reads these and flips the sign itself. The printed
        // table above shows the negative of this -- the shift of the RESULT -- so that its
        // sign matches which of syst+/syst- the source feeds. Same numbers, opposite sign.
        h_deltas[k]->Write("h_delta_" + shift_names[k]);   // signed nominal - variation
        h_vars[k]  ->Write("h_var_"   + shift_names[k]);   // the variation's own EEC curve
    }
    for (size_t k = 0; k < h_gen_curves.size(); ++k)
        h_gen_curves[k]->Write(gen_curve_names[k]);        // particle-level MC, unit area
    fout.Close();

    std::cout << "\nWrote " << fout_name << std::endl;

    // ---- What of all this is actually a result --------------------------------------
    // The results area accumulates one folder per flag combination, and most of them are
    // INPUTS to the band, not results. Print the short list, so the deliverables can be
    // found without reading the code or guessing from folder names.
    {
        const TString out = outFolder;
        printf("\n");
        printf("================================ FINAL OUTPUTS ================================\n");
        printf("All in %s\n\n", out.Data());
        printf("  %-46s %s\n", TString("final_" + label + "_with_systematics.root").Data(),
                                "the result: central values, stat, syst+/-");
        printf("  %-46s %s\n", TString("data_vs_gen_" + label + ".pdf").Data(),
                                "the measurement vs particle-level MC");
        printf("  %-46s %s\n", TString("systematics_curves_" + label + ".pdf").Data(),
                                "per-source breakdown");
        bool has_sf_cmp = false;
        for (const TString &n : shift_names) if (n == "sfupart_off") has_sf_cmp = true;
        if (has_sf_cmp)
            printf("  %-46s %s\n", TString("sfupart_before_after_" + label + ".pdf").Data(),
                                    "the UParT SF on its own (comparison)");
        bool has_jp_cmp = false, has_jp_up = false;
        for (const TString &n : shift_names) {
            if (n == "sfupart_jphf")      has_jp_cmp = true;
            if (n == "sfupart_jpsyst_up") has_jp_up  = true;
        }
        if (has_jp_cmp && has_jp_up)
            printf("  %-46s %s\n", TString("jphf_treatments_" + label + ".pdf").Data(),
                                    "JP HF: curve swap vs +/- syst (comparison)");
        printf("\nEvery OTHER unfolding_* folder is an INPUT to the band -- one systematic\n");
        printf("variation each -- not a result. The variation list above names them.\n");
        printf("===============================================================================\n");
    }
}

// ============================================================================
// The tracking-efficiency uncertainty: EEC vs yields, side by side
// ============================================================================
//
// The same 3% track drop, propagated through the same chain twice -- once with the
// (pt1*pt2)^n weight on (the EEC) and once with it off (the dN/dr yield). Both productions
// use the SAME seed, and TrkEffSyst hashes (seed, entry, track index), so the two runs drop
// exactly the same tracks. The difference between the two bands is therefore the EEC weight
// and nothing else, which is the only reason this comparison means anything.
//
// Reads the two final files apply_weights_and_systematics() already wrote:
//     unfolding_both_pythia_upartv2/final_..._with_systematics.root         (EEC)
//     unfolding_both_pythia_noeecw_upartv2/final_..._noeecw_..._root       (yields)
// so it needs no unfolding of its own. Run the macro for BOTH observables first.
//
// TOP:    the two unfolded distributions, each unit-area normalised (as the macro writes
//         them) so their shapes are directly comparable, each with its tracking band.
// BOTTOM: the same two bands as a ratio to their own nominal -- this is where the SIZE of
//         the uncertainty is read off, which is the actual question.
//
//   root -l -b -q -e '.L apply_weights_and_systematics.C' \
//                 -e 'plot_tracking_eec_vs_yield("both","pythia")'
// ============================================================================
void plot_tracking_eec_vs_yield(TString SAMPLE = "both", TString GENERATOR = "pythia",
                                int test_mode = 2, bool unfoldBayes = false,
                                TString TF_GENERATOR = "pythia")
{
    auto finalFile = [&](bool eecw_off) {
        return resultFolder(SAMPLE, GENERATOR, unfoldBayes, TF_GENERATOR, false, "nominal",
                            "nominal", eecw_off)
             + "final_" + resultLabel(SAMPLE, GENERATOR, test_mode, unfoldBayes, TF_GENERATOR,
                                      false, "nominal", "nominal", eecw_off)
             + "_with_systematics.root";
    };

    TFile *fe = TFile::Open(finalFile(false));   // EEC
    TFile *fy = TFile::Open(finalFile(true));    // yields
    if (!fe || fe->IsZombie() || !fy || fy->IsZombie()) {
        std::cerr << "ERROR: need BOTH final files. Missing one of:\n  " << finalFile(false)
                  << "\n  " << finalFile(true)
                  << "\n   -> run apply_weights_and_systematics(...) with EEC_WEIGHT_OFF false "
                  << "and again with true." << std::endl;
        return;
    }
    // h_syst_<name> is the SYMMETRISED |nominal - variation| that entered the band.
    TH1D *ne = (TH1D *) fe->Get("h_nominal");
    TH1D *se = (TH1D *) fe->Get("h_syst_tracking_eff");
    TH1D *ny = (TH1D *) fy->Get("h_nominal");
    TH1D *sy = (TH1D *) fy->Get("h_syst_tracking_eff");
    if (!ne || !se || !ny || !sy) {
        std::cerr << "ERROR: h_nominal / h_syst_tracking_eff missing -- was 'tracking_eff' "
                  << "booked and its unfolding present in both runs?" << std::endl;
        return;
    }
    const int nb = ne->GetNbinsX();
    if (ny->GetNbinsX() != nb) {
        std::cerr << "ERROR: the two results have different binning ("
                  << nb << " vs " << ny->GetNbinsX() << ")" << std::endl;
        return;
    }

    // Same idiom as the JP HF comparison above: two bands, one solid and one hatched, so
    // whichever is smaller stays visible where one contains the other.
    const Color_t col_eec = (Color_t) TColor::GetColor("#4F8F52");   // EEC
    const Color_t col_yld = (Color_t) TColor::GetColor("#8C6BB1");   // yields
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0); gStyle->SetLegendFillColor(0);

    const Float_t font_scale  = 1200. / 800.;
    const Style_t font_code   = 43;
    const Float_t label_size  = 15. * font_scale;
    const Float_t title_size  = 15. * font_scale;
    const Float_t legend_size = 14. * font_scale;

    TCanvas *c = new TCanvas("c_trk_eec_vs_yield", "", 800, 850);
    TPad *p_top = new TPad("k_top", "", 0., 0.40, 1., 1.);
    TPad *p_bot = new TPad("k_bot", "", 0., 0.,   1., 0.40);
    for (TPad *q : {p_top, p_bot}) { q->SetTicks(1, 0); q->SetFillColor(0); }
    p_top->SetMargin(0.13, 0.05, 0.0,  0.08);
    p_bot->SetMargin(0.13, 0.05, 0.18, 0.0);
    c->cd(); p_top->Draw(); p_bot->Draw();

    TGraphAsymmErrors *ge = new TGraphAsymmErrors(nb), *gy = new TGraphAsymmErrors(nb);
    TGraphAsymmErrors *re = new TGraphAsymmErrors(nb), *ry = new TGraphAsymmErrors(nb);
    double ymax = 0., rmax = 0.;
    for (int i = 1; i <= nb; ++i) {
        const double x = ne->GetBinCenter(i), hw = 0.5 * ne->GetBinWidth(i);
        const double ve = ne->GetBinContent(i), vy = ny->GetBinContent(i);
        const double ee = se->GetBinContent(i), ey = sy->GetBinContent(i);
        ge->SetPoint(i - 1, x, ve); ge->SetPointError(i - 1, hw, hw, ee, ee);
        gy->SetPoint(i - 1, x, vy); gy->SetPointError(i - 1, hw, hw, ey, ey);
        const double fe_ = (ve != 0.) ? ee / std::fabs(ve) : 0.;
        const double fy_ = (vy != 0.) ? ey / std::fabs(vy) : 0.;
        re->SetPoint(i - 1, x, 1.); re->SetPointError(i - 1, hw, hw, fe_, fe_);
        ry->SetPoint(i - 1, x, 1.); ry->SetPointError(i - 1, hw, hw, fy_, fy_);
        ymax = std::max(ymax, std::max(ve + ee, vy + ey));
        rmax = std::max(rmax, std::max(fe_, fy_));
    }
    for (TGraphAsymmErrors *g : {ge, re}) {
        g->SetFillColorAlpha(col_eec, 0.45); g->SetFillStyle(1001);
        g->SetLineColor(col_eec); g->SetLineWidth(1);
    }
    for (TGraphAsymmErrors *g : {gy, ry}) {
        g->SetFillColor(col_yld); g->SetFillStyle(3354);
        g->SetLineColor(col_yld); g->SetLineWidth(2);
    }

    // ---- top: both distributions with their tracking bands ----
    p_top->cd();
    TH1D *fr = (TH1D *) ne->Clone("h_trkcmp_frame"); fr->SetDirectory(nullptr); fr->Reset();
    fr->SetTitle("");
    fr->GetYaxis()->SetRangeUser(0., ymax * 1.45);
    fr->GetYaxis()->SetTitle("normalised distribution");
    fr->GetYaxis()->CenterTitle(true);
    fr->GetYaxis()->SetTitleFont(font_code); fr->GetYaxis()->SetTitleSize(title_size);
    fr->GetYaxis()->SetTitleOffset(1.6);
    fr->GetYaxis()->SetLabelFont(font_code); fr->GetYaxis()->SetLabelSize(label_size);
    fr->GetYaxis()->ChangeLabel(1, -1, 0.);
    fr->GetXaxis()->SetTitleSize(0); fr->GetXaxis()->SetLabelSize(0);
    fr->Draw("AXIS");
    ge->Draw("2 same"); gy->Draw("2 same");
    gPad->RedrawAxis();

    TLegend *lg = new TLegend(0.45, 0.66, 0.93, 0.90);
    lg->SetTextFont(font_code); lg->SetTextSize(legend_size);
    lg->AddEntry(ge, "EEC, tracking band", "f");
    lg->AddEntry(gy, "Yield dN/d#Delta r, tracking band", "f");
    lg->Draw();

    // ---- bottom: the two bands as a ratio to their own nominal ----
    p_bot->cd();
    TH1D *rf = (TH1D *) ne->Clone("h_trkcmp_rframe"); rf->SetDirectory(nullptr); rf->Reset();
    rf->SetTitle("");
    rf->GetYaxis()->SetRangeUser(1. - 1.9 * rmax, 1. + 1.9 * rmax);
    rf->GetYaxis()->SetTitle("Ratio to nominal");
    rf->GetXaxis()->SetTitle("#Delta r");
    rf->GetXaxis()->CenterTitle(true); rf->GetYaxis()->CenterTitle(true);
    rf->GetYaxis()->SetNdivisions(505);
    for (TAxis *ax : { rf->GetXaxis(), rf->GetYaxis() }) {
        ax->SetTitleFont(font_code); ax->SetTitleSize(title_size);
        ax->SetLabelFont(font_code); ax->SetLabelSize(label_size);
    }
    rf->GetXaxis()->SetTitleOffset(1.1); rf->GetYaxis()->SetTitleOffset(1.6);
    rf->Draw("AXIS");
    re->Draw("2 same"); ry->Draw("2 same");
    TLine *l1 = new TLine(rf->GetXaxis()->GetXmin(), 1., rf->GetXaxis()->GetXmax(), 1.);
    l1->SetLineColor(kGray + 1); l1->SetLineStyle(2); l1->Draw();
    p_bot->RedrawAxis();

    TLegend *lg2 = new TLegend(0.16, 0.82, 0.88, 0.98);
    lg2->SetNColumns(2);
    lg2->SetTextFont(font_code); lg2->SetTextSize(legend_size);
    lg2->AddEntry(re, "EEC", "f");
    lg2->AddEntry(ry, "Yield", "f");
    lg2->Draw();

    const TString stem = resultFolder(SAMPLE, GENERATOR, unfoldBayes, TF_GENERATOR)
                       + "tracking_eec_vs_yield_" + SAMPLE + "_" + GENERATOR;
    c->Print(stem + ".pdf");
    c->Print(stem + ".png");
    std::cout << "\nWrote " << stem << ".{pdf,png}" << std::endl;

    printf("\n  Tracking-efficiency uncertainty, EEC vs yields (symmetrised, %% of nominal)\n");
    printf("  %3s %13s | %9s %8s | %9s %8s | %8s\n",
           "bin", "dr range", "EEC", "rel%", "yield", "rel%", "yld/EEC");
    printf("  -----------------+--------------------+--------------------+---------\n");
    for (int i = 1; i <= nb; ++i) {
        const double ve = ne->GetBinContent(i), vy = ny->GetBinContent(i);
        const double ee = se->GetBinContent(i), ey = sy->GetBinContent(i);
        const double pe = ve ? 100. * ee / std::fabs(ve) : 0.;
        const double py = vy ? 100. * ey / std::fabs(vy) : 0.;
        printf("  %3d [%5.2f,%5.2f] | %9.5f %7.2f%% | %9.5f %7.2f%% | %8.2f\n",
               i, ne->GetXaxis()->GetBinLowEdge(i), ne->GetXaxis()->GetBinUpEdge(i),
               ee, pe, ey, py, pe ? py / pe : 0.);
    }
    fe->Close(); fy->Close();
}


// ============================================================================
// Herwig vs Pythia8 summary plots -- moved in from plot_systematics_summary.C on
// 2026-09-25 to keep the number of files down. Unchanged, except that the colour blend uses
// blendWithWhite() from the top of this file (it was an identical copy, blendW()).
//
// Relative shift of every Herwig variation with respect to the Pythia8 nominal.
//
//   each correction / step swapped to Herwig ON ITS OWN   (thin lines)
//   the whole workflow run with Herwig                    (thick red)
//
// All curves are (variation - nominal) / nominal, so a flat zero means the swap does not
// move the normalised EEC at all.
//
// Inputs are the two final_*_with_systematics.root files written by
// apply_weights_and_systematics() above -- run it for the Pythia and the Herwig nominal first.
//
//   root -l -b -q -e '.L apply_weights_and_systematics.C' \
//                 -e 'plot_systematics_summary("both")'
// ============================================================================

// The two final files both plots read: the Pythia nominal and the all-Herwig run.
// test_mode 2 (data) -- a generator comparison of anything else is not a measurement.
inline TString finalBandFile(const TString &sample, const TString &generator,
                             const TString &tfGenerator, bool unfoldBayes)
{
    return resultFolder(sample, generator, unfoldBayes, tfGenerator)
         + "final_" + resultLabel(sample, generator, 2, unfoldBayes, tfGenerator)
         + "_with_systematics.root";
}

void plot_systematics_summary(TString SAMPLE = "both", bool unfoldBayes = false)
{
  // Both plots compare two full productions, so they belong with the production they were
  // made from: methodDir keeps the matrix-inversion pair out of the Bayesian one's way
  // instead of overwriting it under the same name.
  const TString R   = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/"
                    + methodDir(unfoldBayes);

  const TString f_pyt = finalBandFile(SAMPLE, "pythia", "pythia", unfoldBayes);
  const TString f_her = finalBandFile(SAMPLE, "herwig", "herwig", unfoldBayes);

  TFile *fp = TFile::Open(f_pyt);
  TFile *fh = TFile::Open(f_her);
  if (!fp || fp->IsZombie()) { std::cerr << "MISSING " << f_pyt << std::endl; return; }
  if (!fh || fh->IsZombie()) { std::cerr << "MISSING " << f_her << std::endl; return; }

  TH1D *nom    = (TH1D *) fp->Get("h_nominal");
  TH1D *nom_hw = (TH1D *) fh->Get("h_nominal");
  if (!nom || !nom_hw) { std::cerr << "ERROR: h_nominal not found" << std::endl; return; }

  // h_delta_* is (nominal - variation); flip the sign so every curve reads as the shift
  // the Herwig swap produces.
  struct Src { const char *hist; const char *label; const char *hex; int style; };
  std::vector<Src> srcs = {
    { "h_delta_template_fit",    "template fit",        "#4C72B0", 2 },
    { "h_delta_unfolding_model", "unfolding MC",        "#4F8F52", 3 },
    { "h_delta_eff_2sv_btag",    "2SV+UParT efficiency","#8C6BB1", 4 },
    { "h_delta_eec_weight",      "EEC weight",          "#C97430", 5 },
  };

  const int nb = nom->GetNbinsX();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TCanvas *c = new TCanvas("c_syst_summary", "", 800, 600);
  c->SetLeftMargin(0.13); c->SetRightMargin(0.04);
  c->SetTopMargin(0.09);  c->SetBottomMargin(0.14);
  c->SetTicks(1, 1);

  TH1D *frame = (TH1D *) nom->Clone("frame_summary");
  frame->SetDirectory(nullptr);
  frame->Reset();

  // Coherent: the whole workflow run with Herwig.
  TH1D *coh = (TH1D *) nom->Clone("h_rel_all_herwig");
  coh->SetDirectory(nullptr);
  double rmax = 0.;
  for (int i = 1; i <= nb; ++i) {
    const double v = nom->GetBinContent(i);
    const double r = (v != 0.) ? (nom_hw->GetBinContent(i) - v) / v : 0.;
    coh->SetBinContent(i, r); coh->SetBinError(i, 0.);
    rmax = std::max(rmax, std::fabs(r));
  }

  std::vector<TH1D *> curves;
  for (const Src &s : srcs) {
    TH1D *d = (TH1D *) fp->Get(s.hist);
    if (!d) { std::cerr << "WARNING: " << s.hist << " missing" << std::endl; continue; }
    TH1D *rel = (TH1D *) nom->Clone(Form("h_rel_%s", s.label));
    rel->SetDirectory(nullptr);
    for (int i = 1; i <= nb; ++i) {
      const double v = nom->GetBinContent(i);
      const double r = (v != 0.) ? -d->GetBinContent(i) / v : 0.; // flip sign
      rel->SetBinContent(i, r); rel->SetBinError(i, 0.);
      rmax = std::max(rmax, std::fabs(r));
    }
    rel->SetLineColor((Color_t) TColor::GetColor(s.hex));
    rel->SetLineWidth(2);
    rel->SetLineStyle(s.style);
    curves.push_back(rel);
  }

  frame->GetYaxis()->SetRangeUser(-1.45 * rmax, 1.75 * rmax);
  frame->GetXaxis()->SetTitle("gen #Delta r_{BB}");
  frame->GetYaxis()->SetTitle("(Herwig #minus Pythia8) / Pythia8");
  frame->GetXaxis()->SetTitleFont(43); frame->GetXaxis()->SetTitleSize(23);
  frame->GetYaxis()->SetTitleFont(43); frame->GetYaxis()->SetTitleSize(23);
  frame->GetXaxis()->SetLabelFont(43); frame->GetXaxis()->SetLabelSize(19);
  frame->GetYaxis()->SetLabelFont(43); frame->GetYaxis()->SetLabelSize(19);
  frame->GetXaxis()->SetTitleOffset(1.3);
  frame->GetYaxis()->SetTitleOffset(1.6);
  frame->GetYaxis()->SetNdivisions(505);
  frame->Draw("AXIS");

  TLine *l0 = new TLine(frame->GetXaxis()->GetXmin(), 0., frame->GetXaxis()->GetXmax(), 0.);
  l0->SetLineColor(kGray + 1); l0->SetLineStyle(2);
  l0->Draw();

  for (TH1D *h : curves) h->Draw("HIST same");

  coh->SetLineColor((Color_t) TColor::GetColor("#C44E52"));
  coh->SetLineWidth(4);
  coh->Draw("HIST same");

  TLegend *leg = new TLegend(0.16, 0.70, 0.62, 0.90);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->SetTextFont(43); leg->SetTextSize(19);
  leg->AddEntry(coh, "whole workflow (Herwig)", "l");
  for (size_t k = 0; k < curves.size() && k < srcs.size(); ++k)
    leg->AddEntry(curves[k], srcs[k].label, "l");
  leg->Draw();
  c->RedrawAxis();

  const TString stem = R + "systematics_summary_" + SAMPLE;
  c->Print(stem + ".pdf");
  c->Print(stem + ".png");

  printf("\n %-24s", "bin");
  for (int i = 1; i <= nb; ++i) printf(" %6d", i);
  printf("\n %-24s", "whole workflow [%]");
  for (int i = 1; i <= nb; ++i) printf(" %6.2f", 100 * coh->GetBinContent(i));
  for (size_t k = 0; k < curves.size(); ++k) {
    printf("\n %-24s", srcs[k].label);
    for (int i = 1; i <= nb; ++i) printf(" %6.2f", 100 * curves[k]->GetBinContent(i));
  }
  printf("\n\nWrote %s.{pdf,png}\n", stem.Data());
}


// =====================================================================================
// The final result from the two full analyses -- everything Pythia8, everything Herwig --
// with the generator systematic drawn as a symmetric band, +/- |Herwig - Pythia8|,
// around the Pythia8 nominal.
//
//   root -l -b -q -e '.L apply_weights_and_systematics.C' \
//                 -e 'plot_final_generator_band("both")'

void plot_final_generator_band(TString SAMPLE = "both", bool unfoldBayes = false)
{
  // Both plots compare two full productions, so they belong with the production they were
  // made from: methodDir keeps the matrix-inversion pair out of the Bayesian one's way
  // instead of overwriting it under the same name.
  const TString R   = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/"
                    + methodDir(unfoldBayes);

  const TString f_pyt = finalBandFile(SAMPLE, "pythia", "pythia", unfoldBayes);
  const TString f_her = finalBandFile(SAMPLE, "herwig", "herwig", unfoldBayes);

  TFile *fp = TFile::Open(f_pyt);
  TFile *fh = TFile::Open(f_her);
  if (!fp || fp->IsZombie()) { std::cerr << "MISSING " << f_pyt << std::endl; return; }
  if (!fh || fh->IsZombie()) { std::cerr << "MISSING " << f_her << std::endl; return; }

  TH1D *py = (TH1D *) fp->Get("h_nominal");
  TH1D *hw = (TH1D *) fh->Get("h_nominal");
  if (!py || !hw) { std::cerr << "ERROR: h_nominal missing" << std::endl; return; }
  py = (TH1D *) py->Clone("h_full_pythia"); py->SetDirectory(nullptr);
  hw = (TH1D *) hw->Clone("h_full_herwig"); hw->SetDirectory(nullptr);

  const int nb = py->GetNbinsX();

  // ---- identical style to the systematics-curves plot in apply_weights_and_systematics.C
  const Color_t col_nom  = (Color_t) TColor::GetColor("#C44E52"); // measurement red
  const Color_t col_alt  = (Color_t) TColor::GetColor("#4C72B0"); // blue
  const Color_t col_band = blendWithWhite("#C44E52", 0.30);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);

  const Float_t font_scale  = 1200. / 800.;
  const Style_t font_code   = 43;
  const Float_t label_size  = 15. * font_scale;
  const Float_t title_size  = 15. * font_scale;
  const Float_t legend_size = 15. * font_scale;

  TCanvas *c = new TCanvas("c_gen_band", "", 800, 800);
  TPad *p_top = new TPad("p_top_gb", "", 0., 0.3, 1., 1.);
  TPad *p_bot = new TPad("p_bot_gb", "", 0., 0.,  1., 0.3);
  for (TPad *p : {p_top, p_bot}) {
    p->SetTicks(1, 0);
    p->SetFillColor(0);
  }
  p_top->SetMargin(0.1, 0.1, 0.0,  0.1);   // left, right, bottom, top
  p_bot->SetMargin(0.1, 0.1, 0.23, 0.0);
  c->cd();
  p_top->Draw(); p_bot->Draw();

  // Band: +/- |Herwig - Pythia8| around the Pythia8 result, and its relative version.
  TGraphErrors *g_band     = new TGraphErrors(nb);
  TGraphErrors *g_band_rel = new TGraphErrors(nb);
  g_band->SetName("g_generator_band");
  g_band_rel->SetName("g_generator_band_rel");
  double ymax = 0., rmax = 0.;
  for (int i = 1; i <= nb; ++i) {
    const double v = py->GetBinContent(i);
    const double d = std::fabs(hw->GetBinContent(i) - v);
    g_band->SetPoint(i - 1, py->GetBinCenter(i), v);
    g_band->SetPointError(i - 1, 0.5 * py->GetBinWidth(i), d);
    g_band_rel->SetPoint(i - 1, py->GetBinCenter(i), 0.);
    g_band_rel->SetPointError(i - 1, 0.5 * py->GetBinWidth(i), (v != 0.) ? d / v : 0.);
    ymax = std::max(ymax, v + py->GetBinError(i) + d);
    ymax = std::max(ymax, hw->GetBinContent(i) + hw->GetBinError(i));
    if (v != 0.) rmax = std::max(rmax, d / v + hw->GetBinError(i) / v);
  }
  g_band->SetFillColor(col_band);     g_band->SetLineWidth(0);
  g_band_rel->SetFillColor(col_band); g_band_rel->SetLineWidth(0);

  // ---- top ----
  p_top->cd();
  py->SetTitle("");
  py->GetYaxis()->SetRangeUser(0., ymax * 1.75);
  py->GetYaxis()->SetTitle("EEC(#Delta r)");
  py->GetYaxis()->CenterTitle(true);
  py->GetYaxis()->SetTitleFont(font_code); py->GetYaxis()->SetTitleSize(title_size);
  py->GetYaxis()->SetTitleOffset(1.5);
  py->GetYaxis()->SetLabelFont(font_code); py->GetYaxis()->SetLabelSize(label_size);
  py->GetXaxis()->SetTitleSize(0);
  py->GetXaxis()->SetLabelSize(0);
  py->SetLineColor(col_nom); py->SetMarkerColor(col_nom);
  py->SetMarkerStyle(kFullCircle); py->SetMarkerSize(1); py->SetLineWidth(2);
  // Herwig as a dashed line, like the variation curves on the systematics plots.
  hw->SetLineColor(col_alt); hw->SetMarkerColor(col_alt);
  hw->SetLineStyle(2); hw->SetLineWidth(3);

  py->Draw("AXIS");
  g_band->Draw("2 same");
  hw->Draw("HIST same");
  py->Draw("PE X0 same");

  TLegend *leg = new TLegend(0.15, 0.55, 0.55, 0.80);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetMargin(0.15);
  leg->SetTextFont(font_code);
  leg->SetTextSize(legend_size);
  leg->AddEntry(py,     "full analysis, Pythia8", "pe1");
  leg->AddEntry(hw,     "full analysis, Herwig",  "l");
  leg->AddEntry(g_band, "generator uncertainty",  "f");
  leg->Draw();

  TLatex cms_bl;
  cms_bl.SetNDC();
  cms_bl.SetTextFont(62); cms_bl.SetTextSize(0.045); cms_bl.DrawLatex(0.10, 0.945, "CMS");
  cms_bl.SetTextFont(52); cms_bl.SetTextSize(0.037); cms_bl.DrawLatex(0.205, 0.945, "Work In Progress");
  cms_bl.SetTextFont(42); cms_bl.SetTextSize(0.037); cms_bl.DrawLatex(0.66, 0.945, "pp #sqrt{s} = 5.36 TeV");

  // ---- bottom: same "difference relative to the nominal" axis as the syst plot ----
  p_bot->cd();
  TH1D *h_frame = (TH1D *) py->Clone("h_gb_frame");
  h_frame->SetDirectory(nullptr);
  h_frame->Reset();
  if (rmax <= 0.) rmax = 0.05;
  h_frame->SetTitle("");
  h_frame->GetYaxis()->SetRangeUser(-1.6 * rmax, 1.6 * rmax);
  h_frame->GetXaxis()->SetTitle("#Delta r");
  h_frame->GetYaxis()->SetTitle("Herwig / Pythia8 - 1");
  h_frame->GetXaxis()->CenterTitle(true);
  h_frame->GetYaxis()->CenterTitle(true);
  h_frame->GetYaxis()->SetNdivisions(505);
  for (TAxis *ax : { h_frame->GetXaxis(), h_frame->GetYaxis() }) {
    ax->SetTitleFont(font_code);
    ax->SetTitleSize(title_size);
    ax->SetLabelFont(font_code);
    ax->SetLabelSize(label_size);
  }
  h_frame->GetXaxis()->SetTitleOffset(1.1);
  h_frame->GetYaxis()->SetTitleOffset(1.5);
  h_frame->Draw("AXIS");
  g_band_rel->Draw("2 same");

  TH1D *h_rel = (TH1D *) hw->Clone("h_rel_hw");
  h_rel->SetDirectory(nullptr);
  for (int i = 1; i <= nb; ++i) {
    const double v = py->GetBinContent(i);
    h_rel->SetBinContent(i, (v != 0.) ? hw->GetBinContent(i) / v - 1. : 0.);
    h_rel->SetBinError(i,   (v != 0.) ? hw->GetBinError(i) / v : 0.);
  }
  h_rel->SetLineColor(col_alt); h_rel->SetMarkerColor(col_alt);
  h_rel->SetLineStyle(2); h_rel->SetLineWidth(3);
  h_rel->Draw("HIST same");

  TLine *l0 = new TLine(h_frame->GetXaxis()->GetXmin(), 0.,
                        h_frame->GetXaxis()->GetXmax(), 0.);
  l0->SetLineStyle(2);
  l0->Draw();
  p_bot->RedrawAxis();

  const TString stem = R + "final_generator_band_" + SAMPLE;
  c->Print(stem + ".pdf");
  c->Print(stem + ".png");

  printf("\n bin | Pythia8 |  Herwig | band (+/-) |  rel\n");
  for (int i = 1; i <= nb; ++i) {
    const double v = py->GetBinContent(i), d = std::fabs(hw->GetBinContent(i) - v);
    printf(" %3d | %7.4f | %7.4f | %10.4f | %5.2f%%\n", i, v, hw->GetBinContent(i), d,
           v != 0. ? 100 * d / v : 0.);
  }
  printf("\nWrote %s.{pdf,png}\n", stem.Data());
}


// ============================================================================
// Unfolding method comparison: matrix inversion vs Bayesian (added 2026-10-01)
// ============================================================================
// Reads the two final band files apply_weights_and_systematics() writes for the nominal
// (one per method) and compares them. Both must have been built from the same variation
// list, each variation unfolded with that method. NITER is the iteration count the Bayesian
// productions were made with -- the path does not record it, so it is passed here for the
// labels. Chosen 2026-10-01 by the convergence criterion of apply_unfolding_2d.C scan_niter
// (one more iteration moves the result by < 0.5 sigma_stat): 8 for dr, 11 for B without JEC;
// 7 for dr, 11 for B with it.
//
//   plot 1  <obs>_unfolding_MI_vs_bayes                 the two results with their bands;
//                                                        lower pad Bayes/MI vs the MI band
//   plot 2  <obs>_unfolding_MI_vs_bayes_uncertainties   each source's signed shift and the
//                                                        stat error, in % of that method's
//                                                        own nominal, MI vs Bayes
//
//   root -l -b -q 'apply_weights_and_systematics.C+' -e 'plot_unfolding_method_comparison("dr")'
// or, interpreted:
//   root -l -b -q -e '.L apply_weights_and_systematics.C' -e 'plot_unfolding_method_comparison("B")'
void plot_unfolding_method_comparison(TString OBSERVABLE = "dr", int NITER = -1,
                                      TString SAMPLE = "both")
{
  if (NITER < 0) NITER = (OBSERVABLE == "dr") ? 14 : 29;   // refolding p >= 0.95 (data + MC stat), JEC inputs
  const bool    eecOff = nominalEecWeightOff(OBSERVABLE);
  const TString sfNom  = nominalSfupart(OBSERVABLE);
  auto bandFile = [&](bool bayes) -> TString {
    return resultFolder(SAMPLE, "pythia", bayes, "pythia", false, "nominal", sfNom, eecOff, OBSERVABLE)
         + "final_" + resultLabel(SAMPLE, "pythia", 2, bayes, "pythia", false, "nominal", sfNom,
                                  eecOff, OBSERVABLE)
         + "_with_systematics.root";
  };
  const TString f_mi = bandFile(false), f_by = bandFile(true);
  TFile *fm = TFile::Open(f_mi);
  TFile *fb = TFile::Open(f_by);
  if (!fm || fm->IsZombie()) { std::cerr << "MISSING " << f_mi << std::endl; return; }
  if (!fb || fb->IsZombie()) { std::cerr << "MISSING " << f_by << std::endl; return; }

  auto get = [](TFile *f, const char *n, const char *tag) -> TH1D * {
    TH1D *h = dynamic_cast<TH1D *>(f->Get(n));
    if (!h) return nullptr;
    h = (TH1D *) h->Clone(Form("%s_%s", n, tag));
    h->SetDirectory(nullptr);
    return h;
  };
  TH1D *nm = get(fm, "h_nominal", "mi"),         *nb = get(fb, "h_nominal", "by");
  TH1D *um = get(fm, "h_syst_total_up", "mi"),   *ub = get(fb, "h_syst_total_up", "by");
  TH1D *dm = get(fm, "h_syst_total_down", "mi"), *db = get(fb, "h_syst_total_down", "by");
  if (!nm || !nb || !um || !ub || !dm || !db) {
    std::cerr << "ERROR: h_nominal / h_syst_total_up / _down missing in a band file" << std::endl;
    return;
  }
  const int nbins = nm->GetNbinsX();
  if (nb->GetNbinsX() != nbins) { std::cerr << "ERROR: binning differs" << std::endl; return; }

  // The sources, as booked in the rho-era lists. A source missing from either file is skipped.
  struct Src { TString name, title; };
  const std::vector<Src> srcs = {
    { "template_fit",    "MC template modeling (Herwig fit)" },
    { "unfolding_model", "Detector response (Herwig unfolding)" },
    { "herwig_full",     "MC generator (full Herwig chain)" },
    { "tracking_eff",    "Tracking efficiency" },
    { "mistag_0B_up",   "Light-jet mistag (0B #times 2 / #times 0)" },
    // { "sv_sf_up",       "SV b scale factor (#pm1#sigma)" },   // the one-pair version, off since 2026-10-13
    { "sv_sf_stat_up",   "SV b SF: stat" },
    { "sv_sf_nomuon_up", "SV b SF: no muon" },
    { "sv_sf_bb_up",     "SV b SF: bb model" },
    { "jp_hf_up",        "JP HF swap" },
    { "qq_rate_up",      "qq rate" },
  };

  // ---- printout ----
  printf("\n  Unfolding method comparison, %s (%s): MI vs Bayesian\n", OBSERVABLE.Data(),
         eecOff ? "yields" : "EEC weighted");
  printf("  %3s | %9s %9s %6s | %7s %7s | %7s %7s | %7s %7s\n", "bin", "MI", "Bayes", "By/MI",
         "stat MI", "stat By", "syst+MI", "syst+By", "syst-MI", "syst-By");
  for (int i = 1; i <= nbins; ++i) {
    const double vm = nm->GetBinContent(i), vb = nb->GetBinContent(i);
    auto pc = [](double e, double v) { return v != 0. ? 100. * e / std::fabs(v) : 0.; };
    printf("  %3d | %9.4f %9.4f %6.3f | %6.1f%% %6.1f%% | %6.1f%% %6.1f%% | %6.1f%% %6.1f%%\n",
           i, vm, vb, vm != 0. ? vb / vm : 0.,
           pc(nm->GetBinError(i), vm), pc(nb->GetBinError(i), vb),
           pc(um->GetBinContent(i), vm), pc(ub->GetBinContent(i), vb),
           pc(dm->GetBinContent(i), vm), pc(db->GetBinContent(i), vb));
  }

  // ---- common style ----
  const ObsDef  obs       = obsByName(OBSERVABLE);
  const TString obs_sym   = obsSymbol(OBSERVABLE);
  const TString obs_axis  = (OBSERVABLE == "dr") ? "#Delta r" : obs.axis;
  const TString obs_title = eecOff ? ("dN/d" + obs_sym) : ("EEC(" + obs_sym + ")");
  const Color_t col_mi  = (Color_t) TColor::GetColor("#C44E52");   // red: the nominal result
  const Color_t col_by  = (Color_t) TColor::GetColor("#8C6BB1");   // purple: the alternative
  const Color_t box_mi  = blendWithWhite("#C44E52", 0.30);
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  gStyle->SetHatchesLineWidth(2);
  const Float_t font_scale = 1200. / 800.;
  const Style_t font_code  = 43;
  const Float_t lab = 15. * font_scale, tit = 15. * font_scale, leg_s = 14. * font_scale;
  auto styleAx = [&](TH1 *f) {
    for (TAxis *ax : { f->GetXaxis(), f->GetYaxis() }) {
      ax->CenterTitle(true);
      ax->SetTitleFont(font_code); ax->SetTitleSize(tit);
      ax->SetLabelFont(font_code); ax->SetLabelSize(lab);
    }
    if (OBSERVABLE != "dr") f->GetXaxis()->SetRange(1, f->GetNbinsX());
  };
  auto cmsLabel = [](double x, double y, double s) {
    TLatex t; t.SetNDC();
    t.SetTextFont(62); t.SetTextSize(s);         t.DrawLatex(x, y, "CMS");
    t.SetTextFont(52); t.SetTextSize(0.82 * s);  t.DrawLatex(x + 2.35 * s, y, "Work In Progress");
  };

  const TString outdir = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/unfolding_method_comparison"
                       + svsfTag() + trackBdtSfTag() + eraTag() + "/";
  gSystem->mkdir(outdir, true);
  const TString stem = outdir + OBSERVABLE + "_unfolding_MI_vs_bayes";

  // ================= plot 1: the two results =================
  {
    TCanvas *c = new TCanvas("c_meth", "", 800, 850);
    TPad *pt = new TPad("pt_meth", "", 0., 0.32, 1., 1.);
    TPad *pb = new TPad("pb_meth", "", 0., 0.,   1., 0.32);
    for (TPad *p : {pt, pb}) { p->SetTicks(1, 1); p->SetFillColor(0); }
    pt->SetMargin(0.13, 0.05, 0.0, 0.08);
    pb->SetMargin(0.13, 0.05, 0.26, 0.0);
    c->cd(); pt->Draw(); pb->Draw();

    // Bands: MI solid, Bayes hatched on top, so the smaller stays visible inside the larger.
    TGraphAsymmErrors *gm = new TGraphAsymmErrors(nbins), *gb = new TGraphAsymmErrors(nbins);
    TGraphErrors *pm = new TGraphErrors(nbins), *pbp = new TGraphErrors(nbins);
    double ymax = 0.;
    for (int i = 1; i <= nbins; ++i) {
      const double x = nm->GetBinCenter(i), hw = 0.5 * nm->GetBinWidth(i), sh = 0.15 * nm->GetBinWidth(i);
      gm->SetPoint(i - 1, x, nm->GetBinContent(i));
      gm->SetPointError(i - 1, hw, hw, dm->GetBinContent(i), um->GetBinContent(i));
      gb->SetPoint(i - 1, x, nb->GetBinContent(i));
      gb->SetPointError(i - 1, hw, hw, db->GetBinContent(i), ub->GetBinContent(i));
      // Points side by side within the bin so the two stat bars do not overlap.
      pm->SetPoint(i - 1, x - sh, nm->GetBinContent(i));  pm->SetPointError(i - 1, 0., nm->GetBinError(i));
      pbp->SetPoint(i - 1, x + sh, nb->GetBinContent(i)); pbp->SetPointError(i - 1, 0., nb->GetBinError(i));
      ymax = std::max({ymax, nm->GetBinContent(i) + um->GetBinContent(i),
                       nb->GetBinContent(i) + ub->GetBinContent(i)});
    }
    gm->SetFillColor(box_mi); gm->SetFillStyle(1001); gm->SetLineWidth(0);
    gb->SetFillColor(col_by); gb->SetFillStyle(3354); gb->SetLineColor(col_by); gb->SetLineWidth(0);
    pm->SetMarkerStyle(20);  pm->SetMarkerColor(col_mi);  pm->SetLineColor(col_mi);  pm->SetMarkerSize(1.1);
    pbp->SetMarkerStyle(21); pbp->SetMarkerColor(col_by); pbp->SetLineColor(col_by); pbp->SetMarkerSize(1.0);

    pt->cd();
    TH1D *fr = (TH1D *) nm->Clone("fr_meth"); fr->Reset(); fr->SetDirectory(nullptr);
    fr->SetMinimum(0.001); fr->SetMaximum(1.55 * ymax);
    fr->GetYaxis()->SetTitle(obs_title);
    styleAx(fr); fr->GetYaxis()->SetTitleOffset(1.5);
    fr->Draw("AXIS");
    gm->Draw("2 same"); gb->Draw("2 same");
    pm->Draw("P same"); pbp->Draw("P same");
    cmsLabel(0.13, 0.935, 0.05);
    TLatex rt; rt.SetNDC(); rt.SetTextFont(42); rt.SetTextSize(0.04); rt.SetTextAlign(31);
    rt.DrawLatex(0.95, 0.935, "pp #sqrt{s} = 5.36 TeV");
    TLegend *lg = new TLegend(0.50, 0.62, 0.94, 0.89);
    lg->SetTextFont(font_code); lg->SetTextSize(leg_s); lg->SetFillStyle(0); lg->SetBorderSize(0);
    lg->AddEntry(pm, "Matrix inversion (stat)", "pe");
    lg->AddEntry(gm, "Matrix inversion, syst", "f");
    lg->AddEntry(pbp, Form("Bayesian, %d it. (stat)", NITER), "pe");
    lg->AddEntry(gb, "Bayesian, syst", "f");
    lg->Draw();
    pt->RedrawAxis();

    // Lower: Bayes / MI central values, with the MI band (relative) around 1.
    pb->cd();
    TH1D *rf = (TH1D *) nm->Clone("rf_meth"); rf->Reset(); rf->SetDirectory(nullptr);
    TGraphAsymmErrors *rm = new TGraphAsymmErrors(nbins);
    TGraphErrors *rr = new TGraphErrors(nbins);
    double rdev = 0.;
    for (int i = 1; i <= nbins; ++i) {
      const double vm = nm->GetBinContent(i), vb = nb->GetBinContent(i);
      const double x = nm->GetBinCenter(i), hw = 0.5 * nm->GetBinWidth(i);
      rm->SetPoint(i - 1, x, 1.);
      rm->SetPointError(i - 1, hw, hw, vm ? dm->GetBinContent(i) / vm : 0., vm ? um->GetBinContent(i) / vm : 0.);
      const double r = vm ? vb / vm : 0.;
      rr->SetPoint(i - 1, x, r);
      rr->SetPointError(i - 1, 0., vm ? nb->GetBinError(i) / vm : 0.);
      rdev = std::max({rdev, std::fabs(r - 1.) + (vm ? nb->GetBinError(i) / vm : 0.),
                       vm ? um->GetBinContent(i) / vm : 0., vm ? dm->GetBinContent(i) / vm : 0.});
    }
    const double rh = std::min(0.9, std::max(0.15, 1.15 * rdev));
    rf->SetMinimum(1. - rh); rf->SetMaximum(1. + rh);
    rf->GetYaxis()->SetTitle("Bayes / MI");
    rf->GetXaxis()->SetTitle(obs_axis);
    styleAx(rf); rf->GetYaxis()->SetTitleOffset(1.5); rf->GetXaxis()->SetTitleOffset(1.1);
    rf->GetYaxis()->SetNdivisions(505);
    rf->Draw("AXIS");
    rm->SetFillColor(box_mi); rm->SetLineWidth(0); rm->Draw("2 same");
    TLine one(obs.bins[0], 1., obs.bins[obs.nbins], 1.); one.SetLineStyle(2); one.Draw();
    rr->SetMarkerStyle(21); rr->SetMarkerColor(col_by); rr->SetLineColor(col_by); rr->Draw("P same");
    pb->RedrawAxis();

    c->Print(stem + ".pdf"); c->Print(stem + ".png");
  }

  // ================= plot 2: uncertainty breakdown =================
  {
    TCanvas *c = new TCanvas("c_meth_unc", "", 1600, 800);
    c->Divide(4, 2, 0.002, 0.002);
    // Each panel: x = observable, y = shift of the result in % of that method's nominal.
    auto relShift = [&](TFile *f, TH1D *nom, const TString &name, const char *tag) -> TH1D * {
      TH1D *d = dynamic_cast<TH1D *>(f->Get("h_delta_" + name));
      if (!d) return nullptr;
      TH1D *r = (TH1D *) nom->Clone("rel_" + name + "_" + tag); r->Reset(); r->SetDirectory(nullptr);
      for (int i = 1; i <= nbins; ++i) {
        const double v = nom->GetBinContent(i);
        r->SetBinContent(i, v ? -100. * d->GetBinContent(i) / v : 0.);   // h_delta = nom - var
      }
      return r;
    };
    auto relErr = [&](TH1D *nom, TH1D *err, double sign, const char *nm_) -> TH1D * {
      TH1D *r = (TH1D *) nom->Clone(nm_); r->Reset(); r->SetDirectory(nullptr);
      for (int i = 1; i <= nbins; ++i) {
        const double v = nom->GetBinContent(i);
        const double e = err ? err->GetBinContent(i) : nom->GetBinError(i);
        r->SetBinContent(i, v ? sign * 100. * e / std::fabs(v) : 0.);
      }
      return r;
    };
    auto styleLine = [](TH1D *h, Color_t col, Style_t ls) {
      h->SetLineColor(col); h->SetLineWidth(3); h->SetLineStyle(ls); h->SetFillStyle(0); h->SetMarkerSize(0);
    };

    struct Panel { TString title; std::vector<TH1D *> mi, by; };
    std::vector<Panel> panels;
    panels.push_back({ "Statistical",
                       { relErr(nm, nullptr, +1., "st_mi_p"), relErr(nm, nullptr, -1., "st_mi_m") },
                       { relErr(nb, nullptr, +1., "st_by_p"), relErr(nb, nullptr, -1., "st_by_m") } });
    for (const Src &s : srcs) {
      Panel p; p.title = s.title;
      std::vector<TString> names = { s.name };
      if (s.name == "mistag_0B_up") names.push_back("mistag_0B_down");
      if (s.name == "sv_sf_up")     names.push_back("sv_sf_down");
      // the per-contribution pairs (2026-10-13): every "<x>_up" has its "<x>_down"
      if (s.name.EndsWith("_up") && s.name != "sv_sf_up" && s.name != "mistag_0B_up") {
        TString dn = s.name; dn.Replace(dn.Length() - 3, 3, "_down"); names.push_back(dn);
      }
      for (const TString &n : names) {
        if (TH1D *h = relShift(fm, nm, n, "mi")) p.mi.push_back(h);
        if (TH1D *h = relShift(fb, nb, n, "by")) p.by.push_back(h);
      }
      if (!p.mi.empty() || !p.by.empty()) panels.push_back(p);
    }
    panels.push_back({ "Total systematic (+ / #minus)",
                       { relErr(nm, um, +1., "tot_mi_p"), relErr(nm, dm, -1., "tot_mi_m") },
                       { relErr(nb, ub, +1., "tot_by_p"), relErr(nb, db, -1., "tot_by_m") } });

    // One shared y range across panels, so their sizes compare at a glance.
    double amax = 0.;
    for (Panel &p : panels)
      for (auto *v : {&p.mi, &p.by}) for (TH1D *h : *v)
        for (int i = 1; i <= nbins; ++i) amax = std::max(amax, std::fabs(h->GetBinContent(i)));
    amax = std::min(100., 1.2 * amax);

    // Up to 7 panels; the last pad of the 4x2 grid holds the legend.
    const size_t n_pan = std::min<size_t>(panels.size(), 7);
    for (size_t k = 0; k < n_pan; ++k) {
      c->cd(k + 1);
      gPad->SetMargin(0.17, 0.04, 0.15, 0.10); gPad->SetTicks(1, 1);
      TH1D *f = (TH1D *) nm->Clone(Form("fr_unc_%zu", k)); f->Reset(); f->SetDirectory(nullptr);
      f->SetMinimum(-amax); f->SetMaximum(amax);
      f->GetYaxis()->SetTitle("shift of the result [%]");
      f->GetXaxis()->SetTitle(obs_axis);
      styleAx(f);
      f->GetXaxis()->SetLabelSize(0.75 * lab); f->GetYaxis()->SetLabelSize(0.75 * lab);
      f->GetXaxis()->SetTitleSize(0.8 * tit);  f->GetYaxis()->SetTitleSize(0.8 * tit);
      f->GetYaxis()->SetTitleOffset(2.0); f->GetXaxis()->SetTitleOffset(1.3);
      f->Draw("AXIS");
      TLine z(obs.bins[0], 0., obs.bins[obs.nbins], 0.); z.SetLineStyle(2); z.DrawClone();
      for (TH1D *h : panels[k].mi) { styleLine(h, col_mi, kSolid);  h->Draw("hist same"); }
      for (TH1D *h : panels[k].by) { styleLine(h, col_by, kDashed); h->Draw("hist same"); }
      TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.055);
      t.DrawLatex(0.20, 0.92, panels[k].title);
    }
    // Last pad: legend + labels.
    c->cd(8);
    cmsLabel(0.08, 0.82, 0.08);
    TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.06);
    t.DrawLatex(0.08, 0.72, "pp #sqrt{s} = 5.36 TeV");
    t.DrawLatex(0.08, 0.64, obs_title + (eecOff ? ", yields" : ", EEC weighted"));
    TH1D *lm = (TH1D *) nm->Clone("leg_mi"); styleLine(lm, col_mi, kSolid);
    TH1D *lb = (TH1D *) nm->Clone("leg_by"); styleLine(lb, col_by, kDashed);
    TLegend *lg = new TLegend(0.06, 0.25, 0.95, 0.55);
    lg->SetTextFont(42); lg->SetTextSize(0.065); lg->SetBorderSize(0); lg->SetFillStyle(0);
    lg->AddEntry(lm, "Matrix inversion", "l");
    lg->AddEntry(lb, Form("Bayesian, %d iterations", NITER), "l");
    lg->Draw();
    t.SetTextSize(0.05);
    t.DrawLatex(0.08, 0.15, "Each in % of its own nominal");
    c->Print(stem + "_uncertainties.pdf"); c->Print(stem + "_uncertainties.png");
  }
  printf("\nWrote %s{,_uncertainties}.{pdf,png}\n", stem.Data());
}
