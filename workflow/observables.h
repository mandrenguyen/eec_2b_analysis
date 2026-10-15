#ifndef OBSERVABLES_H
#define OBSERVABLES_H

// observables.h -- the measured observables, and the axis abstraction that lets the chain
// fill more than one of them in a single pass over the data.
//
// WHAT AN "OBSERVABLE" IS HERE
// ----------------------------
// The chain measures a quantity built from the two aggregated B hadrons of a 2b jet, binned
// against jet pT, with m_2B as the template-fit variable:
//
//     templates    TH3D(m_2B, X, pT)
//     response     2D (X_reco, pT_reco) <- (X_gen, pT_gen)
//
// X is the observable. For the whole history of this analysis X was dR and nothing else, so
// "the second axis" and "dR" were the same thing and dR was written into every histogram
// name, loop bound and binning lookup by hand.
//
// X is now a parameter. ObsDef carries everything that differs between observables -- the
// binning, the overflow fold, and the suffix its histogram names carry -- so adding a third
// observable is a binning vector plus one ObsDef, not another copy of the chain.
//
// WHY THE VALUE IS NOT IN HERE
// ----------------------------
// ObsDef deliberately does NOT carry a function that computes X from the two B four-vectors.
// dR is computed by tTree::calc_dr, which returns Float_t and takes Float_t arguments, so it
// truncates to single precision at two points. Reimplementing that here to satisfy a uniform
// interface would change dR in its last bits, which would change every existing result for
// no reason. The caller computes the value -- three lines per observable, at the point where
// the four-vectors already exist -- and passes it in. The genericity that was actually worth
// having is in the booking, filling and writing of ~40 histograms per observable.
//
// ADDING AN OBSERVABLE
//   1. binning vector + fill cap in binning_histos_small.h
//   2. an ObsDef for it below, and a line in obsByName()
//   3. compute its value next to dr/B in create_files_for_template_fit.cpp, and add it
//      to BOTH obsList() and obsValues() there (a size check refuses to run otherwise)
//   4. its tag in result_paths.h: isKnownObservable, observableTag, obsProdTag, and the two
//      nominal* conventions (does it carry the UParT SF? is it an EEC or a yield?)
// Everything downstream of that iterates the ObsDef list.
//
// ⚠️ An observable only exists in MC produced AFTER the code that fills it. obsProdTag()
// says which production carries which observable, and a run pointed at an older one fails
// on a missing histogram rather than silently reading the wrong axis.
//
// WHAT IS MEASURED RIGHT NOW (2026-09-22)
// ---------------------------------------
// Two observables: dr (the EEC axis) and B (the momentum balance). The three b/non-b
// momentum fractions -- fnb, fb and lnfb -- are COMMENTED OUT, here and in
// binning_histos_small.h, result_paths.h and create_files_for_template_fit.cpp. They are
// left in place rather than deleted so that re-enabling one is four uncommented blocks.
//
// ⚠️ THE MOMENTUM BALANCE IS CALLED B, NOT z (renamed 2026-09-22 at Zoe's request, so the
// code matches the plots). The rename reached the histogram-name suffix ("_B") and the
// production tag, so MC written by the earlier "_z" code cannot be read by this code --
// step 1 has to be re-run. See obsProdTag() in result_paths.h.

#include "TString.h"
#include "binning_histos_small.h"   // guarded, so including this header twice is safe

#include <algorithm>
#include <cmath>

// ---- B: momentum balance between the two aggregated B hadrons -------------------------
// ONE definition, shared by the chain and by momentum_balance_mc_study(). Duplicating it
// is exactly the failure mode result_paths.h was written to end: two copies of a mapping
// drift, and the disagreement shows up as a plausible-looking number rather than an error.
//
// The namespace keeps its name: MomBalance says what the quantity is, and the rename was
// z -> B on the OBSERVABLE's name, not on the thing that computes it.
namespace MomBalance {

// [0.5, 1) by construction. Symmetric in its arguments, so the order the two B hadrons come
// in does not matter -- which is why reco (SV order) and gen (best-EEC-pair order) can be
// compared directly without a matching convention.
inline double value(double pt1, double pt2)
{
    if (!(pt1 > 0.) || !(pt2 > 0.)) return -1.;   // also catches the -999 ntuple sentinels
    return std::max(pt1, pt2) / (pt1 + pt2);
}

} // namespace MomBalance

// ============================================================================================
// DISABLED 2026-09-22: the three b/non-b momentum fractions.
//
// fnb, lnfb and fb are commented out -- the measurement is dr (EEC) and the momentum balance
// B only. Their value functions are kept verbatim below, together with the reasoning about
// what counts as "b" and "non-b" on each side, which is the part that would be expensive to
// reconstruct. Re-enabling one means uncommenting here, in binning_histos_small.h, in
// result_paths.h (isKnownObservable / observableTag / obsProdTag / the nominal* conventions)
// and in create_files_for_template_fit.cpp (obsList and obsValues, which size-check against
// each other), then re-running step 1.
// ============================================================================================
/*
// ---- fnb: how much of the jet is NOT in the two B hadrons ------------------------------
// Requested as 1/(1 + pT_b/pT_nonb), which is algebraically pT_nonb/(pT_b + pT_nonb) -- and
// that is the form used here, because it is still defined when pT_nonb = 0 (a jet every one
// of whose tracks went into the two Bs is fnb = 0, not a division by zero).
//
// Same functional family as z, which is also 1/(1 + ratio): z = pT_lead/(pT1 + pT2) is
// 1/(1 + pT_sublead/pT_lead). Both map an unbounded ratio onto [0, 1].
//
// WHAT COUNTS AS "b" AND "non-b" (decided with Zoe, 2026-09-17):
//   pT_b     scalar sum of the two aggregated B hadrons' pT -- reco_sv[0].Pt() +
//            reco_sv[1].Pt() at reco, the best-EEC gen pair at gen. Scalar, not the vector
//            sum, to match how z combines the same two objects.
//   pT_nonb  the tracks of the jet that did NOT go into those two B hadrons:
//            RECO  jet tracks with pT > 1 whose BDT score fails the b-track cut. Nothing
//                  else is left over -- the aggregation absorbs every no-SV track that
//                  passes the cut into the nearer vertex -- so this is exactly the
//                  complement of the B reconstruction.
//            GEN   jet refTracks with pT > 1 and refTrkSta < 100, i.e. not from a B decay.
//                  The exact analogue: PartialBsAggregation builds the gen Bs from status
//                  >= 100 with the same pT cut.
// No b-tagging BDT cut is applied on the non-b side: selecting the non-b component with a
// b-likeness score would bias it towards b-like tracks, and gen has no BDT to match anyway.
namespace NonBFraction {

inline double value(double pt_b, double pt_nonb)
{
    // Same sentinel convention as MomBalance: a negative return means "unusable for this
    // observable in this jet", and every fill site skips it.
    if (!(pt_b > 0.))     return -1.;   // no B pair (also catches the -999 ntuple sentinels)
    if (!(pt_nonb >= 0.)) return -1.;   // NaN, or an uninitialised sum

    // THE OBSERVABLE IS 1/(1 + pT_b/pT_nonb). What is returned is pT_nonb/(pT_b + pT_nonb),
    // which is the same number for every pT_nonb > 0 and is also right at pT_nonb = 0, where
    // the literal form divides by zero: a jet whose every selected track went into the two
    // Bs has an infinite ratio and so the observable is 0, which is what this returns.
    // Written this way for that one edge case only -- the axis label on every plot, and the
    // definition everywhere else, is 1/(1 + pT_b/pT_nonb).
    return pt_nonb / (pt_b + pt_nonb);
}

} // namespace NonBFraction

// ---- lnfb: ln(1 + pT_b/(pT_b + pT_nonb)) ------------------------------------------------
// THE DEFINITION IS THE FORMULA IN THE NAME OF THIS BLOCK, and the code below is that
// formula written out: the b momentum over the total, then log1p of it. Nothing is routed
// through a 1/(1 + ratio) rewrite (confirmed with Zoe, 2026-09-21) -- the earlier phrasing
// "f_b = 1/(1 + pT_nonb/pT_b)" was an algebraic aside that made a direct quantity look
// indirect, so it is gone from the code, the comments and the axis labels alike.
//
// ⚠️ The VALUE is unchanged by that rewrite. pT_b/(pT_b + pT_nonb) is what this function
// always evaluated, so every number produced before 2026-09-21 is still correct and no MC
// production had to be redone -- only the labels, which the plot macros read from ObsDef at
// draw time rather than from the stored histogram titles.
//
// This is the complement of fnb above (their sum is 1), so it is the same split of the jet,
// monotonically transformed. It carries no new information -- what changes is where the
// resolution goes.
//
// Range: the argument of the log is in (0, 1], so the value is in (0, ln 2] = (0, 0.6931].
// A jet whose every selected track went into the two Bs sits exactly at ln 2, and the bulk of
// the distribution crowds against that edge.
namespace LogBFraction {

inline double value(double pt_b, double pt_nonb)
{
    // Same sentinel convention as the other two: negative means "unusable in this jet".
    if (!(pt_b > 0.))     return -1.;
    if (!(pt_nonb >= 0.)) return -1.;
    // ln(1 + pT_b/(pT_b + pT_nonb)), straight from the two momenta. The denominator cannot be
    // zero here: pt_b > 0 and pt_nonb >= 0 are both checked above.
    return std::log1p(pt_b / (pt_b + pt_nonb));   // log1p, not log(1+x): accurate for a small ratio
}

} // namespace LogBFraction


// ---- fb: pT_b/(pT_b + pT_nonb), the b momentum fraction --------------------------------
// Requested 2026-09-21 as "pTb / pTtot" on a linear axis: the plain fraction, no log.
//
// ⚠️ THIS IS EXACTLY 1 - fnb, and on the uniform [0, 1] binning both observables use, bin k
// of fb is bin (nbins+1-k) of fnb -- the response matrix mirrors too. It is therefore NOT new
// physics, and a disagreement between an fb result and the corresponding fnb one is a bug in
// one of them, not a measurement. It exists as its own observable so that the axis a plot
// carries is the quantity that was asked for, rather than a mirrored relabelling of another
// one, and so it can be unfolded and banded by the same chain as everything else.
//
// It is also the argument of the log in LogBFraction above: lnfb = ln(1 + fb). Same three
// numbers, three ways -- fnb, fb, lnfb -- filled in one pass because they cost one extra set
// of histograms each and let the binnings be judged against each other.
//
// Range [0, 1]. fb = 1 is a jet whose every selected track went into the two B hadrons, which
// does happen, so the top bin is real and closed.
namespace BFraction {

inline double value(double pt_b, double pt_nonb)
{
    // Same sentinel convention as the others: negative means "unusable in this jet".
    if (!(pt_b > 0.))     return -1.;
    if (!(pt_nonb >= 0.)) return -1.;
    // Straight from the two momenta, like LogBFraction -- no 1/(1 + ratio) rewrite.
    return pt_b / (pt_b + pt_nonb);
}

} // namespace BFraction
*/



// ---- The axis abstraction --------------------------------------------------------------
struct ObsDef {
    TString         name;      // "dr" | "B" -- for log messages only
    TString         suffix;    // histogram-name suffix: "" for dr, "_B" for B
    TString         axis;      // axis title, for the histograms' own labels
    Int_t           nbins;
    const Double_t *bins;
    Double_t        maxFill;   // fold values at or above the top edge into the last bin
    Double_t        minReco;   // reco-level lower cut; < 0 means no cut
    // Is the TOP bin open-ended? dR folds every value above dr_max into its last bin, so
    // that bin really does mean "0.40 and above" and is quoted as "-> infinity". B is a
    // ratio pT_lead/(pT1+pT2) and cannot exceed 1, so its last bin is the closed [0.9, 1.0]
    // and quoting infinity would misstate the physics. Used by the plot labels.
    bool            upperOpen;

    // Overflow fold. Values below the first edge are NOT folded -- for dr that spike is cut
    // by minReco, and for B the range is closed at 0.5 by construction.
    double fill(double x) const { return (x >= bins[nbins]) ? maxFill : x; }

    // The reco-level domain cut. dr <= 0.005 is a reconstruction artefact (two SVs that are
    // really one), excluded from templates and response alike so both cover the same domain.
    // B has no such spike, so minReco is negative and this is always true for it.
    bool passReco(double x) const { return (minReco < 0.) || (x > minReco); }

    // Histogram name for this observable. dr's suffix is "", so dr keeps the names every
    // existing reader already knows -- template_fit.cpp, apply_unfolding_2d.C and every
    // result file on disk are unaffected by this refactor. That is verified, not assumed:
    // see the A/B check in README.md (repo root).
    TString n(const char *base) const { return TString(base) + suffix; }
};

// dr: the original observable. suffix "" -- its names must not change.
inline ObsDef obsDr()
{
    return ObsDef{"dr", "", "#DeltaR", dr_bins, dr_binsVector, dr_max_fill, 0.005, /*upperOpen=*/true};
}

// B: the momentum balance. No lower cut -- the distribution is closed at 0.5.
//
// ⚠️ The suffix is "_B". It was "_z" until 2026-09-22, and the two are NOT interchangeable
// when reading a file: an MC production written before the rename has its balance histograms
// named "..._z" and this code will not find them.
//
// Axis title notation: b1 is the HARDER of the two B hadrons and b2 the softer, so the
// label reads the same quantity MomBalance::value computes, max/(sum). The value itself is
// symmetric in its two arguments, so no ordering is imposed anywhere in the code -- b1/b2
// is a labelling convention for the plot, not an index into the SV or gen-pair arrays.
inline ObsDef obsB()
{
    return ObsDef{"B", "_B", "p_{T}^{b1}/(p_{T}^{b1} + p_{T}^{b2})",
                  B_bins, B_binsVector, B_max_fill, -1., /*upperOpen=*/false};
}

// ---- DISABLED 2026-09-22: the fraction observables' ObsDefs ----------------------------
// fnb, fb and lnfb are commented out along with their binnings and value functions above.
// Note that the comments below still say "z" where they mean the momentum balance -- they
// predate the rename and are left as they were; anything uncommented here needs B_*
// spellings, not z_*. The only edit made to the code itself is that the trailing
// /*upperOpen=*/ argument comments are gone: their */ would have closed this block early.
/*
// fnb: the non-b momentum fraction. Closed on [0, 1] by construction, like z -- fnb = 1
// would need the two B hadrons to carry no pT at all -- so the top bin is a real bin and
// nothing is folded into it from above.
inline ObsDef obsFnb()
{
    // Axis label is the requested form, 1/(1 + pT_b/pT_nonb), not the algebraically identical
    // one the code evaluates -- the plots say what the observable IS.
    return ObsDef{"fnb", "_fnb", "1/(1 + p_{T}^{b}/p_{T}^{non-b})",
                  fnb_bins, fnb_binsVector, fnb_max_fill, -1., false};
}

// fb: pT_b/(pT_b + pT_nonb), linear. Closed on [0, 1], so nothing is folded from above except
// the exact fb = 1 jet that fb_max_fill catches.
inline ObsDef obsFb()
{
    return ObsDef{"fb", "_fb", "p_{T}^{b}/(p_{T}^{b} + p_{T}^{non-b})",
                  fb_bins, fb_binsVector, fb_max_fill, -1., false};
}

// lnfb: ln(1 + pT_b/(pT_b + pT_nonb)). Closed at ln 2 = 0.6931 by construction, so upperOpen
// is false and nothing is folded from above -- the fold guard exists only for a value sitting
// exactly on the top edge, which is a jet with no non-b tracks at all.
//
// The axis title is the formula itself, not a named fraction. It is long, but it is what the
// observable IS, and a plot that says "f_b" sends the reader looking for a definition that
// lives in a header they do not have open.
inline ObsDef obsLnFb()
{
    return ObsDef{"lnfb", "_lnfb",
                  "ln(1 + p_{T}^{b}/(p_{T}^{b} + p_{T}^{non-b}))",
                  lnfb_bins, lnfb_binsVector, lnfb_max_fill, -1., false};
}
*/

// ---- The "current observable", for the template fit ------------------------------------
// template_fit.cpp and its drawing helpers in Help_Functions.h / Draw_EEC.h were written
// when the second axis could only be dR, so the binning and the axis labels are reached from
// a dozen places that have no argument to thread an ObsDef through. This global is set ONCE
// by the template_fit() driver before anything is booked, exactly like sDirname /
// sDirname_www already are, and read by CheckInputBinning() and the label helpers.
//
// It defaults to dR, so any code path that does NOT set it behaves exactly as before.
inline ObsDef &gFitObs()
{
    static ObsDef o = obsDr();
    return o;
}

inline void setFitObservable(const ObsDef &o) { gFitObs() = o; }

// The SHORT symbol for an observable, for a y-axis title like "dN/dB" where ObsDef::axis
// (the full formula) does not fit. One definition, used by apply_unfolding_2d.C and
// apply_weights_and_systematics.C, so the two never label the same plot differently.
// dr keeps the "#Delta r" spelling its plots have always had.
inline TString obsSymbol(const TString &name)
{
    if (name == "dr")  return "#Delta r";
    // 2026-10-12: Zoe renamed the balance "A" -- in the PLOTS only. The internal name, the
    // "_B" file/histogram/folder tags and every OBSERVABLE="B" argument stay as they are, so
    // the existing productions and results are still found.
    // if (name == "B")   return "B";
    if (name == "B")   return "A";
    // DISABLED 2026-09-22, with the fraction observables themselves:
    //   if (name == "fnb")  return "f_{nb}";   // the axis says 1/(1 + pT_b/pT_nonb)
    //   if (name == "fb")   return "p_{T}^{b}/p_{T}^{b+nb}";
    //   if (name == "lnfb") return "ln(1 + p_{T}^{b}/p_{T}^{b+nb})";
    return name;
}

// "dr" | "B" -> the ObsDef, or a zero-nbins ObsDef for anything else so the caller
// can reject it rather than silently fitting the wrong axis.
inline ObsDef obsByName(const TString &name)
{
    if (name == "dr")  return obsDr();
    if (name == "B")   return obsB();
    // DISABLED 2026-09-22 -- uncomment with the matching ObsDefs above:
    //   if (name == "fnb")  return obsFnb();
    //   if (name == "fb")   return obsFb();
    //   if (name == "lnfb") return obsLnFb();
    return ObsDef{"", "", "", 0, nullptr, 0., -1., false};
}

#endif // OBSERVABLES_H
