// Include guard: this header defines global variables, so a translation unit that reached
// it twice (directly and through observables.h) would fail to compile on the redefinitions.
// Guarding can only prevent that error, never change what any macro computes.
#ifndef BINNING_HISTOS_SMALL_H
#define BINNING_HISTOS_SMALL_H

// needed headers
#include "TH2D.h"
#include "TAxis.h"
#include "TString.h"
#include <cstdlib>


// Define binnings

//dr
const Int_t dr_binsVectorSize = 10;
Int_t bins_dr = dr_binsVectorSize - 1;
Int_t dr_bins = bins_dr;
Double_t dr_binsVector[dr_binsVectorSize] = {
  0,    
  0.05,
  0.10,
  0.15,
  0.20,
  0.25,
  0.30,
  0.35,
  0.40,
  0.45
};


Double_t dr_min = dr_binsVector[0];
Double_t dr_max = dr_binsVector[bins_dr];


//Values for histogram filling
Double_t dr_shiftbin = 0.00001;
Double_t dr_max_fill = 0.44;

// B -- momentum balance between the two aggregated B hadrons,
//
//     B = pT_lead / (pT1 + pT2)        in [0.5, 1)
//
// the SECOND measured observable, filled on its own axis in parallel with dr (same jets,
// same weight, different x). B = 0.5 is a balanced pair, B -> 1 a very asymmetric one.
//
// ⚠️ NAMED B, NOT z (renamed 2026-09-22 at Zoe's request, so the code says what the plots
// say). The rename went all the way to the histogram-name suffix "_B" and the production
// tag, so files written by the earlier "_z" code are NOT readable by this code -- see
// obsProdTag() in result_paths.h.
//
// Five uniform bins of 0.1 (requested 2026-09-15), replacing the equal-occupancy
// {0.500, 0.580, 0.660, 0.760, 1.000}.
//
// ⚠️ This is a round-number binning, not a migration-driven one, and the TOP BIN IS WEAK.
// momentum_balance_mc_study() scores uniform-5 at purity/stability 0.251/0.220 in
// 0.9 < B < 1.0, against 0.354/0.373 for the worst dr bin -- so that bin is reconstructing
// more migrated-in jets than its own. Cause: B_reco is pulled towards the bulk (~0.67) from
// both ends, and at B_gen ~ 0.99 the bias reaches -0.175, larger than the 0.100 resolution
// itself. The other four bins are fine (0.62/0.58 down to 0.38/0.44).
//
// Read the 0.9-1.0 bin with that in mind, and check it against the test_mode=0 closure
// before believing it. Re-run momentum_balance_mc_study() (in
// plot_purity_efficiency_response.cpp) before changing these again.
//
// 2026-10-21, tried and dropped (Zoe: back to the 5 bins, Bayesian 122 iterations):
//  - FOUR uniform bins of 0.125 -- full A chain re-run, Bayesian niter 125 (no speed-up), view in
//    results/workflow_rhojecveto_svsf_A4bins_20261021, files in results/backup_A4bins_20261021.
//  - TEN uniform bins of 0.05 -- nominal-only test (condor/run_A10bins_test.sh,
//    results/A10bins_test_20261021/out10): no iteration up to 200 reaches refolding p >= 0.95,
//    matrix inversion oscillates (one bin negative).
// The leading singular values of P(reco|gen) are the same in every binning (1.0, 0.41, 0.23,
// 0.15): the detector resolves ~3-4 A shape components; the binning cannot fix the convergence,
// only a better reconstruction (track aggregation) can.
// const Int_t B_binsVectorSize = 5;
// Double_t B_binsVector[B_binsVectorSize] = {0.500, 0.625, 0.750, 0.875, 1.000};
const Int_t B_binsVectorSize = 6;
Int_t B_bins = B_binsVectorSize - 1;
Int_t bins_B = B_bins;
Double_t B_binsVector[B_binsVectorSize] = {
    0.500,
    0.600,
    0.700,
    0.800,
    0.900,
    1.000
};
Double_t B_min = B_binsVector[0];
Double_t B_max = B_binsVector[B_bins];

// The upper edge is open (B = 1 needs one B hadron with zero pT), so fold it into the last
// bin rather than losing it to ROOT's overflow -- the same treatment dr gets via dr_max_fill.
Double_t B_max_fill = 0.999999;

// ============================================================================================
// DISABLED 2026-09-22: fnb, fb and lnfb are commented out.
//
// The measurement is now the momentum balance B and the EEC (dr) only, by request. The three
// b/non-b momentum fractions below are kept verbatim -- binning, warnings and the measured
// occupancies that justify them -- so that re-enabling one is uncommenting this block plus its
// counterparts in observables.h, result_paths.h and create_files_for_template_fit.cpp, and not
// a reconstruction from the readme.
//
// ⚠️ Their productions on disk ("_upartv2_3obs", "_upartv2_fb", "_upartv2_lnfblin") predate
// the z -> B rename, so re-enabling them ALSO means re-running step 1: the histogram names in
// those files still carry the old "_z" spelling for the balance axis.
// ============================================================================================
/*
// fnb -- the fraction of the jet's charged pT that is NOT in the two aggregated B hadrons,
//
//     fnb = 1/(1 + pT_b/pT_nonb) = pT_nonb / (pT_b + pT_nonb)     in [0, 1]
//
// the THIRD measured observable, filled on its own axis in parallel with dr and z (same
// jets, same weight, different x). fnb = 0 is a jet whose every selected track went into
// the two Bs; large fnb is a jet with a lot of non-b activity around them.
// See NonBFraction in observables.h for what counts as b and non-b on each side.
//
// Five uniform bins of 0.2 (requested 2026-09-18), replacing the shape-driven
// {0, .05, .10, .20, .30, .50, 1.0}. Round numbers over the full [0, 1] range, the same
// choice z made when its equal-occupancy binning was replaced by uniform 0.1 bins.
//
// What that costs, measured on 400k events of Pythia8 QCD (share of the weight per bin):
//
//                gen     reco    data (raw)
//     [0.0,0.2]   81%     47%       36%
//     [0.2,0.4]   14%     33%       30%
//     [0.4,0.6]    3%     11%       19%
//     [0.6,0.8]    2%      9%       13%
//     [0.8,1.0]  0.2%    0.6%        2%
//
// ⚠️ The GEN distribution is very lopsided on this binning -- four fifths of it is in the
// first bin, and the top bin holds a few per mille. That is the price of round numbers on an
// observable whose particle-level mean is 0.122: the bins are widest exactly where the
// distribution is. Read the top two bins knowing they are built on very little truth.
//
// ⚠️ And reco does not sit where gen does (means 0.264 vs 0.122; see NonBFraction in
// observables.h for why), so the response is strongly off-diagonal on ANY binning here. Five
// bins is the well-conditioned end of the choice -- ten 0.1-wide bins were tried on paper and
// the top three would have held under 1% of the truth each, which unregularised matrix
// inversion does not handle gracefully.
//
// The top edge is 1.0, which fnb cannot exceed, so nothing is folded.
// Still owed: a purity/stability study per bin, the equivalent of momentum_balance_mc_study().
const Int_t fnb_binsVectorSize = 6;
Int_t fnb_bins = fnb_binsVectorSize - 1;
Int_t bins_fnb = fnb_bins;
Double_t fnb_binsVector[fnb_binsVectorSize] = {
    0.0,
    0.2,
    0.4,
    0.6,
    0.8,
    1.0
};
Double_t fnb_min = fnb_binsVector[0];
Double_t fnb_max = fnb_binsVector[fnb_bins];

// fnb = 1 exactly (two B hadrons carrying no pT) cannot happen, so this guard never fires in
// practice; it is here so that a value at the top edge lands in the last bin rather than in
// ROOT's overflow, the same contract every other observable has.
Double_t fnb_max_fill = 0.999999;

// fb -- pT_b/(pT_b + pT_nonb), the b momentum fraction, on a LINEAR axis (requested
// 2026-09-21). Range [0, 1]; see BFraction in observables.h.
//
// Uniform 0.2 bins, the same five edges fnb uses -- deliberately, because fb = 1 - fnb and
// matching edges make the two exact mirrors of each other (bin k of fb is bin 6-k of fnb).
// That is the cross-check this observable brings: the two unfolded results must be each
// other's reverse, and if they are not, one of them is wrong.
//
// ⚠️ Expect the WEAK bins at the LOW end here, mirroring fnb's unusable top bin: particle
// level puts ~0.2% of its weight below fb = 0.2, so the first bin is recovered almost
// entirely from migration. Same warning as lnfb, same cause -- the reco definition (tracks
// the BDT rejected) and the gen one (tracks not from a B decay) agree only as far as the BDT
// separates them, and reco sits at mean 0.736 against gen 0.878 on this axis.
const Int_t fb_binsVectorSize = 6;
Int_t fb_bins = fb_binsVectorSize - 1;
Int_t bins_fb = fb_bins;
Double_t fb_binsVector[fb_binsVectorSize] = {
    0.0,
    0.2,
    0.4,
    0.6,
    0.8,
    1.0
};
Double_t fb_min = fb_binsVector[0];
Double_t fb_max = fb_binsVector[fb_bins];

// fb = 1 exactly IS reachable -- a jet with no non-b tracks at all -- so unlike fnb's guard
// this one really fires. It puts that jet in the last bin rather than in ROOT's overflow.
Double_t fb_max_fill = 0.999999;

// lnfb -- ln(1 + pT_b/(pT_b + pT_nonb)), computed straight from the two momenta.
// Range (0, ln 2 = 0.69315]; see LogBFraction in observables.h. It is the complement of fnb
// (the two fractions sum to 1), but it is NOT evaluated through fnb or through any
// 1/(1 + ratio) rewrite -- the formula in this line is the code.
//
// UNIFORM, by request (2026-09-21): five equal-width bins of ln2/5 = 0.138629 over the whole
// range (0, ln2]. This REPLACED a shape-driven set of edges, {0, .45, .55, .62, .66, .69315},
// and the trade it makes is deliberate and expensive -- read the warning before quoting a
// result off this axis.
//
// ⚠️ THE DISTRIBUTION IS NOT UNIFORM ON THIS AXIS. Measured on 400k events of Pythia8 QCD,
// the weight per 0.0347-wide slice piles against the top edge: gen has 41% in the last slice
// alone, 72% in the top two, and 4% below 0.45. Under these edges that puts the great bulk of
// the gen distribution in the top bin [0.55452, 0.69315] and leaves the first two bins --
// everything below 0.27726 -- nearly empty at particle level. The old shape-driven edges gave
// gen 4/10/14/30/41% and reco 18/22/32/19/9% instead.
//
// What that costs is the failure fnb's uniform 0.2 binning already showed: unregularised
// matrix inversion asked to recover an almost-empty truth bin from a reco bin filled mostly
// by migration returns an error of order 100%. Expect the LOW bins here to be the unusable
// ones (the mirror of fnb, where it was the top bin), and read the split test bin by bin
// before quoting any of them.
//
// The first edge is 0, not 0.45: ObsDef folds the HIGH side only, so anything below the
// first edge would go to ROOT's underflow and be lost -- and that is 18% of reco.
//
// ⚠️ The log does NOT spread the crowded end. The ratio pT_b/(pT_b + pT_nonb) sits near 1 for
// most jets, and there d ln(1+x)/dx = 1/(1+x) is about 0.5 -- so the log is close to a linear
// rescaling of that ratio and inherits its shape. Uniform bins in ln(1 + pT_b/(pT_b+pT_nonb))
// are therefore nearly uniform bins in the ratio itself, which is why the log does not rescue
// this binning. If the goal was to open up the b-dominated end, a ratio that is unbounded
// there -- ln(pT_b/pT_nonb), say -- is what does it. That is a DIFFERENT observable, not a
// relabelling of this one, and it would need its own production.
//
// Same reco/gen offset as fnb, in the other direction: mean 0.545 reco against 0.627 gen,
// because the BDT moves pT from the b side to the non-b side. Still owed: a purity/stability
// study per bin.
//
// ⚠️ THE BINNING IS BAKED INTO THE MC PRODUCTION -- create_files_for_template_fit.cpp books
// the templates and the response from these edges. Changing them means re-running step 1, not
// re-running the unfolding. That is why the uniform set got its own production tag,
// obsProdTag("lnfb") = "_upartv2_lnfblin": the shape-driven results already on disk stay
// reproducible from the "_upartv2_lnfb" production they were made from.
const Int_t lnfb_binsVectorSize = 6;
Int_t lnfb_bins = lnfb_binsVectorSize - 1;
Int_t bins_lnfb = lnfb_bins;
// ln2/5 = 0.1386294..., written out rather than computed so the edges are readable and the
// top one is the same literal 0.69315 the fold guard below sits just under.
Double_t lnfb_binsVector[lnfb_binsVectorSize] = {
    0.00000,
    0.13863,
    0.27726,
    0.41589,
    0.55452,
    0.69315
};
Double_t lnfb_min = lnfb_binsVector[0];
Double_t lnfb_max = lnfb_binsVector[lnfb_bins];

// The top edge is 0.69315, a hair ABOVE ln 2 = 0.6931472, so a jet with no non-b tracks at
// all -- which sits exactly at ln 2, and is the one value that can reach the top -- falls
// inside the last bin rather than on its edge. The fold below therefore never fires; it is
// here so every observable honours the same contract.
Double_t lnfb_max_fill = 0.693149;
*/


//mB
const Int_t mb_binsVectorSize = 11;
Int_t mb_bins = mb_binsVectorSize - 1;
Double_t mb_binsVector[mb_binsVectorSize] = {
    0., 
    1.,
    2.,
    3.,
    4.,
    5.,
    6.,
    7.,
    8,
    9,
    10
};
Double_t mb_min = mb_binsVector[0];
Double_t mb_max = mb_binsVector[mb_bins];
Int_t bins_mb = 10;
Double_t mb_max_fill = 9.9;

//EEC
const Int_t eec_binsVectorSize = 12;
Int_t bins_eec = eec_binsVectorSize - 1;
Int_t eec_bins = bins_eec;
Double_t eec_binsVector[eec_binsVectorSize] = {
    0., 
    100.,
    200.,
    300.,
    400.,
    500.,
    600.,
    700.,
    800.,
    900.,
    1000., 
    1e+8
};
Double_t eec_max = eec_binsVector[bins_eec];
Double_t eec_min = eec_binsVector[0];
Double_t eec_step = (eec_max-eec_min)/bins_eec;
Double_t eec_max_fill = 1e+8 -1;

//Jet Pt
const Int_t jtpt_binsVectorSize = 3;
Int_t jtpt_bins = jtpt_binsVectorSize - 1;
Double_t jtpt_binsVector[jtpt_binsVectorSize] = {
    80., 
    100., 
    120.
};
Double_t jtpt_min = jtpt_binsVector[0];
Double_t jtpt_max = jtpt_binsVector[jtpt_bins];

// The last pT bin is open-ended: "the last edge and above". Jets above jtpt_max are filled
// at jtpt_max_fill so they land in that bin rather than ROOT's overflow bin, which the pT
// slices in template_fit.cpp skip. Derived from the last edge.
Double_t jtpt_max_fill = jtpt_max - 0.1;
inline Double_t jtpt_fill(Double_t pt){ return (pt >= jtpt_max) ? jtpt_max_fill : pt; }


// bins_pt / pt_min / pt_max alias the jet-pT binning above and must never disagree with it:
// code loops to bins_pt then indexes jtpt_binsVector[ibin_pt]. Derived here, so adding an
// edge to jtpt_binsVector is the only edit needed to add a pT bin.
Int_t bins_pt = jtpt_bins;
Double_t pt_min = jtpt_min;
Double_t pt_max = jtpt_max;

// ---- 3D unfolding in (dR, A, jet pT), added 2026-10-27 ------------------------------------
// NON-SQUARE response: the side bins exist at GEN only -- pT underflow 60-80, pT overflow
// >= 120, dR overflow >= 0.45 -- so jets that migrate into the reco range from outside it are
// inside the matrix instead of being purity-corrected away. RooUnfold 2.0 ignores UseOverflow
// for 2D/3D, so the side bins are ordinary bins and the gen values are clamped into them
// (u3d_gen_*_fill), never left to ROOT's under/overflow.
// Reco has hard edges, pT 80-120 and dR < 0.45 (plus the dr > 0.005 cut): a reco jet outside
// them is a fake for the 3D purity. A is [0.5, 1] by construction on both sides, B_binsVector.
// Reco pT and dR are twice as fine as gen, so n_reco = 4x18x5 = 360 > n_gen = 4x10x5 = 200.
// ⚠️ The measured pT bin is 100-120 GeV here, NOT the open ">= 100" of the 2D chain.
// const Int_t u3d_reco_pt_size = 5;
// Double_t u3d_reco_pt[u3d_reco_pt_size] = {80., 90., 100., 110., 120.};
// const Int_t u3d_gen_pt_size = 5;
// Double_t u3d_gen_pt[u3d_gen_pt_size] = {60., 80., 100., 120., 1000.};   // 60-80 uf, 120+ of
// 2026-10-29 (Zoe): the nominal gen bin is 100-200, the overflow 200-1000; reco reaches 200 GeV.
const Int_t u3d_reco_pt_size = 6;
Double_t u3d_reco_pt[u3d_reco_pt_size] = {80., 90., 100., 125., 150., 200.};
const Int_t u3d_gen_pt_size = 5;
Double_t u3d_gen_pt[u3d_gen_pt_size] = {60., 80., 100., 200., 1000.};   // 60-80 uf, 200+ of
const Int_t u3d_reco_dr_size = 19;
Double_t u3d_reco_dr[u3d_reco_dr_size] = {
    0.,    0.025, 0.05,  0.075, 0.10,  0.125, 0.15,  0.175, 0.20,  0.225,
    0.25,  0.275, 0.30,  0.325, 0.35,  0.375, 0.40,  0.425, 0.45
};
const Int_t u3d_gen_dr_size = 11;
Double_t u3d_gen_dr[u3d_gen_dr_size] = {
    0., 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
    1.0    // 0.45+ overflow
};
inline Double_t u3d_gen_pt_fill(Double_t pt)
{
    const Double_t top = u3d_gen_pt[u3d_gen_pt_size - 1];
    return (pt >= top) ? top - 1. : pt;
}
inline Double_t u3d_gen_dr_fill(Double_t dr)
{
    const Double_t top = u3d_gen_dr[u3d_gen_dr_size - 1];
    return (dr >= top) ? top - 0.001 : dr;
}
inline Double_t u3d_A_fill(Double_t A) { return (A >= B_max) ? B_max_fill : A; }
inline bool u3d_reco_in(Double_t dr, Double_t pt)
{
    return pt >= u3d_reco_pt[0] && pt < u3d_reco_pt[u3d_reco_pt_size - 1]
        && dr < u3d_reco_dr[u3d_reco_dr_size - 1];
}

//Get the dimension label
Int_t mb_dim = 0;
Int_t dr_dim = 1;
Int_t eec_dim = 2;
Int_t pt_dim = 3;

//Recover the weighted dr distribution from a 2D distribution where dr and eec axes are separated
void recover_eec_distr(TH1D* &h_1D, TH2D* &h, Double_t last_eec = 1000){

    //save the eec along the bins (evaluate at the center of the bin)
    Float_t eec_step = (eec_max-eec_min)/eec_bins;
    Float_t eec = 0;
    Int_t eec_bins_tot = eec_bins;
    
    //save bin entries for Fill()
    Float_t mB, dr, x;

    //For each dr fill in each eec bin the corresponding bin content times the eec weight (at the center of the bin)
    for(Int_t ibin_dr = 1; ibin_dr <= bins_dr; ibin_dr++){

        dr = h->GetXaxis()->GetBinLowEdge(ibin_dr)+dr_shiftbin; //shift by a small amount to make sure the correct bin is filled
        Float_t dr_error = 0;

        //dr = h->GetXaxis()->GetBinCenter(ibin_dr);
        //Float_t dr_error = 0;   
        for(Int_t ibin_eec = 1; ibin_eec <= eec_bins_tot; ibin_eec++){

            eec = h->GetYaxis()->GetBinLowEdge(ibin_eec)+eec_step/2;
            //If last eec bin fill with a higher weight
            if(eec_bins_tot == ibin_eec) eec = last_eec;
            
            //Save the bin errors (added in quadrature)
            Float_t dr_err_temp = h->GetBinError(ibin_dr, ibin_eec)*eec;
            dr_error += pow(dr_err_temp, 2);

            h_1D->Fill(dr, h->GetBinContent(ibin_dr, ibin_eec)*eec);

        }
        //Approximate the bin error
        h_1D->SetBinError(ibin_dr, std::sqrt(dr_error));

    }

}

// ---- Jet-pT range labels (2026-10-29) ----
// JET_PT_MAX (env) is the hard upper jet-pT cut of create_files_for_template_fit (KinematicConfig::ptHigh):
// with it the open last bin (lo .. jtpt_max, everything above folded in) really is lo .. JET_PT_MAX.
inline double jetPtMaxCut() { const char *e = std::getenv("JET_PT_MAX"); return e ? std::atof(e) : 0.; }
// Gen-level jet-pT binning of the 2D response, efficiencies and gen corrections. Without JET_PT_MAX it is
// jtpt_binsVector (square response, as before). With it: 80 | 100 | JET_PT_MAX | 1000, the last bin the
// gen overflow (everything above 1000 folded in) -- the reco cut does not cut the gen jet (Zoe, 2026-10-29).
inline int jtptGenNbins() { return jetPtMaxCut() > 0. ? 3 : jtpt_bins; }
inline const Double_t *jtptGenEdges()
{
    static Double_t e[4];
    if (jetPtMaxCut() <= 0.) return jtpt_binsVector;
    e[0] = jtpt_binsVector[0]; e[1] = jtpt_binsVector[1]; e[2] = jetPtMaxCut(); e[3] = 1000.;
    return e;
}
inline Double_t jtptGenFill(Double_t pt)
{
    if (jetPtMaxCut() <= 0.) return jtpt_fill(pt);
    return (pt >= 1000.) ? 999.9 : pt;
}

inline TString jetPtRangeLabel(double lo, double hi, const char *var = "p_{T}^{jet}")
{
    if (hi >= jtpt_max - 1e-6) {   // the open-ended last bin
        const double m = jetPtMaxCut();
        return m > 0. ? TString::Format("%g < %s < %g GeV", lo, var, m) : TString::Format("%s > %g GeV", var, lo);
    }
    return TString::Format("%g < %s < %g GeV", lo, var, hi);
}

// ---- DR_MIN: lower dR cut (2026-10-30) ----
// Zoe: "run the whole analysis but without the first dr bin (start at 0.05)", cut at reco AND gen, and
// as a jet selection, so A loses those jets too. DR_MIN (env) must be an edge of dr_binsVector; the
// dR axis then starts there (obsDr() in observables.h drops the drFirstEdge() edges below it) and the
// reco domain cut dr > 0.005 becomes dr > DR_MIN. Unset = no cut, the 9 bins from 0 as before.
// From dr_binsVectorSize, not bins_dr: template_fit.cpp overwrites bins_dr with the fitted observable's count.
inline double drMinCut() { const char *e = std::getenv("DR_MIN"); return e ? std::atof(e) : 0.; }
inline int drFirstEdge()
{
    int k = 0;
    while (k < dr_binsVectorSize - 2 && dr_binsVector[k] < drMinCut() - 1e-9) ++k;
    return k;
}
// The gen side of the cut (the reco side is ObsDef::passReco of obsDr()). Always true when unset.
inline bool passDrMinGen(double dr_gen) { return drFirstEdge() == 0 || dr_gen > dr_binsVector[drFirstEdge()]; }

#endif // BINNING_HISTOS_SMALL_H

