#include "binning_histos_small.h"
#include "result_paths.h"
#include "observables.h"   // ObsDef + gFitObs(): which observable is on the second axis   // the ONE definition of every result path and variation tag
#include <algorithm>
#include <map>
#include <vector>
#include "TRandom3.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixDEigen.h"

namespace ROCColor {
  Color_t blue()   { return TColor::GetColor("#4C72B0"); }
  Color_t red()    { return TColor::GetColor("#C44E52"); }
  Color_t green()  { return TColor::GetColor("#4F8F52"); }
  Color_t purple() { return TColor::GetColor("#8C6BB1"); }
  Color_t orange() { return TColor::GetColor("#C97430"); }
  // Neutral highlight for the ratio-pad gap bands: a hue used by none of the five curves above,
  // so a filled region never reads as one of the plotted quantities. Drawn at two alphas.
  Color_t teal()   { return TColor::GetColor("#1F8A8A"); }
}

template <typename T>
T *getOrWarn(TFile *f, const TString &name)
{
    T *obj = dynamic_cast<T *>(f->Get(name.Data()));
    if (!obj)
        std::cerr << "ERROR: '" << name << "' (" << T::Class_Name() << ") not found in "
                  << f->GetName() << std::endl;
    return obj;
}

TFile *openOrWarn(const TString &name)
{
    TFile *f = TFile::Open(name);
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open " << name << std::endl;
        return nullptr;
    }
    return f;
}

std::vector<TFile *> openAllOrWarn(const std::vector<TString> &names)
{
    std::vector<TFile *> files;
    for (const TString &n : names) {
        TFile *f = openOrWarn(n);
        if (!f) {
            for (TFile *g : files) g->Close();
            files.clear();
            return files;
        }
        files.push_back(f);
    }
    return files;
}

// ---- Summing inputs across files ----------------------------------------------------
// One list covers both levels of summing: the per-block outputs of run_agg_ntuple_chunks.sh,
// and -- with sample "both" -- the bjet and qcd samples stacked on top of that. Same
// operation as merge_RMatrixTH2D_bjet_qcd_merged.C, done in memory so no pre-merged file
// is needed.
//
// ONLY count histograms and the response may be summed like this. The pre-divided ratios
// (h_full_purity_tf, h_full_efficiency_tf, ...) must NOT be: adding N of them gives N x the
// true value. They are recomputed from the summed counts by ratioFromCounts() below.
//
// 2026-10-24 (Zoe): "both" is a WEIGHTED qcd + bjet combination, as in template_fit.cpp. The
// b-enriched sample's weights are normalised like an inclusive sample (its 2B yield is ~14.5x
// qcd's), so every bjet file is scaled to qcd's 2B yield (k = S_qcd / S_bjet) and the two are
// mixed with statistics weights (w_bjet = Neff_bjet / (Neff_qcd + Neff_bjet), w_qcd = 1 - w_bjet).
// gMixWeight holds the factor per sample directory (QCD, bJet, QCDHerwig, bJetHerwig); a file
// in none of them, or any single-sample run, gets 1. Set by setMixWeights() below.
std::map<TString, double> gMixWeight;
double fileMixWeight(const TFile *f)
{
    const TString path = f->GetName();
    for (const auto &kv : gMixWeight)
        if (path.Contains("/" + kv.first + "/agg_ntuple_chunks/")) return kv.second;
    return 1.;
}
template <typename T>
T *sumOverFiles(const std::vector<TFile *> &files, const TString &name)
{
    T *out = nullptr;
    for (TFile *f : files) {
        T *h = getOrWarn<T>(f, name);
        if (!h) { delete out; return nullptr; }
        const double w = fileMixWeight(f);
        if (!out) {
            out = (T *) h->Clone(name + "_sum");
            out->SetDirectory(nullptr); // outlives the input files
            if (w != 1.) out->Scale(w);
        } else {
            // out->Add(h);   // until 2026-10-24: raw weights
            out->Add(h, w);
        }
    }
    if (!out) std::cerr << "ERROR: no input files to sum '" << name << "' from" << std::endl;
    return out;
}

RooUnfoldResponse *reweightPrior(const RooUnfoldResponse *r, const TH1 *prior)
{
    TH1 *meas  = (TH1 *) r->Hmeasured()->Clone("prior_meas");
    TH1 *truth = (TH1 *) r->Htruth()->Clone("prior_truth");
    TH2 *resp  = (TH2 *) r->Hresponse()->Clone("prior_resp");
    meas->SetDirectory(nullptr); truth->SetDirectory(nullptr); resp->SetDirectory(nullptr);
    const int nt = r->GetNbinsTruth(), nm = r->GetNbinsMeasured();
    double sumT = 0., sumP = 0.;
    for (int j = 0; j < nt; ++j) { sumT += RooUnfoldResponse::GetBinContent(truth, j); sumP += RooUnfoldResponse::GetBinContent(prior, j); }
    for (int j = 0; j < nt; ++j) {
        const int gj = RooUnfoldResponse::GetBin(truth, j);
        const double T = truth->GetBinContent(gj);
        const double w = (T > 0. && sumP > 0.) ? RooUnfoldResponse::GetBinContent(prior, j) * (sumT / sumP) / T : 1.;
        for (int i = 0; i < nm; ++i) {
            const double old = resp->GetBinContent(i + 1, j + 1);
            const int gi = RooUnfoldResponse::GetBin(meas, i);
            meas->SetBinContent(gi, meas->GetBinContent(gi) + old * (w - 1.));
            resp->SetBinContent(i + 1, j + 1, old * w);
            resp->SetBinError(i + 1, j + 1, resp->GetBinError(i + 1, j + 1) * w);
        }
        truth->SetBinContent(gj, T * w);
        truth->SetBinError(gj, truth->GetBinError(gj) * w);
    }
    return new RooUnfoldResponse(meas, truth, resp, TString(r->GetName()) + "_prior", r->GetTitle());
}

// A response with every histogram scaled by w (the weighted combination above). Built the
// same way as reweightPrior(): (measured, truth, response), fakes = measured - matched.
RooUnfoldResponse *scaledResponse(const RooUnfoldResponse *r, double w)
{
    TH1 *meas  = (TH1 *) r->Hmeasured()->Clone("mix_meas");  meas->SetDirectory(nullptr);  meas->Scale(w);
    TH1 *truth = (TH1 *) r->Htruth()->Clone("mix_truth");    truth->SetDirectory(nullptr); truth->Scale(w);
    TH2 *resp  = (TH2 *) r->Hresponse()->Clone("mix_resp");  resp->SetDirectory(nullptr);  resp->Scale(w);
    return new RooUnfoldResponse(meas, truth, resp, r->GetName(), r->GetTitle());
}

RooUnfoldResponse *sumResponseOverFiles(const std::vector<TFile *> &files, const TString &name)
{
    RooUnfoldResponse *out = nullptr;
    for (TFile *f : files) {
        RooUnfoldResponse *r = getOrWarn<RooUnfoldResponse>(f, name);
        if (!r) { delete out; return nullptr; }
        const double w = fileMixWeight(f);
        if (w != 1.) r = scaledResponse(r, w);
        if (!out) {
            out = new RooUnfoldResponse(*r);
            out->SetName(r->GetName());
            out->SetTitle(r->GetTitle());
        } else {
            out->Add(*r);
        }
    }
    if (!out) std::cerr << "ERROR: no input files to sum '" << name << "' from" << std::endl;
    return out;
}

// ---- Sample selection ---------------------------------------------------------------
// sample "qcd" | "bjet" | "both", generator "pythia" | "herwig" -> the agg_ntuple_chunks
// directories written by run_agg_ntuple_chunks.sh for those flags.
std::vector<TString> samplesIn(const TString &sample)
{
    if (sample == "both") return {"qcd", "bjet"};
    return {sample};
}

TString sampleSubdir(const TString &sample, const TString &generator)
{
    const bool herwig = (generator == "herwig");
    if (sample == "qcd")  return herwig ? "QCDHerwig"  : "QCD";
    if (sample == "bjet") return herwig ? "bJetHerwig" : "bJet";
    std::cerr << "ERROR: unknown sample '" << sample << "' (use qcd | bjet | both)" << std::endl;
    return "";
}

// Which of template_fit.cpp's variations supplies the signal fraction, and every path tag,
// now live in result_paths.h -- shared with apply_weights_and_systematics.C so the writer
// and the reader cannot disagree. varNames[] in Help_Functions.h are the same strings.

// The template fit that supplies the signal fraction, the one generator-dependent input
// that is NOT part of the unfolding MC. Returns "" when that fit does not exist, so both
// the driver and apply_unfolding() can refuse before doing any work.
TString templateFitFile(const TString &sample, const TString &tfGenerator,
                        bool track_eff_unc = false,
                        const TString &tfVariation = "nominal",
                        bool eec_weight_off = false,
                        const TString &observable = "dr")
{
    // Written by template_fit.cpp(SAMPLE, GENERATOR), one directory per combination.
    // Its SAMPLE is qcd | both only -- there is no bjet-only fit, since h3D_0b exists in
    // the qcd sample alone -- so a bjet unfolding takes the qcd fit.
    // const TString tf_sample = (sample == "both") ? "both" : "qcd";
    // Diagnostic (2026-10-24): TF_SAMPLE = qcd | both takes the fit of that sample whatever the
    // response sample is, to separate the two places the bjet sample enters. Unset = as before.
    // The results path carries it as tfSampleTag() (result_paths.h), so a TF_SAMPLE run cannot
    // land on the nominal result.
    const char *tfs = std::getenv("TF_SAMPLE");
    const TString tf_sample = (tfs && *tfs) ? TString(tfs) : ((sample == "both") ? TString("both") : TString("qcd"));
    // With track_eff_unc the fit itself is redone on the 3%-track-drop templates
    // (template_fit.cpp(SAMPLE, GENERATOR, true)), so the signal fraction moves with the
    // tracking efficiency instead of being held at its nominal value. The data being
    // fitted is the same either way.
    const TString trk_tag = track_eff_unc ? "_trkdrop030" : "";
    // The yield run fits the UNWEIGHTED templates -- template_fit.cpp(SAMPLE, GENERATOR,
    // false, true) -- so the signal fraction is the one a yield measurement needs rather
    // than the EEC-weighted one. Its own directory, never the EEC fit's.
    const TString eecw_tag = eec_weight_off ? "_noeecw" : "";
    // The observable tag template_fit() stamps on its own output directory. It must match,
    // or the B unfolding would read the dR fit's h_sig_fraction_fit -- which has the right
    // NAME and the wrong binning, the quietest possible way to be wrong.
    const TString obs_tag = (observable == "dr") ? "" : ("_" + observable);
    // TF_METHOD=chi2scan (2026-10-27): template_fit.cpp's chi2-scan fit, same file layout in
    // "<dir>_chi2scan". Unset = the RooFit directory, as before. tfMethodTag() tags the result.
    const TString method_dir = (tfMethodTag() == "_tfchi2scan") ? "_chi2scan" : "";
    return "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/TemplateFit_Run3/TemplateFits_"
         + tf_sample + "_" + tfGenerator + trk_tag + svsfTag() + trackBdtSfTag()
         + eecw_tag + obs_tag + eraTag() + method_dir + "/"
         + tfVariation + "_Run3_TemplateFits_histos_3d_80_inf.root";

    /* ---- disabled (kept for reference): Afnan's fits, Pythia8 only ----
    // These were the nominal until template_fit.cpp gained its own flags. Mixing them with
    // a Herwig fit produced here would confound generator with producer and b-tag WP
    // naming (btagWP712 vs btagWP0712), so both generators now come from the same macro.
    "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/results/TemplateFit_Run3/TemplateFits_btagWP712_qcdbjet_upartv2/nominal_Run3_TemplateFits_histos_3d_80_inf.root" // qcd + bjet
    "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/results/TemplateFit_Run3/TemplateFits_btagWP712_qcd_upartv2/nominal_Run3_TemplateFits_histos_3d_80_inf.root"     // qcd
    "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/TemplateFits_Run3_minHLT60_LinearBin_upartv2/nominal_Run3_TemplateFits_histos_3d_80_inf.root"            // Zoe, older
    ---- */
}

// MC systematic-variation tag, as create_files_for_template_fit.cpp writes it into the
// per-block INPUT filenames. A different namespace from the results-path tags -- this one
// names MC productions, those name result folders -- but it has to be the same literal, so
// it delegates to trkTag() in result_paths.h rather than repeating "_trkdrop030".
// (That literal tracks TrkEffSyst::kDropFraction in tracking_efficiency_syst.h.)
// It sits between the MCGEN part of the name and the run script's OUT_TAG, so a variation
// run looks for <varTag>_upartv2.
// obsProdTag() -- which MC production carries an observable -- lives in result_paths.h,
// shared with template_fit.cpp, which has to open the same blocks.

TString mcVarTag(bool track_eff_unc, bool eec_weight_off = false)
{
    // eecWeightTag() is "_noeecw": the yield production, written by the same macro with the
    // EEC weight switched off. Same delegation rule as above -- one literal, in
    // result_paths.h, shared by the MC-production names and the result folders.
    // svsfTag(): the SV b scale-factor production when SV_SF is set (2026-10-12), in the
    // position create_files_for_template_fit.cpp writes it (trk, svsf, noeecw).
    return trkTag(track_eff_unc) + svsfTag() + trackBdtSfTag()
         + eecWeightTag(eec_weight_off);
}

// One per-block file, or the top-level merged one when block < 0. prefix is "RMatrix_" or
// "" (the MCGEN templates), mid is "" or "MCGEN".
TString aggChunkFile(const TString &base, const TString &subdir, int block,
                     const TString &prefix, const TString &tag, const TString &mid,
                     const TString &btagTag, const TString &outSuffix)
{
    TString path = base + "/" + subdir + "/agg_ntuple_chunks/";
    if (block >= 0) path += Form("block_%04d/", block);
    return path + prefix + "Run3_" + btagTag + "_template_for_fit_histos_3D_" + tag + "_f"
                + mid + outSuffix + ".root";
}

// Every per-block file that exists, for one or both samples. Block counts differ per sample
// (10 qcd/pythia, 9 bjet/pythia, 8 qcd/herwig, 9 bjet/herwig), so take what is on disk
// rather than a fixed range -- same rule as the glob in run_agg_ntuple_chunks.sh.
std::vector<TString> aggChunkFiles(const TString &base, const TString &sample,
                                   const TString &generator, const TString &prefix,
                                   const TString &mid, const TString &btagTag,
                                   const TString &outSuffix, int max_blocks = 50)
{
    std::vector<TString> files;
    for (const TString &s : samplesIn(sample)) {
        const TString subdir = sampleSubdir(s, generator);
        if (subdir.Length() == 0) return {};
        int found = 0;
        for (int b = 0; b < max_blocks; ++b) {
            TString f = aggChunkFile(base, subdir, b, prefix, s, mid, btagTag, outSuffix);
            if (gSystem->AccessPathName(f) == kFALSE) { files.push_back(f); ++found; }
        }
        std::cout << "   " << subdir << ": " << found << " block files" << std::endl;
        if (found == 0)
            std::cerr << "ERROR: no block files under " << base << "/" << subdir
                      << "/agg_ntuple_chunks/ -- run run_agg_ntuple_chunks.sh for this sample"
                      << std::endl;
    }
    return files;
}


// The weighted qcd + bjet combination (see gMixWeight above), from the SAME numbers as
// template_fit.cpp: the 2B MC templates (h3D_bb, all bins) of the per-block MCGEN files, whose
// hadd the fit reads. S = sum of weights, Q = sum of squared weights, Neff = S^2 / Q.
// Adds the two factors of this generator to gMixWeight; sample != "both" adds nothing.
bool setMixWeights(const TString &base, const TString &sample, const TString &generator,
                   const TString &btagTag, const TString &outSuffix, const ObsDef &obs)
{
    if (sample != "both") return true;
    double S[2] = {0., 0.}, Q[2] = {0., 0.};
    const TString smp[2] = {"qcd", "bjet"};
    for (int s = 0; s < 2; ++s) {
        for (int b = 0; b < 50; ++b) {
            const TString f = aggChunkFile(base, sampleSubdir(smp[s], generator), b, "", smp[s], "MCGEN", btagTag, outSuffix);
            if (gSystem->AccessPathName(f)) continue;
            TFile *in = TFile::Open(f);
            TH1 *h = in ? dynamic_cast<TH1 *>(in->Get(obs.n("h3D_bb"))) : nullptr;
            if (!h) { std::cerr << "ERROR: no " << obs.n("h3D_bb") << " in " << f << std::endl; return false; }
            // Cells are 0 .. GetNcells()-1; "<=" read the last one twice (ROOT clamps the index).
            // for (int i = 0; i <= h->GetNcells(); ++i) { S[s] += h->GetBinContent(i); Q[s] += std::pow(h->GetBinError(i), 2); }   // until 2026-10-26
            for (int i = 0; i < h->GetNcells(); ++i) { S[s] += h->GetBinContent(i); Q[s] += std::pow(h->GetBinError(i), 2); }
            in->Close();
        }
    }
    if (S[0] <= 0. || S[1] <= 0. || Q[0] <= 0. || Q[1] <= 0.) {
        std::cerr << "ERROR: weighted combination: empty 2B templates (qcd " << S[0] << ", bjet " << S[1] << ")" << std::endl;
        return false;
    }
    const double k = S[0] / S[1], nq = S[0] * S[0] / Q[0], nb = S[1] * S[1] / Q[1];
    const double wb = nb / (nq + nb), wq = 1. - wb;
    gMixWeight[sampleSubdir("qcd", generator)]  = wq;
    gMixWeight[sampleSubdir("bjet", generator)] = wb * k;
    printf("qcd+bjet weighted combination (%s): k = %.5f (bjet 2B / qcd 2B = %.2f), Neff qcd %.0f, bjet %.0f"
           " -> w_qcd %.4f, w_bjet %.4f (bjet scaled by %.5f)\n", generator.Data(), k, 1. / k, nq, nb, wq, wb, wb * k);
    return true;
}

struct RefoldGof {
    double chi2   = 0.;
    int    ndf    = 0;
    double pvalue = -1.;
    double chi2ndf() const { return ndf > 0 ? chi2 / ndf : -1.; }
};

// h_mc_meas (optional, added 2026-10-01): the response's own measured MC distribution. When
// given, the refolded prediction carries the MC's relative statistical error in each bin, added
// in quadrature to the data error. Not in CMS-AN-2024-133 (chi2 of refolded vs input data),
// and not used for the choice since 2026-10-19 (Zoe: wrong); kept as a printed reference.
// Without it (the default) only the data error enters.
RefoldGof refoldChi2(const TH2D *h_refolded, const TH2D *h_input,
                     int ix_min, int ix_max, int iy_min, int iy_max,
                     const TH2 *h_mc_meas = nullptr)
{
    RefoldGof g;
    if (!h_refolded || !h_input) return g;
    for (int iy = iy_min; iy <= iy_max; ++iy) {
        for (int ix = ix_min; ix <= ix_max; ++ix) {
            const double d  = h_input->GetBinContent(ix, iy);
            const double ed = h_input->GetBinError(ix, iy);
            const double r  = h_refolded->GetBinContent(ix, iy);
            if (d == 0. || ed <= 0.) continue;
            double e2 = ed * ed;
            if (h_mc_meas) {
                const double m = h_mc_meas->GetBinContent(ix, iy);
                if (m > 0.) e2 += std::pow(r * h_mc_meas->GetBinError(ix, iy) / m, 2);
            }
            const double pull = (r - d) / std::sqrt(e2);
            g.chi2 += pull * pull;
            ++g.ndf;
        }
    }
    // ndf = number of bins compared. Refolding is not a fit to the input distribution,
    // so there are no fitted parameters to subtract.
    if (g.ndf > 0) g.pvalue = TMath::Prob(g.chi2, g.ndf);
    return g;
}

RefoldGof refoldChi2Toys(const RooUnfoldResponse *response, const TH2D *h_refolded,
                         const TH2D *h_input, int niter,
                         int ix_min, int ix_max, int iy_min, int iy_max,
                         int n_toys = 1000, UInt_t seed = 20261019,
                         RefoldGof *g_diag_r = nullptr, RefoldGof *g_cov_dr = nullptr)
{
    RefoldGof g;
    if (!response || !h_refolded || !h_input) return g;
    std::vector<std::pair<int, int>> bins;
    for (int iy = iy_min; iy <= iy_max; ++iy)
        for (int ix = ix_min; ix <= ix_max; ++ix)
            if (h_input->GetBinContent(ix, iy) != 0. && h_input->GetBinError(ix, iy) > 0.)
                bins.emplace_back(ix, iy);
    const int n = (int) bins.size();
    if (n == 0 || n_toys < n + 3) return g;

    const Bool_t add_dir = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    TRandom3 rng(seed);
    TH2D *h_toy = (TH2D *) h_input->Clone("h_refold_toy");
    TVectorD sum(n);
    TMatrixD sum2(n, n);
    TVectorD sum_r(n), sum_r2(n);
    TMatrixD sum_rr(n, n);
    for (int t = 0; t < n_toys; ++t) {
        for (const auto &b : bins)
            h_toy->SetBinContent(b.first, b.second, h_input->GetBinContent(b.first, b.second)
                                 + rng.Gaus(0., h_input->GetBinError(b.first, b.second)));
        RooUnfoldBayes unf(response, h_toy, niter);
        unf.SetVerbose(-1);
        TH1 *u = unf.Hreco(RooUnfold::kNoError);
        TH1 *r = response->ApplyToTruth(u, "h_refold_toy_r");
        TVectorD res(n);
        for (int k = 0; k < n; ++k) {
            const double rk = r->GetBinContent(bins[k].first, bins[k].second);
            res(k) = rk - h_toy->GetBinContent(bins[k].first, bins[k].second);
            sum_r(k) += rk;
            sum_r2(k) += rk * rk;
        }
        for (int a = 0; a < n; ++a) {
            const double ra = r->GetBinContent(bins[a].first, bins[a].second);
            for (int c = 0; c < n; ++c)
                sum_rr(a, c) += ra * r->GetBinContent(bins[c].first, bins[c].second);
        }
        sum += res;
        for (int a = 0; a < n; ++a)
            for (int c = 0; c < n; ++c) sum2(a, c) += res(a) * res(c);
        delete u;
        delete r;
    }
    delete h_toy;
    TH1::AddDirectory(add_dir);

    if (g_diag_r) {
        *g_diag_r = RefoldGof();
        for (int k = 0; k < n; ++k) {
            const double var_r = (sum_r2(k) - sum_r(k) * sum_r(k) / n_toys) / (n_toys - 1);
            const double ed = h_input->GetBinError(bins[k].first, bins[k].second);
            const double pull = (h_refolded->GetBinContent(bins[k].first, bins[k].second)
                               - h_input->GetBinContent(bins[k].first, bins[k].second))
                              / std::sqrt(ed * ed + std::max(var_r, 0.));
            g_diag_r->chi2 += pull * pull;
            ++g_diag_r->ndf;
        }
        g_diag_r->pvalue = TMath::Prob(g_diag_r->chi2, g_diag_r->ndf);
    }
    if (g_cov_dr) {
        *g_cov_dr = RefoldGof();
        TMatrixDSym W(n);
        TVectorD delta(n);
        for (int a = 0; a < n; ++a) {
            delta(a) = h_refolded->GetBinContent(bins[a].first, bins[a].second)
                     - h_input->GetBinContent(bins[a].first, bins[a].second);
            for (int c = 0; c < n; ++c)
                W(a, c) = (sum_rr(a, c) - sum_r(a) * sum_r(c) / n_toys) / (n_toys - 1);
            W(a, a) += std::pow(h_input->GetBinError(bins[a].first, bins[a].second), 2);
        }
        double det = 0.;
        W.Invert(&det);
        if (det > 0.) {
            g_cov_dr->chi2 = W.Similarity(delta);
            g_cov_dr->ndf  = n;
            g_cov_dr->pvalue = TMath::Prob(g_cov_dr->chi2, g_cov_dr->ndf);
        }
    }

    TMatrixDSym V(n);
    for (int a = 0; a < n; ++a)
        for (int c = 0; c < n; ++c)
            V(a, c) = (sum2(a, c) - sum(a) * sum(c) / n_toys) / (n_toys - 1);
    TMatrixDSymEigen eig(V);
    const TVectorD lam = eig.GetEigenValues();
    const TMatrixD vec = eig.GetEigenVectors();
    if (lam(0) <= 0.) return g;
    for (int k = 0; k < n; ++k) {
        if (lam(k) <= 1e-10 * lam(0)) continue;
        double proj = 0.;
        for (int a = 0; a < n; ++a)
            proj += vec(a, k) * (h_refolded->GetBinContent(bins[a].first, bins[a].second)
                               - h_input->GetBinContent(bins[a].first, bins[a].second));
        g.chi2 += proj * proj / lam(k);
        ++g.ndf;
    }
    g.chi2 *= double(n_toys - g.ndf - 2) / (n_toys - 1);
    if (g.ndf > 0) g.pvalue = TMath::Prob(g.chi2, g.ndf);
    return g;
}

bool normalizeToUnitArea(TH1D *h)
{
    if (!h) return false;
    const Double_t integral = h->Integral();
    if (integral <= 0.) {
        std::cerr << "WARNING: " << h->GetName() << " has integral " << integral
                  << ", skipping normalisation" << std::endl;
        return false;
    }
    h->Scale(1. / integral, "width");
    return true;
}


void apply_unfolding(TString &label, TString &folder, bool btag, Int_t n, TString pT_selection,
                     int test_mode, bool unfoldBayes, bool scan_niter,
                     TString sample = "qcd", TString generator = "pythia",
                     TString tfGenerator = "pythia", bool track_eff_unc = false,
                     TString tfVariation = "nominal", TString sfupartVariation = "nominal",
                     bool eec_weight_off = false, TString observable = "dr",
                     int niter_fixed = 7)
{
    // Every histogram name below carries this observable's suffix -- "" for dr, so the dr
    // path reads exactly the histograms it always read.
    const ObsDef obs = obsByName(observable);

    // Unfolding options. test_mode and unfoldBayes are passed in as arguments:
    //   0 = FULL-MC closure : h3D_bb, full-sample corrections, truth h_full_efficiency_denominator_tf.
    //       Same events in and out, so a technical (not statistically independent) closure.
    //   1 = SPLIT test      : h3D_pseudodata_bb (even half), odd-half corrections,
    //       truth h_pseudodata_truth_tf (even half) -> real, independent closure.
    //   2 = DATA            : h3D_data + template-fit signal fraction, full-sample corrections.
    //   unfoldBayes: true = Bayesian, false = matrix inversion.
    //   scan_niter:  true = scan niter 1..100 and pick the optimal one from the refolding
    //       goodness-of-fit test; false = a single unfolding at niter_fixed. Bayesian only.
    const bool split_test       = (test_mode == 1);  // -> pseudo (odd-half) corrections + even-half truth
    const bool is_data          = (test_mode == 2);
    const bool multiply_sigfrac = is_data;           // only real data needs the bb template fit
    // Matching purity = (reco_pass && gen_pass) / (reco_pass), both already restricted to
    // true-2b jets. It is NOT a bb-purity, so h3D_pseudodata_bb still needs it: that input
    // contains jets whose gen jet fails gen_pass, and the response has no truth row for them.
    bool apply_purity = true;


    // What the observable is called on every axis. The yield run is the same machinery
    // measuring dN/dr, so calling it "EEC" on the plots would be simply wrong.
    // Name the OBSERVABLE, not always dr: "dN/dB" and "EEC(B)" for the momentum balance.
    // obs.axis is the full label (p_{T}^{b1}/(...)), too long for a y-axis, so use a short
    // symbol for the quantity and put the full label on the x axis.
    const TString obs_sym   = obsSymbol(observable);
    const TString obs_title = eec_weight_off ? ("dN/d" + obs_sym) : ("EEC(" + obs_sym + ")");
    // The x-axis title, used by every plot below. dr spells it "#Delta r" rather than taking
    // obs.axis ("#DeltaR"), so that all of this chain's plots -- these and the ones
    // apply_weights_and_systematics.C draws next to them -- label the axis identically.
    const TString obs_axis  = (observable == "dr") ? "#Delta r" : obs.axis;

    // An axis with no explicit range is padded by ROOT when it draws it: B lives on [0.5, 1.0]
    // and comes out drawn to 1.1, empty frame where the observable cannot go (B = 1 needs one
    // B hadron with zero pT). Asking for the bin range explicitly is what stops it. Call this on
    // whichever histogram is drawn FIRST in a pad -- that one owns the pad's axes.
    // dr already draws as exactly [0, 0.45], so it is left alone and its plots are unchanged.
    // Same fix, same reason, as fixXRange() in apply_weights_and_systematics.C.
    auto fixXRange = [&](TH1 *h) {
        if (observable != "dr" && h) h->GetXaxis()->SetRange(1, h->GetNbinsX());
    };

    const Color_t blue = ROCColor::blue();
    const Color_t red = ROCColor::red();
    const Color_t green = ROCColor::green();
    const Color_t purple = ROCColor::purple();
    const Color_t orange = ROCColor::orange();
    const Color_t teal = ROCColor::teal();

    // ---- MC inputs follow the sample/generator flags ---------------------------------
    // The response, the corrections and the MC truth all come from the same file list, as
    // they must: one set of events, one pT binning, or the full-MC closure cannot be 1.
    // sample "both" appends the bjet blocks to the qcd ones and everything is summed.
    // Switching generator switches all of it -- that is the unfolding-model systematic.
    const TString mc_base   = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024";
    const TString btag_tag  = "btagWP0712"; // follows BTAG_WP in the run scripts (%04d of 1000*WP)
    // OUT_TAG in the run scripts, preceded by the systematic-variation tag the macro
    // writes. track_eff_unc = false is the nominal and gives exactly "_upartv2" as before;
    // true reads the 3%-tracking-efficiency production instead. The response, the
    // corrections, the truth and (in modes 0/1) the input all come from that same list,
    // so the whole unfolding moves with the variation -- which is the point.
    // ⚠️ The B histograms exist ONLY in the "_upartv2_B" production, run from this working
    // copy after the 2026-09-22 rename. The nominal "_upartv2" blocks hold dR only, and the
    // older "_zfirst" / "_3obs" / "_lnfblin" / "_fb" blocks carry the balance axis under its
    // OLD name "_z" -- a B run pointed at any of them fails on a missing h3D_data_B.
    const TString out_tag   = mcVarTag(track_eff_unc, eec_weight_off)
                            + obsProdTag(observable);

    std::cout << "Sample: " << sample << ", unfolding generator: " << generator
              << ", template-fit generator: " << tfGenerator << std::endl;

    std::vector<TString> filenames_response =
        aggChunkFiles(mc_base, sample, generator, "RMatrix_", "", btag_tag, out_tag);
    if (filenames_response.empty()) return;
    // The weighted qcd + bjet combination, for every sum below (2026-10-24).
    gMixWeight.clear();
    if (!setMixWeights(mc_base, sample, generator, btag_tag, out_tag, obs)) return;

    /* ---- disabled (kept for reference): response-matrix-only variation ----
    // Took ONLY the migration matrix from the other generator and left the corrections,
    // the truth and the input alone, isolating the response (the closure is then not
    // expected to be 1 -- that difference IS the uncertainty). Driven by a
    // responseGenerator argument, with a _resp<gen> tag on the folder and filenames.
    std::vector<TString> filenames_corrections =
        aggChunkFiles(mc_base, sample, generator, "RMatrix_", "", btag_tag, out_tag);
    std::vector<TString> filenames_response = filenames_corrections;
    if (responseGenerator != generator)
        filenames_response =
            aggChunkFiles(mc_base, sample, responseGenerator, "RMatrix_", "", btag_tag, out_tag);
    // ... and the response was then read from its own open-file list, everything else
    // from the corrections list.
    ---- */

    /* ---- disabled (kept for reference): pre-merged single response files ----
    // Afnan's merged files, Pythia8 only. The bjet+qcd one was produced by
    // merge_RMatrixTH2D_bjet_qcd_merged.C; sample="both" now does that sum in memory.
    "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/MergedResult_btagWP712_MattProd/RMatrix_Run3_btagWP712_template_for_fit_histos_3D_qcd_f_80_9999_2_merged.root"
    "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/bJet/agg_ntuple_chunks/MergedResult_btagWP712_MattProd/RMatrix_Run3_btagWP712_template_for_fit_histos_3D_bjet_f_80_9999_2_merged.root"
    "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/bJet/agg_ntuple_chunks/MergedResult_btagWP712_MattProd/RMatrix_Run3_btagWP712_template_for_fit_histos_3D_bjet_qcd_merged.root"
    ---- */

    //--  Signal fraction from the template fit. This is the OTHER generator-dependent
    // input and it is deliberately its own flag: tfGenerator, independent of the
    // unfolding generator above. The sample still picks qcd vs qcd+bjet templates.
    TString filename_template_fit = templateFitFile(sample, tfGenerator, track_eff_unc,
                                                    tfVariation, eec_weight_off, observable);

    if (multiply_sigfrac) {
        if (filename_template_fit.Length() == 0) {
            std::cerr << "ERROR: no template fit available for tfGenerator '" << tfGenerator
                      << "' -- produce it and set the path in apply_unfolding_2d.C" << std::endl;
            return;
        }
        std::cout << "Using template file (" << tfGenerator << "): "
                  << filename_template_fit << std::endl;
    }

    // Dataset ntuples to unfold: data or MC (test mode decision).
    // Modes 0/1 unfold MC (h3D_bb / h3D_pseudodata_bb) from the MCGEN files of the same
    // sample as the response, so they are summed the same way.
    // Mode 2 unfolds real data, which is one file and independent of the flags.
    std::vector<TString> filenames_data;
    if (is_data) {
        // The data file has to match the OBSERVABLE, not just the selection: h3D_data is
        // filled with eec * weight_tree, so unfolding the EEC-weighted file through an
        // unweighted response would hand back the EEC. The yield run takes its own data
        // production instead (make_hardprobes_condor_scripts.sh, EEC_WEIGHT_OFF=true).
        //
        // ⚠️ It also has to match the PRODUCTION: only the "_upartv2_B" data carries the
        // B histograms at all. The dR-era files below predate the observable axis, and the
        // pre-rename ones name the balance axis "_z", so a B run pointed at either fails on
        // a missing h3D_data_B rather than reading the wrong thing -- which is the right way
        // round, but it has to be pointed correctly here.
        // dataTemplateFile() (result_paths.h) is the ONE definition of which data file an
        // observable is measured from, shared with template_fit.cpp and plot_raw_yields.C --
        // the fit, the unfolding and the raw-yield plot must all read the same jets.
        filenames_data.push_back(dataTemplateFile(eec_weight_off, observable));
        // filenames_data.push_back("/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/HardProbes/agg_template_chunks/Run3_btagWP0712_template_for_fit_histos_3D_data_fMCGEN_upartv2.root"); // Zoe
    } else {
        filenames_data = aggChunkFiles(mc_base, sample, generator, "", "MCGEN", btag_tag, out_tag);
        if (filenames_data.empty()) return;
    }

    std::cout << "Getting data from " << filenames_data.size() << " file(s), first: "
              << filenames_data.front() << std::endl;
    //Select central pT bin
    int ibin_pt = 2;
    //Print options
    std::cout << "Options:"
              << "\n\tunfoldBayes:" << unfoldBayes
              << "\n\tibin_pt:" << ibin_pt
              << "\n\ttest_mode:" << test_mode
              << "\n\tmultiply_sigfrac:" << multiply_sigfrac
              << "\n\tsplit_test:" << split_test
              << "\n\tapply_purity:" << apply_purity
              << "\n\tscan_niter:" << scan_niter
              << std::endl;



    label += unfoldBayes ? "_bayesian" : "_MI";
    if(test_mode==0)  label += "_full_closure";
    if(test_mode==1)  label += "_split_test";
    if(test_mode==2)  label += "_data";
    if(!apply_purity) label += "_nopurity";
    // A scan writes its result at the LAST pass (niter_max), not at a chosen niter. Its result
    // file and the per-pass correction-stage / debug plots carry "_iterative" (2026-10-26), so a
    // scan can no longer overwrite the fixed-niter run's in the same folder. The scan-only
    // outputs (iteration PDF/PNGs, refolding and convergence scans) keep `label`: they never
    // clashed, and stage_final_svsf.sh reads them by that name.
    const TString out_label = scan_niter ? label + "_iterative" : label;
    TString fout_name;
    // fout_name = folder + "histos_" + label + "_after_unfolding_2D.root";   // until 2026-10-26
    fout_name = folder + "histos_" + out_label + "_after_unfolding_2D.root";

    // ---------------- Plotting setup ------------
    gSystem->Load("libRooUnfold.so");
    gStyle->SetErrorX(0.5);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    TGaxis::SetMaxDigits(3);   // "x10^-3" on the axis instead of 0.000123 on every label

    const Float_t  font_scale = 1200. / 800.;
    const Style_t  font_code  = 43;
    const Float_t  label_size = 15. * font_scale;
    const Float_t  title_size = 15. * font_scale;
    const Float_t  legend_size = 15. * font_scale;
    // One offset for every axis title, so they all sit the same distance from their numbers.
    const Float_t  title_offset = 3.0; // 1.0 not enough
    // ... except on the X axis of a RATIO pad. With precision-43 (pixel) fonts the offset is
    // measured in units of the PAD height, and the ratio pads here are 0.3 of the canvas, so
    // 3.0 puts the x title about a full pad below the axis -- off the bottom of the canvas,
    // drawn nowhere. That is why the correction-stages and bottomline plots came out with a
    // bare x axis: the title was set all along, it just landed outside the picture.
    // The Y titles keep title_offset, where the same scaling works out fine.
    // 1.3 still clipped the descender of the "#Delta" against the canvas edge; 1.0 clears the
    // tick labels and leaves the glyphs whole in both pads' bottom margins (0.23 and 0.18).
    const Float_t  ratio_x_title_offset = 1.0;

    // ---- Grab response matrix + corrections
    std::cout << "Getting response + corrections from " << filenames_response.size()
              << " file(s), first: " << filenames_response.front() << std::endl;
    std::vector<TFile *> fin_unfolding = openAllOrWarn(filenames_response);
    if (fin_unfolding.empty()) return;

    // hadd (and any summing merger) sums the per-block RATIO histograms across blocks,
    // giving N_blocks x the true value (e.g. hEECweightEff ~ 10 x 1.67 ~ 17). So DON'T read
    // the pre-divided ratios; recompute every correction here from the COUNT histograms,
    // which sum correctly. That applies to summing bjet onto qcd as well.
    // opt "b" = binomial (num is a subset of den), "" = normal errors.
    // TString, not const char*: the callers pass obs.n(...), which returns a TString.
    auto ratioFromCounts = [&](const TString &num, const TString &den, const char* opt) -> TH2D* {
        TH2D* hn = sumOverFiles<TH2D>(fin_unfolding, num);
        TH2D* hd = sumOverFiles<TH2D>(fin_unfolding, den);
        if (!hn || !hd) return nullptr;
        TH2D* r = (TH2D*) hn->Clone(Form("ratio_%s_over_%s", num.Data(), den.Data()));
        r->SetDirectory(0);
        r->Divide(hn, hd, 1., 1., opt);
        return r;
    };

    // Apply an MC-derived correction to the data as a pure bin-by-bin RESCALING: the
    // content and the error in each bin are multiplied by the SAME factor, so the data's
    // RELATIVE error is unchanged.
    //
    // This is deliberately NOT TH2::Multiply/Divide. Those propagate the correction's own
    // statistical error into the result, which mixes the MC's limited statistics into the
    // DATA's statistical error. The convention here is to treat every MC-derived correction
    // as EXACT and to carry its statistical precision separately, as its own systematic --
    // otherwise the same uncertainty is neither cleanly attributable nor separable.
    //
    // Runs over underflow and overflow too, so out-of-range bins stay consistent with the
    // in-range ones rather than being silently left uncorrected.
    //
    // A zero correction cannot be rescaled: the bin is emptied and counted, and the count is
    // reported. TH2::Divide would also produce 0 there, but silently.
    //
    // Returns false when it could not apply the correction (missing histogram, other binning),
    // and every caller then stops: until 2026-10-26 it printed an error and the run went on,
    // writing a result without that correction.
    //
    // Same bin count AND same edges on both axes. dR (9 bins) and A (5) already differ in count,
    // but two binnings with equal counts and different edges would pass a count-only check.
    auto sameBinning = [](const TH2 *a, const TH2 *b) -> bool {
        const TAxis *ax[2] = {a->GetXaxis(), a->GetYaxis()}, *bx[2] = {b->GetXaxis(), b->GetYaxis()};
        for (int k = 0; k < 2; ++k) {
            if (ax[k]->GetNbins() != bx[k]->GetNbins()) return false;
            for (int i = 1; i <= ax[k]->GetNbins() + 1; ++i)   // nbins + 1 = the upper edge
                if (std::fabs(ax[k]->GetBinLowEdge(i) - bx[k]->GetBinLowEdge(i)) > 1e-9) return false;
        }
        return true;
    };
    // auto applyCorrection = [](TH2D *h, const TH2D *corr, bool divide, const char *what) {   // until 2026-10-26
    auto applyCorrection = [&](TH2D *h, const TH2D *corr, bool divide, const char *what) -> bool {
        if (!h || !corr) {
            std::cerr << "ERROR: " << what << " is missing -- NOT applied, stopping" << std::endl;
            return false;
        }
        // if (corr->GetNbinsX() != h->GetNbinsX() || corr->GetNbinsY() != h->GetNbinsY()) {   // until 2026-10-26
        if (!sameBinning(h, corr)) {
            std::cerr << "ERROR: " << what << " has " << corr->GetNbinsX() << "x"
                      << corr->GetNbinsY() << " bins but the data has " << h->GetNbinsX()
                      << "x" << h->GetNbinsY() << " (or the edges differ) -- NOT applied, stopping"
                      << std::endl;
            return false;
        }
        int n_zero = 0;
        for (int ix = 0; ix <= h->GetNbinsX() + 1; ++ix) {
            for (int iy = 0; iy <= h->GetNbinsY() + 1; ++iy) {
                const double c = corr->GetBinContent(ix, iy);
                if (c == 0.) {
                    if (h->GetBinContent(ix, iy) != 0.) ++n_zero;
                    h->SetBinContent(ix, iy, 0.);
                    h->SetBinError(ix, iy, 0.);
                    continue;
                }
                const double f = divide ? (1. / c) : c;
                h->SetBinContent(ix, iy, h->GetBinContent(ix, iy) * f);
                h->SetBinError  (ix, iy, h->GetBinError(ix, iy)   * f);   // same factor
            }
        }
        std::cout << "\t---->" << (divide ? "Dividing by " : "Multiplying by ") << what
                  << " (bin-by-bin rescale; MC treated as exact)";
        if (n_zero) std::cout << " -- WARNING: " << n_zero
                              << " non-empty bin(s) had a zero correction and were emptied";
        std::cout << std::endl;
        return true;
    };

    // ----------- Grab data -----------
    std::vector<TFile *> fin_data = openAllOrWarn(filenames_data);
    if (fin_data.empty()) return;
    TString histname = (test_mode == 0) ? obs.n("h3D_bb")
                     : (test_mode == 1) ? obs.n("h3D_pseudodata_bb")
                     :                     obs.n("h3D_data");
    std::cout << "Using input histogram: " << histname << std::endl;
    TH3D *h_data_reco_3D_in = sumOverFiles<TH3D>(fin_data, histname);
    if (!h_data_reco_3D_in) {
        std::cerr << "       re-run create_files_for_template_fit.cpp to produce it." << std::endl;
        return;
    }
    TH3D *h_data_reco_3D = (TH3D*) h_data_reco_3D_in->Clone("h_data_reco_3D"); //att
    TH2D *h_data_reco = (TH2D*)h_data_reco_3D->Project3D("zy");
    TH2D *h_data_after_fit = (TH2D*) h_data_reco->Clone("h_data_after_fit");



    // Multiply histograms by signal fraction

    if (multiply_sigfrac) {
        std::cout << "\t---->Multiplying by signal fraction" << std::endl;
        // Grab signal fraction from template fit
        TString fname_fit = filename_template_fit;
        std::cout << "Getting signal fraction from " << fname_fit << std::endl;
        TFile *fin_fit = openOrWarn(fname_fit);
        if (!fin_fit) return;
        TH2D *h_sig_fraction = getOrWarn<TH2D>(fin_fit, "h_sig_fraction_fit");
        if (!h_sig_fraction) return;
        // A fit of another observable has the right name and the wrong binning; TH2::Multiply
        // would print an error and leave the data unscaled. Refuse instead (2026-10-26).
        if (!sameBinning(h_data_after_fit, h_sig_fraction)) {
            std::cerr << "ERROR: h_sig_fraction_fit (" << h_sig_fraction->GetNbinsX() << "x"
                      << h_sig_fraction->GetNbinsY() << ") does not have the data's binning ("
                      << h_data_after_fit->GetNbinsX() << "x" << h_data_after_fit->GetNbinsY()
                      << ") in " << fname_fit << std::endl;
            return;
        }
        // h_data_after_fit->Multiply(h_sig_fraction);   // until 2026-10-26: return value unchecked
        if (!h_data_after_fit->Multiply(h_sig_fraction)) {
            std::cerr << "ERROR: multiplying by the signal fraction failed" << std::endl;
            return;
        }

    }
    else {
        std::cout << "\t---->Not multiplying by signal fraction" << std::endl;
    }

    // ---- UParT b-tag efficiency scale factor, on the RECO-level data ------------------
    // Applied HERE, before unfolding, and not to the unfolded result: the SF is measured in
    // reco R_BB bins, so it belongs on the distribution that is still in reco R_BB. The
    // response matrix then carries it through the migration like any other feature of the
    // data. Applying it afterwards would put a reco-binned correction on particle-level
    // bins and silently ignore that migration.
    //
    // It is DIVIDED in: the chain corrects the data with the efficiency measured in MC, and
    // the true efficiency is eps_data = SF * eps_MC, so data/eps_MC has to be further
    // divided by SF. Data only -- there is no data in the closure modes, and the SF must
    // never touch MC.
    if (is_data && sfupartVariation == "off") {
        std::cout << "\t---->UParT SF NOT applied (SFUPART_VARIATION=off): this is the BEFORE-SF "
                  << "result, for comparison only" << std::endl;
    }
    else if (is_data) {
        const TString sf_hist = sfupartHist(sfupartVariation, observable);
        std::cout << "\t---->Dividing reco data by the UParT SF (" << sfupartVariation
                  << " -> " << sf_hist << ")" << std::endl;
        TFile *fin_sf = openOrWarn(sfupartFile(observable));
        if (!fin_sf) return;
        TH1D *h_sf = getOrWarn<TH1D>(fin_sf, sf_hist);
        if (!h_sf) return;

        // For the variations that propagate a QUOTED UNCERTAINTY rather than swapping in a
        // different curve (jpsyst_up/jpsyst_dn), the per-bin error to add sits in its own
        // histogram in the same file. Null for every other variation, and the arithmetic
        // below then reduces to exactly what it was.
        TH1D *h_sf_syst = nullptr;
        const TString sf_syst_hist = sfupartSystHist(sfupartVariation, observable);
        if (sf_syst_hist.Length()) {
            h_sf_syst = getOrWarn<TH1D>(fin_sf, sf_syst_hist);
            if (!h_sf_syst) return;
            // Indexed with the SAME bin number as the SF, so a mismatch would silently apply
            // the wrong bin's uncertainty. Refuse instead.
            if (h_sf_syst->GetNbinsX() != h_sf->GetNbinsX()) {
                std::cerr << "ERROR: " << sf_syst_hist << " has " << h_sf_syst->GetNbinsX()
                          << " bins, " << sf_hist << " has " << h_sf->GetNbinsX() << std::endl;
                return;
            }
        }

        const int nx = h_data_after_fit->GetNbinsX();   // x is dr, y is pt
        const int ny = h_data_after_fit->GetNbinsY();

        // The SF has FEWER dr bins than the analysis: its last bin spans 0.35-0.45, which
        // covers the analysis's last TWO bins. So the analysis bin is matched to the SF bin
        // by BIN CENTRE rather than by index -- bins 8 and 9 both land in SF bin 8 on their
        // own, with no special case to get wrong if the binning changes again.
        auto sfBinFor = [&](int ix) -> int {
            // Underflow/overflow have no meaningful centre: clamp them to the nearest real
            // SF bin, so every bin of the data gets corrected by something sensible.
            if (ix < 1)  return 1;
            if (ix > nx) return h_sf->GetNbinsX();
            const double c = h_data_after_fit->GetXaxis()->GetBinCenter(ix);
            int b = h_sf->GetXaxis()->FindBin(c);
            if (b < 1)                     b = 1;
            if (b > h_sf->GetNbinsX())     b = h_sf->GetNbinsX();
            return b;
        };

        // Bin by bin: read the content, multiply by the INVERSE of the SF for that dr bin,
        // write it back.
        //
        // INVERSE, i.e. divided in, because SF = eff_b_data / eff_b_mc (verified against
        // h_eff_b_data_* / h_eff_b_mc_* in the same file) and the chain corrects the data
        // with the efficiency measured in MC -- Divide(h_full_efficiency) and
        // Divide(h_svbtag_eff) above. The true efficiency is eps_data = SF * eps_MC, so
        // data/eps_MC still has to be divided by SF. Physically: SF < 1 means data tags
        // LESS efficiently than MC, so the yield has to be scaled UP, which is 1/SF > 1.
        // Multiplying would move it the wrong way. Flip this one flag if that reading is
        // ever wrong -- everything else follows from it.
        const bool kDivideBySF = true;
        //
        // No jet pT dependence is assumed: the same dr-dependent SF is applied to every pT
        // bin. The loop below runs over all of y for that reason.
        //
        // The loops run 0..n+1, so the UNDERFLOW and OVERFLOW bins are corrected too. Both
        // are empty in the current data (dr underflow = 0 against 7.9e8 in range), so this
        // is a no-op today -- but it keeps the out-of-range bins consistent with the rest
        // if they are ever filled, rather than leaving them silently uncorrected.
        // statup/statdn shift the central SF by +/- its own statistical error, coherently in
        // every dr bin. That error is the calibration's statistical precision; it is booked
        // as its OWN systematic rather than being folded into the data's statistical error,
        // so that the data error bar stays a data error bar.
        const double stat_shift = (sfupartVariation == "statup") ?  1.
                                : (sfupartVariation == "statdn") ? -1. : 0.;
        // The JP HF CALIBRATION uncertainty, applied the same coherent way: +/- the quoted
        // per-bin systematic in every dr bin at once. This is the second of the two ways to
        // use the JP HF measurement -- see isKnownSfupartVariation() in result_paths.h for
        // why applying the alternative curve (jpcalib_hf) is not the same thing.
        const double syst_shift = sfupartSystShift(sfupartVariation);

        // The SF of one SF bin, with the stat / JP HF shifts of this variation.
        auto sfValue = [&](int isf) -> double {
            return h_sf->GetBinContent(isf) + stat_shift * h_sf->GetBinError(isf)
                 + (h_sf_syst ? syst_shift * h_sf_syst->GetBinContent(isf) : 0.);
        };
        // The SF for analysis bin ix (added 2026-10-21, with the 4-bin A): the average of the SF
        // bins it overlaps, weighted by the overlap width. An analysis bin inside ONE SF bin (every
        // dr bin, and the old 5-bin A, which IS the SF binning) gets exactly that bin's SF, as
        // sfBinFor() gave. A bin straddling two SF bins (4-bin A: 0.625-0.75 spans SF 0.6-0.7 and
        // 0.7-0.8) gets their mean instead of the one its centre falls in -- by centre, SF bin
        // 0.7-0.8 was used by no A bin at all. Under/overflow keep the clamped sfBinFor().
        // SFUPART_SMOOTH=line (2026-10-30, A only, result_paths.h): the per-bin values of THIS
        // variation fitted by a straight line in A (weights 1/stat error^2 of the SF bins), and every
        // analysis bin gets the line at its centre (= the line's average over the bin).
        const bool smoothLine = (sfupartSmoothTag(observable) != "");
        double lineA = 0., lineB = 0.;
        if (smoothLine) {
            double S = 0., Sx = 0., Sy = 0., Sxx = 0., Sxy = 0.;
            for (int k = 1; k <= h_sf->GetNbinsX(); ++k) {
                const double e = h_sf->GetBinError(k), w = e > 0. ? 1. / (e * e) : 1.;
                const double x = h_sf->GetXaxis()->GetBinCenter(k), y = sfValue(k);
                S += w; Sx += w * x; Sy += w * y; Sxx += w * x * x; Sxy += w * x * y;
            }
            lineB = (S * Sxy - Sx * Sy) / (S * Sxx - Sx * Sx);
            lineA = (Sy - lineB * Sx) / S;
            printf("\n  UParT SF SMOOTHED (SFUPART_SMOOTH=line): SF(A) = %.5f + %.5f * A; per bin:", lineA, lineB);
            for (int k = 1; k <= h_sf->GetNbinsX(); ++k)
                printf("  %.4f->%.4f", sfValue(k), lineA + lineB * h_sf->GetXaxis()->GetBinCenter(k));
            printf("\n");
        }
        auto sfFor = [&](int ix) -> double {
            if (smoothLine) {
                const double c = (ix < 1 || ix > nx) ? h_sf->GetXaxis()->GetBinCenter(sfBinFor(ix))
                                                     : h_data_after_fit->GetXaxis()->GetBinCenter(ix);
                return lineA + lineB * c;
            }
            if (ix < 1 || ix > nx) return sfValue(sfBinFor(ix));
            const double lo = h_data_after_fit->GetXaxis()->GetBinLowEdge(ix);
            const double hi = h_data_after_fit->GetXaxis()->GetBinUpEdge(ix);
            double sum = 0., wsum = 0.;
            int nover = 0, kover = -1;
            for (int k = 1; k <= h_sf->GetNbinsX(); ++k) {
                const double o = std::min(hi, h_sf->GetXaxis()->GetBinUpEdge(k))
                               - std::max(lo, h_sf->GetXaxis()->GetBinLowEdge(k));
                if (o <= 0.) continue;
                sum += o * sfValue(k); wsum += o; ++nover; kover = k;
            }
            // Inside ONE SF bin: that bin's SF exactly. (o * SF) / o is not always SF to the last
            // bit, which moved the dR MI band by 3e-13 relative (2026-10-22 restore check).
            if (nover == 1) return sfValue(kover);
            return wsum > 0. ? sum / wsum : sfValue(sfBinFor(ix));
        };

        printf("\n  UParT SF applied bin by bin (reco level, pT bin %d, %s%s%s)\n",
               ibin_pt, kDivideBySF ? "data / SF" : "data * SF",
               stat_shift > 0 ? ", SF shifted UP by its stat error"
             : stat_shift < 0 ? ", SF shifted DOWN by its stat error" : "",
               syst_shift > 0 ? ", SF shifted UP by the JP HF syst"
             : syst_shift < 0 ? ", SF shifted DOWN by the JP HF syst" : "");
        printf("  %3s %13s %5s %9s %9s %14s %14s %8s\n",
               "bin", "dr range", "sfbin", "SF", "factor", "before", "after", "ratio");
        printf("  -------------------------------------------------------------------------------------\n");

        for (int ix = 0; ix <= nx + 1; ++ix) {
            const int  isf = sfBinFor(ix);   // the centre's SF bin, for the printout only
            // const double sf = h_sf->GetBinContent(isf)
            //                 + stat_shift * h_sf->GetBinError(isf)
            //                 + (h_sf_syst ? syst_shift * h_sf_syst->GetBinContent(isf) : 0.);
            const double sf = sfFor(ix);
            if (sf == 0.) {
                std::cerr << "ERROR: UParT SF is zero in its bin " << isf << std::endl;
                return;
            }
            const double factor = kDivideBySF ? (1. / sf) : sf;

            const double before = (ix >= 1 && ix <= nx)
                                ? h_data_after_fit->GetBinContent(ix, ibin_pt) : 0.;

            for (int iy = 0; iy <= ny + 1; ++iy) {
                // A pure rescale: the content and the DATA's statistical error are both
                // multiplied by the same factor, so the data's relative error is unchanged
                // and the error bar still means what it says. The SF's own error is NOT
                // added in quadrature here -- it is the statup/statdn systematic instead.
                // Same treatment the MC corrections get further down.
                h_data_after_fit->SetBinContent(ix, iy,
                                                h_data_after_fit->GetBinContent(ix, iy) * factor);
                h_data_after_fit->SetBinError  (ix, iy,
                                                h_data_after_fit->GetBinError(ix, iy)   * factor);
            }

            if (ix < 1 || ix > nx) continue;   // nothing useful to print for under/overflow
            const double after = h_data_after_fit->GetBinContent(ix, ibin_pt);
            printf("  %3d [%5.2f,%5.2f] %5d %9.5f %9.5f %14.6g %14.6g %8.4f%s\n", ix,
                   h_data_after_fit->GetXaxis()->GetBinLowEdge(ix),
                   h_data_after_fit->GetXaxis()->GetBinUpEdge(ix),
                   isf, sf, factor, before, after, (before != 0. ? after / before : 0.),
                   (ix > 1 && isf == sfBinFor(ix - 1)) ? "   <- shares the SF bin above" : "");
        }
        printf("\n");
        fin_sf->Close();
    }



     //Dimension of the response matrix
    // The OBSERVABLE's bin count, not the global dr one (fixed 2026-10-21): bins_dr is dr's 9,
    // so for A (5 bins) every ibin_dr_max loop ran to 9. The chi2 / bias numbers were unharmed
    // (ROOT clamps bins 6-9 to the empty overflow), but GetBinUpEdge(9) = 1.4 on the A axis, so
    // the bottomline ratio pad's unity line and grey band ran past A = 1.
    // int dim = bins_pt*bins_dr;
    // int ibin_dr_max = bins_dr;
    const int bins_obs = obs.nbins;
    int dim = bins_pt*bins_obs;
    int ibin_dr_min = 1;
    int ibin_dr_max = bins_obs;

    std::cout << "dim = " << dim << std::endl;
    std::cout << "bins_pt = " << bins_pt << std::endl;
    std::cout << "bins (" << obs.name << ") = " << bins_obs << std::endl;

    // Note: Result = unfold(raw * purity) * 1 / (efficiency)
    //       fakes are negligible

    


    TH2D *h_full_purity = nullptr;
    TH2D *h_full_efficiency = nullptr;
    TH2D *h_mc_reco = nullptr;
    RooUnfoldResponse *response = nullptr;
    // ⚠️ ASSIGNED BUT NEVER USED (checked 2026-09-22). Nothing downstream reads it, so the
    // two branches below are dead. Left in place because it is the obvious handle for a
    // "before the efficiency correction" comparison -- but note that in the split test it is
    // assigned the ODD half (the response half's matched gen), NOT the even half the
    // pseudodata comes from. If it is ever wired up to a plot, it must be swapped for the
    // even-half equivalent first, or that comparison would be against the matrix's own truth.
    TH2D *h_mc_true_no_eff = nullptr;

    // The reco-level MC comparison must match what the data curve has had done to it:
    // purity-corrected data -> matched reco (numerator); uncorrected data -> all reco (denominator).
    if (split_test){
        std::cout << "\t----> Doing split test" << std::endl;
        // --- 
        // h_full_purity = getOrWarn<TH2D>(fin_unfolding, obs.n("h_full_pseudo_purity_tf")); // wrong ratio due to merged files 
        // h_full_efficiency = getOrWarn<TH2D>(fin_unfolding, obs.n("h_full_pseudo_efficiency_tf")); // wrong ratio due to merged files 
        // --- 

        // The third argument is the Divide() OPTION, not a name: "h_pseudo_purity" has no "b", so
        // these got normal instead of binomial errors (until 2026-10-26). Numerator is a subset of
        // the denominator, same per-jet weight -> "b". The errors are not used downstream
        // (applyCorrection treats corrections as exact), so the result does not change.
        // h_full_purity = ratioFromCounts(obs.n("h_full_pseudo_purity_numerator_tf"), obs.n("h_full_pseudo_purity_denominator_tf"), "h_pseudo_purity");
        // h_full_efficiency = ratioFromCounts(obs.n("h_full_pseudo_efficiency_numerator_tf"), obs.n("h_full_pseudo_efficiency_denominator_tf"), "h_pseudo_efficiency");
        h_full_purity = ratioFromCounts(obs.n("h_full_pseudo_purity_numerator_tf"), obs.n("h_full_pseudo_purity_denominator_tf"), "b");
        h_full_efficiency = ratioFromCounts(obs.n("h_full_pseudo_efficiency_numerator_tf"), obs.n("h_full_pseudo_efficiency_denominator_tf"), "b");
        
        h_mc_reco = sumOverFiles<TH2D>(fin_unfolding, apply_purity ? obs.n("h_full_pseudo_purity_numerator_tf")
                                                                  : obs.n("h_full_pseudo_purity_denominator_tf"));
        response = sumResponseOverFiles(fin_unfolding, obs.n("response_tf_pseudo_full"));
        h_mc_true_no_eff = sumOverFiles<TH2D>(fin_unfolding, obs.n("h_full_pseudo_efficiency_numerator_tf"));
    }
    else {
        // ---  
        // h_full_purity = getOrWarn<TH2D>(fin_unfolding, obs.n("h_full_purity_tf")); // Wrong ratio due to merged files 
        // h_full_efficiency = getOrWarn<TH2D>(fin_unfolding, obs.n("h_full_efficiency_tf")); // wrong ratio due to merged files 
        // --- 

        // Binomial, as in the split test above (until 2026-10-26 the option was "h_purity" /
        // "h_efficiency", i.e. normal errors).
        // h_full_purity     = ratioFromCounts(obs.n("h_full_purity_numerator_tf"), obs.n("h_full_purity_denominator_tf"), "h_purity");
        // h_full_efficiency = ratioFromCounts(obs.n("h_full_efficiency_numerator_tf"), obs.n("h_full_efficiency_denominator_tf"), "h_efficiency");
        h_full_purity     = ratioFromCounts(obs.n("h_full_purity_numerator_tf"), obs.n("h_full_purity_denominator_tf"), "b");
        h_full_efficiency = ratioFromCounts(obs.n("h_full_efficiency_numerator_tf"), obs.n("h_full_efficiency_denominator_tf"), "b");

        h_mc_reco = sumOverFiles<TH2D>(fin_unfolding, apply_purity ? obs.n("h_full_purity_numerator_tf")
                                                                  : obs.n("h_full_purity_denominator_tf"));
        response = sumResponseOverFiles(fin_unfolding, obs.n("response_tf_full"));
        h_mc_true_no_eff = sumOverFiles<TH2D>(fin_unfolding, obs.n("h_full_efficiency_numerator_tf"));
    }
    if (!h_full_purity || !h_full_efficiency || !h_mc_reco || !response) return;

    // ---- Unfolding prior from another generator (UNFOLD_PRIOR, priorTag() in result_paths.h) ----
    // Everything above -- the migrations, purity, efficiency and the corrections after the
    // unfolding -- stays the unfolding generator's (the detector response). Only the truth the
    // Bayesian iterations start from is taken from the prior generator's response truth.
    {
        const char *pe = std::getenv("UNFOLD_PRIOR");
        const TString prior_gen = (pe && *pe) ? TString(pe) : TString("");
        if (prior_gen.Length()) {
            std::vector<TString> filenames_prior =
                aggChunkFiles(mc_base, sample, prior_gen, "RMatrix_", "", btag_tag, out_tag);
            if (!setMixWeights(mc_base, sample, prior_gen, btag_tag, out_tag, obs)) return;   // the prior generator's own mix
            std::vector<TFile *> fin_prior = openAllOrWarn(filenames_prior);
            if (fin_prior.empty()) return;
            RooUnfoldResponse *r_prior = sumResponseOverFiles(fin_prior,
                obs.n(split_test ? "response_tf_pseudo_full" : "response_tf_full"));
            if (!r_prior) return;
            if (r_prior->GetNbinsTruth() != response->GetNbinsTruth()) {
                std::cerr << "ERROR: prior and response truth binnings differ" << std::endl; return; }
            std::cout << "\t---->Unfolding prior: the " << prior_gen << " truth (" << filenames_prior.size()
                      << " files), detector response: " << generator << std::endl;
            response = reweightPrior(response, r_prior->Htruth());
        }
    }



    // ---- Print condition number
    // TDecompSVD *svd= new TDecompSVD(response->Mresponse());  // response is a RooUnfold response object, svd is the singular value decomposition (SVD) matrix. the response->Mresponse() returns the normalized migration matrix
    // 2026-10-29: the response can be non-square (JET_PT_MAX: a gen overflow bin, more truth than reco
    // bins); TDecompSVD needs rows >= columns, and the transpose has the same singular values.
    const TMatrixD &m_resp = response->Mresponse();
    const bool resp_square = (m_resp.GetNrows() == m_resp.GetNcols());
    TDecompSVD *svd = (m_resp.GetNrows() >= m_resp.GetNcols()) ? new TDecompSVD(m_resp)
                                                               : new TDecompSVD(TMatrixD(TMatrixD::kTransposed, m_resp));
    auto singular_values = svd->GetSig(); //this is a vector with the singular values, i.e., the diagonal elements of S. They are ordered from largest to smallest.
    // double cond_number = singular_values.Max() / singular_values.Min();   // until 2026-10-26
    // CMS-SMP-22-007 Sec. 5 (2026-10-26): largest over smallest NONZERO value of the probability
    // matrix -- Mresponse(), whose columns (truth bins) sum to 1 over the detector bins here, since
    // the response has no misses. An empty truth bin gives a zero, which is now skipped instead of
    // making the number infinite. The paper says "eigenvalues"; the matrix is not symmetric, so the
    // eigenvalue version (moduli of the complex eigenvalues) is printed next to the singular-value
    // one, which stays on the plots until the authors say which they meant.
    double sv_max = singular_values.Max(), sv_min = 0.;
    for (int k = 0; k < singular_values.GetNrows(); ++k) {
        const double s = singular_values(k);
        if (s > 1e-12 * sv_max && (sv_min == 0. || s < sv_min)) sv_min = s;
    }
    const double cond_number = sv_min > 0. ? sv_max / sv_min : -1.;
    double ev_max = 0., ev_min = 0.;
    if (resp_square) {   // eigenvalues exist only for a square matrix (2026-10-29)
        TMatrixDEigen eig(response->Mresponse());
        const TVectorD re = eig.GetEigenValuesRe(), im = eig.GetEigenValuesIm();
        for (int k = 0; k < re.GetNrows(); ++k) ev_max = std::max(ev_max, std::hypot(re(k), im(k)));
        for (int k = 0; k < re.GetNrows(); ++k) {
            const double m = std::hypot(re(k), im(k));
            if (m > 1e-12 * ev_max && (ev_min == 0. || m < ev_min)) ev_min = m;
        }
    }
    std::cout << "\t---->Condition number nominal = " << cond_number
              << "  (singular values; eigenvalue moduli: " << (ev_min > 0. ? ev_max / ev_min : -1.) << ")"
              << std::endl;

    // ---- Grab the truth level MC ----
    // Must be the SAME production as the response files: the truth, the response and the
    // purity/efficiency corrections have to come from one set of events with one pT
    // binning, otherwise the full-MC closure (test_mode 0) cannot be exactly 1.
    std::vector<TFile *> &fin_response_truth = fin_unfolding;
    std::cout << "Getting truth from the same " << fin_response_truth.size()
              << " response file(s)" << std::endl;
    TH2D *h_mc_true = nullptr;
    if (test_mode == 1) {
        // ⚠️ THE TRUTH OF THE PSEUDODATA, NOT THE TRUTH THAT FILLED THE MATRIX.
        // h_pseudodata_truth_tf is filled in create_files_for_template_fit.cpp under
        // `ient % 2 == 0` -- the EVEN half, the same events that fill h3D_pseudodata_bb,
        // which is what this mode unfolds. The response and the purity/efficiency ratios are
        // the ODD half (`ient % 2 == 1`, the ps_* histograms). So the closure asked here is
        // "does the odd-half response recover the even half's own gen distribution", which is
        // a real independent test. Comparing against the odd half's gen instead
        // (h_full_pseudo_efficiency_denominator_tf) would be testing the unfolding against
        // the prior it was built from -- it would close better and mean nothing.
        // Both halves carry the same gen_pass gate and the same w_reco weight, so the
        // normalisations are directly comparable; they are ~half the sample each, so a
        // half-swap would NOT show up as a gross scale error. It has to be read off the name.
        h_mc_true = sumOverFiles<TH2D>(fin_response_truth, obs.n("h_pseudodata_truth_tf"));
    } else if (test_mode == 0) {
        // Full-MC closure: the full-sample all-gen distribution the chain recovers by construction.
        h_mc_true = sumOverFiles<TH2D>(fin_response_truth, obs.n("h_full_efficiency_denominator_tf"));
    } else {
        // Data: gen reference. After the combined SV-reco + b-tag correction below, the
        // result is at the "all true 2b" level, so compare against hgenjet_2b_all (the
        // combined-efficiency denominator), not hgenjet_2b_passbtag (2SV+btag level).
        h_mc_true = sumOverFiles<TH2D>(fin_response_truth, obs.n("hgenjet_2b_all"));
    }
    if (!h_mc_true) return;


    //------- Apply purity correction
    TH2D *h_data_purity_corrected = (TH2D *) h_data_after_fit->Clone("h_data_purity_corrected");
    if (apply_purity) {
        // applyCorrection(h_data_purity_corrected, h_full_purity, /*divide=*/false, "purity");   // until 2026-10-26
        if (!applyCorrection(h_data_purity_corrected, h_full_purity, /*divide=*/false, "purity")) return;
    } else {
        std::cout << "\t---->NOT multiplying data by purity" << std::endl;
    }

    //double underflowY = h_data_purity_corrected->Integral(0, h_data_purity_corrected->GetNbinsX()+1, 0, 0);
    //double overflowY  = h_data_purity_corrected->Integral(0, h_data_purity_corrected->GetNbinsX()+1,
                          //h_data_purity_corrected->GetNbinsY()+1, h_data_purity_corrected->GetNbinsY()+1);
    //std::cout << "underflowY=" << underflowY << " overflowY=" << overflowY << std::endl;

    // ---- Scan the Bayesian iteration count, one PDF page per value (only when scan_niter).
    // ---- disabled (kept for reference): iteration scan output in the LLR home dir ----
    // TString home_llr = "/home/llr/cms/zaidan/";
    // label already encodes the option flags (_bayesian / _split_test / _nopurity), so toggling
    // those bools writes to distinct files instead of overwriting the previous run's outputs.
    TString iter_pdf = folder + "bayesian_unfolding_iterations_" + label + ".pdf";  // multi-page PDF of all iterations
    TString png_dir  = folder + "bayesian_unfolding_iterations_" + label + "_png/"; // one PNG per iteration lands here
    if (scan_niter) gSystem->mkdir(png_dir, kTRUE);                     // create the PNG folder if it does not exist
    // 2026-11-04: a shorter scan left the pages of a longer one behind (dR 51..200 of 10-28 next to 1..50), and they were staged
    if (scan_niter) gSystem->Exec("rm -f " + png_dir + "iteration_*.png");
    const int niter_min = 0;
    // 1..100 since 2026-10-01 (was 1..50): B does not reach a refolding p of 0.95 by 50.
    // const int niter_max = scan_niter ? 99 : 0;   // scan niter = 1..niter_max+1; one pass when not scanning
    // const int niter_max = scan_niter ? 199 : 0;   // 1..200 since 2026-10-19 (Zoe)
    // SCAN_NITER_MAX=N (env, 2026-10-22) scans 1..N instead, e.g. 150 for the talk plots.
    // Unset = 200, as before. Only the range changes: the choice is the FIRST p > 0.05
    // (CMS-SMP-22-007, since 2026-10-26; the first p >= 0.95 before).
    const char *scan_max_env = std::getenv("SCAN_NITER_MAX");
    // const int niter_max = scan_niter ? ((scan_max_env && std::atoi(scan_max_env) > 0) ? std::atoi(scan_max_env) - 1 : 199) : 0;
    // 2026-10-28 (Zoe): dR scans 1..50 by default (its niter is 16), A keeps 1..200 (126); env still wins.
    const int scan_default = (obs.name == "dr") ? 50 : 200;
    const int niter_max = scan_niter ? ((scan_max_env && std::atoi(scan_max_env) > 0) ? std::atoi(scan_max_env) - 1 : scan_default - 1) : 0;
    const Int_t prev_ignore = gErrorIgnoreLevel;
    if (scan_niter) gErrorIgnoreLevel = kError;  // silence "Replacing existing histogram" spam from re-cloning each pass

    // Refolding goodness-of-fit vs iteration (see refoldChi2 above). Filled on every pass;
    // only used for the regularisation plot when scan_niter is on.
    // The optimal iteration is the FIRST one whose refolding p-value clears this threshold:
    // earlier iterations are still non-closing, later ones only overfit the input further.
    const double pvalue_threshold = 0.05;
    std::vector<double> v_niter, v_pvalue_all, v_pvalue_pt, v_chi2ndf_all, v_chi2ndf_pt;
    // The same test with the response MC statistics included (see refoldChi2).
    std::vector<double> v_pvalue_all_mc, v_chi2ndf_all_mc;
    std::vector<double> v_pvalue_all_cov, v_chi2ndf_all_cov;
    std::vector<double> v_pvalue_all_r, v_chi2ndf_all_r;
    std::vector<double> v_pvalue_all_rcov, v_chi2ndf_all_rcov;
    const TH2 *h_mc_meas_gof = dynamic_cast<const TH2 *>(response->Hmeasured());
    // The unfolded result in the measured pT bin at every iteration (scan only), for the
    // bias-vs-variance choice of niter after the loop (see "Bias vs statistical error" below).
    std::vector<TH1D *> v_unf_pt;
    std::vector<TH1D *> v_fc_pt;   // the fully corrected result (SFs, efficiencies) per iteration, for the iteration test
    // chi2 bottom-line test per iteration (2026-10-29, Zoe, AN Fig. 68): unfolded data vs gen MC, reco data vs reco MC
    std::vector<double> v_bl_unf, v_bl_reco;
    TMatrixD cov_unf_now;   // covariance of the unfolded result at the current iteration (scan only)

    for (int iter = niter_min; iter <= niter_max; ++iter) {

    // Single non-scan pass keeps the nominal 7 iterations; the scan walks 1..100.
    // scan_niter walks 1..100; otherwise use niter_fixed. The old hardcoded 7 is the default,
    // so an unchanged call is unchanged -- but for a new observable the right value has to be
    // found with scan_niter=true first. For the balance the scan says 21-26; 7 refolds at
    // chi2/ndf = 16 (p = 2e-29), i.e. badly under-iterated and still close to the prior.
    //
    // ⚠️ A scan_niter run WRITES THE LAST PASS (niter_max + 1, 200 unless SCAN_NITER_MAX), not
    // the optimum it reports -- to the "_iterative" file (out_label), never the fixed run's. Use
    // the scan to find the number, then re-run with scan_niter=false and niter_fixed set.
    const Int_t niter_now = scan_niter ? iter + 1 : niter_fixed;

    // ---- Unfold
    // Say WHICH unfolding: an "niter" in the log of a matrix inversion reads as a Bayesian
    // run at that iteration count, and the two results are not the same measurement.
    std::cout << (unfoldBayes ? Form("\t---->Unfolding, Bayesian (niter = %d)", (int) niter_now)
                              : "\t---->Unfolding, matrix inversion (no regularisation)")
              << std::endl;
    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
    TH2D *h_data_unfolded = nullptr;
    // ---- disabled (kept for reference): the two covariance matrices ----
    // Nothing downstream ever read them, and RooUnfold hands back matrices whose dimension
    // is not dim = bins_pt*bins_dr, so each assignment printed
    //     Error in <operator=(const TMatrixT &)>: matrices not compatible
    // and left the matrix as it was. Harmless, but it is an "Error" in every log, and matrix
    // inversion is now the default -- a real inversion problem has to be visible in that log.
    // The per-bin errors that ARE used come from Hreco(kCovariance) just below, which
    // computes the same covariance internally.
    // TMatrixD covariance_matrix_before_unfolding(dim,dim);
    // TMatrixD covariance_matrix_after_unfolding(dim,dim);
    if (unfoldBayes) {
        RooUnfoldBayes unfold(response, h_data_purity_corrected, niter_now);
        // Clone before `unfold` leaves scope: some RooUnfold versions hand back a cached
        // histogram that the unfolder owns and deletes with itself.
        h_data_unfolded = (TH2D *) unfold.Hreco(errorTreatment)->Clone("h_data_unfolded");
        if (scan_niter) {   // for the chi2 bottom-line test (resize first: TMatrixD::operator= needs equal shapes)
            const TMatrixD e = unfold.Ereco(RooUnfold::kCovariance);
            cov_unf_now.ResizeTo(e); cov_unf_now = e;
        }
        // covariance_matrix_before_unfolding = unfold.GetMeasuredCov();
        // covariance_matrix_after_unfolding = unfold.Ereco();
    }
    /* ---- disabled 2026-10-27 (kept for reference): matrix inversion ----
    // Moving to a 3D (jet pT, dR, A) Bayesian unfolding with a NON-SQUARE response that keeps
    // the under/overflow bins. RooUnfoldInvert needs a square, invertible matrix, so it has no
    // place there; unfoldBayes = false is refused in apply_unfolding_2d() instead.
    else {
        RooUnfoldInvert unfold(response, h_data_purity_corrected);
        h_data_unfolded = (TH2D *) unfold.Hreco(errorTreatment)->Clone("h_data_unfolded");
        // covariance_matrix_before_unfolding = unfold.GetMeasuredCov();
        // covariance_matrix_after_unfolding = unfold.Ereco();
    }
    ---- end disabled ---- */

    // ---- Fold back
    std::cout << "\t---->Refolding" << std::endl;
    TH2D *h_data_refolded = (TH2D *) response->ApplyToTruth(h_data_unfolded, "h_data_refolded");

    // ---- Refolding goodness of fit, against the distribution that went INTO the unfolding
    // (h_data_purity_corrected -- the same object handed to RooUnfoldBayes above), not the
    // raw h_data_reco: the refolded distribution lives at the purity-corrected level too.
    // "all" = the whole 2D reco space that was unfolded; "pt" = the pT bin actually plotted.
    const RefoldGof gof_all = refoldChi2(h_data_refolded, h_data_purity_corrected,
                                         ibin_dr_min, ibin_dr_max, 1, bins_pt);
    const RefoldGof gof_pt  = refoldChi2(h_data_refolded, h_data_purity_corrected,
                                         ibin_dr_min, ibin_dr_max, ibin_pt, ibin_pt);
    // "niter = " only where there are iterations; matrix inversion would otherwise report the
    // untouched niter_fixed and read as a Bayesian run at that value. Built as its own
    // TString, NOT a nested Form(): Form() returns a pointer into one static buffer, so a
    // Form() inside a Form() has the outer call overwrite what the inner one returned.
    const TString reg_label = unfoldBayes ? TString::Format("niter = %3d", (int) niter_now)
                                          : TString("inversion");
    std::cout << Form("\t---->Refolding GoF  %-11s | 2D: chi2/ndf = %8.2f/%3d = %6.3f, "
                      "p = %8.3e | pT bin %d: chi2/ndf = %8.2f/%3d = %6.3f, p = %8.3e",
                      reg_label.Data(),
                      gof_all.chi2, gof_all.ndf, gof_all.chi2ndf(), gof_all.pvalue,
                      ibin_pt, gof_pt.chi2, gof_pt.ndf, gof_pt.chi2ndf(), gof_pt.pvalue)
              << std::endl;
    v_niter.push_back(niter_now);
    v_pvalue_all.push_back(gof_all.pvalue);
    v_pvalue_pt.push_back(gof_pt.pvalue);
    v_chi2ndf_all.push_back(gof_all.chi2ndf());
    v_chi2ndf_pt.push_back(gof_pt.chi2ndf());
    {
        const RefoldGof gof_mc = refoldChi2(h_data_refolded, h_data_purity_corrected,
                                            ibin_dr_min, ibin_dr_max, 1, bins_pt, h_mc_meas_gof);
        v_pvalue_all_mc.push_back(gof_mc.pvalue);
        v_chi2ndf_all_mc.push_back(gof_mc.chi2ndf());
        std::cout << Form("\t---->Refolding GoF, data (+) MC stat | 2D: chi2/ndf = %8.2f/%3d = %6.3f, p = %8.3e",
                          gof_mc.chi2, gof_mc.ndf, gof_mc.chi2ndf(), gof_mc.pvalue) << std::endl;
    }
    if (unfoldBayes) {
        RefoldGof gof_r, gof_rcov;
        const RefoldGof gof_cov = refoldChi2Toys(response, h_data_refolded, h_data_purity_corrected,
                                                 niter_now, ibin_dr_min, ibin_dr_max, 1, bins_pt,
                                                 1000, 20261019, &gof_r, &gof_rcov);
        v_pvalue_all_rcov.push_back(gof_rcov.pvalue);
        v_chi2ndf_all_rcov.push_back(gof_rcov.chi2ndf());
        std::cout << Form("\t---->Refolding GoF, data (+) refolded stat, full cov | 2D: chi2/ndf = %8.2f/%3d = %6.3f, p = %8.3e",
                          gof_rcov.chi2, gof_rcov.ndf, gof_rcov.chi2ndf(), gof_rcov.pvalue) << std::endl;
        v_pvalue_all_cov.push_back(gof_cov.pvalue);
        v_chi2ndf_all_cov.push_back(gof_cov.chi2ndf());
        v_pvalue_all_r.push_back(gof_r.pvalue);
        v_chi2ndf_all_r.push_back(gof_r.chi2ndf());
        std::cout << Form("\t---->Refolding GoF, data (+) refolded stat | 2D: chi2/ndf = %8.2f/%3d = %6.3f, p = %8.3e",
                          gof_r.chi2, gof_r.ndf, gof_r.chi2ndf(), gof_r.pvalue) << std::endl;
        std::cout << Form("\t---->Refolding GoF, cov(r - d) toys, data stat | 2D: chi2/ndf = %8.2f/%3d = %6.3f, p = %8.3e",
                          gof_cov.chi2, gof_cov.ndf, gof_cov.chi2ndf(), gof_cov.pvalue) << std::endl;
    }
    if (scan_niter) {
        // x is dr, y is pt: the measured pT bin, with the covariance-based errors of Hreco.
        TH1D *u = h_data_unfolded->ProjectionX(Form("h_unf_pt_it%d", (int) niter_now),
                                               ibin_pt, ibin_pt, "e");
        u->SetDirectory(nullptr);
        v_unf_pt.push_back(u);
    }

    // ---- Apply efficiency correction
    TH2D *h_data_efficiency_corrected = (TH2D *) h_data_unfolded->Clone("h_data_efficiency_corrected");
    // applyCorrection(h_data_efficiency_corrected, h_full_efficiency, /*divide=*/true,
    //                 "reconstruction efficiency");   // until 2026-10-26: result written anyway on failure
    if (!applyCorrection(h_data_efficiency_corrected, h_full_efficiency, /*divide=*/true,
                         "reconstruction efficiency")) { gErrorIgnoreLevel = prev_ignore; return; }

    // ---- Final corrections
    TH2D *h_data_fully_corrected = (TH2D *) h_data_efficiency_corrected->Clone("h_data_fully_corrected");

    // ---- Combined SV-reconstruction + b-tag efficiency correction (data only) ----
    // Peels BOTH the 2-SV reconstruction and the b-tag selection off the truth level,
    // taking the result from "2b jets reconstructed with 2 SV and tagged" to "all true 2b".
    // Orthogonal to the reconstruction-kinematic efficiency (h_full_efficiency) applied above.
    // Skipped in closure modes (0/1), whose truth targets live at the 2SV+btag level.
    if (test_mode == 2) {
        TH2D *h_svbtag_eff = ratioFromCounts(obs.n("hgenjet_2b_reco_btag"), obs.n("hgenjet_2b_all"), "b");
        // if (h_svbtag_eff) {   // until 2026-10-26: a missing correction was skipped and the result written
        //     applyCorrection(h_data_fully_corrected, h_svbtag_eff, /*divide=*/true,
        //                     "combined SV-reco + b-tag efficiency");
        // }
        // systematic hook: scale h_svbtag_eff by the CMS b-tag SF map here, then vary SF +/-.
        if (!applyCorrection(h_data_fully_corrected, h_svbtag_eff, /*divide=*/true,
                             "combined SV-reco + b-tag efficiency")) { gErrorIgnoreLevel = prev_ignore; return; }
        // EEC-weight correction: convert the reco-EEC-weighted result to gen-EEC-weighted.
        // r_eec = sum(eec_gen)/sum(eec_reco) over reconstructed 2b jets -> MULTIPLY.
        //
        // THE EEC MEASUREMENT ONLY. Its two inputs are the same jets filled with w_reco =
        // weight * eec_reco and w_gen = weight * eec_gen, so with the weight off
        // (EecWeight::value() returns 1) they are the SAME histogram and the ratio is exactly
        // 1 -- verified bit-for-bit on one of the (now disabled) fraction productions, all
        // 10 cells. Applying it to a
        // yield was therefore harmless but wrong to show: the correction-stages plot drew a
        // "+ EEC weight" stage sitting exactly on the one before it, and the log announced a
        // correction that corrected nothing. A yield has no EEC weight to correct back to.
        TH2D *h_eec_weight_eff = eec_weight_off
                               ? nullptr
                               : ratioFromCounts(obs.n("hgenjet_2b_reco_btag"),
                                                 obs.n("hgenjet_2b_passbtag"), "");
        // if (h_eec_weight_eff) {   // until 2026-10-26: a missing EEC-weight correction was skipped silently
        //     applyCorrection(h_data_fully_corrected, h_eec_weight_eff, /*divide=*/false,
        //                     "EEC-weight (reco->gen) correction");
        // } else if (eec_weight_off) {
        if (!eec_weight_off) {
            if (!applyCorrection(h_data_fully_corrected, h_eec_weight_eff, /*divide=*/false,
                                 "EEC-weight (reco->gen) correction")) { gErrorIgnoreLevel = prev_ignore; return; }
        } else {
            std::cout << "\t---->No EEC-weight correction: the weight is off, this is a yield"
                      << std::endl;
        }
    }
    // ---- chi2 bottom-line test (2026-10-29, Zoe; AN Fig. 68) in the measured pT bin ----
    // chi2_unf  = (u - s g)^T V^-1 (u - s g): u = fully corrected data, V = the unfolding covariance
    //             propagated through the per-bin corrections (c_i c_j V_ij, c = corrected / unfolded)
    //             + the gen-MC stat on the diagonal; g = gen MC scaled to the data integral (s).
    // chi2_reco = sum (r - s' m)^2 / (sr^2 + s'^2 sm^2): r = reco data (purity-corrected, the unfolding
    //             input), m = matched reco MC scaled to r. Shapes only: both MC are scaled to the data.
    // The test passes when chi2_unf <= chi2_reco (the unfolding does not create (dis)agreement).
    if (scan_niter && h_mc_true && h_mc_reco && cov_unf_now.GetNrows() > 0) {
        const int nb = ibin_dr_max - ibin_dr_min + 1, nxu = h_data_unfolded->GetNbinsX();
        double su = 0., sg = 0., sr = 0., sm = 0.;
        for (int i = ibin_dr_min; i <= ibin_dr_max; ++i) {
            su += h_data_fully_corrected->GetBinContent(i, ibin_pt); sg += h_mc_true->GetBinContent(i, ibin_pt);
            sr += h_data_purity_corrected->GetBinContent(i, ibin_pt); sm += h_mc_reco->GetBinContent(i, ibin_pt);
        }
        const double s_g = sg > 0. ? su / sg : 0., s_m = sm > 0. ? sr / sm : 0.;
        TMatrixD V(nb, nb); TVectorD d(nb);
        for (int a = 0; a < nb; ++a) {
            const int i = ibin_dr_min + a;
            const double ua = h_data_unfolded->GetBinContent(i, ibin_pt);
            const double ca = ua != 0. ? h_data_fully_corrected->GetBinContent(i, ibin_pt) / ua : 0.;
            d[a] = h_data_fully_corrected->GetBinContent(i, ibin_pt) - s_g * h_mc_true->GetBinContent(i, ibin_pt);
            for (int b = 0; b < nb; ++b) {
                const int j = ibin_dr_min + b;
                const double ub = h_data_unfolded->GetBinContent(j, ibin_pt);
                const double cb = ub != 0. ? h_data_fully_corrected->GetBinContent(j, ibin_pt) / ub : 0.;
                V(a, b) = ca * cb * cov_unf_now((i - 1) + nxu * (ibin_pt - 1), (j - 1) + nxu * (ibin_pt - 1));
            }
            V(a, a) += std::pow(s_g * h_mc_true->GetBinError(i, ibin_pt), 2);
        }
        double chi2_unf = -1.;
        TDecompSVD svd(V);
        if (svd.Decompose()) { TVectorD x = d; Bool_t ok = kFALSE; x = svd.Solve(x, ok); if (ok) chi2_unf = d * x; }
        double chi2_reco = 0.;
        for (int i = ibin_dr_min; i <= ibin_dr_max; ++i) {
            const double r = h_data_purity_corrected->GetBinContent(i, ibin_pt), m = h_mc_reco->GetBinContent(i, ibin_pt);
            const double v = std::pow(h_data_purity_corrected->GetBinError(i, ibin_pt), 2) + std::pow(s_m * h_mc_reco->GetBinError(i, ibin_pt), 2);
            if (v > 0.) chi2_reco += std::pow(r - s_m * m, 2) / v;
        }
        v_bl_unf.push_back(chi2_unf); v_bl_reco.push_back(chi2_reco);
        std::cout << Form("\t---->Bottom-line chi2 (pT bin %d, %d bins): unfolded vs gen %.2f, reco vs reco MC %.2f, ratio %.3f",
                          ibin_pt, nb, chi2_unf, chi2_reco, chi2_reco > 0. ? chi2_unf / chi2_reco : -1.) << std::endl;
    }

    if (scan_niter) {   // 2026-10-28 (Zoe): the iteration test on the nominal result, SFs applied
        TH1D *u = h_data_fully_corrected->ProjectionX(Form("h_fc_pt_it%d", (int) niter_now), ibin_pt, ibin_pt, "e");
        u->SetDirectory(nullptr);
        v_fc_pt.push_back(u);
    }


    // ---- Correction-stage plot: raw EEC (top) + incremental ratios (bottom), data only ----
    // Top pad: EEC(dr) after unfolding with each correction stacked on (absolute).
    // Bottom pad: each correction isolated as the fractional effect (ratio - 1); 0 = no effect.
    // Same main/ratio layout, palette and fonts as the bottomline plot.
    if (test_mode == 2) {
        // normalise_corr_stages = true  -> every curve to unit area; the ratio pad shows the
        //   SHAPE effect of each correction. Writes ..._norm.{pdf,png}.
        // false -> absolute EEC; the ratio pad shows the true SIZE of each correction, and the
        //   gen MC is scaled to the data integral as a shape reference. Writes unsuffixed files.
        // Currently OFF, because normalising divides out exactly what this plot is meant to show.
        // To restore it set this true AND uncomment the normalizeToUnitArea loop below.
        const bool normalise_corr_stages = false;
        std::cout << "\t---->Making correction-stage EEC + ratio plot ("
                  << (normalise_corr_stages ? "normalised" : "absolute") << ")" << std::endl;
        TH2D *h_btag_c = ratioFromCounts(obs.n("hgenjet_2b_passbtag"),  obs.n("hgenjet_2b"),          "b"); // b-tag only, after 2SV
        TH2D *h_svb_c  = ratioFromCounts(obs.n("hgenjet_2b_reco_btag"), obs.n("hgenjet_2b_all"),      "b"); // 2SV + b-tag combined
        // Null for a yield run, exactly as above: no EEC weight, no EEC-weight stage.
        TH2D *h_eec_c  = eec_weight_off ? nullptr
                       : ratioFromCounts(obs.n("hgenjet_2b_reco_btag"), obs.n("hgenjet_2b_passbtag"), "");  // EEC weight

        // 2SV alone. No count pair gives it directly (hgenjet_2b / hgenjet_2b_all mixes in the
        // reco/gen EEC-weight ratio, different weights), so factorise the combined efficiency:
        //   eff(2SV) = eff(2SV + b-tag) / eff(b-tag) = h_svb_c / h_btag_c.
        // The errors treat the two ratios as independent when they are not; fine for a diagnostic.
        TH2D *h_2sv_c = nullptr;
        if (h_svb_c && h_btag_c) {
            h_2sv_c = (TH2D*) h_svb_c->Clone("ratio_2sv_only");
            h_2sv_c->SetDirectory(0);
            h_2sv_c->Divide(h_btag_c);
        }

        // Stage order: after unfolding -> + 2SV -> + 2SV + b-tag -> + EEC weight. The combined
        // stage is the same histogram whichever cut is peeled first (h_svb_c covers both), so
        // reordering only changes which single correction is applied on its own in between.
        TH2D *h2_base   = (TH2D*) h_data_efficiency_corrected->Clone("h2_cmp_base");
        TH2D *h2_2sv    = (TH2D*) h_data_efficiency_corrected->Clone("h2_cmp_2sv");
        TH2D *h2_svbtag = (TH2D*) h_data_efficiency_corrected->Clone("h2_cmp_svbtag");
        TH2D *h2_full   = (TH2D*) h_data_efficiency_corrected->Clone("h2_cmp_full");
        if (h_2sv_c) h2_2sv->Divide(h_2sv_c);
        if (h_svb_c)  h2_svbtag->Divide(h_svb_c);
        if (h_svb_c)  h2_full->Divide(h_svb_c);
        if (h_eec_c)  h2_full->Multiply(h_eec_c);

        TH1D *p_base   = h2_base  ->ProjectionX(Form("p_cmp_base_pt%d",   ibin_pt), ibin_pt, ibin_pt);
        TH1D *p_2sv    = h2_2sv   ->ProjectionX(Form("p_cmp_2sv_pt%d",    ibin_pt), ibin_pt, ibin_pt);
        TH1D *p_svbtag = h2_svbtag->ProjectionX(Form("p_cmp_svbtag_pt%d", ibin_pt), ibin_pt, ibin_pt);
        TH1D *p_full   = h2_full  ->ProjectionX(Form("p_cmp_full_pt%d",   ibin_pt), ibin_pt, ibin_pt);
        // Gen MC at the level the fully corrected result targets -- for test_mode 2 that is
        // hgenjet_2b_all, "all true 2b". This is the curve p_full is meant to land on, so the
        // plot shows not just what each correction does but whether they arrive somewhere sane.
        TH1D *p_true   = h_mc_true->ProjectionX(Form("p_cmp_true_pt%d",   ibin_pt), ibin_pt, ibin_pt);

        int b1 = 1, b2 = p_base->GetNbinsX();
        while (b1 <= p_base->GetNbinsX() && p_base->GetBinContent(b1) <= 0.) ++b1;
        while (b2 >= 1                   && p_base->GetBinContent(b2) <= 0.) --b2;

        // Whatever the mode, do this BEFORE the ratios below so they inherit it.
        for (TH1D* h : {p_base, p_2sv, p_svbtag, p_full, p_true}) h->GetXaxis()->SetRange(b1, b2);
        if (normalise_corr_stages) {
            // Unit-area normalisation of every stage -- commented out, see normalise_corr_stages above.
            //for (TH1D* h : {p_base, p_2sv, p_svbtag, p_full, p_true}) normalizeToUnitArea(h);
        } else {
            // The gen MC cannot be overlaid absolutely: it carries an MC event weight with no
            // luminosity scaling, so its integral has no relation to the data's. Put it on the
            // data's scale instead -- shape reference only, and the legend says so.
            const double I_data_cmp = p_full->Integral(), I_true_cmp = p_true->Integral();
            if (I_data_cmp > 0. && I_true_cmp > 0.) p_true->Scale(I_data_cmp / I_true_cmp);
        }

        // ratio of each stage to the PREVIOUS one, isolating that single correction:
        //   2SV = 2sv/after-unfolding, +b-tag = svbtag/2sv, +EEC = full/svbtag.  1 = no change.
        TH1D *r_2sv  = (TH1D*) p_2sv   ->Clone(Form("r_2sv_pt%d",  ibin_pt));  r_2sv ->Divide(p_base);
        TH1D *r_btag = (TH1D*) p_svbtag->Clone(Form("r_btag_pt%d", ibin_pt));  r_btag->Divide(p_2sv);
        TH1D *r_eec  = (TH1D*) p_full  ->Clone(Form("r_eec_pt%d",  ibin_pt));  r_eec ->Divide(p_svbtag);

        // ROCColor palette from the bottomline plot; each step keeps its colour top & bottom.
        // Red is the measurement (p_full, the fully corrected result) and blue the gen MC it is
        // compared to, as everywhere else; the intermediate stages take the remaining colours,
        // in the order they are applied (first added correction purple, second orange).
        p_base  ->SetLineColor(green);  p_base  ->SetMarkerColor(green);  p_base  ->SetMarkerStyle(kOpenCircle);
        p_2sv   ->SetLineColor(purple); p_2sv   ->SetMarkerColor(purple); p_2sv   ->SetMarkerStyle(kFullCircle);
        p_svbtag->SetLineColor(orange); p_svbtag->SetMarkerColor(orange); p_svbtag->SetMarkerStyle(kFullSquare);
        p_full  ->SetLineColor(red);    p_full  ->SetMarkerColor(red);    p_full  ->SetMarkerStyle(kFullTriangleUp);
        p_true  ->SetLineColor(blue);   p_true  ->SetMarkerColor(blue);   p_true  ->SetMarkerStyle(kOpenTriangleUp);
        r_2sv ->SetLineColor(purple); r_2sv ->SetMarkerColor(purple); r_2sv ->SetMarkerStyle(kFullCircle);
        r_btag->SetLineColor(orange); r_btag->SetMarkerColor(orange); r_btag->SetMarkerStyle(kFullSquare);
        r_eec ->SetLineColor(red);    r_eec ->SetMarkerColor(red);    r_eec ->SetMarkerStyle(kFullTriangleUp);
        for (TH1D* h : {p_base, p_2sv, p_svbtag, p_full, p_true, r_2sv, r_btag, r_eec}) {
            h->SetStats(0); h->SetLineWidth(2); h->SetMarkerSize(1.2); h->GetXaxis()->SetRange(b1, b2);
        }

        double ymax = 0.;
        for (TH1D* h : {p_base, p_2sv, p_svbtag, p_full, p_true})
            for (int i = b1; i <= b2; ++i) ymax = std::max(ymax, h->GetBinContent(i) + h->GetBinError(i));
        double rmin = 1., rmax = 1.;
        for (TH1D* h : {r_2sv, r_btag, r_eec})
            for (int i = b1; i <= b2; ++i) { double v = h->GetBinContent(i); if (v <= 0.) continue; rmin = std::min(rmin, v); rmax = std::max(rmax, v); }
        double rpad = 0.15 * (rmax - rmin) + 1e-6;

        Float_t ptlo = h_data_reco->GetYaxis()->GetBinLowEdge(ibin_pt);
        Float_t pthi = h_data_reco->GetYaxis()->GetBinUpEdge(ibin_pt);

        // canvas + two pads (same main/ratio split as the bottomline plot)
        TCanvas *c_cmp = new TCanvas(Form("c_corr_stages_pt%d", ibin_pt), "", 800, 800);
        TPad *pad_main  = new TPad(Form("pad_cmp_main_%d",  ibin_pt), "", 0., 0.3, 1., 1.);
        TPad *pad_ratio = new TPad(Form("pad_cmp_ratio_%d", ibin_pt), "", 0., 0.,  1., 0.3);
        for (TPad *p : {pad_main, pad_ratio}) { p->SetTicks(1, 0); p->SetFillColor(0); }
        // pad_main ->SetMargin(0.13, 0.05, 0.00, 0.08);
        // pad_ratio->SetMargin(0.13, 0.05, 0.32, 0.00);
        pad_main->SetMargin(0.1, 0.1, 0.0, 0.1);
        // Bottom margin 0.30, not 0.23: the x title has to fit UNDER the tick labels, and B's
        // label (p_{T}^{b1}/(p_{T}^{b1} + p_{T}^{b2})) is two storeys tall -- at 0.23 its
        // subscripts were cut off by the canvas edge. #Delta r is unaffected either way.
        pad_ratio->SetMargin(0.1, 0.1, 0.30, 0.0);


        c_cmp->cd(); pad_main->Draw(); pad_ratio->Draw();

        // ----- top: raw EEC(dr) with each correction -----
        pad_main->cd();
        p_base->SetTitle("");
        p_base->GetYaxis()->SetRangeUser(0., ymax * 1.6);
        p_base->GetYaxis()->SetTitle(normalise_corr_stages ? "normalised " + obs_title
                                                          : obs_title);
        p_base->GetYaxis()->CenterTitle(true);
        p_base->GetYaxis()->SetTitleFont(font_code); p_base->GetYaxis()->SetTitleSize(title_size); p_base->GetYaxis()->SetTitleOffset(1.5);
        p_base->GetYaxis()->SetLabelFont(font_code); p_base->GetYaxis()->SetLabelSize(label_size);
        // The y range starts at 0 and this pad's bottom edge is shared with the ratio pad, so the
        // "0" label sits exactly on the join and gets clipped in half. Hide it (label size 0).
        p_base->GetYaxis()->ChangeLabel(1, -1, 0.);
        p_base->GetXaxis()->SetLabelSize(0); p_base->GetXaxis()->SetTitleSize(0);   // x belongs to the ratio pad
        // Every stage, "after unfolding" included, as points (X0 = no horizontal / bin-width
        // error bars), so the starting point reads as one of the stages rather than as a frame.
        p_base->SetLineWidth(2);
        // Gen MC as an outline with its error band, like the bottomline plot, so it reads as the
        // reference rather than as another correction stage.
        fixXRange(p_base);
        p_base->Draw("PE X0"); p_true->Draw("HIST E SAME");
        p_2sv->Draw("PE X0 SAME"); p_svbtag->Draw("PE X0 SAME");
        // p_full only exists as a separate stage when there is an EEC weight to correct: in a
        // yield run h_eec_c is null, so p_full IS p_svbtag and drawing it again would put a
        // second marker on top of the first and claim a stage that did nothing.
        if (h_eec_c) p_full->Draw("PE X0 SAME");

        TLegend *lg = new TLegend(0.17, 0.60, 0.60, 0.85);
        lg->SetFillStyle(0); lg->SetBorderSize(0); lg->SetMargin(0.15);
        lg->SetTextFont(font_code); lg->SetTextSize(legend_size * 0.8);
        // lg->SetHeader(pthi == jtpt_max ? Form("p_{T}^{jet} > %.0f GeV", ptlo)
        //                                : Form("%.0f < p_{T}^{jet} < %.0f GeV", ptlo, pthi));
        lg->SetHeader(jetPtRangeLabel(ptlo, pthi));
        lg->AddEntry(p_base,   "after unfolding",            "pe1");
        lg->AddEntry(p_2sv,    "+ 2SV eff.",                 "pe1");
        lg->AddEntry(p_svbtag, "+ 2SV + b-tag eff.",         "pe1");
        if (h_eec_c) lg->AddEntry(p_full,   "+ 2SV + b-tag + EEC weight", "pe1");
        lg->AddEntry(p_true,   normalise_corr_stages ? "Gen MC" : "Gen MC (scaled to data)", "l");
        lg->Draw();

        // CMS Internal above the top axis (outside the frame)
        TLatex cms; cms.SetNDC();
        cms.SetTextFont(62); cms.SetTextSize(0.042); cms.DrawLatex(0.13, 0.945, "CMS");
        cms.SetTextFont(52); cms.SetTextSize(0.034); cms.DrawLatex(0.235, 0.945, "Work In Progress");
        cms.SetTextFont(42); cms.SetTextSize(0.034); cms.DrawLatex(0.66, 0.945, "pp #sqrt{s} = 5.36 TeV");
        pad_main->RedrawAxis();

        // ----- bottom: incremental ratios (ratio - 1), colours matched to the added stage -----
        pad_ratio->cd();
        // r_2sv is drawn first, so it owns this pad's axes: style them on it and nowhere else.
        r_2sv->SetTitle("");
        r_2sv->GetYaxis()->SetRangeUser(rmin - rpad, rmax + rpad);
        // Short enough to fit the ratio pad's height: "ratio to prev. stage" runs past the pad
        // edge and gets clipped mid-word.
        r_2sv->GetYaxis()->SetTitle("ratio to prev.");
        r_2sv->GetYaxis()->CenterTitle(true);
        r_2sv->GetYaxis()->SetTitleFont(font_code); r_2sv->GetYaxis()->SetTitleSize(title_size); r_2sv->GetYaxis()->SetTitleOffset(1.5);
        r_2sv->GetYaxis()->SetLabelFont(font_code); r_2sv->GetYaxis()->SetLabelSize(label_size);
        r_2sv->GetYaxis()->SetNdivisions(505);
        r_2sv->GetXaxis()->SetTitle(obs_axis);
        r_2sv->GetXaxis()->CenterTitle(true);

        // Precision-43 offsets scale off the (short) ratio-pad height, so this is
        // ratio_x_title_offset and NOT title_offset -- see where they are defined.
        r_2sv->GetXaxis()->SetTitleFont(font_code); r_2sv->GetXaxis()->SetTitleSize(title_size); r_2sv->GetXaxis()->SetTitleOffset(ratio_x_title_offset);
        r_2sv->GetXaxis()->SetLabelFont(font_code); r_2sv->GetXaxis()->SetLabelSize(label_size);
        fixXRange(r_2sv);
        r_2sv->Draw("PE X0"); r_btag->Draw("PE X0 SAME");
        if (h_eec_c) r_eec->Draw("PE X0 SAME");

        // Unity reference. A cloned histogram drawn with "HIST L" joins BIN CENTRES, so it stopped
        // half a bin short of each end of the axis; a TLine spans the plotted range edge to edge.
        TLine *ref1 = new TLine(r_2sv->GetXaxis()->GetBinLowEdge(b1), 1.,
                                r_2sv->GetXaxis()->GetBinUpEdge(b2),  1.);
        ref1->SetLineColor(kGray+2); ref1->SetLineStyle(2); ref1->SetLineWidth(1);
        ref1->Draw("same");

        // labels for the three isolated corrections (colours match the added stage above)
        TLegend *lg_r = new TLegend(0.15, 0.84, 0.95, 0.99);
        lg_r->SetNColumns(3);
        lg_r->SetFillStyle(0); lg_r->SetBorderSize(0); lg_r->SetMargin(0.12);
        lg_r->SetTextFont(font_code); lg_r->SetTextSize(legend_size * 0.8);
        lg_r->AddEntry(r_2sv,  "2SV",  "pe1");
        lg_r->AddEntry(r_btag, "+b-tag", "pe1");
        if (h_eec_c) lg_r->AddEntry(r_eec,  "+EEC weight", "pe1");
        lg_r->Draw();
        pad_ratio->RedrawAxis();

        // Suffix so the normalised and absolute versions do not overwrite each other.
        // TString cmp_stem = folder + "correction_stages_" + label + Form("_pt%d", ibin_pt)   // until 2026-10-26
        TString cmp_stem = folder + "correction_stages_" + out_label + Form("_pt%d", ibin_pt)
                         + (normalise_corr_stages ? "_norm" : "");
        c_cmp->Print(cmp_stem + ".pdf");
        c_cmp->Print(cmp_stem + ".png");
    }


    // ---- DEBUG: reco SV pt vs gen B pt, to explain the size of the EEC-weight correction ----
    if (test_mode == 2) {
        TH2D *h_svpt = sumOverFiles<TH2D>(fin_unfolding, "h_svpt_vs_bpt");
        TH1D *h_ptr  = sumOverFiles<TH1D>(fin_unfolding, "h_ptratio");
        if (h_svpt) {
            std::cout << "\t---->Making DEBUG reco-vs-gen pt plot" << std::endl;
            gStyle->SetPalette(kViridis);
            TCanvas *c_dbg = new TCanvas("c_debug_svpt", "", 850, 700);
            c_dbg->SetLeftMargin(0.13); c_dbg->SetRightMargin(0.18); c_dbg->SetTicks(1, 1);
            h_svpt->SetStats(0);
            h_svpt->SetTitle(";gen partial-B p_{T} [GeV];reco SV p_{T} [GeV];weighted entries");
            h_svpt->GetXaxis()->SetTitleFont(font_code); h_svpt->GetXaxis()->SetTitleSize(title_size); h_svpt->GetXaxis()->SetTitleOffset(1.2);
            h_svpt->GetXaxis()->SetLabelFont(font_code); h_svpt->GetXaxis()->SetLabelSize(label_size);
            h_svpt->GetYaxis()->SetTitleFont(font_code); h_svpt->GetYaxis()->SetTitleSize(title_size); h_svpt->GetYaxis()->SetTitleOffset(1.5);
            h_svpt->GetYaxis()->SetLabelFont(font_code); h_svpt->GetYaxis()->SetLabelSize(label_size);
            h_svpt->GetZaxis()->SetTitleFont(font_code); h_svpt->GetZaxis()->SetTitleSize(title_size); h_svpt->GetZaxis()->SetTitleOffset(1.4);
            h_svpt->GetZaxis()->SetLabelFont(font_code); h_svpt->GetZaxis()->SetLabelSize(label_size);
            h_svpt->Draw("COLZ");
            // mean reco/gen pt on an opaque white box (readable over the colour map)
            if (h_ptr) {
                TLegend *lg_d = new TLegend(0.16, 0.80, 0.50, 0.87);
                lg_d->SetFillColor(kWhite); lg_d->SetFillStyle(1001); lg_d->SetBorderSize(1);
                lg_d->SetTextFont(font_code); lg_d->SetTextSize(legend_size * 0.75);
                lg_d->AddEntry((TObject*)0, Form("#LTp_{T}^{reco}/p_{T}^{gen}#GT = %.2f", h_ptr->GetMean()), "");
                lg_d->Draw();
            }
            // CMS Internal above the top axis
            TLatex cmsd; cmsd.SetNDC();
            cmsd.SetTextFont(62); cmsd.SetTextSize(0.042); cmsd.DrawLatex(0.13, 0.945, "CMS");
            cmsd.SetTextFont(52); cmsd.SetTextSize(0.034); cmsd.DrawLatex(0.235, 0.945, "Work In Progress");
            // c_dbg->Print(folder + "debug_svpt_vs_bpt_" + label + ".pdf");   // until 2026-10-26
            // c_dbg->Print(folder + "debug_svpt_vs_bpt_" + label + ".png");
            c_dbg->Print(folder + "debug_svpt_vs_bpt_" + out_label + ".pdf");
            c_dbg->Print(folder + "debug_svpt_vs_bpt_" + out_label + ".png");
        }
    }


    // ---- Graphical bottomline test

    std::cout << "Performing graphical bottomline test" << std::endl;
    TH1D *h_mc_reco_2D = h_mc_reco->ProjectionX("h_mc_reco_2D", ibin_pt, ibin_pt);
    TH1D *h_mc_true_2D = h_mc_true->ProjectionX("h_mc_true_2D", ibin_pt, ibin_pt);
    TH1D *h_data_purity_corrected_2D = h_data_purity_corrected->ProjectionX("h_data_purity_corrected_2D", ibin_pt, ibin_pt);
    TH1D *h_data_fully_corrected_2D = h_data_fully_corrected->ProjectionX("h_data_fully_corrected_2D", ibin_pt, ibin_pt);
    TH1D *h_data_refolded_2D = h_data_refolded->ProjectionX("h_data_refolded_2D", ibin_pt, ibin_pt);
    TH1D *h_data_after_fit_2D = h_data_after_fit->ProjectionX("h_data_after_fit_2D", ibin_pt, ibin_pt);
    TH1D *h_data_unfolded_2D = h_data_unfolded->ProjectionX("h_data_unfolded_2D", ibin_pt, ibin_pt);

    // A Divide() across incompatible binnings returns silently and leaves the ratio empty.
    if (h_mc_true_2D->GetNbinsX() != h_data_fully_corrected_2D->GetNbinsX()) {
        std::cerr << "ERROR: truth has " << h_mc_true_2D->GetNbinsX() << " dr bins but the unfolded result has "
                  << h_data_fully_corrected_2D->GetNbinsX() << std::endl;
        return;
    }

    if (true) {

        double ymax = 0.;
        for (auto h : {
                    h_mc_reco_2D,
                    h_mc_true_2D,
                    h_data_purity_corrected_2D,
                    h_data_fully_corrected_2D,
                    h_data_refolded_2D,
                    h_data_after_fit_2D,
                    h_data_unfolded_2D
                    }) {
                        if (!h) continue;
                        h->GetXaxis()->SetRange(ibin_dr_min, ibin_dr_max);
                        normalizeToUnitArea(h);
                    }
        // The y range has to clear the error bars of the curves that are actually drawn, not
        // just their central values, and not the ones left out of the plot.
        for (auto h : {h_mc_reco_2D, h_mc_true_2D, h_data_purity_corrected_2D, h_data_fully_corrected_2D}) {
            for (int i = ibin_dr_min; i <= ibin_dr_max; ++i)
                ymax = std::max(ymax, h->GetBinContent(i) + h->GetBinError(i));
        }

        Float_t pt_min_plot = h_data_reco->GetYaxis()->GetBinLowEdge(ibin_pt);
        Float_t pt_max_plot = h_data_reco->GetYaxis()->GetBinUpEdge(ibin_pt);

        TLegend *leg = new TLegend(0.15, 0.55, 0.55, 0.80);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetMargin(0.15);
        leg->SetTextFont(font_code);
        leg->SetTextSize(legend_size);
        // The last pT bin is open-ended (jets above jtpt_max are folded into it by
        // jtpt_fill()), so quote it as a threshold rather than a closed range.
        // leg->SetHeader(pt_max_plot == jtpt_max
        //                ? Form("p_{T}^{jet} > %.0f GeV", pt_min_plot)
        //                : Form("%.0f < p_{T}^{jet} < %.0f GeV", pt_min_plot, pt_max_plot));
        leg->SetHeader(jetPtRangeLabel(pt_min_plot, pt_max_plot));

        // Name every curve by what it is and, in the split test, WHICH HALF it came from.
        // The halves are not interchangeable and the plot has to say which is which:
        //   even half (ient % 2 == 0) = the pseudodata, and its own gen truth
        //   odd  half (ient % 2 == 1) = the response matrix and the purity/efficiency ratios
        // The closure being tested is "unfolded even half" against "even half truth" -- the
        // gen of the PSEUDODATA, not the gen that filled the matrix. A legend that just said
        // "Gen MC" left that unsaid, which is the one thing a reader of this plot has to know.
        TString lbl_reco     = split_test ? "Reco pseudodata (even half)"  : "Reco data";
        TString lbl_mc_reco  = split_test ? (apply_purity ? "Reco MC, matched (odd half)"
                                                          : "Reco MC, all reco (odd half)")
                                          : "Reco MC";
        TString lbl_unfolded = split_test ? "Unfolded pseudodata (even half)"  : "Unfolded data";
        TString lbl_mc_true  = split_test ? "Gen MC, pseudodata truth (even half)"  : "Gen MC";

        // Only the first histogram drawn owns the pad's axes; style them there and nowhere else.
        h_data_purity_corrected_2D->SetTitle("");//Data, " + label + " " + " response matrix");
        h_data_purity_corrected_2D->SetStats(0);
        h_data_purity_corrected_2D->SetMarkerColor(purple);
        h_data_purity_corrected_2D->SetLineColor(purple);
        h_data_purity_corrected_2D->SetMarkerStyle(kFullCircle);
        h_data_purity_corrected_2D->SetMarkerSize(1);
        // 1.75x headroom keeps the legend clear of the EEC peak.
        h_data_purity_corrected_2D->GetYaxis()->SetRangeUser(0., ymax*1.75);
        // TLatex (#Delta), not TMathText (\Delta\mbox{r}): TMathText axis titles do not get
        // placed reliably — the x-axis title was silently dropped altogether.
        h_data_purity_corrected_2D->GetYaxis()->SetTitle(obs_title);
        h_data_purity_corrected_2D->GetYaxis()->CenterTitle(true);
        h_data_purity_corrected_2D->GetYaxis()->SetTitleFont(font_code);
        h_data_purity_corrected_2D->GetYaxis()->SetTitleSize(title_size);
        h_data_purity_corrected_2D->GetYaxis()->SetTitleOffset(1.5);
        h_data_purity_corrected_2D->GetYaxis()->SetLabelFont(font_code);
        h_data_purity_corrected_2D->GetYaxis()->SetLabelSize(label_size);
        h_data_purity_corrected_2D->GetXaxis()->SetLabelSize(0); // x labels belong to the ratio pad
        h_data_purity_corrected_2D->GetXaxis()->SetTitleSize(0);
        leg->AddEntry(h_data_purity_corrected_2D, lbl_reco, "pe1");


        //h_data_after_fit_2D->SetMarkerColor(orange);
        //h_data_after_fit_2D->SetLineColor(orange);
        //h_data_after_fit_2D->SetMarkerStyle(kFullTriangleUp);
        //h_data_after_fit_2D->SetMarkerSize(1);
        //leg->AddEntry(h_data_after_fit_2D, "Before purity correction", "pe1");

        //h_data_unfolded_2D->SetMarkerColor(green);
        //h_data_unfolded_2D->SetLineColor(green);
        //h_data_unfolded_2D->SetMarkerStyle(kFullTriangleUp);
        //h_data_unfolded_2D->SetMarkerSize(1);
        //leg->AddEntry(h_data_unfolded_2D, "After unfolding (no efficiency correction)", "pe1");

        h_mc_reco_2D->SetStats(0);
        h_mc_reco_2D->SetMarkerColor(green);
        h_mc_reco_2D->SetLineColor(green);
        h_mc_reco_2D->SetMarkerStyle(kFullTriangleUp);
        h_mc_reco_2D->SetMarkerSize(1);
        leg->AddEntry(h_mc_reco_2D, lbl_mc_reco, "pe1");

        // Red: this is the measurement. Reserved for it across every plot in the analysis.
        h_data_fully_corrected_2D->SetMarkerColor(red);
        h_data_fully_corrected_2D->SetLineColor(red);
        h_data_fully_corrected_2D->SetMarkerStyle(kOpenCross);
        h_data_fully_corrected_2D->SetMarkerSize(1);
        leg->AddEntry(h_data_fully_corrected_2D, lbl_unfolded, "pe1");

        // Blue: the particle-level MC the measurement is compared to (the "theory" curve).
        h_mc_true_2D->SetMarkerColor(blue);
        h_mc_true_2D->SetLineColor(blue);
        h_mc_true_2D->SetMarkerStyle(kOpenTriangleUp);
        h_mc_true_2D->SetMarkerSize(1);
        leg->AddEntry(h_mc_true_2D, lbl_mc_true, "lp");   // line + open-triangle marker

        // Reco-level data, like h_data_purity_corrected_2D above, so it shares that curve's
        // purple and is told apart by the marker; blue now belongs to the gen MC.
        h_data_refolded_2D->SetMarkerColor(purple);
        h_data_refolded_2D->SetLineColor(purple);
        h_data_refolded_2D->SetMarkerStyle(kFullCross);
        h_data_refolded_2D->SetMarkerSize(1);


      
        
        TCanvas *c_unfold = new TCanvas("c_unfold", "", 800, 800);

        TPad *pad_main  = new TPad("pad_main",  "", 0., 0.3, 1., 1.);   // placeholder coords;
        TPad *pad_ratio = new TPad("pad_ratio", "", 0., 0.,  1., 0.3);  // applySquareLayout owns them
        
    
        
        
        for (TPad *p : {pad_main, pad_ratio}) {
            p->SetTicks(1, 0);   // mirrored ticks on all four sides
            p->SetFillColor(0);
        
        }
        //pad_ratio->SetLogx();
        //pad_main->SetLogx();

        pad_main->SetMargin(0.1, 0.1, 0.0, 0.1);
        // Bottom margin 0.30, not 0.23: the x title has to fit UNDER the tick labels, and B's
        // label (p_{T}^{b1}/(p_{T}^{b1} + p_{T}^{b2})) is two storeys tall -- at 0.23 its
        // subscripts were cut off by the canvas edge. #Delta r is unaffected either way.
        pad_ratio->SetMargin(0.1, 0.1, 0.30, 0.0);

    
        // Attach the pads to the canvas before drawing into them, so gPad and the pad's
        // primitive list agree at every Draw().
        c_unfold->cd();
        pad_main->Draw();
        pad_ratio->Draw();

        pad_main->cd();
        h_data_purity_corrected_2D->GetXaxis()->SetRange(ibin_dr_min, ibin_dr_max);
        fixXRange(h_data_purity_corrected_2D);
        h_data_purity_corrected_2D->Draw("PE X0");

        //h_data_unfolded_2D->Draw("pe1 same");
        h_mc_reco_2D->Draw("PE X0 same");
        h_data_fully_corrected_2D->Draw("PE X0 same");
        //h_data_after_fit_2D->Draw("pe1 same");
        h_mc_true_2D->Draw("hist E same");

        leg->Draw();
        TLatex *test_info_text = new TLatex;
        
        test_info_text->SetNDC();
        test_info_text->SetTextFont(font_code);
        test_info_text->SetTextSize(label_size);
        test_info_text->DrawLatex(0.57, 0.75,
            test_mode == 0 ? "MC full-sample closure" :
            test_mode == 1 ? "MC split closure test"  :
                             "Data Run 3");
        test_info_text->DrawLatex(0.57, 0.70, "2D unfolding");
        if(unfoldBayes) {
            test_info_text->DrawLatex(0.57, 0.65, "Bayesian unfolding");
            test_info_text->DrawLatex(0.57, 0.60, Form("N iter = %d", niter_now));
            // Refolding GoF of THIS iteration, so every scan page carries its own number.
            // Data only: in the closure modes the input is (mode 0) or comes from (mode 1) the
            // same MC as the response, so this number is a technical check, not a measured GoF,
            // and quoting it on the plot invites reading it as one. It is still printed to stdout
            // for both modes, and still drives the iteration scan.
            if (is_data) {
                test_info_text->DrawLatex(0.57, 0.55,
                    Form("refold #chi^{2}/ndf = %.2f", gof_all.chi2ndf()));
                test_info_text->DrawLatex(0.57, 0.50, Form("p = %.3f", gof_all.pvalue));
                test_info_text->DrawLatex(0.57, 0.45, Form("Response matrix CN = %.1f", cond_number));
            }
            else {
                test_info_text->DrawLatex(0.57, 0.55, Form("Response matrix CN = %.1f", cond_number));
            }
        }
        else {
            test_info_text->DrawLatex(0.57, 0.65, "Matrix inversion unfolding");
            test_info_text->DrawLatex(0.57, 0.60, Form("Response matrix CN = %.1f", cond_number));
        }
        //drawHeader();

        // CMS Internal above the top axis (outside the frame)
        TLatex cms_bl; cms_bl.SetNDC();
        cms_bl.SetTextFont(62); cms_bl.SetTextSize(0.045); cms_bl.DrawLatex(0.10, 0.945, "CMS"); // cms_bl.SetTextSize(0.042); 
        cms_bl.SetTextFont(52); cms_bl.SetTextSize(0.037); cms_bl.DrawLatex(0.205, 0.945, "Work In Progress"); // cms_bl.SetTextSize(0.034); 
        cms_bl.SetTextFont(42); cms_bl.SetTextSize(0.037); cms_bl.DrawLatex(0.66, 0.945, "pp #sqrt{s} = 5.36 TeV"); // cms_bl.SetTextSize(0.034);

        pad_main->RedrawAxis();

        
        TH1D *h_data_mc_reco_ratio = (TH1D *) h_data_purity_corrected_2D->Clone("h_data_mc_reco_ratio");
        
        h_data_mc_reco_ratio->SetTitle("");
        h_data_mc_reco_ratio->SetStats(0);
        h_data_mc_reco_ratio->Divide(h_mc_reco_2D);
        h_data_mc_reco_ratio->SetMarkerStyle(kFullCircle);
        h_data_mc_reco_ratio->SetMarkerColor(purple);
        h_data_mc_reco_ratio->SetLineColor(purple);
        h_data_mc_reco_ratio->SetMarkerSize(1.1);

        TH1D *h_data_mc_true_ratio = (TH1D *) h_data_fully_corrected_2D->Clone("h_data_mc_true_ratio");
        h_data_mc_true_ratio->Divide(h_mc_true_2D);
        h_data_mc_true_ratio->SetMarkerStyle(kOpenCross);
        h_data_mc_true_ratio->SetMarkerColor(red);
        h_data_mc_true_ratio->SetLineColor(red);
        h_data_mc_true_ratio->SetMarkerSize(1.3);

        TH1D *h_mc_gen_reco_ratio = (TH1D *) h_mc_true_2D->Clone("h_mc_gen_reco_ratio");


        h_mc_gen_reco_ratio->Divide(h_mc_reco_2D);
        // Each ratio carries the colour of its NUMERATOR in the main pad: this one is the gen MC,
        // so blue, matching the gen MC curve above.
        h_mc_gen_reco_ratio->SetLineColor(blue);
        h_mc_gen_reco_ratio->SetLineWidth(1);
        h_mc_gen_reco_ratio->SetMarkerSize(0);
                h_mc_gen_reco_ratio->GetYaxis()->SetTitle("ratio");
                h_mc_gen_reco_ratio->GetXaxis()->SetTitle(obs_axis);
                        h_mc_gen_reco_ratio->GetYaxis()->CenterTitle(true);
                        h_mc_gen_reco_ratio->GetXaxis()->CenterTitle(true);
        h_mc_gen_reco_ratio->GetXaxis()->SetTitleSize(0.2); 
        h_mc_gen_reco_ratio->GetXaxis()->SetTitleOffset(0.28);

        h_mc_gen_reco_ratio->GetXaxis()->SetLabelSize(0.06);
        h_mc_gen_reco_ratio->GetYaxis()->SetLabelSize(0.1);
        h_mc_gen_reco_ratio->GetYaxis()->SetNdivisions(505);
        h_mc_gen_reco_ratio->GetYaxis()->SetTitleOffset(0.24);
        h_mc_gen_reco_ratio->GetYaxis()->SetTitleSize(0.2);
        h_mc_gen_reco_ratio->GetYaxis()->SetLimits(0.5, 2);
        h_mc_gen_reco_ratio->SetLineStyle(2);
        h_mc_gen_reco_ratio->GetYaxis()->SetRangeUser(0.5, 2);

        // Match the ratio-pad axes to the main pad: pixel (precision-43) fonts with the same
        // title size, label size and title offset, so every axis name is the same size and the
        // same distance from its numbers. Overrides the fractional sizes set just above.
        // h_mc_gen_reco_ratio is drawn FIRST in the ratio pad, so it owns that pad's axes and
        // these are the settings that render -- the ones on h_data_mc_reco_ratio below do not.
        for (TAxis *ax : { h_mc_gen_reco_ratio->GetXaxis(), h_mc_gen_reco_ratio->GetYaxis() }) {
            ax->SetTitleFont(font_code);
            ax->SetTitleSize(title_size);
            ax->SetLabelFont(font_code);
            ax->SetLabelSize(label_size);
        }
        // X gets the ratio-pad offset, Y keeps the canvas-wide one: at title_offset the x title
        // is drawn a full pad height below the axis, i.e. off the canvas, which is why this
        // plot had no x title at all.
        h_mc_gen_reco_ratio->GetYaxis()->SetTitleOffset(title_offset);
        h_mc_gen_reco_ratio->GetXaxis()->SetTitleOffset(ratio_x_title_offset);


        // Orange, not red: its numerator is the unfolded data, but red is already drawn in this pad
        // as the unfolded/gen-MC points, and a dashed red line next to red crosses reads as the same
        // quantity. Green and purple are taken by the reco MC and the reco data, so orange it is.
        TH1D *h_unfolded_reco_ratio = (TH1D *) h_data_fully_corrected_2D->Clone("h_unfolded_reco_ratio");
        h_unfolded_reco_ratio->Divide(h_data_purity_corrected_2D);
        h_unfolded_reco_ratio->SetLineColor(orange);
        h_unfolded_reco_ratio->SetLineWidth(1);
        h_unfolded_reco_ratio->SetMarkerSize(0);
        h_unfolded_reco_ratio->SetLineStyle(2);
        

        h_data_mc_reco_ratio->GetYaxis()->SetRangeUser(0.5, 2);
        h_data_mc_reco_ratio->GetYaxis()->SetTitle("ratio");
        h_data_mc_reco_ratio->GetYaxis()->CenterTitle(true);
        h_data_mc_reco_ratio->GetYaxis()->SetTitleFont(font_code);
        h_data_mc_reco_ratio->GetYaxis()->SetTitleSize(title_size);
        h_data_mc_reco_ratio->GetYaxis()->SetTitleOffset(0.05);
        h_data_mc_reco_ratio->GetYaxis()->SetLabelFont(font_code);
        h_data_mc_reco_ratio->GetYaxis()->SetLabelSize(label_size);
        h_data_mc_reco_ratio->GetYaxis()->SetNdivisions(505);
        h_data_mc_reco_ratio->GetXaxis()->SetTitle(obs_axis);
        h_data_mc_reco_ratio->GetXaxis()->CenterTitle(true);
        h_data_mc_reco_ratio->GetXaxis()->SetTitleFont(font_code);
        h_data_mc_reco_ratio->GetXaxis()->SetTitleSize(title_size);
        
        // With precision-3 fonts this offset scales off the (short) pad height: 3.2 pushed the
        // title clean off the canvas. 1.35 lands it just under the tick labels.
        h_data_mc_reco_ratio->GetXaxis()->SetTitleOffset(0.05);
        h_data_mc_reco_ratio->GetXaxis()->SetLabelFont(12);
        h_data_mc_reco_ratio->GetXaxis()->SetLabelSize(16);
        h_data_mc_reco_ratio->GetXaxis();
       

        // -- Add style also for h_mc_gen_reco_ratio (drawn first in ratio)
        h_mc_gen_reco_ratio->GetYaxis()->SetRangeUser(0.5, 2);
        h_mc_gen_reco_ratio->GetYaxis()->SetTitle("ratio");
        h_mc_gen_reco_ratio->GetYaxis()->CenterTitle(true);
        h_mc_gen_reco_ratio->GetYaxis()->SetTitleFont(font_code);
        h_mc_gen_reco_ratio->GetYaxis()->SetTitleSize(title_size);
        h_mc_gen_reco_ratio->GetYaxis()->SetTitleOffset(1.5);
        h_mc_gen_reco_ratio->GetYaxis()->SetLabelFont(font_code);
        h_mc_gen_reco_ratio->GetYaxis()->SetLabelSize(label_size);
        h_mc_gen_reco_ratio->GetYaxis()->SetNdivisions(505);
        h_mc_gen_reco_ratio->GetXaxis()->SetTitle(obs_axis);
        h_mc_gen_reco_ratio->GetXaxis()->CenterTitle(true);
        h_mc_gen_reco_ratio->GetXaxis()->SetTitleFont(font_code);
        h_mc_gen_reco_ratio->GetXaxis()->SetTitleSize(title_size);

       

        

        // The unity line and the gap bands span the axis actually plotted, not the full dr
        // range: the two disagree as soon as ibin_dr_min/ibin_dr_max restrict the view.
        const Double_t x_lo = h_mc_reco_2D->GetXaxis()->GetBinLowEdge(ibin_dr_min);
        const Double_t x_hi = h_mc_reco_2D->GetXaxis()->GetBinUpEdge(ibin_dr_max);

        pad_ratio->cd();
        //h_data_mc_reco_ratio->Draw("axis");   // axes first, so the bands land behind everything

      
        TLine *line = new TLine(x_lo, 1., x_hi, 1.);
        line->SetLineWidth(1.);
        line->SetLineStyle(kSolid);   // solid, so the unity line reads apart from the dashed blue/orange curves
        line->SetLineColor(kGray+2);

        // "Before" band: for each dr bin, fill the full bin width between the gen MC/reco MC line
        // (blue) and the unfolded/reco line (orange) -- the size of the correction. Stepped closed
        // polygon: forward along blue at the bin edges, back along orange, so each bin is filled
        // edge to edge. Drawn behind, translucent light grey.
        const int nb_band = ibin_dr_max - ibin_dr_min + 1;
        TGraph *band_before_after = new TGraph(4 * nb_band);
        int ip = 0;
        for (int i = 0; i < nb_band; ++i) {              // top edge: along blue, left -> right
            const int bin = ibin_dr_min + i;
            const double xlo = h_mc_gen_reco_ratio->GetXaxis()->GetBinLowEdge(bin);
            const double xhi = h_mc_gen_reco_ratio->GetXaxis()->GetBinUpEdge(bin);
            const double yr  = h_mc_gen_reco_ratio->GetBinContent(bin);
            band_before_after->SetPoint(ip++, xlo, yr);
            band_before_after->SetPoint(ip++, xhi, yr);
        }
        for (int i = nb_band - 1; i >= 0; --i) {         // bottom edge: along orange, right -> left
            const int bin = ibin_dr_min + i;
            const double xlo = h_unfolded_reco_ratio->GetXaxis()->GetBinLowEdge(bin);
            const double xhi = h_unfolded_reco_ratio->GetXaxis()->GetBinUpEdge(bin);
            const double yg  = h_unfolded_reco_ratio->GetBinContent(bin);
            band_before_after->SetPoint(ip++, xhi, yg);
            band_before_after->SetPoint(ip++, xlo, yg);
        }
        // Solid light grey, NOT SetFillColorAlpha: the PDF backend does not render alpha
        // transparency (it shows on screen in XQuartz but disappears in the PDF), whereas a
        // solid fill renders identically everywhere. Curves are redrawn on top, so opaque is fine.
        band_before_after->SetFillColor(TColor::GetColor(224, 224, 224));
        band_before_after->SetFillStyle(1001);
        band_before_after->SetLineWidth(0);

        fixXRange(h_mc_gen_reco_ratio);
        h_mc_gen_reco_ratio->Draw("hist");     // axes + blue line
        band_before_after->Draw("f same");     // band behind the curves
        h_mc_gen_reco_ratio->Draw("hist same");// redraw blue crisp over the band
        line->Draw("same");

        h_unfolded_reco_ratio->Draw("hist same");
        h_data_mc_reco_ratio->Draw("PE X0 same");
        h_data_mc_true_ratio->Draw("PE X0 same");

        // Curves grouped by family (points = closure, lines = size of correction), then the two
        // band swatches, so the legend reads top-down as the story: correction, closure, and the
        // two gaps that connect them.
        TLegend *leg_ratio = new TLegend(0.16, 0.66, 0.96, 0.985);
        leg_ratio->SetFillStyle(0);
        leg_ratio->SetBorderSize(0);
        leg_ratio->SetTextFont(font_code);
        leg_ratio->SetTextSize(legend_size - 2.);
        leg_ratio->SetNColumns(2);
        leg_ratio->AddEntry(h_data_mc_reco_ratio, "reco data/reco MC", "pe1");
        leg_ratio->AddEntry(h_mc_gen_reco_ratio, "gen MC/reco MC", "l");
        leg_ratio->AddEntry(h_data_mc_true_ratio, "unfolded data/gen MC", "pe1");
        leg_ratio->AddEntry(h_unfolded_reco_ratio, "unfolded/reco", "l");
        leg_ratio->Draw();
        pad_ratio->RedrawAxis();
        c_unfold->cd();
        c_unfold->Update();
        
        h_data_mc_reco_ratio->GetXaxis();
        c_unfold->Modified();
        c_unfold->Update();
       
        if (scan_niter) {
            // One page per iteration: "(" opens the PDF on the first page, ")" closes it on the last.
            TString iter_page = iter_pdf;
            if (iter == niter_min)      iter_page += "(";
            else if (iter == niter_max) iter_page += ")";
            c_unfold->Print(iter_page, "pdf");

            // Also drop a standalone PNG per iteration into png_dir (zero-padded so they sort 001..100).
            c_unfold->Print(Form("%siteration_%03d.png", png_dir.Data(), iter + 1), "png");
        } else {
            // Single run: the usual one-off bottomline plot next to the other outputs.
            TString plot_stem = folder + "unfolding_plot_" + label + "_bottomline_test_eec_2D";
            c_unfold->Print(plot_stem + ".pdf");
            c_unfold->Print(plot_stem + ".png");
        }

        // Only delete the canvas when another pass will recreate it (avoids the name clash during
        // the scan). On the final/only pass keep it, so the window stays open to view.
        if (iter < niter_max) delete c_unfold;


    }

    // Write the corrected histograms once, from the final iteration's result.
    if (iter == niter_max) {
    std::cout << "Creating file " << fout_name << std::endl;
    TFile *fout = new TFile(fout_name, "recreate");

    h_data_fully_corrected->SetName("h_data_unfolded_corrected");
    h_data_fully_corrected->Write();
    h_data_after_fit->SetName("h_data");
    h_data_after_fit->Write();
    h_mc_reco->SetName("h_mc_reco");
    h_mc_reco->Write();

    //Write 2D histograms for plotting
    h_mc_reco_2D->Write();
    // The particle-level gen reference -- hgenjet_2b_all in data mode. Re-enabled so
    // apply_weights_and_systematics.C can draw it against the final result: it is the only
    // curve in this file that lives at gen level, and without it the theory comparison
    // would have to re-open and re-sum every per-block MC file itself. Already restricted
    // to [ibin_dr_min, ibin_dr_max] and unit-normalised by the bottomline block above,
    // exactly like h_data_fully_corrected_2D next to it.
    h_mc_true_2D->Write();
    h_data_purity_corrected_2D->Write();
    h_data_fully_corrected_2D->Write();
    h_data_refolded_2D->Write();

    fout->Close();
    delete fout;
    }

    }   // end loop over iteration counts
    gErrorIgnoreLevel = prev_ignore;

    // ---- Regularisation: pick the number of iterations from the refolding test ------------
    // p-value of refolded vs input against the iteration number. It rises as the nonclosure
    // shrinks, then saturates. The optimal iteration is the smallest niter with
    // p >= pvalue_threshold; going further only feeds fluctuations back through the prior.
    if (unfoldBayes && v_niter.size() > 1) {

        auto firstAbove = [&](const std::vector<double> &p) -> int {
            for (size_t i = 0; i < p.size(); ++i)
                if (p[i] >= pvalue_threshold) return (int) v_niter[i];
            return -1;   // never reaches the threshold within the scanned range
        };
        const int best_all = firstAbove(v_pvalue_all);
        const int best_pt  = firstAbove(v_pvalue_pt);

        // Secondary diagnostic: with a well-conditioned response the p-value can clear the
        // threshold at the first iteration, and chi2/ndf still shows where the refolding starts
        // over-describing the input (chi2/ndf well below 1).
        auto closestToOne = [&](const std::vector<double> &c) -> int {
            int best = -1; double dmin = 1e30;
            for (size_t i = 0; i < c.size(); ++i) {
                if (c[i] <= 0.) continue;
                const double d = std::fabs(c[i] - 1.);
                if (d < dmin) { dmin = d; best = (int) v_niter[i]; }
            }
            return best;
        };
        const int c2one_all = closestToOne(v_chi2ndf_all);

        std::cout << "\n================ Refolding goodness-of-fit scan ================\n"
                  << "  input   : h_data_purity_corrected (reconstructed level)\n"
                  << "  refolded: response->ApplyToTruth(unfolded)\n"
                  << Form("  criterion: smallest niter with p >= %.3f\n", pvalue_threshold);
        std::cout << Form("  %6s  %12s %10s   %12s %10s   %14s %10s   %14s %10s   %14s %10s   %14s %10s\n",
                          "niter", "chi2/ndf 2D", "p 2D", "chi2/ndf pT", "p pT",
                          "chi2/ndf 2D+MC", "p 2D+MC", "chi2/ndf cov", "p cov",
                          "chi2/ndf 2D+r", "p 2D+r", "chi2/ndf 2D+Vr", "p 2D+Vr");
        for (size_t i = 0; i < v_niter.size(); ++i)
            std::cout << Form("  %6.0f  %12.3f %10.3e   %12.3f %10.3e   %14.3f %10.3e   %14.3f %10.3e   %14.3f %10.3e   %14.3f %10.3e\n",
                              v_niter[i], v_chi2ndf_all[i], v_pvalue_all[i],
                              v_chi2ndf_pt[i], v_pvalue_pt[i],
                              v_chi2ndf_all_mc[i], v_pvalue_all_mc[i],
                              v_chi2ndf_all_cov[i], v_pvalue_all_cov[i],
                              v_chi2ndf_all_r[i], v_pvalue_all_r[i],
                              v_chi2ndf_all_rcov[i], v_pvalue_all_rcov[i]);
        // CMS-AN-2024-133 rule: first iteration with p >= 0.95, data (+) refolded stat, uncorrelated.
        auto firstAboveT = [&](const std::vector<double> &p, double t) -> int {
            for (size_t i = 0; i < p.size(); ++i) if (p[i] >= t) return (int) v_niter[i];
            return -1;
        };
        // SUPERSEDED 2026-10-26 by CMS-SMP-22-007 (below). Was THE CHOICE (Zoe, 2026-10-19),
        // CMS-AN-2024-133 Sec. 5.6/6: the first iteration whose
        // refolded distribution reaches p >= 0.95 against the input, per-bin error
        // sigma_d (+) sigma_r -- sigma_r the spread of the refolded result in data-only toys
        // (refoldChi2Toys), taken as uncorrelated with d. No response MC stat (until 10-19 it was
        // added -- Zoe: wrong). The analysis driver reads the OPTIMAL line below.
        // const double an_threshold = 0.95;
        // const int best_an = firstAboveT(v_pvalue_all_mc, an_threshold);
        // const int best_an = firstAboveT(v_pvalue_all, an_threshold);
        // const int best_an = firstAboveT(v_pvalue_all_r, an_threshold);   // 10-19 .. 10-26
        //
        // THE CHOICE since 2026-10-26 (Zoe): CMS-SMP-22-007 (JHEP 05 (2024) 116), Sec. 5 -- stop at
        // the first iteration whose forward-folded distribution has p ABOVE 0.05 against the
        // measured one, using the detector-level DATA statistical covariance only. Ours is
        // diagonal: one entry per jet per observable, and every correction before the unfolding
        // acts per bin, so the data bins are uncorrelated -- that is v_pvalue_all (refoldChi2
        // without the MC term), over the full 2D (observable x jet pT) space that is unfolded.
        // The AN-2024-133 rule is still printed below as a cross-check.
        // const double an_threshold = 0.05;
        // int best_an = -1;
        // for (size_t i = 0; i < v_pvalue_all.size(); ++i)
        //     if (v_pvalue_all[i] > an_threshold) { best_an = (int) v_niter[i]; break; }   // CMS-SMP-22-007, 10-26 only
        //
        // THE CHOICE since 2026-10-26 evening (Zoe): CMS-AN-2024-133 (v8, Baldenegro, Cunqueiro,
        // Mijuskovic, Vladimirov) Sec. 5.6 / 6 -- the first iteration whose refolded distribution
        // reaches p >= 0.95 against the input. The AN does not say which errors enter the chi2;
        // data statistics only (v_pvalue_all), as in CMS-SMP-22-007 from the same group. One niter
        // for every variation for now (the AN re-chooses it per systematic, Table 3 -- not done).
        // The CMS-SMP-22-007 rule (first p > 0.05) is printed below as a cross-check.
        const double an_threshold = 0.95;
        const int best_an = firstAboveT(v_pvalue_all, an_threshold);
        int best_smp = -1;
        for (size_t i = 0; i < v_pvalue_all.size(); ++i)
            if (v_pvalue_all[i] > 0.05) { best_smp = (int) v_niter[i]; break; }
        std::cout << Form("  ----> first p >= 0.95, data (+) refolded stat, full cov = %d\n", firstAboveT(v_pvalue_all_rcov, 0.95))
                  << Form("  ----> first p >= 0.95, data (+) refolded stat = %d\n", firstAboveT(v_pvalue_all_r, 0.95))
                  << Form("  ----> first p >= 0.95, cov(r - d) toys, data stat = %d\n", firstAboveT(v_pvalue_all_cov, 0.95))
                  << Form("  ----> first p >= 0.05, cov(r - d) toys, data stat = %d\n", firstAboveT(v_pvalue_all_cov, 0.05))
                  << Form("  ----> first p >= 0.95, data stat only       = %d\n", firstAboveT(v_pvalue_all, 0.95))
                  << Form("  ----> first p >= 0.05, data (+) MC stat     = %d\n", firstAboveT(v_pvalue_all_mc, 0.05))
                  << Form("  ----> first p >= 0.95, data (+) MC stat     = %d\n", firstAboveT(v_pvalue_all_mc, 0.95))
                  // << Form("  ----> OPTIMAL niter (AN-2024-133: first p >= %.2f, data (+) refolded stat) = %d\n",
                  //         an_threshold, best_an);   // until 2026-10-26
                  // << Form("  ----> cross-check, AN-2024-133 (first p >= 0.95, data (+) refolded stat) = %d\n",
                  //         firstAboveT(v_pvalue_all_r, 0.95))
                  // << Form("  ----> OPTIMAL niter (CMS-SMP-22-007: first p > %.2f, data stat) = %d\n",
                  //         an_threshold, best_an);   // 10-26 afternoon only
                  << Form("  ----> cross-check, CMS-SMP-22-007 (first p > 0.05, data stat) = %d\n", best_smp)
                  << Form("  ----> OPTIMAL niter (AN-2024-133: first p >= %.2f, data stat) = %d\n",
                          an_threshold, best_an);
        std::cout << Form("  ----> optimal niter (full 2D reco space) = %d\n", best_all)
                  << Form("  ----> optimal niter (pT bin %d only)      = %d\n", ibin_pt, best_pt)
                  << Form("  ----> cross-check, chi2/ndf closest to 1 = %d\n", c2one_all)
                  << "================================================================\n"
                  << std::endl;

        // ---- Bias vs statistical error (added 2026-10-01) ------------------------------
        // The refolding test above compares against the data's own statistical errors, which
        // are tiny here (~1e6 jets): it keeps rejecting until the unfolding has nearly stopped
        // regularising, i.e. it pushes niter towards matrix inversion. This instead balances
        // the two things an extra iteration trades:
        //   bias(n)  = u_n - u_MI   the pull of the prior that is still left -- matrix
        //                           inversion is the unregularised (unbiased) limit
        //   stat(n)  = sigma_n      the statistical error, which grows with n
        // in the measured pT bin, each distribution unit-area normalised (the result is a
        // shape), all relative to u_MI. The optimum minimises the rms over bins of
        // sqrt(bias^2 + stat^2). Cross-check: dchange = rms of (u_n - u_{n-1}) / sigma_n, how
        // much one more iteration moves the result in units of its own stat error.
        // Kept outside the block below so the convergence plot further down can draw them.
        // ---- chi2 bottom-line test vs iteration (2026-10-29, Zoe; AN Fig. 68) ----
        if (scan_niter && v_bl_unf.size() == v_niter.size() && !v_bl_unf.empty()) {
            const int nn = (int) v_niter.size();
            TGraph *g_unf = new TGraph(nn), *g_reco = new TGraph(nn), *g_rat = new TGraph(nn);
            double ymx = 0., rmx = 0., ymn = 1e30, rmn = 1e30;
            for (int k = 0; k < nn; ++k) {
                for (double v : {v_bl_unf[k], v_bl_reco[k]}) if (v > 0.) ymn = std::min(ymn, v);
                g_unf->SetPoint(k, v_niter[k], v_bl_unf[k]); g_reco->SetPoint(k, v_niter[k], v_bl_reco[k]);
                const double rr = v_bl_reco[k] > 0. ? v_bl_unf[k] / v_bl_reco[k] : 0.;
                g_rat->SetPoint(k, v_niter[k], rr);
                if (rr > 0.) rmn = std::min(rmn, rr);
                ymx = std::max(ymx, std::max(v_bl_unf[k], v_bl_reco[k])); rmx = std::max(rmx, rr);
            }
            const char *ref_env_bl = std::getenv("NITER_REF");
            const int n_ref_bl = (ref_env_bl && std::atoi(ref_env_bl) > 0) ? std::atoi(ref_env_bl) : ((obs.name == "dr") ? 16 : 126);
            // const double xmax_bl = v_niter.back() + 0.5, ytop = 1.45 * ymx, rtop = std::max(1.5, 1.15 * rmx);
            // 2026-10-29 (Zoe): log y on the chi2 panel; a decade of room above for the legend
            // const double xmax_bl = v_niter.back() + 0.5, ytop = 4. * ymx, rtop = std::max(1.5, 1.15 * rmx);
            // const double ybot = (ymn < 1e30) ? 0.2 * ymn : 0.1;
            // Then (Zoe): whole decades, labelled 10^n -- from the decade below the smallest chi2 to one
            // decade above the largest (room for the legend).
            const double xmax_bl = v_niter.back() + 0.5, rtop = std::max(1.5, 1.15 * rmx);
            const double ybot = (ymn < 1e30) ? std::pow(10., std::floor(std::log10(ymn))) : 0.1;
            const double ytop = std::pow(10., std::ceil(std::log10(3. * ymx)));   // the decade above 3x the largest
            // the ratio panel in log y too (Zoe), whole decades around 1
            const double rbot = (rmn < 1e30) ? std::min(0.1, std::pow(10., std::floor(std::log10(rmn)))) : 0.1;
            const double rtop_log = std::max(10., std::pow(10., std::ceil(std::log10(3. * rmx))));
            TCanvas *c_bl = new TCanvas("c_bottomline_chi2", "", 900, 900);
            TPad *ptb = new TPad("pbl_top", "", 0., 0.35, 1., 1.), *pbb = new TPad("pbl_bot", "", 0., 0., 1., 0.35);
            ptb->SetBottomMargin(0.02); ptb->SetLeftMargin(0.14); ptb->SetRightMargin(0.04); ptb->SetTopMargin(0.10);
            pbb->SetTopMargin(0.02);    pbb->SetBottomMargin(0.30); pbb->SetLeftMargin(0.14); pbb->SetRightMargin(0.04);
            pbb->SetLogy();
            ptb->SetLogy();
            c_bl->cd(); ptb->Draw(); pbb->Draw(); ptb->cd();
            // TH1F *fr = ptb->DrawFrame(0.5, 0., xmax_bl, ytop);
            TH1F *fr = ptb->DrawFrame(0.5, ybot, xmax_bl, ytop);
            fr->GetYaxis()->SetTitle("#chi^{2}"); fr->GetYaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetLabelSize(0.045);
            fr->GetYaxis()->SetTitleOffset(1.1); fr->GetXaxis()->SetLabelSize(0);
            // fr->GetYaxis()->SetMoreLogLabels(); fr->GetYaxis()->SetNoExponent();   // 20, 30, ... 500: Zoe wants 10^n
            g_unf->SetMarkerStyle(20); g_unf->SetMarkerColor(ROCColor::red()); g_unf->SetLineColor(ROCColor::red()); g_unf->SetLineWidth(2);
            g_reco->SetMarkerStyle(21); g_reco->SetMarkerColor(ROCColor::blue()); g_reco->SetLineColor(ROCColor::blue()); g_reco->SetLineWidth(2);
            g_unf->Draw("PL same"); g_reco->Draw("PL same");
            // TLine *lref = new TLine(n_ref_bl, ybot, n_ref_bl, ytop);
            // stops just under the legend (NDC 0.62 of a frame spanning 0.02-0.90; log y)
            const double yline = std::pow(10., std::log10(ybot) + (0.62 - 0.02) / 0.88 * (std::log10(ytop) - std::log10(ybot)));
            TLine *lref = new TLine(n_ref_bl, ybot, n_ref_bl, yline);
            lref->SetLineStyle(2); lref->SetLineColor(TColor::GetColor("#f29fc5")); lref->SetLineWidth(2); lref->Draw();
            // white legend box: the chosen-niter line runs under it (A, N = 126)
            TLegend *lg = new TLegend(0.45, 0.64, 0.94, 0.86); lg->SetBorderSize(0); lg->SetFillStyle(1001); lg->SetFillColor(kWhite); lg->SetTextSize(0.042);
            lg->AddEntry(g_unf, "#chi^{2}(unfolded data, gen MC)", "pl");
            lg->AddEntry(g_reco, "#chi^{2}(reco data, reco MC)", "pl");
            lg->AddEntry(lref, Form("chosen N_{iter} = %d", n_ref_bl), "l");
            lg->Draw();
            TLatex tl; tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.042);
            // below the curves (the legend owns the top right)
            tl.DrawLatex(0.18, 0.14, test_mode == 2 ? "Data" : (split_test ? "Split test" : "Full-MC closure"));
            tl.DrawLatex(0.18, 0.08, Form("%d %s bins, %s, MC scaled to data", ibin_dr_max - ibin_dr_min + 1,
                                          obs.name == "dr" ? "#DeltaR" : "A",
                                          jetPtRangeLabel(jtpt_binsVector[ibin_pt - 1], jtpt_binsVector[ibin_pt]).Data()));
            TLatex cmsl; cmsl.SetNDC();
            cmsl.SetTextFont(62); cmsl.SetTextSize(0.055); cmsl.DrawLatex(0.14, 0.92, "CMS");
            cmsl.SetTextFont(52); cmsl.SetTextSize(0.045); cmsl.DrawLatex(0.25, 0.92, "Work In Progress");
            cmsl.SetTextFont(42); cmsl.SetTextSize(0.045); cmsl.DrawLatex(0.64, 0.92, "pp #sqrt{s} = 5.36 TeV");
            pbb->cd();
            // TH1F *fb = pbb->DrawFrame(0.5, 0., xmax_bl, rtop);
            TH1F *fb = pbb->DrawFrame(0.5, rbot, xmax_bl, rtop_log);
            fb->GetYaxis()->SetTitle("#chi^{2}_{unf} / #chi^{2}_{reco}"); fb->GetYaxis()->CenterTitle(true);
            fb->GetYaxis()->SetTitleSize(0.09); fb->GetYaxis()->SetTitleOffset(0.65); fb->GetYaxis()->SetLabelSize(0.08); fb->GetYaxis()->SetNdivisions(505);
            fb->GetXaxis()->SetTitle("number of iterations"); fb->GetXaxis()->SetTitleSize(0.10); fb->GetXaxis()->SetLabelSize(0.085);
            g_rat->SetMarkerStyle(20); g_rat->SetMarkerColor(kBlack); g_rat->SetLineWidth(2); g_rat->Draw("PL same");
            TLine *l1b = new TLine(0.5, 1., xmax_bl, 1.); l1b->SetLineStyle(2); l1b->Draw();
            TLine *lrefb = new TLine(n_ref_bl, rbot, n_ref_bl, rtop_log);
            lrefb->SetLineStyle(2); lrefb->SetLineColor(TColor::GetColor("#f29fc5")); lrefb->SetLineWidth(2); lrefb->Draw();
            const TString bl_stem = folder + "bottomline_chi2_vs_iteration_" + out_label;
            c_bl->Print(bl_stem + ".pdf"); c_bl->Print(bl_stem + ".png");
            std::cout << "\t---->chi2 bottom-line test vs iteration: " << bl_stem << ".png" << std::endl;
        }

        // ---- Iteration bottom-line test (2026-10-28, Zoe): the unfolded distribution in the measured
        // pT bin at every iteration 1..N_ref, and its ratio to iteration N_ref (the chosen niter).
        // N_ref = NITER_REF (env), default 16 for dR and 126 for A (the AN-2024-133 rule on 10-28).
        // // After the unfolding, before the efficiency corrections: those are the same factors for every
        // // iteration, so the ratios are the ones of the final result.
        // Since the same evening: the nominal result of each iteration, every correction and SF applied
        // (h_data_fully_corrected). A shows at most ~16 curves.
        // if (scan_niter && v_unf_pt.size() == v_niter.size()) {   // first version: before the corrections
        // Then (Zoe, "with the nominal result (SF applied)"): the fully corrected result of each iteration.
        if (scan_niter && v_fc_pt.size() == v_niter.size()) {
            const std::vector<TH1D *> &v_it = v_fc_pt;
            const char *ref_env = std::getenv("NITER_REF");
            const int n_ref = (ref_env && std::atoi(ref_env) > 0) ? std::atoi(ref_env) : ((obs.name == "dr") ? 16 : 126);
            int k_ref = -1;
            for (size_t k = 0; k < v_niter.size(); ++k) if ((int) v_niter[k] == n_ref) k_ref = (int) k;
            if (k_ref < 0) {
                std::cout << "\t---->iteration test: niter " << n_ref << " not in the scan, skipped" << std::endl;
            } else {
                const int step = std::max(1, (int) std::ceil(n_ref / 16.));
                TCanvas *c_it = new TCanvas("c_iter_test", "", 900, 1000);
                TPad *pt = new TPad("pit_top", "", 0., 0.35, 1., 1.), *pb = new TPad("pit_bot", "", 0., 0., 1., 0.35);
                pt->SetBottomMargin(0.02); pt->SetLeftMargin(0.14); pt->SetRightMargin(0.04); pt->SetTopMargin(0.13);
                pb->SetTopMargin(0.02);    pb->SetBottomMargin(0.30); pb->SetLeftMargin(0.14); pb->SetRightMargin(0.04);
                c_it->cd(); pt->Draw(); pb->Draw();
                TH1D *ref = (TH1D *) v_it[k_ref]->Clone("h_iter_ref"); ref->SetDirectory(nullptr);
                ref->Scale(1., "width");
                std::vector<int> ks;
                for (int k = 0; k <= k_ref; ++k) if (((int) v_niter[k] - 1) % step == 0 || k == k_ref) ks.push_back(k);
                gStyle->SetPalette(kViridis);   // the palette this function already sets for its 2D plots
                TLegend *lg = new TLegend(0.55, 0.60, 0.95, 0.86); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.026);
                lg->SetNColumns(2);
                double ymax = 0., rlo = 1., rhi = 1.;
                std::vector<TH1D *> hs, rs;
                for (int k : ks) {
                    TH1D *h = (TH1D *) v_it[k]->Clone(Form("h_iter_%d", (int) v_niter[k])); h->SetDirectory(nullptr);
                    h->Scale(1., "width");
                    TH1D *r = (TH1D *) h->Clone(Form("r_iter_%d", (int) v_niter[k])); r->SetDirectory(nullptr);
                    r->Divide(ref);
                    // the ratio's error is the iteration's own relative error (the two are strongly correlated)
                    for (int i = 1; i <= r->GetNbinsX(); ++i)
                        r->SetBinError(i, h->GetBinContent(i) > 0. ? r->GetBinContent(i) * h->GetBinError(i) / h->GetBinContent(i) : 0.);
                    // the chosen niter in pastel pink (Zoe), the others along the palette
                    // const int col = (k == k_ref) ? kBlack : TColor::GetColorPalette((int) (254. * k / std::max(1, k_ref)));
                    const int col = (k == k_ref) ? TColor::GetColor("#f29fc5") : TColor::GetColorPalette((int) (254. * k / std::max(1, k_ref)));
                    for (TH1D *x : {h, r}) { x->SetLineColor(col); x->SetMarkerColor(col); x->SetLineWidth(k == k_ref ? 3 : 2); x->SetMarkerStyle(k == k_ref ? 20 : 1); }
                    for (int i = ibin_dr_min; i <= ibin_dr_max; ++i) {
                        ymax = std::max(ymax, h->GetBinContent(i));
                        if (r->GetBinContent(i) > 0.) { rlo = std::min(rlo, r->GetBinContent(i)); rhi = std::max(rhi, r->GetBinContent(i)); }
                    }
                    hs.push_back(h); rs.push_back(r);
                    lg->AddEntry(h, Form("N_{iter} = %d%s", (int) v_niter[k], k == k_ref ? " (chosen)" : ""), "l");
                }
                // closures (2026-10-28, Zoe): the truth too, dashed, and its ratio to N_ref below
                TH1D *ht = nullptr, *rt = nullptr;
                if (test_mode != 2 && h_mc_true) {
                    ht = h_mc_true->ProjectionX("h_iter_truth", ibin_pt, ibin_pt, "e"); ht->SetDirectory(nullptr);
                    ht->Scale(1., "width");
                    rt = (TH1D *) ht->Clone("r_iter_truth"); rt->SetDirectory(nullptr); rt->Divide(ref);
                    for (TH1D *x : {ht, rt}) { x->SetLineColor(kBlack); x->SetLineStyle(2); x->SetLineWidth(3); x->SetMarkerSize(0); }
                    for (int i = ibin_dr_min; i <= ibin_dr_max; ++i) {
                        ymax = std::max(ymax, ht->GetBinContent(i));
                        if (rt->GetBinContent(i) > 0.) { rlo = std::min(rlo, rt->GetBinContent(i)); rhi = std::max(rhi, rt->GetBinContent(i)); }
                    }
                    lg->AddEntry(ht, split_test ? "truth (even half)" : "truth (gen MC)", "l");
                }
                pt->cd();
                hs.front()->SetMaximum(1.6 * ymax); hs.front()->SetMinimum(0.);   // room for the legend
                hs.front()->GetXaxis()->SetRange(ibin_dr_min, ibin_dr_max);
                // hs.front()->GetYaxis()->SetTitle(Form("fully corrected dN/d%s (nominal, SFs applied)", obs.axis.Data()));
                hs.front()->GetYaxis()->SetTitle(test_mode == 2 ? Form("fully corrected dN/d%s (nominal, SFs applied)", obs.axis.Data())
                                                                : Form("unfolded + corrected dN/d%s (MC)", obs.axis.Data()));
                hs.front()->GetYaxis()->SetTitleSize(0.040); hs.front()->GetXaxis()->SetLabelSize(0);
                hs.front()->Draw("hist");
                for (TH1D *h : hs) h->Draw("hist same");
                hs.back()->Draw("PE same");
                if (ht) ht->Draw("hist same");
                lg->Draw();
                TLatex tl; tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.035);
                TString pt_txt = jetPtRangeLabel(jtpt_binsVector[ibin_pt - 1], jtpt_binsVector[ibin_pt]);
                // tl.DrawLatex(0.17, 0.81, "Data, Bayesian unfolding");
                tl.DrawLatex(0.17, 0.81, test_mode == 2 ? TString("Data, Bayesian unfolding")
                                         : split_test   ? TString("Split test, Bayesian unfolding")
                                                        : TString::Format("Full-MC closure (%s), Bayesian", generator.Data()));
                tl.DrawLatex(0.17, 0.76, pt_txt);
                if (step > 1) tl.DrawLatex(0.17, 0.71, Form("every %d iterations shown", step));
                TLatex cmsl; cmsl.SetNDC();
                cmsl.SetTextFont(62); cmsl.SetTextSize(0.050); cmsl.DrawLatex(0.14, 0.94, "CMS");
                cmsl.SetTextFont(52); cmsl.SetTextSize(0.040); cmsl.DrawLatex(0.25, 0.94, "Work In Progress");
                cmsl.SetTextFont(42); cmsl.SetTextSize(0.040); cmsl.DrawLatex(0.66, 0.94, "pp #sqrt{s} = 5.36 TeV");
                pb->cd();
                TH1D *rf = rs.front();
                const double pad = 0.1 * std::max(rhi - 1., 1. - rlo) + 0.01;
                rf->SetMinimum(std::max(0., rlo - pad)); rf->SetMaximum(rhi + pad);
                rf->GetXaxis()->SetRange(ibin_dr_min, ibin_dr_max);
                rf->GetYaxis()->SetTitle(Form("N_{iter} / N_{iter} = %d", n_ref)); rf->GetYaxis()->CenterTitle(true);
                rf->GetYaxis()->SetTitleSize(0.075); rf->GetYaxis()->SetTitleOffset(0.8); rf->GetYaxis()->SetLabelSize(0.065);
                rf->GetYaxis()->SetNdivisions(505);
                rf->GetXaxis()->SetTitle(obs.axis); rf->GetXaxis()->SetTitleSize(0.09); rf->GetXaxis()->SetLabelSize(0.075);
                rf->Draw("hist");
                for (TH1D *r : rs) r->Draw("hist same");
                rs.back()->Draw("E1 same");
                if (rt) rt->Draw("hist same");
                TLine *l1 = new TLine(rf->GetXaxis()->GetBinLowEdge(ibin_dr_min), 1., rf->GetXaxis()->GetBinUpEdge(ibin_dr_max), 1.);
                l1->SetLineStyle(2); l1->Draw();
                const TString it_stem = folder + "iteration_test_1_to_" + Form("%d", n_ref) + "_" + out_label;
                c_it->Print(it_stem + ".pdf"); c_it->Print(it_stem + ".png");
                std::cout << "\t---->iteration test (1.." << n_ref << ", ratio to " << n_ref << "): " << it_stem << ".png" << std::endl;
            }
        }

        std::vector<double> v_cv_n, v_cv_bias, v_cv_stat, v_cv_dn, v_cv_dch;
        const double dchange_threshold = 0.5;
        int best_conv = -1, best_tot = -1;
        double mi_stat_pct = 0.;
        // ---- disabled 2026-10-27 (kept for reference): bias vs MI needs the matrix inversion ----
        // Every quantity here is measured against u_mi = RooUnfoldInvert, which is gone (see the
        // Unfold step). With this block off v_cv_* stay empty and the convergence plot below
        // is skipped by its own !v_cv_dch.empty() guard.
        /*
        if (v_unf_pt.size() == v_niter.size()) {
            RooUnfoldInvert inv(response, h_data_purity_corrected);
            TH2D *h_mi = (TH2D *) inv.Hreco(RooUnfold::kCovariance)->Clone("h_unf_mi_scan");
            TH1D *u_mi = h_mi->ProjectionX("h_unf_pt_mi", ibin_pt, ibin_pt, "e");
            u_mi->SetDirectory(nullptr);
            auto unitArea = [&](TH1D *h) {
                const double s = h->Integral(ibin_dr_min, ibin_dr_max);
                if (s > 0.) h->Scale(1. / s);
            };
            unitArea(u_mi);
            for (TH1D *u : v_unf_pt) unitArea(u);

            // SUPERSEDED (10-19 by AN-2024-133, 10-26 by CMS-SMP-22-007), kept as a cross-check.
            // Was THE CHOICE (Zoe, 2026-10-01): the first niter at which one more iteration moves
            // the result by less than half its statistical error (dchange < 0.5) -- converged
            // at the level of the measurement's precision. Gives 8 for dr, 11 for B on the
            // _rho inputs; 7 for dr, 11 for B with the JEC (_rhojec). The analysis driver reads
            // this line from the scan log and unfolds every Bayesian variation with it. The bias-vs-stat minimum is printed alongside: it sits much higher
            // (41 dr, >50 B) because Bayesian approaches the MI result very slowly, so the
            // distance to MI at the chosen niter is the prior dependence still left in the
            // result -- read it from the "rms bias" column.
            double min_tot = 1e30;
            std::cout << "\n================ Bias vs statistical error (pT bin " << ibin_pt
                      << ", unit area, % of MI) ================\n";
            std::cout << Form("  %6s  %9s %9s %9s   %9s\n", "niter", "rms bias", "rms stat",
                              "rms tot", "dchange");
            for (size_t k = 0; k < v_unf_pt.size(); ++k) {
                double sb = 0., ss = 0., st = 0., sc = 0.; int n = 0;
                for (int i = ibin_dr_min; i <= ibin_dr_max; ++i) {
                    const double m = u_mi->GetBinContent(i);
                    if (m <= 0.) continue;
                    const double b = (v_unf_pt[k]->GetBinContent(i) - m) / m;
                    const double s = v_unf_pt[k]->GetBinError(i) / m;
                    sb += b * b; ss += s * s; st += b * b + s * s;
                    if (k > 0 && v_unf_pt[k]->GetBinError(i) > 0.) {
                        const double c = (v_unf_pt[k]->GetBinContent(i) - v_unf_pt[k - 1]->GetBinContent(i))
                                       / v_unf_pt[k]->GetBinError(i);
                        sc += c * c;
                    }
                    ++n;
                }
                if (n == 0) continue;
                const double tot = std::sqrt(st / n);
                if (tot < min_tot) { min_tot = tot; best_tot = (int) v_niter[k]; }
                if (best_conv < 0 && k > 0 && std::sqrt(sc / n) < dchange_threshold)
                    best_conv = (int) v_niter[k];
                v_cv_n.push_back(v_niter[k]);
                v_cv_bias.push_back(100. * std::sqrt(sb / n));
                v_cv_stat.push_back(100. * std::sqrt(ss / n));
                if (k > 0) { v_cv_dn.push_back(v_niter[k]); v_cv_dch.push_back(std::sqrt(sc / n)); }
                std::cout << Form("  %6.0f  %8.2f%% %8.2f%% %8.2f%%   %9.3f\n", v_niter[k],
                                  100. * std::sqrt(sb / n), 100. * std::sqrt(ss / n),
                                  100. * tot, k > 0 ? std::sqrt(sc / n) : 0.);
            }
            double ss_mi = 0.; int n_mi = 0;
            for (int i = ibin_dr_min; i <= ibin_dr_max; ++i) {
                const double m = u_mi->GetBinContent(i);
                if (m > 0.) { ss_mi += std::pow(u_mi->GetBinError(i) / m, 2); ++n_mi; }
            }
            mi_stat_pct = n_mi ? 100. * std::sqrt(ss_mi / n_mi) : 0.;
            std::cout << Form("  matrix inversion: rms stat = %.2f%%\n", mi_stat_pct)
                      << Form("  ----> cross-check, converged: dchange < %.1f          = %d\n",
                              dchange_threshold, best_conv)
                      << Form("  ----> cross-check, min rms sqrt(bias^2 + stat^2)        = %d\n", best_tot)
                      << "================================================================\n"
                      << std::endl;
        }
        ---- end disabled ---- */

        const int np = (int) v_niter.size();
        TGraph *g_p_all    = new TGraph(np, v_niter.data(), v_pvalue_all.data());
        TGraph *g_p_pt     = new TGraph(np, v_niter.data(), v_pvalue_pt.data());
        TGraph *g_c2_all   = new TGraph(np, v_niter.data(), v_chi2ndf_all.data());
        TGraph *g_c2_pt    = new TGraph(np, v_niter.data(), v_chi2ndf_pt.data());
        // The test that picks the iteration is data (+) refolded stat (purple); data only (blue)
        // and data (+) MC stat (grey) are drawn for reference.
        TGraph *g_p_mc     = new TGraph(np, v_niter.data(), v_pvalue_all_mc.data());
        TGraph *g_c2_mc    = new TGraph(np, v_niter.data(), v_chi2ndf_all_mc.data());
        TGraph *g_p_cov    = new TGraph(np, v_niter.data(), v_pvalue_all_cov.data());
        TGraph *g_c2_cov   = new TGraph(np, v_niter.data(), v_chi2ndf_all_cov.data());
        for (TGraph *g : {g_p_cov, g_c2_cov}) {
            g->SetLineColor(green); g->SetMarkerColor(green); g->SetMarkerStyle(kFullSquare);
            g->SetLineWidth(2); g->SetMarkerSize(0.9);
        }
        TGraph *g_p_r      = new TGraph(np, v_niter.data(), v_pvalue_all_r.data());
        TGraph *g_c2_r     = new TGraph(np, v_niter.data(), v_chi2ndf_all_r.data());
        for (TGraph *g : {g_p_r, g_c2_r}) {
            g->SetLineColor(purple); g->SetMarkerColor(purple); g->SetMarkerStyle(kFullTriangleUp);
            g->SetLineWidth(2); g->SetMarkerSize(0.9);
        }
        TGraph *g_p_rcov   = new TGraph(np, v_niter.data(), v_pvalue_all_rcov.data());
        TGraph *g_c2_rcov  = new TGraph(np, v_niter.data(), v_chi2ndf_all_rcov.data());
        for (TGraph *g : {g_p_rcov, g_c2_rcov}) {
            g->SetLineColor(orange); g->SetMarkerColor(orange); g->SetMarkerStyle(kOpenTriangleUp);
            g->SetLineWidth(2); g->SetMarkerSize(0.9);
        }
        g_p_all ->SetLineColor(blue);   g_p_all ->SetMarkerColor(blue);   g_p_all ->SetMarkerStyle(kFullCircle);
        g_p_pt  ->SetLineColor(orange); g_p_pt  ->SetMarkerColor(orange); g_p_pt  ->SetMarkerStyle(kOpenSquare);
        g_c2_all->SetLineColor(blue);   g_c2_all->SetMarkerColor(blue);   g_c2_all->SetMarkerStyle(kFullCircle);
        g_c2_pt ->SetLineColor(orange); g_c2_pt ->SetMarkerColor(orange); g_c2_pt ->SetMarkerStyle(kOpenSquare);
        for (TGraph *g : {g_p_all, g_p_pt, g_c2_all, g_c2_pt}) { g->SetLineWidth(2); g->SetMarkerSize(0.9); }
        for (TGraph *g : {g_p_mc, g_c2_mc}) {
            g->SetLineColor(blue); g->SetMarkerColor(blue); g->SetMarkerStyle(kFullCircle);
            g->SetLineWidth(2); g->SetMarkerSize(0.9);
        }
        // for (TGraph *g : {g_p_all, g_c2_all}) {
        for (TGraph *g : {g_p_mc, g_c2_mc}) {
            g->SetLineColor(kGray + 1); g->SetMarkerColor(kGray + 1); g->SetMarkerStyle(kOpenCircle);
        }

        const double x_lo_it = v_niter.front() - 0.5;
        const double x_hi_it = v_niter.back()  + 0.5;

        TCanvas *c_reg = new TCanvas("c_regularisation", "", 800, 800);
        TPad *pad_p  = new TPad("pad_reg_p",  "", 0., 0.42, 1., 1.);
        TPad *pad_c2 = new TPad("pad_reg_c2", "", 0., 0.,   1., 0.42);
        for (TPad *p : {pad_p, pad_c2}) { p->SetTicks(1, 1); p->SetFillColor(0); }
        pad_p ->SetMargin(0.13, 0.05, 0.00, 0.08);
        pad_c2->SetMargin(0.13, 0.05, 0.22, 0.00);
        c_reg->cd(); pad_p->Draw(); pad_c2->Draw();

        // ---- top: p-value vs iteration, full 2D reco space only (the space that was unfolded,
        // and the one the selected iteration is taken from). Log y: the first iterations sit
        // many orders of magnitude below the threshold, and the approach to it has to be readable.
        pad_p->cd();
        //pad_p->SetLogy();
        double p_min = 1.;
        for (double p : v_pvalue_all) if (p > 0.) p_min = std::min(p_min, p);
        // One whole decade below the smallest p, so (a) the first iterations stay inside the
        // frame -- with 20 iterations the first p can be ~1e-17 -- and (b) the lowest labelled
        // decade sits one step above the frame edge instead of being clipped by it.
        const double y_lo_p = std::max(std::pow(10., std::floor(std::log10(p_min)) - 1.), 1e-30);
        // p is bounded by 1; the band above it (to 1.45) holds the labels and the legend, which
        // otherwise sit on one of the two rising curves wherever they are put.
        const double y_hi_p = 1.45;
        TH1F *fr_p = pad_p->DrawFrame(x_lo_it, y_lo_p, x_hi_it, y_hi_p);
        fr_p->GetYaxis()->SetTitle("p-value");
        fr_p->GetYaxis()->CenterTitle(true);
        fr_p->GetYaxis()->SetTitleFont(font_code); fr_p->GetYaxis()->SetTitleSize(title_size); fr_p->GetYaxis()->SetTitleOffset(1.5);
        fr_p->GetYaxis()->SetLabelFont(font_code); fr_p->GetYaxis()->SetLabelSize(label_size);
        // Label only every few decades: over ~10 decades one label per power of ten is unreadable.
        // On a log axis the primary-division count is how many decades ROOT puts between labels.
        fr_p->GetYaxis()->SetNdivisions(506);
        // The two pads share an edge, and this pad's LOWEST y label sits exactly on it, so it
        // gets clipped in half (and would collide with the bottom pad's highest label anyway).
        // Hide it: size 0. labNum = 1 is the first (lowest) label.
        fr_p->GetYaxis()->ChangeLabel(1, -1, 0.);
        fr_p->GetXaxis()->SetLabelSize(0); fr_p->GetXaxis()->SetTitleSize(0);  // x belongs to the bottom pad
        // data (+) MC stat (grey) and cov(r - d) (green) off the plot since 2026-10-21 (Zoe);
        // still computed, printed in the log table and written to the .root below.
        g_p_all->Draw("PL SAME");
        // g_p_mc ->Draw("PL SAME");
        // g_p_cov->Draw("PL SAME");
        g_p_r  ->Draw("PL SAME");
        // g_p_rcov->Draw("PL SAME");   // V_d + V_r full cov.: not drawn since 2026-10-22 (Zoe); still printed + in the .root

        TLine *l_thr = new TLine(x_lo_it, an_threshold, x_hi_it, an_threshold);
        l_thr->SetLineColor(kGray + 2); l_thr->SetLineStyle(2); l_thr->Draw("same");
        // Vertical marker at the chosen iteration. TLine takes user coordinates, which on a
        // log pad are the values themselves (not their logs), so the frame range is what to use.
        if (best_an > 0) {
            TLine *l_best = new TLine(best_an, y_lo_p, best_an, y_hi_p);
            l_best->SetLineColor(red); l_best->SetLineStyle(7); l_best->SetLineWidth(2);
            l_best->Draw("same");
        }
        // Label the threshold line where it sits rather than adding a legend. User coordinates
        // so the label tracks the line; it goes above the line on linear y, below on log y.
        TLatex tx_thr;
        tx_thr.SetTextFont(font_code); tx_thr.SetTextSize(label_size); tx_thr.SetTextColor(kGray + 2);
        // Mid-axis, just under the line: clear of both rises in the dr and the B scans.
        // tx_thr.DrawLatex(x_lo_it + 0.50 * (x_hi_it - x_lo_it),
        //                  an_threshold - 0.08 * (y_hi_p - y_lo_p),
        //                  Form("p = %.2f", an_threshold));   // until 2026-10-26 (line at 0.95)
        // Just ABOVE the line for a low threshold (the p = 0.05 rule, under it is outside the frame),
        // just under it for a high one (p = 0.95, the AN rule).
        // tx_thr.DrawLatex(x_lo_it + 0.50 * (x_hi_it - x_lo_it),
        //                  an_threshold + 0.03 * (y_hi_p - y_lo_p),
        //                  Form("p = %.2f", an_threshold));
        tx_thr.DrawLatex(x_lo_it + 0.50 * (x_hi_it - x_lo_it),
                         an_threshold + (an_threshold < 0.5 ? 0.03 : -0.08) * (y_hi_p - y_lo_p),
                         Form("p = %.2f", an_threshold));

        TLatex tx_reg; tx_reg.SetNDC();
        tx_reg.SetTextFont(font_code); tx_reg.SetTextSize(label_size);
        // Right of the rise, in the empty band under the saturated curve: on the left the text
        // sits on top of the part of the curve the plot is about.
        // In the headroom band above p = 1 (see y_hi_p), clear of both curves.
        tx_reg.DrawLatex(0.17, 0.86,
            test_mode == 0 ? "MC full-sample closure" :
            test_mode == 1 ? "MC split closure test"  : "Data Run 3");
        tx_reg.DrawLatex(0.17, 0.80, "D'Agostini, early stopping");
        // tx_reg.DrawLatex(0.17, 0.74, best_an > 0 ? Form("chosen N iter = %d (first p #geq %.2f)", best_an, an_threshold)
        // tx_reg.DrawLatex(0.17, 0.74, best_an > 0 ? Form("chosen N iter = %d (first p > %.2f)", best_an, an_threshold)   // 10-26 afternoon
        tx_reg.DrawLatex(0.17, 0.74, best_an > 0 ? Form("chosen N iter = %d (first p #geq %.2f)", best_an, an_threshold)
                                                 : "no iteration reaches the threshold");
        // TLegend *lg_reg = new TLegend(0.54, 0.58, 0.90, 0.89);   // five curves, until 2026-10-21
        // TLegend *lg_reg = new TLegend(0.54, 0.70, 0.90, 0.89);   // three curves, until 2026-10-22
        TLegend *lg_reg = new TLegend(0.54, 0.76, 0.90, 0.89);
        lg_reg->SetBorderSize(0); lg_reg->SetFillStyle(0);
        lg_reg->SetTextFont(font_code); lg_reg->SetTextSize(label_size);
        // lg_reg->AddEntry(g_p_mc,  "data #oplus MC stat. (used)", "pl");
        // lg_reg->AddEntry(g_p_all, "data stat. only", "pl");
        // lg_reg->AddEntry(g_p_all, "data stat.", "pl");   // until 2026-10-26
        lg_reg->AddEntry(g_p_all, "data stat. (used)", "pl");
        // lg_reg->AddEntry(g_p_mc,  "data #oplus MC stat.", "pl");
        // lg_reg->AddEntry(g_p_cov, "cov(r#minusd), data stat.", "pl");
        // lg_reg->AddEntry(g_p_r,   "data #oplus refolded stat. (used)", "pl");   // until 2026-10-26
        lg_reg->AddEntry(g_p_r,   "data #oplus refolded stat.", "pl");
        // lg_reg->AddEntry(g_p_rcov, "V_{d} + V_{r}, full cov.", "pl");
        lg_reg->Draw();
        TLatex cms_reg; cms_reg.SetNDC();
        cms_reg.SetTextFont(62); cms_reg.SetTextSize(0.042); cms_reg.DrawLatex(0.13, 0.945, "CMS");
        cms_reg.SetTextFont(52); cms_reg.SetTextSize(0.034); cms_reg.DrawLatex(0.205, 0.945, "Work In Progress");
        cms_reg.SetTextFont(42); cms_reg.SetTextSize(0.034); cms_reg.DrawLatex(0.66, 0.945, "pp #sqrt{s} = 5.36 TeV");
        pad_p->RedrawAxis();

        // ---- bottom: the same information as chi2/ndf, where 1 is the target
        pad_c2->cd();
        pad_c2->SetLogy();
        double c2_max = 1., c2_min = 1.;
        for (double c : v_chi2ndf_all) if (c > 0.) { c2_max = std::max(c2_max, c); c2_min = std::min(c2_min, c); }
        // for (double c : v_chi2ndf_all_mc) if (c > 0.) { c2_max = std::max(c2_max, c); c2_min = std::min(c2_min, c); }
        // for (double c : v_chi2ndf_all_cov) if (c > 0.) { c2_max = std::max(c2_max, c); c2_min = std::min(c2_min, c); }
        for (double c : v_chi2ndf_all_r) if (c > 0.) { c2_max = std::max(c2_max, c); c2_min = std::min(c2_min, c); }
        // for (double c : v_chi2ndf_all_rcov) if (c > 0.) { c2_max = std::max(c2_max, c); c2_min = std::min(c2_min, c); }
        TH1F *fr_c2 = pad_c2->DrawFrame(x_lo_it, c2_min * 0.5, x_hi_it, c2_max * 2.);
        fr_c2->GetYaxis()->SetTitle("#chi^{2}/ndf");
        fr_c2->GetYaxis()->CenterTitle(true);
        fr_c2->GetYaxis()->SetTitleFont(font_code); fr_c2->GetYaxis()->SetTitleSize(title_size); fr_c2->GetYaxis()->SetTitleOffset(1.5);
        fr_c2->GetYaxis()->SetLabelFont(font_code); fr_c2->GetYaxis()->SetLabelSize(label_size);
        // Same shared edge from the other side: hide this pad's HIGHEST y label (labNum = -1).
        fr_c2->GetYaxis()->ChangeLabel(-1, -1, 0.);
        fr_c2->GetXaxis()->SetTitle("number of iterations");
        fr_c2->GetXaxis()->CenterTitle(true);
        // Offset in precision-43 fonts scales off the (short) pad height: anything near 2 puts
        // the title off the canvas entirely. 1.3 lands it just below the tick labels.
        fr_c2->GetXaxis()->SetTitleFont(font_code); fr_c2->GetXaxis()->SetTitleSize(title_size); fr_c2->GetXaxis()->SetTitleOffset(1.3);
        fr_c2->GetXaxis()->SetLabelFont(font_code); fr_c2->GetXaxis()->SetLabelSize(label_size);
        g_c2_all->Draw("PL SAME");
        // g_c2_mc ->Draw("PL SAME");
        // g_c2_cov->Draw("PL SAME");
        g_c2_r  ->Draw("PL SAME");
        // g_c2_rcov->Draw("PL SAME");
        TLine *l_one = new TLine(x_lo_it, 1., x_hi_it, 1.);
        l_one->SetLineColor(kGray + 2); l_one->SetLineStyle(2); l_one->Draw("same");
        pad_c2->RedrawAxis();

        // Next to the per-iteration scan outputs, which this plot summarises.
        TString reg_stem = folder + "refolding_pvalue_vs_iteration_" + label;
        c_reg->Print(reg_stem + ".pdf");
        c_reg->Print(reg_stem + ".png");

        // Keep the scan itself, so the choice can be re-plotted without re-unfolding.
        TFile *f_reg = new TFile(reg_stem + ".root", "recreate");
        g_p_all ->Write("g_pvalue_vs_niter_2D");
        g_p_pt  ->Write(Form("g_pvalue_vs_niter_pt%d", ibin_pt));
        g_c2_all->Write("g_chi2ndf_vs_niter_2D");
        g_c2_pt ->Write(Form("g_chi2ndf_vs_niter_pt%d", ibin_pt));
        g_p_mc  ->Write("g_pvalue_vs_niter_2D_datamcstat");
        g_c2_mc ->Write("g_chi2ndf_vs_niter_2D_datamcstat");
        g_p_cov ->Write("g_pvalue_vs_niter_2D_covtoys");
        g_c2_cov->Write("g_chi2ndf_vs_niter_2D_covtoys");
        g_p_r   ->Write("g_pvalue_vs_niter_2D_datarefoldedstat");
        g_c2_r  ->Write("g_chi2ndf_vs_niter_2D_datarefoldedstat");
        g_p_rcov ->Write("g_pvalue_vs_niter_2D_datarefoldedcov");
        g_c2_rcov->Write("g_chi2ndf_vs_niter_2D_datarefoldedcov");
        f_reg->Close();
        delete f_reg;
        std::cout << "Wrote regularisation scan to " << reg_stem << ".{pdf,png,root}" << std::endl;

        // ---- The convergence criterion -- the iteration count that IS used (2026-10-01) -----
        // Top: how much one more iteration moves the result, rms over bins in units of its own
        // stat error, with the 0.5 threshold and the chosen count. Bottom: what the choice
        // trades -- the rms stat error grows with niter while the distance to the unregularised
        // matrix-inversion result (the prior dependence still left) shrinks.
        if (!v_cv_dch.empty()) {
            TGraph *g_dch  = new TGraph((int) v_cv_dn.size(), v_cv_dn.data(), v_cv_dch.data());
            TGraph *g_bias = new TGraph((int) v_cv_n.size(),  v_cv_n.data(),  v_cv_bias.data());
            TGraph *g_stat = new TGraph((int) v_cv_n.size(),  v_cv_n.data(),  v_cv_stat.data());
            g_dch ->SetLineColor(blue);   g_dch ->SetMarkerColor(blue);   g_dch ->SetMarkerStyle(kFullCircle);
            g_stat->SetLineColor(blue);   g_stat->SetMarkerColor(blue);   g_stat->SetMarkerStyle(kFullCircle);
            g_bias->SetLineColor(orange); g_bias->SetMarkerColor(orange); g_bias->SetMarkerStyle(kOpenSquare);
            for (TGraph *g : {g_dch, g_bias, g_stat}) { g->SetLineWidth(2); g->SetMarkerSize(0.9); }

            TCanvas *c_cv = new TCanvas("c_convergence", "", 800, 800);
            TPad *pc_t = new TPad("pad_cv_t", "", 0., 0.42, 1., 1.);
            TPad *pc_b = new TPad("pad_cv_b", "", 0., 0.,   1., 0.42);
            for (TPad *p : {pc_t, pc_b}) { p->SetTicks(1, 1); p->SetFillColor(0); }
            pc_t->SetMargin(0.13, 0.05, 0.00, 0.08);
            pc_b->SetMargin(0.13, 0.05, 0.22, 0.00);
            c_cv->cd(); pc_t->Draw(); pc_b->Draw();

            pc_t->cd(); pc_t->SetLogy();
            double d_min = 1e30, d_max = 0.;
            for (double d : v_cv_dch) if (d > 0.) { d_min = std::min(d_min, d); d_max = std::max(d_max, d); }
            const double y_lo = d_min * 0.5, y_hi = d_max * 3.;
            TH1F *fr = pc_t->DrawFrame(x_lo_it, y_lo, x_hi_it, y_hi);
            fr->GetYaxis()->SetTitle("change per iteration [#sigma_{stat}]");
            fr->GetYaxis()->CenterTitle(true);
            fr->GetYaxis()->SetTitleFont(font_code); fr->GetYaxis()->SetTitleSize(title_size); fr->GetYaxis()->SetTitleOffset(1.5);
            fr->GetYaxis()->SetLabelFont(font_code); fr->GetYaxis()->SetLabelSize(label_size);
            fr->GetYaxis()->ChangeLabel(1, -1, 0.);
            fr->GetXaxis()->SetLabelSize(0); fr->GetXaxis()->SetTitleSize(0);
            g_dch->Draw("PL SAME");
            TLine *l_t = new TLine(x_lo_it, dchange_threshold, x_hi_it, dchange_threshold);
            l_t->SetLineColor(kGray + 2); l_t->SetLineStyle(2); l_t->Draw("same");
            if (best_conv > 0) {
                TLine *l_b = new TLine(best_conv, y_lo, best_conv, y_hi);
                l_b->SetLineColor(kGray + 2); l_b->SetLineStyle(3); l_b->SetLineWidth(2); l_b->Draw("same");
            }
            if (best_an > 0) {
                TLine *l_c = new TLine(best_an, y_lo, best_an, y_hi);
                l_c->SetLineColor(red); l_c->SetLineStyle(7); l_c->SetLineWidth(2); l_c->Draw("same");
            }
            TLatex tx; tx.SetNDC(); tx.SetTextFont(font_code); tx.SetTextSize(label_size);
            tx.DrawLatex(0.45, 0.80, test_mode == 0 ? "MC full-sample closure" :
                                     test_mode == 1 ? "MC split closure test"  : "Data Run 3");
            tx.DrawLatex(0.45, 0.74, Form("D'Agostini, converged at %.1f #sigma_{stat}", dchange_threshold));
            tx.DrawLatex(0.45, 0.68, best_conv > 0 ? Form("converged at N iter = %d (grey)", best_conv)
                                                   : "not converged in the scanned range");
            // if (best_an > 0) tx.DrawLatex(0.45, 0.62, Form("chosen N iter = %d (red, refolding p #geq 0.95)", best_an));   // until 2026-10-26
            // if (best_an > 0) tx.DrawLatex(0.45, 0.62, Form("chosen N iter = %d (red, refolding p > %.2f)", best_an, an_threshold));   // 10-26 afternoon
            if (best_an > 0) tx.DrawLatex(0.45, 0.62, Form("chosen N iter = %d (red, refolding p #geq %.2f)", best_an, an_threshold));
            TLatex cms; cms.SetNDC();
            cms.SetTextFont(62); cms.SetTextSize(0.042); cms.DrawLatex(0.13, 0.945, "CMS");
            cms.SetTextFont(52); cms.SetTextSize(0.034); cms.DrawLatex(0.205, 0.945, "Work In Progress");
            cms.SetTextFont(42); cms.SetTextSize(0.034); cms.DrawLatex(0.66, 0.945, "pp #sqrt{s} = 5.36 TeV");
            pc_t->RedrawAxis();

            pc_b->cd();
            double b_max = mi_stat_pct;
            for (size_t i = 0; i < v_cv_n.size(); ++i) b_max = std::max({b_max, v_cv_bias[i], v_cv_stat[i]});
            TH1F *fb = pc_b->DrawFrame(x_lo_it, 0., x_hi_it, 1.15 * b_max);
            fb->GetYaxis()->SetTitle("rms over bins [%]");
            fb->GetYaxis()->CenterTitle(true);
            fb->GetYaxis()->SetTitleFont(font_code); fb->GetYaxis()->SetTitleSize(title_size); fb->GetYaxis()->SetTitleOffset(1.5);
            fb->GetYaxis()->SetLabelFont(font_code); fb->GetYaxis()->SetLabelSize(label_size);
            fb->GetYaxis()->ChangeLabel(-1, -1, 0.);
            fb->GetXaxis()->SetTitle("number of iterations");
            fb->GetXaxis()->CenterTitle(true);
            fb->GetXaxis()->SetTitleFont(font_code); fb->GetXaxis()->SetTitleSize(title_size); fb->GetXaxis()->SetTitleOffset(1.3);
            fb->GetXaxis()->SetLabelFont(font_code); fb->GetXaxis()->SetLabelSize(label_size);
            g_stat->Draw("PL SAME");
            g_bias->Draw("PL SAME");
            TLine *l_mi = new TLine(x_lo_it, mi_stat_pct, x_hi_it, mi_stat_pct);
            l_mi->SetLineColor(kGray + 2); l_mi->SetLineStyle(2); l_mi->Draw("same");
            if (best_conv > 0) {
                TLine *l_b2 = new TLine(best_an > 0 ? best_an : best_conv, 0., best_an > 0 ? best_an : best_conv, 1.15 * b_max);
                l_b2->SetLineColor(red); l_b2->SetLineStyle(7); l_b2->SetLineWidth(2); l_b2->Draw("same");
            }
            TLegend *lg = new TLegend(0.50, 0.60, 0.93, 0.95);
            lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextFont(font_code); lg->SetTextSize(label_size);
            lg->AddEntry(g_stat, "stat. error (Bayesian)", "pl");
            lg->AddEntry(g_bias, "distance to matrix inversion", "pl");
            lg->AddEntry(l_mi,   "stat. error (matrix inversion)", "l");
            lg->Draw();
            pc_b->RedrawAxis();

            const TString cv_stem = folder + "convergence_vs_iteration_" + label;
            c_cv->Print(cv_stem + ".pdf");
            c_cv->Print(cv_stem + ".png");
            std::cout << "Wrote convergence scan to " << cv_stem << ".{pdf,png}" << std::endl;
        }
    }

    // gApplication -> Terminate(0);
}

// SAMPLE:      qcd | bjet | both   ("both" sums the bjet response, corrections and truth
//              onto the qcd ones -- counts and RooUnfoldResponse are added, the
//              purity/efficiency ratios are recomputed from the summed counts)
// There are two independent generator choices, one per half of the chain:
//
// UNFOLDING_GENERATOR: pythia | herwig
//              Everything the unfolding uses: the migration matrix, the purity and
//              efficiency corrections, the truth reference, h_mc_reco, the hgenjet_2b_*
//              corrections, and in modes 0/1 the MC being unfolded. Read from the
//              per-block output of run_agg_ntuple_chunks.sh, so no pre-merged file.
// TF_GENERATOR: pythia | herwig
//              The template fit that supplies the signal fraction (data mode only).
//              Independent of UNFOLDING_GENERATOR: only a Pythia fit exists today, so
//              a Herwig unfolding still takes its signal fraction from the Pythia fit.
//              Setting it to herwig stops with a message until that fit is produced.
// test_mode:   0 = full-MC closure, 1 = split test, 2 = data
// unfoldBayes: true = Bayesian, false = matrix inversion.
//              DEFAULT false: matrix inversion is the unfolding of BOTH observables. It has
//              no regularisation and so no iteration count to justify per observable -- the
//              Bayesian run needed niter = 7 for dr and 21-26 for B, two numbers tuned on two
//              refolding scans, and a band built from variations unfolded at a regularisation
//              strength of their own is not a band. Pass true to get the Bayesian result back.
//              SINCE 2026-10-27: matrix inversion is commented out (non-square 3D response
//              with under/overflow), the default is true and false is refused.
// scan_niter:  true = scan niter = 1..100 and pick the optimal one from the refolding
//              goodness-of-fit test; false = one unfolding at the nominal niter = 4.
//              Bayesian only -- matrix inversion has no iterations to scan, so it is refused
//              with unfoldBayes = false rather than silently ignored.
// TRACK_EFF_UNC: false (nominal) | true
//              Which MC production to unfold with. true reads the files written with 3%
//              of the reconstructed tracks thrown away during the B reconstruction
//              (TRACK_EFF_UNC=true in run_agg_ntuple_chunks.sh) -- the tracking-efficiency
//              systematic. It varies BOTH halves of the MC: the response/corrections AND
//              the template fit that supplies the signal fraction, so the propagation is
//              complete -- which means template_fit.cpp(SAMPLE, TF_GENERATOR, true) has to
//              have been run first. Only the MC changes; the data being unfolded and the
//              data being fitted are the same files either way.
//              Results go to their own folder, so the nominal is never overwritten.
// TF_VARIATION: nominal | var0B_2 | var0B_0
//              Which of template_fit.cpp's variations supplies the signal fraction (data
//              mode only). var0B_2 / var0B_0 are the LIGHT-JET MISTAG variations: the 0B
//              template -- the jets with no gen b hadron that passed the b tag, i.e.
//              mistagged light and charm -- has its contribution to the background PDF
//              doubled or removed, and the fit is redone. Both fits already exist in every
//              template-fit directory (template_fit.cpp writes one file per variation), so
//              this costs one more unfolding and no refit. Results go to their own folder.
//              Unlike TRACK_EFF_UNC this varies ONLY the signal fraction: the response and
//              the corrections are the nominal ones, which is right -- the mistag rate
//              affects what fraction of the data is signal, not how signal migrates.
//
// e.g.  root -l -b -q 'apply_unfolding_2d.C("both","pythia")'   // nominal
//       root -l -b -q 'apply_unfolding_2d.C("both","herwig")'   // Herwig unfolding MC
//       root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",true)'
//                                                             // tracking-efficiency variation
//       root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",false,"var0B_2")'
//                                                             // light-jet mistag up
// SFUPART_VARIATION: nominal | off | jpcalib_hf | qqrate_up | qqrate_down | statup | statdn
//                    | jpsyst_up | jpsyst_dn
//              Which UParT b-tag efficiency scale factor to divide into the RECO-level data
//              (data mode only). The SF is measured in reco R_BB bins, so it is applied
//              before unfolding and the response matrix carries it through the migration.
//              nominal is the central SF -- it is a CORRECTION, applied in every run, not an
//              opt-in. The rest are its variations: jpcalib_hf swaps in the JP heavy-flavour
//              calibration curve, qqrate_up/down move the qq rate by +/-25%, statup/statdn
//              shift the central SF by +/- its statistical error, and jpsyst_up/jpsyst_dn
//              shift it by +/- the quoted JP HF systematic.
//              Each lands in its own results folder.
//
//              THE TWO JP HF TREATMENTS ARE DIFFERENT AND BOTH ARE KEPT:
//                jpcalib_hf            = the blue curve as measured. Its shift flips sign
//                                        across dr, so it reshapes and survives the unit-area
//                                        normalisation.
//                jpsyst_up/jpsyst_dn   = +/- |jphf - central| coherently in all dr bins. Same
//                                        magnitudes, one sign; closer to a normalisation, so
//                                        much of it cancels in the normalisation.
//              Same input numbers used two ways -- never put both in the band at once.
void apply_unfolding_2d(TString SAMPLE = "both", TString UNFOLDING_GENERATOR = "pythia",
                        // int test_mode = 2, bool unfoldBayes = false, bool scan_niter = false,   // until 2026-10-27
                        int test_mode = 2, bool unfoldBayes = true, bool scan_niter = false,
                        // "pythia", the nominal. It was "herwig" until 2026-09-24, so the
                        // documented nominal call apply_unfolding_2d.C("both","pythia") in
                        // fact wrote the _tfherwig VARIATION, while the reader
                        // (apply_weights_and_systematics.C) defaults to "pythia".
                        TString TF_GENERATOR = "pythia", bool TRACK_EFF_UNC = false,
                        TString TF_VARIATION = "nominal",
                        TString SFUPART_VARIATION = "nominal",
                        bool EEC_WEIGHT_OFF = false,
                        TString OBSERVABLE = "dr",
                        int NITER = 7){

    // Kept so the old positional name still reads naturally below.
    TString GENERATOR = UNFOLDING_GENERATOR;

    // A typo must not fall through to the default sample: check before reading anything.
    if (SAMPLE != "qcd" && SAMPLE != "bjet" && SAMPLE != "both") {
        std::cerr << "ERROR: unknown SAMPLE '" << SAMPLE << "' (use qcd | bjet | both)" << std::endl;
        return;
    }
    if (GENERATOR != "pythia" && GENERATOR != "herwig") {
        std::cerr << "ERROR: unknown UNFOLDING_GENERATOR '" << GENERATOR
                  << "' (use pythia | herwig)" << std::endl;
        return;
    }
    if (TF_GENERATOR != "pythia" && TF_GENERATOR != "herwig") {
        std::cerr << "ERROR: unknown TF_GENERATOR '" << TF_GENERATOR
                  << "' (use pythia | herwig)" << std::endl;
        return;
    }
    if (!isKnownObservable(OBSERVABLE)) {
        std::cerr << "ERROR: unknown OBSERVABLE '" << OBSERVABLE << "' (use dr | B)" << std::endl;
        return;
    }
    // Matrix inversion has nothing to regularise, so there is no iteration count to scan.
    // Refuse rather than ignore the flag: a scan run that quietly did one inversion would
    // still write a "regularisation" plot, and that plot would mean nothing.
    // if (scan_niter && !unfoldBayes) {   // until 2026-10-27
    //     std::cerr << "ERROR: scan_niter needs unfoldBayes = true -- matrix inversion has no "
    //               << "iterations to scan." << std::endl;
    //     return;
    // }
    // Matrix inversion is commented out (2026-10-27, see the Unfold step in apply_unfolding):
    // the 3D (jet pT, dR, A) response with under/overflow is non-square. Refuse instead of
    // writing a "_MI" result that is not one.
    if (!unfoldBayes) {
        std::cerr << "ERROR: unfoldBayes = false -- matrix inversion is disabled (non-square 3D "
                  << "response). Use unfoldBayes = true." << std::endl;
        return;
    }
    // The UParT scale factor is measured in reco dR bins and has no B equivalent, so there
    // is nothing to apply for B. This is a KNOWN GAP, not a silent skip: the dr result
    // carries a correction (0.933 -> 1.149 across its range) that a B result does not.
    // Refuse an explicit SF request rather than pretend to honour it.
    /* ---- disabled 2026-10-19 (kept for reference): A now has its own UParT SF ----
    if (OBSERVABLE == "B" && SFUPART_VARIATION != "nominal" && SFUPART_VARIATION != "off") {
        std::cerr << "ERROR: SFUPART_VARIATION '" << SFUPART_VARIATION << "' is meaningless "
                  << "for OBSERVABLE 'B': the UParT SF is binned in reco dR and has no B "
                  << "equivalent. Use \"off\"." << std::endl;
        return;
    }
    ---- end disabled ---- */
    // What the nominal run of this observable is -- both rules live in result_paths.h, so the
    // reader (apply_weights_and_systematics.C) builds the same paths from the same two
    // functions instead of reproducing the forcing below from memory.
    if (SFUPART_VARIATION == "nominal" && nominalSfupart(OBSERVABLE) == "off") {
        std::cout << "NOTE: OBSERVABLE '" << OBSERVABLE << "' -- no UParT SF is applied (none "
                  << "exists in " << OBSERVABLE << "). Forcing SFUPART_VARIATION to \"off\" so "
                  << "the results path says so." << std::endl;
        SFUPART_VARIATION = "off";
    }
    // B is a yield measurement: the EEC weight is off, so the result is dN/dB. Forced, not
    // assumed, because a weighted-B production does exist -- and a B result unfolded with the
    // weight on would land in its own folder, be read by nothing, and leave the B band
    // looking like a missing file.
    if (!EEC_WEIGHT_OFF && nominalEecWeightOff(OBSERVABLE)) {
        std::cout << "NOTE: OBSERVABLE '" << OBSERVABLE << "' is measured WITHOUT the EEC "
                  << "weight (dN/d" << OBSERVABLE << "). Forcing EEC_WEIGHT_OFF to true -- see "
                  << "nominalEecWeightOff() in result_paths.h." << std::endl;
        EEC_WEIGHT_OFF = true;
    }
    if (TF_VARIATION != "nominal" && TF_VARIATION != "var0B_2" && TF_VARIATION != "var0B_0") {
        std::cerr << "ERROR: unknown TF_VARIATION '" << TF_VARIATION
                  << "' (use nominal | var0B_2 | var0B_0)" << std::endl;
        return;
    }
    // A varied signal fraction only ever multiplies real data, so in the closure modes it
    // would change nothing while still writing to a separate folder -- which reads as a
    // result and is not one.
    if (TF_VARIATION != "nominal" && test_mode != 2) {
        std::cerr << "ERROR: TF_VARIATION '" << TF_VARIATION << "' needs test_mode 2 (data): "
                  << "the signal fraction is only applied to data" << std::endl;
        return;
    }
    // Same for TF_SAMPLE: it only picks the signal fraction, so a closure run would write a
    // tagged folder that is just the nominal closure.
    // TF_METHOD is allowed in the closure modes (2026-10-28): it does not touch them, but its chain
    // runs at its own niter, so its closures need their own folder next to its results.
    if (tfSampleTag().Length() && test_mode != 2) {
        std::cerr << "ERROR: TF_SAMPLE needs test_mode 2 (data): the signal fraction is only "
                  << "applied to data" << std::endl;
        return;
    }
    if (sfupartHist(SFUPART_VARIATION, OBSERVABLE).Length() == 0) {
        std::cerr << "ERROR: unknown SFUPART_VARIATION '" << SFUPART_VARIATION
                  << "' (use nominal | off | jpcalib_hf | qqrate_up | qqrate_down | statup | "
                  << "statdn | jpsyst_up | jpsyst_dn)" << std::endl;
        return;
    }
    // Same reasoning as TF_VARIATION: the SF is only ever applied to data.
    if (SFUPART_VARIATION != "nominal" && SFUPART_VARIATION != "off" && test_mode != 2) {
        std::cerr << "ERROR: SFUPART_VARIATION '" << SFUPART_VARIATION << "' needs test_mode 2 (data): "
                  << "the UParT SF is only applied to data" << std::endl;
        return;
    }
    // Refuse before creating anything, as above: a missing SF file would otherwise surface
    // only in the log, after the results folder exists.
    if (test_mode == 2 && SFUPART_VARIATION != "off") {
        if (gSystem->AccessPathName(sfupartFile(OBSERVABLE))) {
            std::cerr << "ERROR: missing UParT SF file " << sfupartFile(OBSERVABLE) << std::endl;
            return;
        }
    }

    // Refuse before creating anything: past this point stdout goes to the log file, so a
    // failure inside apply_unfolding() would leave an empty results folder and a silent
    // terminal.
    if (test_mode == 2) {
        const TString tf_file = templateFitFile(SAMPLE, TF_GENERATOR, TRACK_EFF_UNC, TF_VARIATION,
                                                EEC_WEIGHT_OFF, OBSERVABLE);
        if (tf_file.Length() == 0 || gSystem->AccessPathName(tf_file)) {
            std::cerr << "ERROR: no template fit for TF_GENERATOR '" << TF_GENERATOR
                      << "'" << (TRACK_EFF_UNC ? " with TRACK_EFF_UNC" : "")
                      << ", variation '" << TF_VARIATION << "'"
                      << ": " << tf_file << "\n       run template_fit.cpp+(\""
                      << ((SAMPLE == "both") ? "both" : "qcd") << "\",\"" << TF_GENERATOR
                      << "\"," << (TRACK_EFF_UNC ? "true" : "false") << ","
                      << (EEC_WEIGHT_OFF ? "true" : "false") << ",\"" << OBSERVABLE
                      << "\") first" << std::endl;
            return;
        }
    }

    // Tag only what actually differs from the nominal, so nominal names are unchanged.
    // Both built by result_paths.h, the same functions apply_weights_and_systematics.C uses
    // to FIND them. "" tags mean nominal, so nominal paths are exactly what they always
    // were; a variation lands in its own directory next to the nominal, never over it.
    TString dataset = SAMPLE + "_" + GENERATOR
                    + variationTag(TF_GENERATOR, TRACK_EFF_UNC, TF_VARIATION, SFUPART_VARIATION,
                                   EEC_WEIGHT_OFF, OBSERVABLE);

    TString folder = resultFolder(SAMPLE, GENERATOR, unfoldBayes, TF_GENERATOR, TRACK_EFF_UNC,
                                  TF_VARIATION, SFUPART_VARIATION, EEC_WEIGHT_OFF, OBSERVABLE);

    // An unrecognised variation name poisons the tag rather than silently becoming the
    // nominal (see result_paths.h). Catch it here, before a directory is created.
    if (folder.Contains(kUnknownTag())) {
        std::cerr << "ERROR: refusing to run -- one of the variation names is not recognised."
                  << std::endl;
        return;
    }

    gSystem->mkdir(folder, kTRUE);
    TString pT_selection = "80_inf";
    bool btag = true;
    Int_t n = 1;

    // ---- Log file ----
    TDatime now;
    TString timestamp = Form(
        "%04d-%02d-%02d_%02d-%02d-%02d",
        now.GetYear(),
        now.GetMonth(),
        now.GetDay(),
        now.GetHour(),
        now.GetMinute(),
        now.GetSecond()
    );
    TString label;
    label += unfoldBayes ? "_bayesian" : "_MI";
    if(test_mode==0)  label += "_full_closure";
    if(test_mode==1)  label += "_split_test";
    if(test_mode==2)  label += "_data";
    if(scan_niter)    label += "_iterative";
    TString logfile = folder +  "unfolding_" + dataset + label + "_" +  timestamp + ".log";
    gSystem->RedirectOutput(logfile);  // "a" = append
	apply_unfolding(dataset, folder, btag, n, pT_selection, test_mode, unfoldBayes, scan_niter,
	                SAMPLE, GENERATOR, TF_GENERATOR, TRACK_EFF_UNC, TF_VARIATION, SFUPART_VARIATION,
	                EEC_WEIGHT_OFF, OBSERVABLE, NITER);

    // Restore terminal output
    gSystem->RedirectOutput(nullptr);

    // Close the loop between writer and reader: resultFile() is the path
    // apply_weights_and_systematics.C will look for. Check here, at write time, that it
    // actually resolves -- otherwise a mismatch only shows up much later as a variation
    // that is "missing" or, worse, as one that silently resolved to something else.
    TString expected = resultFile(SAMPLE, GENERATOR, test_mode, unfoldBayes,
                                  TF_GENERATOR, TRACK_EFF_UNC, TF_VARIATION,
                                  SFUPART_VARIATION, EEC_WEIGHT_OFF, OBSERVABLE);
    // A scan writes the "_iterative" file instead (out_label in apply_unfolding), which the
    // systematics macro never reads -- check that one, so the scan does not report a missing result.
    if (scan_niter) expected.ReplaceAll("_after_unfolding_2D.root", "_iterative_after_unfolding_2D.root");
    if (gSystem->AccessPathName(expected)) {
        std::cerr << "WARNING: the run finished but its result is not where the systematics "
                  << "macro will look for it:\n         " << expected
                  << "\n         Check " << logfile << " -- the unfolding may have returned "
                  << "early, or the naming in result_paths.h has drifted from what "
                  << "apply_unfolding() writes." << std::endl;
    } else {
        std::cout << "Result: " << expected << std::endl;
    }

}

// ============================================================================
// 3D Bayesian unfolding in (dR, A, jet pT) -- MC closures only (added 2026-10-27)
// ============================================================================
//
// Reads the Response3D objects of create_files_for_template_fit.cpp (the *_3d names): a
// NON-SQUARE response, 450 reco bins (pT 80|90|100|125|150|200, dR 0.025, A 0.1) x 200 gen bins
// (pT 60-80 | 80-100 | 100-200 | 200+, dR 0.05 plus a 0.45+ overflow, A 0.1) since 2026-10-29
// (before: 360 reco bins, pT 80-120 in 10 GeV x gen 60|80|100|120|1000). The side bins
// are at gen only -- see u3d_* in binning_histos_small.h.
//
// MC only for now: there is no 3D signal fraction (the template fit is per 2D bin) and no 3D
// data histogram yet.
//   test_mode 0  full-MC closure: the whole reco sample through the whole-sample response.
//                Closes exactly at every niter (the prior IS the truth) -- a plumbing check.
//   test_mode 1  split test: even-half reco (h_pseudodata_reco_3d) through the odd-half
//                response, purity and efficiency, against the even-half truth.
// Chain as in apply_unfolding(): x purity (reco bins) -> Bayes -> / efficiency (gen bins).
// Both corrections rescale content and error by the same factor (MC treated as exact).
//
// EEC_WEIGHT_OFF picks the production: false = the EEC-weighted one (dr's), true = the
// unweighted _noeecw one (A's). Which weight the 3D result uses is still open (2026-10-27).
// NITER > 0 unfolds at that count. NITER = 0 scans 1..SCAN_NITER_MAX (env, default 200) and
// takes the AN-2024-133 rule -- the first refolding p >= 0.95, input stat only -- over the
// whole 3D reco space. The scan runs without error propagation; errors at the chosen niter only.
//
// Writes $P/results/unfolding_3d/<sample>_<gen><production>/: the 3D histograms, the scan, and
// a PDF with the dR and A projections in gen pT 80-100 and 100-200 (dR < 0.45 only: the 0.45+
// bin is a side bin) against the truth.
//
// e.g.  root -l -b -q -e '.L apply_unfolding_2d.C' -e 'apply_unfolding_3d("both","pythia",1)'
void apply_unfolding_3d(TString SAMPLE = "both", TString GENERATOR = "pythia", int test_mode = 1,
                        bool EEC_WEIGHT_OFF = false, int NITER = 0)
{
    if (test_mode != 0 && test_mode != 1) {
        std::cerr << "ERROR: apply_unfolding_3d is MC only (test_mode 0 | 1): there is no 3D "
                  << "signal fraction or 3D data histogram yet." << std::endl;
        return;
    }
    if (SAMPLE != "qcd" && SAMPLE != "bjet" && SAMPLE != "both") {
        std::cerr << "ERROR: unknown SAMPLE '" << SAMPLE << "' (use qcd | bjet | both)" << std::endl;
        return;
    }
    gSystem->Load("libRooUnfold.so");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    const bool split = (test_mode == 1);

    // ---- Inputs: the same block files and qcd+bjet mix as apply_unfolding() ----
    const TString mc_base  = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024";
    const TString btag_tag = "btagWP0712";
    const TString obs_name = EEC_WEIGHT_OFF ? "B" : "dr";   // only picks the production
    const TString out_tag  = mcVarTag(false, EEC_WEIGHT_OFF) + obsProdTag(obs_name);
    std::vector<TString> names = aggChunkFiles(mc_base, SAMPLE, GENERATOR, "RMatrix_", "", btag_tag, out_tag);
    if (names.empty()) return;
    gMixWeight.clear();
    if (!setMixWeights(mc_base, SAMPLE, GENERATOR, btag_tag, out_tag, EEC_WEIGHT_OFF ? obsB() : obsDr())) return;
    std::vector<TFile *> fin = openAllOrWarn(names);
    if (fin.empty()) return;

    const TString ps = split ? "_pseudo" : "";
    TH3D *pur_num = sumOverFiles<TH3D>(fin, "h_full" + ps + "_purity_numerator_3d");
    TH3D *pur_den = sumOverFiles<TH3D>(fin, "h_full" + ps + "_purity_denominator_3d");
    TH3D *eff_num = sumOverFiles<TH3D>(fin, "h_full" + ps + "_efficiency_numerator_3d");
    TH3D *eff_den = sumOverFiles<TH3D>(fin, "h_full" + ps + "_efficiency_denominator_3d");
    TH3D *h_in    = sumOverFiles<TH3D>(fin, split ? "h_pseudodata_reco_3d"  : "h_full_purity_denominator_3d");
    TH3D *h_truth = sumOverFiles<TH3D>(fin, split ? "h_pseudodata_truth_3d" : "h_full_efficiency_denominator_3d");
    RooUnfoldResponse *response = sumResponseOverFiles(fin, split ? "response_pseudo_full_3d" : "response_full_3d");
    if (!pur_num || !pur_den || !eff_num || !eff_den || !h_in || !h_truth || !response) {
        std::cerr << "ERROR: 3D objects missing -- was step 1 re-run after 2026-10-27?" << std::endl;
        return;
    }
    std::cout << "3D response: " << response->GetNbinsMeasured() << " reco x "
              << response->GetNbinsTruth() << " gen bins" << std::endl;

    // ---- The full (non-square) response, P(detector | particle) (2026-10-28, Zoe) ----
    // RooUnfold flattens a TH3 with x fastest: global bin = dR + n_dR * (A + n_A * pT), so both
    // axes are 4 pT blocks (thick lines), each holding 5 A blocks (thin lines) of dR bins. Rows
    // (gen bins) normalised to 1 as in plot_purity_efficiency_response.cpp, colour scale 0-1.
    {
        const TString rfolder = mc_base + "/results/unfolding_3d/" + SAMPLE + "_" + GENERATOR + out_tag + "/";
        gSystem->mkdir(rfolder, kTRUE);
        TH2 *hR = (TH2 *) response->Hresponse()->Clone("h_response_3d_rownorm");
        hR->SetDirectory(nullptr);
        for (int iy = 1; iy <= hR->GetNbinsY(); ++iy) {
            double row = 0.;
            for (int ix = 1; ix <= hR->GetNbinsX(); ++ix) row += hR->GetBinContent(ix, iy);
            if (row <= 0.) continue;
            for (int ix = 1; ix <= hR->GetNbinsX(); ++ix) {
                hR->SetBinContent(ix, iy, hR->GetBinContent(ix, iy) / row);
                hR->SetBinError  (ix, iy, hR->GetBinError(ix, iy)   / row);
            }
        }
        // same blue -> purple -> red scale as the 2D response plots
        {
            Double_t stops[3] = {0.00, 0.50, 1.00};
            Double_t rr[3] = {0x4C / 255., 0x8C / 255., 0xC4 / 255.};
            Double_t gg[3] = {0x72 / 255., 0x6B / 255., 0x4E / 255.};
            Double_t bb[3] = {0xB0 / 255., 0xB1 / 255., 0x52 / 255.};
            TColor::CreateGradientColorTable(3, stops, rr, gg, bb, 255);
            gStyle->SetNumberContours(255);
        }
        const int ndr_r = u3d_reco_dr_size - 1, ndr_g = u3d_gen_dr_size - 1, nA = B_binsVectorSize - 1;
        const int npt_r = u3d_reco_pt_size - 1, npt_g = u3d_gen_pt_size - 1;
        hR->SetTitle("");
        hR->GetXaxis()->SetTitle(Form("Detector-level bin (%d = %d #DeltaR #times %d A #times %d p_{T})",
                                      hR->GetNbinsX(), ndr_r, nA, npt_r));
        hR->GetYaxis()->SetTitle(Form("Particle-level bin (%d = %d #DeltaR #times %d A #times %d p_{T})",
                                      hR->GetNbinsY(), ndr_g, nA, npt_g));
        hR->GetZaxis()->SetTitle("P(detector | particle)");
        hR->GetXaxis()->SetTitleSize(0.035); hR->GetYaxis()->SetTitleSize(0.035);
        hR->GetYaxis()->SetTitleOffset(1.3); hR->GetZaxis()->SetTitleOffset(1.3);
        // hR->SetMinimum(0.); hR->SetMaximum(1.);   // until 2026-10-28 evening
        // hR->SetMinimum(0.); hR->SetMaximum(0.5);   // Zoe, then: 99% of the cells are below 0.063
        hR->SetMinimum(0.); hR->SetMaximum(0.1);      // Zoe: 0-0.1, the few cells above saturate
        {   // what the cells hold, to choose the scale
            std::vector<double> v;
            for (int ix = 1; ix <= hR->GetNbinsX(); ++ix) for (int iy = 1; iy <= hR->GetNbinsY(); ++iy)
                if (hR->GetBinContent(ix, iy) > 0.) v.push_back(hR->GetBinContent(ix, iy));
            std::sort(v.begin(), v.end());
            if (!v.empty()) {
                double rowmax = 0.; int nrow = 0;
                for (int iy = 1; iy <= hR->GetNbinsY(); ++iy) {
                    double m = 0.; for (int ix = 1; ix <= hR->GetNbinsX(); ++ix) m = std::max(m, hR->GetBinContent(ix, iy));
                    if (m > 0.) { rowmax += m; ++nrow; }
                }
                std::cout << Form("3D response cells: %zu non-empty, max %.3f, 99%% below %.3f, 90%% below %.3f, median %.4f; mean of the row maxima %.3f",
                                  v.size(), v.back(), v[(size_t) (0.99 * (v.size() - 1))], v[(size_t) (0.90 * (v.size() - 1))], v[v.size() / 2], rowmax / std::max(1, nrow)) << std::endl;
            }
        }
        TCanvas *cr = new TCanvas("c_resp3d", "", 1300, 950);
        cr->SetRightMargin(0.14); cr->SetLeftMargin(0.10); cr->SetTopMargin(0.09); cr->SetBottomMargin(0.10);
        cr->SetFrameFillColor(kGray + 1);
        hR->Draw("COLZ");
        const double x0 = hR->GetXaxis()->GetXmin(), x1 = hR->GetXaxis()->GetXmax();
        const double y0 = hR->GetYaxis()->GetXmin(), y1 = hR->GetYaxis()->GetXmax();
        const double dx = (x1 - x0) / hR->GetNbinsX(), dy = (y1 - y0) / hR->GetNbinsY();
        for (int k = 1; k < npt_r * nA; ++k) {   // A blocks thin, pT blocks thick
            TLine *l = new TLine(x0 + k * ndr_r * dx, y0, x0 + k * ndr_r * dx, y1);
            l->SetLineColor(kGray + 3); l->SetLineWidth(k % nA ? 1 : 3); l->SetLineStyle(k % nA ? 3 : 1); l->Draw();
        }
        for (int k = 1; k < npt_g * nA; ++k) {
            TLine *l = new TLine(x0, y0 + k * ndr_g * dy, x1, y0 + k * ndr_g * dy);
            l->SetLineColor(kGray + 3); l->SetLineWidth(k % nA ? 1 : 3); l->SetLineStyle(k % nA ? 3 : 1); l->Draw();
        }
        TLatex tl; tl.SetTextFont(42); tl.SetTextSize(0.022); tl.SetTextAlign(21);
        for (int k = 0; k < npt_r; ++k)   // reco pT of each column block, above the frame
            tl.DrawLatex(x0 + (k + 0.5) * nA * ndr_r * dx, y1 + 2. * dy,
                         Form("reco %g-%g GeV", u3d_reco_pt[k], u3d_reco_pt[k + 1]));
        tl.SetTextAlign(12); tl.SetTextColor(kWhite);
        for (int k = 0; k < npt_g; ++k) {  // gen pT of each row block, inside the frame on the left
            const bool of = (k == npt_g - 1), uf = (k == 0);
            tl.DrawLatex(x0 + 3. * dx, y0 + (k + 0.5) * nA * ndr_g * dy,
                         of ? Form("gen > %g GeV (overflow)", u3d_gen_pt[k])
                            : Form("gen %g-%g GeV%s", u3d_gen_pt[k], u3d_gen_pt[k + 1], uf ? " (underflow)" : ""));
        }
        TLatex th; th.SetNDC(); th.SetTextSize(0.03);
        th.SetTextFont(62); th.DrawLatex(0.10, 0.955, "CMS");
        th.SetTextFont(52); th.DrawLatex(0.16, 0.955, "Work In Progress");
        th.SetTextFont(42); th.SetTextAlign(31);
        th.DrawLatex(0.86, 0.955, Form("3D response (%s, %s), %s", SAMPLE.Data(), GENERATOR.Data(),
                                       EEC_WEIGHT_OFF ? "no EEC weight" : "EEC-weighted"));
        cr->RedrawAxis();
        const TString rstem = rfolder + (split ? "response_matrix_3d_pseudo" : "response_matrix_3d");
        cr->Print(rstem + ".pdf");
        cr->Print(rstem + ".png");
        std::cout << "Wrote " << rstem << ".{pdf,png}" << std::endl;

        // ---- (2026-10-29, Zoe's slides) two closer looks at the same matrix ----
        auto header = [&](const char *what) {
            TLatex t; t.SetNDC(); t.SetTextSize(0.03);
            t.SetTextFont(62); t.DrawLatex(0.10, 0.955, "CMS");
            t.SetTextFont(52); t.DrawLatex(0.16, 0.955, "Work In Progress");
            t.SetTextFont(42); t.SetTextAlign(31); t.DrawLatex(0.86, 0.955, what);
        };
        // (a) zoom on the measured pT: the reco blocks from 100 GeV up x the gen 100-X bin, their A blocks.
        // Since 2026-10-29 generic in the edges (was reco 100-120 = blocks 3-4 x gen 100-120).
        {
            int rb0 = 0; while (rb0 < npt_r && u3d_reco_pt[rb0] < 100. - 1e-6) ++rb0;   // first reco block >= 100 GeV
            const int nrb = npt_r - rb0;                                              // reco blocks in the zoom
            const int rx0 = rb0 * nA * ndr_r, gy0 = 2 * nA * ndr_g;                   // gen block 3 = the nominal one
            TH2 *hz = (TH2 *) hR->Clone("h_response_3d_zoom"); hz->SetDirectory(nullptr);
            hz->GetXaxis()->SetRange(rx0 + 1, rx0 + nrb * nA * ndr_r);
            hz->GetYaxis()->SetRange(gy0 + 1, gy0 + nA * ndr_g);
            TCanvas *cz = new TCanvas("c_resp3d_zoom", "", 1300, 950);
            cz->SetRightMargin(0.14); cz->SetLeftMargin(0.10); cz->SetTopMargin(0.09); cz->SetBottomMargin(0.10);
            cz->SetFrameFillColor(kGray + 1);
            hz->Draw("COLZ");
            const double xa = x0 + rx0 * dx, xb = x0 + (rx0 + nrb * nA * ndr_r) * dx;
            const double ya = y0 + gy0 * dy, yb = y0 + (gy0 + nA * ndr_g) * dy;
            for (int k = 1; k < nrb * nA; ++k) {
                TLine *l = new TLine(xa + k * ndr_r * dx, ya, xa + k * ndr_r * dx, yb);
                l->SetLineColor(kGray + 3); l->SetLineWidth(k % nA ? 1 : 3); l->Draw();
            }
            for (int k = 1; k < nA; ++k) {
                TLine *l = new TLine(xa, ya + k * ndr_g * dy, xb, ya + k * ndr_g * dy);
                l->SetLineColor(kGray + 3); l->SetLineWidth(1); l->Draw();
            }
            hz->GetXaxis()->SetLabelSize(0); hz->GetYaxis()->SetLabelSize(0);   // global bin numbers mean nothing here
            hz->GetXaxis()->SetTickLength(0); hz->GetYaxis()->SetTickLength(0);
            TLatex tz; tz.SetTextFont(42); tz.SetTextSize(0.022); tz.SetTextAlign(21);
            for (int k = 0; k < nrb; ++k)
                tz.DrawLatex(xa + (k + 0.5) * nA * ndr_r * dx, yb + 0.6 * dy, Form("reco %g-%g GeV", u3d_reco_pt[rb0 + k], u3d_reco_pt[rb0 + k + 1]));
            tz.SetTextSize(0.018);
            tz.SetTextSize(nrb > 2 ? 0.013 : 0.018);
            for (int k = 0; k < nrb * nA; ++k)   // A bin of each reco column block
                tz.DrawLatex(xa + (k + 0.5) * ndr_r * dx, ya - 1.6 * dy, Form("A %g-%g", B_binsVector[k % nA], B_binsVector[k % nA + 1]));
            tz.SetTextAlign(32);
            for (int k = 0; k < nA; ++k)       // A bin of each gen row block
                tz.DrawLatex(xa - 1.5 * dx, ya + (k + 0.5) * ndr_g * dy, Form("A %g-%g", B_binsVector[k], B_binsVector[k + 1]));
            hz->GetXaxis()->SetTitle("Detector-level bin (#DeltaR within each A block)");
            hz->GetYaxis()->SetTitle(Form("Particle-level bin, gen %g-%g GeV", u3d_gen_pt[2], u3d_gen_pt[3]));
            hz->GetXaxis()->SetTitleOffset(1.2); hz->GetYaxis()->SetTitleOffset(1.2);
            header(Form("3D response, measured p_{T} (%s)", EEC_WEIGHT_OFF ? "no EEC weight" : "EEC-weighted"));
            cz->RedrawAxis();
            // const TString zs = rfolder + (split ? "response_matrix_3d_pseudo_zoom_pt100to120" : "response_matrix_3d_zoom_pt100to120");
            const TString zs = rfolder + (split ? "response_matrix_3d_pseudo_zoom_nominal_pt" : "response_matrix_3d_zoom_nominal_pt");
            cz->Print(zs + ".pdf"); cz->Print(zs + ".png");
        }
        // (b) projection on (dR, pT): the raw response summed over reco and gen A, then rows normalised
        {
            const TH2 *raw = response->Hresponse();
            TH2D *hp = new TH2D("h_response_3d_proj_dr", "", ndr_r * npt_r, 0., ndr_r * npt_r, ndr_g * npt_g, 0., ndr_g * npt_g);
            hp->SetDirectory(nullptr);
            for (int ix = 1; ix <= raw->GetNbinsX(); ++ix)
                for (int iy = 1; iy <= raw->GetNbinsY(); ++iy) {
                    const int r = ix - 1, g = iy - 1;
                    const int rdr = r % ndr_r, rpt = r / (ndr_r * nA), gdr = g % ndr_g, gpt = g / (ndr_g * nA);
                    hp->AddBinContent(hp->GetBin(rdr + ndr_r * rpt + 1, gdr + ndr_g * gpt + 1), raw->GetBinContent(ix, iy));
                }
            for (int iy = 1; iy <= hp->GetNbinsY(); ++iy) {
                double row = 0.;
                for (int ix = 1; ix <= hp->GetNbinsX(); ++ix) row += hp->GetBinContent(ix, iy);
                if (row > 0.) for (int ix = 1; ix <= hp->GetNbinsX(); ++ix) hp->SetBinContent(ix, iy, hp->GetBinContent(ix, iy) / row);
            }
            hp->SetMinimum(0.); hp->SetMaximum(0.3);
            hp->GetXaxis()->SetTitle(Form("Detector-level bin (%d #DeltaR #times %d p_{T})", ndr_r, npt_r));
            hp->GetYaxis()->SetTitle(Form("Particle-level bin (%d #DeltaR #times %d p_{T})", ndr_g, npt_g));
            hp->GetZaxis()->SetTitle("P(detector | particle)");
            hp->GetXaxis()->SetTitleSize(0.035); hp->GetYaxis()->SetTitleSize(0.035);
            hp->GetYaxis()->SetTitleOffset(1.3); hp->GetZaxis()->SetTitleOffset(1.3);
            TCanvas *cp = new TCanvas("c_resp3d_proj", "", 1300, 950);
            cp->SetRightMargin(0.14); cp->SetLeftMargin(0.10); cp->SetTopMargin(0.09); cp->SetBottomMargin(0.10);
            cp->SetFrameFillColor(kGray + 1);
            hp->Draw("COLZ");
            for (int k = 1; k < npt_r; ++k) { TLine *l = new TLine(k * ndr_r, 0., k * ndr_r, ndr_g * npt_g); l->SetLineColor(kGray + 3); l->SetLineWidth(3); l->Draw(); }
            for (int k = 1; k < npt_g; ++k) { TLine *l = new TLine(0., k * ndr_g, ndr_r * npt_r, k * ndr_g); l->SetLineColor(kGray + 3); l->SetLineWidth(3); l->Draw(); }
            TLatex tp; tp.SetTextFont(42); tp.SetTextSize(0.022); tp.SetTextAlign(21);
            for (int k = 0; k < npt_r; ++k) tp.DrawLatex((k + 0.5) * ndr_r, ndr_g * npt_g + 0.6, Form("reco %g-%g GeV", u3d_reco_pt[k], u3d_reco_pt[k + 1]));
            tp.SetTextAlign(12); tp.SetTextColor(kWhite);
            for (int k = 0; k < npt_g; ++k)
                tp.DrawLatex(0.5, (k + 0.5) * ndr_g, k == npt_g - 1 ? Form("gen > %g GeV", u3d_gen_pt[k]) : Form("gen %g-%g GeV", u3d_gen_pt[k], u3d_gen_pt[k + 1]));
            header(Form("3D response projected on #DeltaR #times p_{T} (%s)", EEC_WEIGHT_OFF ? "no EEC weight" : "EEC-weighted"));
            cp->RedrawAxis();
            const TString ps = rfolder + (split ? "response_matrix_3d_pseudo_proj_dr_pt" : "response_matrix_3d_proj_dr_pt");
            cp->Print(ps + ".pdf"); cp->Print(ps + ".png");
            std::cout << "Wrote " << ps << ".png and the nominal-pT zoom" << std::endl;
        }
    }

    // Bin-by-bin rescaling by num/den (or its inverse), content and error alike.
    auto rescale = [](TH3D *h, const TH3D *num, const TH3D *den, bool divide) {
        for (int i = 0; i < h->GetNcells(); ++i) {
            const double d = den->GetBinContent(i);
            const double r = (d > 0.) ? num->GetBinContent(i) / d : 0.;
            const double f = divide ? ((r > 0.) ? 1. / r : 0.) : r;
            h->SetBinContent(i, h->GetBinContent(i) * f);
            h->SetBinError(i, h->GetBinError(i) * f);
        }
    };
    TH3D *h_pc = (TH3D *) h_in->Clone("h_input_purity_corrected_3d");
    h_pc->SetDirectory(nullptr);
    rescale(h_pc, pur_num, pur_den, false);

    // Refolding chi2 over every reco bin with input, input stat only (AN-2024-133).
    auto refoldP = [&](const TH1 *h_unf, double &chi2, int &ndf) -> double {
        TH1 *r = response->ApplyToTruth(h_unf, "h_refold_3d_tmp");
        chi2 = 0.; ndf = 0;
        for (int ix = 1; ix <= h_pc->GetNbinsX(); ++ix)
            for (int iy = 1; iy <= h_pc->GetNbinsY(); ++iy)
                for (int iz = 1; iz <= h_pc->GetNbinsZ(); ++iz) {
                    const double d = h_pc->GetBinContent(ix, iy, iz), e = h_pc->GetBinError(ix, iy, iz);
                    if (d == 0. || e <= 0.) continue;
                    chi2 += std::pow((r->GetBinContent(ix, iy, iz) - d) / e, 2);
                    ++ndf;
                }
        delete r;
        return ndf > 0 ? TMath::Prob(chi2, ndf) : -1.;
    };
    // Closure in the measured region (gen pT 80-200, dR < 0.45, all A): rms of unfolded/truth - 1.
    const int iz_lo = 2, iz_hi = 3, ix_hi = u3d_gen_dr_size - 2;
    auto closureRms = [&](const TH3D *u) -> double {
        double s = 0.; int n = 0;
        for (int iz = iz_lo; iz <= iz_hi; ++iz)
            for (int iy = 1; iy <= u->GetNbinsY(); ++iy)
                for (int ix = 1; ix <= ix_hi; ++ix) {
                    const double t = h_truth->GetBinContent(ix, iy, iz);
                    if (t <= 0.) continue;
                    s += std::pow(u->GetBinContent(ix, iy, iz) / t - 1., 2);
                    ++n;
                }
        return n ? std::sqrt(s / n) : -1.;
    };

    // ---- niter: fixed, or the scan ----
    const double an_threshold = 0.95;
    int niter = NITER;
    std::vector<double> v_n, v_p, v_rms;
    if (NITER <= 0) {
        const char *env = std::getenv("SCAN_NITER_MAX");
        const int n_max = (env && std::atoi(env) > 0) ? std::atoi(env) : 200;
        std::cout << "\n==== 3D niter scan 1.." << n_max << " (refolding, input stat; closure rms in gen pT 80-200, dR < 0.45) ====\n"
                  << Form("  %5s %12s %5s %10s %12s\n", "niter", "chi2", "ndf", "p", "closure rms");
        for (int n = 1; n <= n_max; ++n) {
            RooUnfoldBayes unf(response, h_pc, n);
            unf.SetVerbose(-1);
            TH3D *u = (TH3D *) unf.Hreco(RooUnfold::kNoError)->Clone("h_unf_3d_scan");
            double chi2 = 0.; int ndf = 0;
            const double p = refoldP(u, chi2, ndf);
            rescale(u, eff_num, eff_den, true);
            const double rms = closureRms(u);
            std::cout << Form("  %5d %12.2f %5d %10.3e %11.2f%%\n", n, chi2, ndf, p, 100. * rms);
            v_n.push_back(n); v_p.push_back(p); v_rms.push_back(100. * rms);
            if (niter <= 0 && p >= an_threshold) niter = n;
            delete u;
        }
        if (niter <= 0) {
            niter = n_max;
            std::cout << "  ----> no niter reaches p >= " << an_threshold << " up to " << n_max
                      << "; using " << n_max << std::endl;
        } else {
            std::cout << Form("  ----> OPTIMAL niter (AN-2024-133: first p >= %.2f, input stat) = %d\n", an_threshold, niter);
        }
    }

    // ---- The unfolding at the chosen niter, with errors ----
    std::cout << "\t---->Unfolding 3D, Bayesian (niter = " << niter << ")" << std::endl;
    RooUnfoldBayes unf(response, h_pc, niter);
    unf.SetVerbose(-1);
    TH3D *h_unf = (TH3D *) unf.Hreco(RooUnfold::kCovariance)->Clone("h_unfolded_before_eff_3d");
    h_unf->SetDirectory(nullptr);
    double chi2 = 0.; int ndf = 0;
    const double p = refoldP(h_unf, chi2, ndf);
    TH3D *h_ref = (TH3D *) response->ApplyToTruth(h_unf, "h_refolded_3d");
    TH3D *h_res = (TH3D *) h_unf->Clone("h_unfolded_3d");
    h_res->SetDirectory(nullptr);
    rescale(h_res, eff_num, eff_den, true);
    std::cout << Form("\t---->Refolding GoF: chi2/ndf = %.2f/%d, p = %.3e; closure rms (gen pT 80-200, dR < 0.45) = %.2f%%",
                      chi2, ndf, p, 100. * closureRms(h_res)) << std::endl;

    // ---- Output ----
    const TString folder = mc_base + "/results/unfolding_3d/" + SAMPLE + "_" + GENERATOR + out_tag + "/";
    gSystem->mkdir(folder, kTRUE);
    const TString stem = folder + (split ? "split_test" : "full_closure") + Form("_niter%d", niter);
    TFile *fout = TFile::Open(stem + ".root", "RECREATE");
    for (TH1 *h : std::initializer_list<TH1 *>{h_in, h_pc, h_unf, h_ref, h_res, h_truth}) h->Write();
    if (!v_n.empty()) {
        TGraph g_p((int) v_n.size(), v_n.data(), v_p.data());     g_p.Write("g_refold_p_vs_niter");
        TGraph g_r((int) v_n.size(), v_n.data(), v_rms.data());   g_r.Write("g_closure_rms_vs_niter");
    }
    fout->Close();

    // Projections in one gen pT bin: dR over A, or A over dR < 0.45.
    auto proj = [&](const TH3D *h, bool dr, int iz, const char *tag) -> TH1D * {
        TH1D *o = dr ? h->ProjectionX(Form("%s_dr_%d", tag, iz), 1, h->GetNbinsY(), iz, iz, "e")
                     : h->ProjectionY(Form("%s_A_%d", tag, iz), 1, ix_hi, iz, iz, "e");
        if (dr) o->GetXaxis()->SetRange(1, ix_hi);
        o->SetDirectory(nullptr);
        return o;
    };
    const Color_t red = ROCColor::red(), blue = ROCColor::blue();
    TCanvas *c = new TCanvas("c_unf3d", "", 1200, 1000);
    c->Divide(2, 2);
    int pad = 0;
    for (int dr = 1; dr >= 0; --dr)
        for (int iz = iz_lo; iz <= iz_hi; ++iz) {
            TPad *pp = (TPad *) c->cd(++pad);
            TPad *top = new TPad(Form("t%d", pad), "", 0., 0.3, 1., 1.);
            TPad *bot = new TPad(Form("b%d", pad), "", 0., 0., 1., 0.3);
            top->SetBottomMargin(0.02); top->SetLeftMargin(0.15);
            bot->SetTopMargin(0.02); bot->SetBottomMargin(0.3); bot->SetLeftMargin(0.15);
            pp->cd(); top->Draw(); bot->Draw();
            TH1D *hu = proj(h_res, dr, iz, "unf"), *ht = proj(h_truth, dr, iz, "truth");
            hu->SetMarkerStyle(kFullCircle); hu->SetMarkerColor(red); hu->SetLineColor(red);
            ht->SetLineColor(blue); ht->SetLineWidth(2);
            top->cd();
            ht->SetMaximum(1.3 * std::max(ht->GetMaximum(), hu->GetMaximum()));
            ht->SetMinimum(0.);
            ht->GetYaxis()->SetTitle(EEC_WEIGHT_OFF ? "weighted jets" : "EEC-weighted jets");
            ht->GetXaxis()->SetLabelSize(0);
            ht->Draw("hist");
            hu->Draw("e1 same");
            TLegend *lg = new TLegend(0.50, 0.62, 0.88, 0.84);
            lg->SetBorderSize(0); lg->SetFillStyle(0);
            lg->AddEntry(hu, Form("unfolded, 3D Bayes N = %d", niter), "pe");
            lg->AddEntry(ht, split ? "truth (even half)" : "truth", "l");
            lg->Draw();
            TLatex tx; tx.SetNDC(); tx.SetTextSize(0.045);
            tx.SetTextFont(62); tx.DrawLatex(0.15, 0.92, "CMS");
            tx.SetTextFont(52); tx.DrawLatex(0.25, 0.92, "Work In Progress");
            tx.SetTextFont(42);
            tx.DrawLatex(0.50, 0.86, Form("%s, gen p_{T} %g-%g GeV", split ? "split test" : "full-MC closure",
                                          u3d_gen_pt[iz - 1], u3d_gen_pt[iz]));
            bot->cd();
            TH1D *r = (TH1D *) hu->Clone(Form("ratio_%d", pad));
            r->Divide(ht);
            r->SetMinimum(0.8); r->SetMaximum(1.2);
            r->GetYaxis()->SetTitle("unf. / truth"); r->GetYaxis()->SetNdivisions(505);
            r->GetYaxis()->SetTitleSize(0.1); r->GetYaxis()->SetLabelSize(0.09); r->GetYaxis()->SetTitleOffset(0.6);
            r->GetXaxis()->SetTitle(dr ? "#DeltaR" : "A = p_{T}^{b1}/(p_{T}^{b1} + p_{T}^{b2})");
            r->GetXaxis()->SetTitleSize(0.11); r->GetXaxis()->SetLabelSize(0.1);
            r->Draw("e1");
            TLine *l1 = new TLine(r->GetXaxis()->GetBinLowEdge(r->GetXaxis()->GetFirst()), 1.,
                                  r->GetXaxis()->GetBinUpEdge(r->GetXaxis()->GetLast()), 1.);
            l1->SetLineStyle(2); l1->Draw();
        }
    c->Print(stem + ".pdf");
    c->Print(stem + ".png");
    std::cout << "Wrote " << stem << ".{root,pdf,png}" << std::endl;
}

// ============================================================================
// The tracking-efficiency variation on the PLAIN dr distribution
// ============================================================================
//
// The band quotes the 3% track drop on the EEC, where every jet enters weighted by
// (pt1*pt2)^n. This draws the same variation on dr with NO EEC weight -- the plain number of
// selected 2b jets per dr bin -- so the weight's contribution to the systematic can be seen
// rather than assumed.
//
// Both histograms are filled side by side in create_files_for_template_fit.cpp and are
// already in every production, nominal and _trkdrop030 alike:
//     h3D_bb     (mB, dr, jtpt) filled with eec * weight_tree   <- what the analysis uses
//     h_count_bb (mB, dr, jtpt) filled with weight_tree only    <- the unweighted counts
// so this needs no new production and no new unfolding.
//
// RECO LEVEL, and it cannot be anything else without new MC. The response, the purity, the
// efficiency and the truth are ALL built with w_reco = weight_tree * eec_reco
// (create_files_for_template_fit.cpp:1169), so there is no unweighted response to unfold an
// unweighted spectrum through. Unfolding the counts with the EEC-weighted response would be
// meaningless. What is shown here is therefore the variation as it enters the chain, before
// unfolding -- the right comparison for "what does the EEC weight do to this systematic",
// and NOT a second estimate of the systematic on the measurement.
//
// It lives in this file because it reuses aggChunkFiles()/sumOverFiles() above, but it is a
// separate entry point, so load the file and call it rather than running the macro:
//
//   root -l -b -q -e '.L apply_unfolding_2d.C' -e 'plot_tracking_eff_dr("both","pythia")'
// ============================================================================
void plot_tracking_eff_dr(TString SAMPLE = "both", TString GENERATOR = "pythia",
                          int ibin_pt = 2)
{
    if (!isKnownSample(SAMPLE) || !isKnownGenerator(GENERATOR)) {
        std::cerr << "ERROR: bad SAMPLE / GENERATOR" << std::endl;
        return;
    }

    const TString mc_base  = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024";
    const TString btag_tag = "btagWP0712";

    // The two productions: identical events, one with 3% of the reco tracks thrown away.
    // "_upartv2" was the pre-_rho era, which is no longer produced (fixed 2026-10-26): the
    // current one is obsProdTag("dr"), the same production apply_unfolding() reads.
    // const TString tag_nom = mcVarTag(false) + "_upartv2";
    // const TString tag_trk = mcVarTag(true)  + "_upartv2";
    const TString tag_nom = mcVarTag(false) + obsProdTag("dr");
    const TString tag_trk = mcVarTag(true)  + obsProdTag("dr");

    std::cout << "Nominal production:" << std::endl;
    std::vector<TString> files_nom_n =
        aggChunkFiles(mc_base, SAMPLE, GENERATOR, "", "MCGEN", btag_tag, tag_nom);
    std::cout << "3% track-drop production:" << std::endl;
    std::vector<TString> files_trk_n =
        aggChunkFiles(mc_base, SAMPLE, GENERATOR, "", "MCGEN", btag_tag, tag_trk);
    if (files_nom_n.empty() || files_trk_n.empty()) {
        std::cerr << "ERROR: need BOTH productions -- run run_agg_ntuple_chunks.sh with "
                  << "TRACK_EFF_UNC=false and again with TRACK_EFF_UNC=true" << std::endl;
        return;
    }

    std::vector<TFile *> f_nom = openAllOrWarn(files_nom_n);
    std::vector<TFile *> f_trk = openAllOrWarn(files_trk_n);
    if (f_nom.empty() || f_trk.empty()) return;

    // h_count_bb = no EEC weight, h3D_bb = the EEC weight. Same events, same selection.
    // Literal names, not obs.n(): this is a dR-only diagnostic, and dR's suffix is "".
    //
    // "both" is the weighted qcd + bjet combination, as in apply_unfolding() (2026-10-26; before,
    // the two samples were summed raw). Each production gets its OWN weights: gMixWeight is keyed
    // by sample directory, not by production, so it is set again before the second pair of sums.
    /* ---- disabled (kept for reference): the raw sums, until 2026-10-26 ----
    TH3D *c_nom = sumOverFiles<TH3D>(f_nom, "h_count_bb");
    TH3D *c_trk = sumOverFiles<TH3D>(f_trk, "h_count_bb");
    TH3D *e_nom = sumOverFiles<TH3D>(f_nom, "h3D_bb");
    TH3D *e_trk = sumOverFiles<TH3D>(f_trk, "h3D_bb");
    ---- end disabled ---- */
    gMixWeight.clear();
    if (!setMixWeights(mc_base, SAMPLE, GENERATOR, btag_tag, tag_nom, obsDr())) return;
    TH3D *c_nom = sumOverFiles<TH3D>(f_nom, "h_count_bb");
    TH3D *e_nom = sumOverFiles<TH3D>(f_nom, "h3D_bb");
    gMixWeight.clear();
    if (!setMixWeights(mc_base, SAMPLE, GENERATOR, btag_tag, tag_trk, obsDr())) return;
    TH3D *c_trk = sumOverFiles<TH3D>(f_trk, "h_count_bb");
    TH3D *e_trk = sumOverFiles<TH3D>(f_trk, "h3D_bb");
    gMixWeight.clear();
    if (!c_nom || !c_trk || !e_nom || !e_trk) return;

    // (mB, dr, jtpt) -> dr, integrating mB over everything and taking one jet-pT bin, which
    // is the same slice the unfolded result is quoted in.
    auto drOf = [&](TH3D *h, const char *name) {
        TH1D *p = (TH1D *) h->ProjectionY(name, 1, h->GetNbinsX(), ibin_pt, ibin_pt);
        p->SetDirectory(nullptr);
        const double I = p->Integral("width");
        if (I > 0.) p->Scale(1. / I, "width");   // unit area, as the EEC result is normalised
        return p;
    };
    TH1D *p_c_nom = drOf(c_nom, "p_count_nom");
    TH1D *p_c_trk = drOf(c_trk, "p_count_trk");
    TH1D *p_e_nom = drOf(e_nom, "p_eec_nom");
    TH1D *p_e_trk = drOf(e_trk, "p_eec_trk");

    const int nb = p_c_nom->GetNbinsX();

    // Ratios: the variation over the nominal, with and without the weight. The ratio of two
    // unit-area spectra, so a pure normalisation change cancels and what is left is shape --
    // exactly what the systematic on the normalised result is made of.
    TH1D *r_count = (TH1D *) p_c_trk->Clone("r_count");
    TH1D *r_eec   = (TH1D *) p_e_trk->Clone("r_eec");
    for (TH1D *h : {r_count, r_eec}) h->SetDirectory(nullptr);
    r_count->Divide(p_c_nom);
    r_eec  ->Divide(p_e_nom);
    for (TH1D *h : {r_count, r_eec})
        for (int i = 1; i <= nb; ++i) h->SetBinError(i, 0.);

    // ---- draw ----
    // Reco-level MC, so neither curve is "the measurement" and neither is particle level:
    // the house red and blue stay unused here. Nominal is the reference (black), the
    // variation and the two ratios take the secondary green and purple.
    const Color_t col_nom   = kBlack;
    const Color_t col_trk   = ROCColor::green();
    const Color_t col_eec_r = ROCColor::purple();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);

    const Float_t font_scale  = 1200. / 800.;
    const Style_t font_code   = 43;
    const Float_t label_size  = 15. * font_scale;
    const Float_t title_size  = 15. * font_scale;
    const Float_t legend_size = 14. * font_scale;

    TCanvas *c = new TCanvas("c_trkeff_dr", "", 800, 850);
    TPad *t_top = new TPad("t_top", "", 0., 0.40, 1., 1.);
    TPad *t_bot = new TPad("t_bot", "", 0., 0.,   1., 0.40);
    for (TPad *q : {t_top, t_bot}) { q->SetTicks(1, 0); q->SetFillColor(0); }
    t_top->SetMargin(0.13, 0.05, 0.0,  0.08);
    t_bot->SetMargin(0.13, 0.05, 0.18, 0.0);
    c->cd();
    t_top->Draw(); t_bot->Draw();

    // ---- top: the unweighted dr spectrum, nominal vs 3% drop ----
    t_top->cd();
    double ymax = 0.;
    for (int i = 1; i <= nb; ++i)
        ymax = std::max(ymax, std::max(p_c_nom->GetBinContent(i), p_c_trk->GetBinContent(i)));

    p_c_nom->SetTitle("");
    p_c_nom->GetYaxis()->SetRangeUser(0., ymax * 1.45);
    p_c_nom->GetYaxis()->SetTitle("Normalised counts");
    p_c_nom->GetYaxis()->CenterTitle(true);
    p_c_nom->GetYaxis()->SetTitleFont(font_code); p_c_nom->GetYaxis()->SetTitleSize(title_size);
    p_c_nom->GetYaxis()->SetTitleOffset(1.6);
    p_c_nom->GetYaxis()->SetLabelFont(font_code); p_c_nom->GetYaxis()->SetLabelSize(label_size);
    // The pads touch (bottom margin 0), so the lowest label sits exactly on the boundary and
    // is drawn half-clipped. Suppress it; the axis is still read from the ones above.
    p_c_nom->GetYaxis()->ChangeLabel(1, -1, 0.);
    p_c_nom->GetXaxis()->SetTitleSize(0);
    p_c_nom->GetXaxis()->SetLabelSize(0);

    p_c_nom->SetLineColor(col_nom); p_c_nom->SetMarkerColor(col_nom);
    p_c_nom->SetMarkerStyle(kFullCircle); p_c_nom->SetMarkerSize(1.0);
    p_c_nom->SetLineWidth(2);
    p_c_trk->SetLineColor(col_trk); p_c_trk->SetMarkerColor(col_trk);
    p_c_trk->SetLineWidth(3); p_c_trk->SetLineStyle(2);

    p_c_nom->Draw("PE X0");
    p_c_trk->Draw("HIST same");
    p_c_nom->Draw("PE X0 same");
    gPad->RedrawAxis();

    TLegend *leg = new TLegend(0.45, 0.66, 0.93, 0.90);
    leg->SetTextFont(font_code); leg->SetTextSize(legend_size);
    leg->AddEntry(p_c_nom, "Nominal", "pe");
    leg->AddEntry(p_c_trk, "3% of reco tracks dropped", "l");
    leg->Draw();

    // ---- bottom: the variation/nominal ratio, with and without the EEC weight ----
    t_bot->cd();
    double rlo = 1., rhi = 1.;
    for (int i = 1; i <= nb; ++i) {
        rlo = std::min(rlo, std::min(r_count->GetBinContent(i), r_eec->GetBinContent(i)));
        rhi = std::max(rhi, std::max(r_count->GetBinContent(i), r_eec->GetBinContent(i)));
    }
    const double rpad = std::max(0.01, 0.55 * (rhi - rlo));

    r_count->SetTitle("");
    r_count->GetYaxis()->SetRangeUser(rlo - rpad, rhi + rpad);
    r_count->GetYaxis()->SetTitle("Dropped / nominal");
    r_count->GetXaxis()->SetTitle("#Delta r");   // dR-only diagnostic; no ObsDef in scope
    r_count->GetXaxis()->CenterTitle(true);
    r_count->GetYaxis()->CenterTitle(true);
    r_count->GetYaxis()->SetNdivisions(505);
    for (TAxis *ax : { r_count->GetXaxis(), r_count->GetYaxis() }) {
        ax->SetTitleFont(font_code); ax->SetTitleSize(title_size);
        ax->SetLabelFont(font_code); ax->SetLabelSize(label_size);
    }
    r_count->GetXaxis()->SetTitleOffset(1.1);
    r_count->GetYaxis()->SetTitleOffset(1.6);

    r_count->SetLineColor(col_trk);   r_count->SetLineWidth(3); r_count->SetMarkerSize(0);
    r_eec  ->SetLineColor(col_eec_r); r_eec  ->SetLineWidth(3); r_eec  ->SetMarkerSize(0);
    r_eec  ->SetLineStyle(2);

    r_count->Draw("HIST");
    r_eec  ->Draw("HIST same");

    TLine *l_one = new TLine(r_count->GetXaxis()->GetXmin(), 1.,
                             r_count->GetXaxis()->GetXmax(), 1.);
    l_one->SetLineColor(kGray + 1); l_one->SetLineStyle(2);
    l_one->Draw();
    t_bot->RedrawAxis();

    TLegend *leg2 = new TLegend(0.16, 0.80, 0.92, 0.97);
    leg2->SetNColumns(2);
    leg2->SetTextFont(font_code); leg2->SetTextSize(legend_size);
    leg2->AddEntry(r_count, "no EEC weight", "l");
    leg2->AddEntry(r_eec,   "with EEC weight", "l");
    leg2->Draw();

    // Nothing here is unfolded -- this is a RECO-level diagnostic -- so it does not belong in
    // the matrix_inversion/ production. It keeps the flat results layout, which is what
    // methodDir(true) gives.
    const TString stem = resultFolder(SAMPLE, GENERATOR, /*unfoldBayes=*/true)
                       + "tracking_eff_dr_"
                       + SAMPLE + "_" + GENERATOR + Form("_pt%d", ibin_pt);
    c->Print(stem + ".pdf");
    c->Print(stem + ".png");
    std::cout << "\nWrote " << stem << ".{pdf,png}" << std::endl;

    printf("\n  3%% track drop on dr, reco level, unit-area normalised (pT bin %d)\n", ibin_pt);
    printf("  %3s %13s | %12s %12s | %12s\n",
           "bin", "dr range", "no weight", "with weight", "ratio of the two");
    printf("  -----------------+---------------------------+-----------------\n");
    for (int i = 1; i <= nb; ++i) {
        const double rc = r_count->GetBinContent(i), re = r_eec->GetBinContent(i);
        printf("  %3d [%5.2f,%5.2f] | %+11.2f%% %+11.2f%% | %12.2f\n", i,
               r_count->GetXaxis()->GetBinLowEdge(i), r_count->GetXaxis()->GetBinUpEdge(i),
               100. * (rc - 1.), 100. * (re - 1.),
               (std::fabs(re - 1.) > 1e-9) ? (rc - 1.) / (re - 1.) : 0.);
    }
    printf("\n  Reco level: there is no unweighted response to unfold through, so this is the\n"
           "  variation as it ENTERS the chain, not a second estimate of the systematic.\n");

    for (TFile *f : f_nom) f->Close();
    for (TFile *f : f_trk) f->Close();
}
