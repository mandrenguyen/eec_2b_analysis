#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TROOT.h"
#include "TColor.h"
#include "TLatex.h"
#include "TLegend.h"
#include "RooUnfold.h"
#include "TLine.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TSystem.h"
#include "binning_histos_small.h"
#include <vector>
#include <cmath>

// Plots go to /data_CMS so this folder keeps only code.
// const char* PLOT_OUTDIR = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results";
// const char* PLOT_OUTDIR_DEFAULT_UNUSED = "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/results/unfolding_upartv2/Response_qcd_upartv2"; // qcd Rmatrix
// Default kept as it was. NOTE it points into a COLLEAGUE'S area, which is not writable
// from this account -- the macro's third argument overrides it, and must be used for any
// run that needs to actually save its plots.
const char* PLOT_OUTDIR = "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/results/unfolding_bjet_upartv2/Response_bjet_upartv2"; // bjet Rmatrix
// const char* PLOT_OUTDIR = "/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/results/unfolding_Rqcdbjet_upartv2/Response_bjetqndqcd_upartv2"; // qcd+bjet merged matrix 


namespace ROCColor {
  Color_t make(Int_t idx, Int_t r, Int_t g, Int_t b) {
      if (!gROOT->GetColor(idx)) new TColor(idx, r / 255., g / 255., b / 255.);
      return (Color_t) idx;
  }
  Color_t blue()   { return make(2401, 0x4C, 0x72, 0xB0); }
  Color_t red()    { return make(2402, 0xC4, 0x4E, 0x52); }
  Color_t green()  { return make(2403, 0x4F, 0x8F, 0x52); }
  Color_t purple() { return make(2404, 0x8C, 0x6B, 0xB1); }
  Color_t orange() { return make(2405, 0xC9, 0x74, 0x30); }
  Color_t teal()   { return make(2406, 0x1F, 0x8A, 0x8A); }
}
Int_t setBluePurpleRedPalette(Int_t ncontours = 255)
{
    const Int_t nstops = 3;
    Double_t stops[nstops] = { 0.00,        0.50,        1.00        };
    Double_t rr[nstops]    = { 0x4C / 255., 0x8C / 255., 0xC4 / 255. };
    Double_t gg[nstops]    = { 0x72 / 255., 0x6B / 255., 0x4E / 255. };
    Double_t bb[nstops]    = { 0xB0 / 255., 0xB1 / 255., 0x52 / 255. };
    const Int_t first = TColor::CreateGradientColorTable(nstops, stops, rr, gg, bb, ncontours);
    gStyle->SetNumberContours(ncontours);
    return first;
}

void drawCorrectionMap(TH2D *h, const char *ztitle, const char *level,
                       double zmin, double zmax)
{
    if (!h) return;

    TH2D *hdraw = (TH2D*) h->Clone(Form("%s_clamped", h->GetName()));
    hdraw->SetDirectory(0);
    for (int ix = 1; ix <= hdraw->GetNbinsX(); ++ix)
        for (int iy = 1; iy <= hdraw->GetNbinsY(); ++iy)
            if (hdraw->GetBinContent(ix, iy) < zmin) hdraw->SetBinContent(ix, iy, zmin);

    hdraw->SetStats(0);
    hdraw->GetXaxis()->SetTitle(Form("%s #Deltar", level));
    hdraw->GetYaxis()->SetTitle(Form("%s jet p_{T} (GeV)", level));
    hdraw->GetZaxis()->SetTitle(ztitle);
    hdraw->GetZaxis()->SetTitleOffset(1.2);
    hdraw->SetMinimum(zmin);
    hdraw->SetMaximum(zmax);

    hdraw->GetXaxis()->SetNdivisions(hdraw->GetNbinsX(), false);   
    hdraw->GetYaxis()->SetNdivisions(hdraw->GetNbinsY(), false);  

    gPad->SetGridx(1);
    gPad->SetGridy(1);

    hdraw->Draw("COLZ");

    gPad->Update();
    hdraw->GetYaxis()->ChangeLabel(-1, -1, -1, -1, -1, -1, "#infty");

    TLatex t;
    t.SetTextAlign(22);          // centred on the cell
    t.SetTextFont(42);
    t.SetTextSize(0.032);
    t.SetTextColor(kWhite);
    t.SetTextAngle(45);
    for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
            const double v = h->GetBinContent(ix, iy);
            if (v == 0.) continue;
            t.DrawLatex(h->GetXaxis()->GetBinCenter(ix),
                        h->GetYaxis()->GetBinCenter(iy),
                        Form("%.2f", v));
        }
    }
    t.SetTextAngle(0);       
  }

void drawCMSHeader(TPad *pad, const char *sublabel = "Work In Progress",
                   const char *rightlabel = "2024 pp (5.36 TeV)")
{
    if (!pad) return;
    pad->cd();
    const double left  = pad->GetLeftMargin();
    const double right = 1.0 - pad->GetRightMargin();
    // +0.02 above the frame; enough to clear a y-axis exponent if one appears.
    const double y     = 1.0 - pad->GetTopMargin() + 0.02;

    // TLatex sizes are a fraction of the PAD HEIGHT, so the same numbers that look
    // right on a 1100px canvas are oversized on the 600-700px ones here.
    const double s_cms = 0.026;          // bold "CMS"
    const double s_sub = 0.020;          // italic sublabel and the beam/energy text

    TLatex tex;
    tex.SetNDC();
    tex.SetTextAlign(11);
     tex.SetTextFont(62); tex.SetTextSize(s_cms); tex.DrawLatex(left, y, "CMS");
    tex.SetTextFont(52); tex.SetTextSize(s_sub); tex.DrawLatex(left + 2.4 * s_cms, y, sublabel);
    tex.SetTextAlign(31);
    tex.SetTextFont(42); tex.SetTextSize(s_sub); tex.DrawLatex(right, y, rightlabel);
}


// Condition number kappa = sigma_max / sigma_min of the response matrix. Completely empty
// bins contribute an exact zero singular value, which would make kappa infinite without
// saying anything about the constrained part of the problem, so they are dropped and the
// SVD runs on the surviving sub-matrix.
double responseConditionNumber(const TH2 *h, const char *tag)
{
    if (!h) return -1.;

    const int nx = h->GetNbinsX();   // detector level
    const int ny = h->GetNbinsY();   // particle level

    std::vector<int> keep_x, keep_y;
    for (int ix = 1; ix <= nx; ++ix) {
        double s = 0.;
        for (int iy = 1; iy <= ny; ++iy) s += std::abs(h->GetBinContent(ix, iy));
        if (s > 0.) keep_x.push_back(ix);
    }
    for (int iy = 1; iy <= ny; ++iy) {
        double s = 0.;
        for (int ix = 1; ix <= nx; ++ix) s += std::abs(h->GetBinContent(ix, iy));
        if (s > 0.) keep_y.push_back(iy);
    }

    const int mx = (int) keep_x.size();
    const int my = (int) keep_y.size();
    std::cout << "\n[" << tag << "] response matrix " << nx << " x " << ny
              << " -> " << mx << " x " << my << " after dropping "
              << (nx - mx) << " empty detector-level and "
              << (ny - my) << " empty particle-level bins" << std::endl;

    if (mx < 1 || my < 1) { std::cout << "  nothing left to decompose" << std::endl; return -1.; }

    // TDecompSVD requires nrows >= ncols; the singular values are the same for the
    // transpose, so orient the matrix that way when the reduced shape is wide.
    const bool transpose = (my < mx);
    const int nrow = transpose ? mx : my;
    const int ncol = transpose ? my : mx;

    TMatrixD m(nrow, ncol);
    for (int i = 0; i < my; ++i)
        for (int j = 0; j < mx; ++j) {
            const double v = h->GetBinContent(keep_x[j], keep_y[i]);
            if (transpose) m(j, i) = v; else m(i, j) = v;
        }

    TDecompSVD svd(m);
    if (!svd.Decompose()) { std::cout << "  SVD failed" << std::endl; return -1.; }

    const TVectorD sig = svd.GetSig();
    const double smax = sig[0];
    const double smin = sig[sig.GetNrows() - 1];

    std::cout << "  sigma_max = " << smax << ", sigma_min = " << smin << std::endl;
    if (smin <= 0.) {
        std::cout << "  kappa = inf (matrix still rank deficient: some non-empty rows/columns"
                     " are linearly dependent)" << std::endl;
        return -1.;
    }
    const double kappa = smax / smin;
    std::cout << "  kappa = " << kappa << "  (log10 = " << std::log10(kappa) << ")" << std::endl;
    return kappa;
}

// infile      : the RMatrix file to read. Default keeps the bjet merged file this macro has
//               always used, so an argument-less call is unchanged.
// observable  : "dr" (default) | "B" -- picks response_tf_full vs response_tf_full_B.
//
// For the momentum balance, the responses live in the per-block "_upartv2_B" production and
// have to be merged first (hadd needs the RooUnfold env, see CLAUDE.md):
//   hadd -f RMatrix_MERGED_qcdbjet_noeecw_upartv2_B.root <QCD and bJet block files>
//   root -l -b -q 'plot_purity_efficiency_response.cpp("<that file>","B")'
//
// ⚠️ The balance was called "z" and its histograms "..._z" before 2026-09-22. A pre-rename
// RMatrix file has no response_tf_full_B in it.
void plot_purity_efficiency_response(const char *infile = nullptr,
                                     const char *observable = "dr",
                                     const char *outdir = nullptr)
{
  // PLOT_OUTDIR is a global the SaveAs calls below read; override it when asked.
  if (outdir) { PLOT_OUTDIR = outdir; gSystem->mkdir(outdir, kTRUE); }
  // CMS aesthetics: no stats box, no title box, ticks on all four sides.
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  setBluePurpleRedPalette();
  gStyle->SetPaintTextFormat("4.3f");


  // btagWP<NNN> follows BTAG_WP in the run scripts.
  // TFile *f = TFile::Open("/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/"    
                     // "RMatrix_Run3_btagWP0712_template_for_fit_histos_3D_qcd_f_upartv2.root","READ");
  

  /// qcd only
  // TFile *f = TFile::Open("/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/MergedResult_btagWP712_MattProd/"
                         // "RMatrix_Run3_btagWP712_template_for_fit_histos_3D_qcd_f_80_9999_2_merged.root","READ"); // QCD matrix 
  /// bjet only 
  TFile *f = infile
      ? TFile::Open(infile, "READ")
      : TFile::Open("/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/bJet/agg_ntuple_chunks/MergedResult_btagWP712_MattProd/"
                    "RMatrix_Run3_btagWP712_template_for_fit_histos_3D_bjet_f_80_9999_2_merged.root","READ"); // bjet matrix
  if (!f || f->IsZombie()) { std::cout << "cannot open input file" << std::endl; return; }
  // "" for dr, so the dr path reads exactly the object it always read.
  const TString obs_sfx = (TString(observable) == "dr") ? "" : ("_" + TString(observable));

  /// Merged Rmatrix: qcd + bjet
  // TFile *f = TFile::Open("/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/bJet/agg_ntuple_chunks/MergedResult_btagWP712_MattProd/"
  //                           "RMatrix_Run3_btagWP712_template_for_fit_histos_3D_bjet_qcd_merged.root", "READ"); // merged Rmatrix (qcd + bjet)


  if(!f || f->IsZombie())
    {std::cout << "Cannot open file" << std::endl;
      return;}

  
    auto ratioFromCounts = [&](const char* num, const char* den, const char* name) -> TH2D* {
        TH2D* hn = (TH2D*) f->Get(num);
        TH2D* hd = (TH2D*) f->Get(den);
        if (!hn || !hd) { std::cout << "Missing " << num << " or " << den << std::endl; return nullptr; }
        TH2D* r = (TH2D*) hn->Clone(name);
        r->SetDirectory(0);
        r->Divide(hn, hd, 1., 1., "b");   // numerator is a subset of the denominator
        return r;
    };

    TH2D *h_purity     = ratioFromCounts("h_full_purity_numerator_tf",
                                         "h_full_purity_denominator_tf",     "h_purity");
    TH2D *h_efficiency = ratioFromCounts("h_full_efficiency_numerator_tf",
                                         "h_full_efficiency_denominator_tf", "h_efficiency");

    if(!h_purity || !h_efficiency)
      {std::cout << "Histograms not found" << std::endl;
	    return;}

    std::cout << "purity     range: " << h_purity->GetMinimum()     << " .. " << h_purity->GetMaximum()     << std::endl;
    std::cout << "efficiency range: " << h_efficiency->GetMinimum() << " .. " << h_efficiency->GetMaximum() << std::endl;

  
    const int n_pt = h_purity->GetNbinsY();

    // ---- 2D maps ---------------------------------------------------------------
    TCanvas *cPurity_jtpt_dr     = new TCanvas("cPurity_jtpt_dr",     "Purity vs jtpt and dR",     800, 700);
    TCanvas *cEfficiency_jtpt_dr = new TCanvas("cEfficiency_jtpt_dr", "Efficiency vs jtpt and dR", 800, 700);

    const double z_lo = 0.75, z_hi = 1.0;

    cPurity_jtpt_dr->cd();
    gPad->SetRightMargin(0.16);
    gPad->SetLeftMargin(0.13);
    gPad->SetTopMargin(0.08);
    drawCorrectionMap(h_purity, "Purity", "reco", z_lo, z_hi);
    drawCMSHeader((TPad*) gPad);

    cEfficiency_jtpt_dr->cd();
    gPad->SetRightMargin(0.16);
    gPad->SetLeftMargin(0.13);
    gPad->SetTopMargin(0.08);
    drawCorrectionMap(h_efficiency, "Efficiency", "gen", z_lo, z_hi);
    drawCMSHeader((TPad*) gPad);

   
    TCanvas *cCorr1D = new TCanvas("cCorr1D", "Efficiency vs dr, nominal pT bin", 800, 700);
    cCorr1D->cd();
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);

    /* ---- disabled (kept for reference): purity was overlaid here too ----
    const Color_t col_purity = ROCColor::purple();
    */ // ---- end disabled block ----
    const Color_t col_eff    = ROCColor::orange();

    TLegend *leg = new TLegend(0.16, 0.13, 0.72, 0.28);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.028);

    // Only the nominal jet-pT bin: the 100-120 GeV slice, ibin_pt = 2 (see
    // TemplateFit_Run3/Draw_EEC.h). Clamped in case the pT binning ever shrinks.
    const int ipt_nom = (n_pt >= 2) ? 2 : n_pt;
    {
        const int ipt = ipt_nom;
        const double ptlo = h_efficiency->GetYaxis()->GetBinLowEdge(ipt);
        const double pthi = h_efficiency->GetYaxis()->GetBinUpEdge(ipt);

        TH1D *e1 = h_efficiency->ProjectionX(Form("p_eff_pt%d", ipt), ipt, ipt);

        e1->SetDirectory(0);
        e1->SetStats(0);
        e1->SetLineWidth(2);
        e1->SetMarkerSize(1.2);
        e1->SetLineStyle(kSolid);
        e1->SetMarkerStyle(kFullCircle);
        e1->SetLineColor(col_eff);
        e1->SetMarkerColor(col_eff);

        e1->SetTitle("");
        e1->GetYaxis()->SetRangeUser(0.9, 1.0);
        e1->GetYaxis()->SetTitle("efficiency");
        e1->GetYaxis()->CenterTitle(true);
        e1->GetYaxis()->SetTitleOffset(1.3);
        e1->GetXaxis()->SetTitle("#Deltar");
        e1->GetXaxis()->CenterTitle(true);
        e1->Draw("PE");

        const TString pt_txt = (ipt == n_pt)
            ? Form("p_{T}^{jet} > %.0f GeV", ptlo)
            : Form("%.0f < p_{T}^{jet} < %.0f GeV", ptlo, pthi);
        leg->AddEntry(e1, Form("Efficiency, %s", pt_txt.Data()), "lp");
    }
    leg->Draw();

    TLine *l_one = new TLine(h_efficiency->GetXaxis()->GetXmin(), 1.0,
                             h_efficiency->GetXaxis()->GetXmax(), 1.0);
    l_one->SetLineColor(kGray + 2);
    l_one->SetLineStyle(2);
    l_one->Draw("same");
    drawCMSHeader((TPad*) gPad);
    gPad->RedrawAxis();

    // ---- Response matrix -------------------------------------------------------
    RooUnfoldResponse *resp = (RooUnfoldResponse*)f->Get("response_tf_full" + obs_sfx);

    if (!resp) {
      std::cout << "Response not found!" << std::endl;
      return;}

    TH2 *hResponse = (TH2*) resp->Hresponse()->Clone("hResponse_norm");
    hResponse->SetDirectory(0);

    // keep the un-normalised counts around: the conditioning is reported for both
    TH2 *hResponseRaw = (TH2*) resp->Hresponse()->Clone("hResponse_raw");
    hResponseRaw->SetDirectory(0);


    for (int iy = 1; iy <= hResponse->GetNbinsY(); ++iy) {
        double row = 0.;
        for (int ix = 1; ix <= hResponse->GetNbinsX(); ++ix) row += hResponse->GetBinContent(ix, iy);
        if (row <= 0.) continue;
        for (int ix = 1; ix <= hResponse->GetNbinsX(); ++ix) {
            hResponse->SetBinContent(ix, iy, hResponse->GetBinContent(ix, iy) / row);
            hResponse->SetBinError  (ix, iy, hResponse->GetBinError(ix, iy)   / row);
        }
    }
/* ---- disabled (kept for reference): the previous global normalisation ----
    hResponse->Scale(1.0 / hResponse->Integral());
*/ // ---- end disabled block ----

    responseConditionNumber(hResponseRaw, "raw counts");
    responseConditionNumber(hResponse,    "row-normalised P(detector | particle)");

    hResponse->SetTitle("");
    hResponse->GetXaxis()->SetTitle("Detector-level bin");
    hResponse->GetYaxis()->SetTitle("Particle-level bin");
    hResponse->GetZaxis()->SetTitle("P(detector | particle)");
    hResponse->GetZaxis()->SetTitleOffset(1.2);

    TCanvas *c = new TCanvas("c","response",800,700);

    
    c->SetFrameFillColor(kGray + 1);
    gPad->SetFrameFillColor(kGray + 1);

    gPad->SetRightMargin(0.16);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.08);
    gPad->SetFixedAspectRatio();
    hResponse->SetStats(0);
    hResponse->Draw("COLZ");

    const int n = hResponse->GetNbinsX();
    const int n_dr_resp = (n_pt > 0) ? n / n_pt : n;

    TLine *diag = new TLine(hResponse->GetXaxis()->GetXmin(), hResponse->GetYaxis()->GetXmin(),
                            hResponse->GetXaxis()->GetXmax(), hResponse->GetYaxis()->GetXmax());
    diag->SetLineColor(kGray + 2);
    diag->SetLineStyle(2);
    diag->SetLineWidth(1);
    diag->Draw("same");

    for (int ib = n_dr_resp; ib < n; ib += n_dr_resp) {
        const double xb = hResponse->GetXaxis()->GetBinUpEdge(ib);
        const double yb = hResponse->GetYaxis()->GetBinUpEdge(ib);
        TLine *lv = new TLine(xb, hResponse->GetYaxis()->GetXmin(), xb, hResponse->GetYaxis()->GetXmax());
        TLine *lh = new TLine(hResponse->GetXaxis()->GetXmin(), yb, hResponse->GetXaxis()->GetXmax(), yb);
        for (TLine* l : {lv, lh}) { l->SetLineColor(kGray + 3); l->SetLineWidth(1); l->Draw("same"); }
    }
    drawCMSHeader((TPad*) gPad);
    gPad->RedrawAxis();

    gSystem->mkdir(PLOT_OUTDIR, kTRUE);

    c->SaveAs(Form("%s/response_matrix.pdf", PLOT_OUTDIR));
    c->SaveAs(Form("%s/response_matrix.png", PLOT_OUTDIR));

    cPurity_jtpt_dr->SaveAs(Form("%s/purity_jtpt_dr.pdf", PLOT_OUTDIR));
    cPurity_jtpt_dr->SaveAs(Form("%s/purity_jtpt_dr.png", PLOT_OUTDIR));
    cEfficiency_jtpt_dr->SaveAs(Form("%s/efficiency_jtpt_dr.pdf", PLOT_OUTDIR));
    cEfficiency_jtpt_dr->SaveAs(Form("%s/efficiency_jtpt_dr.png", PLOT_OUTDIR));
    cCorr1D->SaveAs(Form("%s/efficiency_vs_dr_nominal_pt.pdf", PLOT_OUTDIR));
    cCorr1D->SaveAs(Form("%s/efficiency_vs_dr_nominal_pt.png", PLOT_OUTDIR));
}


// ============================================================================
// momentum_balance_mc_study() -- the purity/stability study that chose the B binning.
// Moved in from momentum_balance_mc_study.C on 2026-09-25 to keep the number of files
// down; unchanged. It belongs next to the purity/efficiency plots: same question (is the
// binning coarse enough for the resolution), asked per bin before the response exists.
//
// MC-only study of a candidate NEW measured observable: the momentum balance between the
// two reconstructed Bs,
//
//     B = pT_lead / (pT1 + pT2)          in [0.5, 1)
//
// where pT1, pT2 are the pTs of the two aggregated B hadrons -- the SAME two objects whose
// pT product is the EEC weight. B = 0.5 is a balanced pair, B -> 1 a very asymmetric one.
//
// WHY THIS MACRO EXISTS
// ---------------------
// Adding B as a parallel measured axis (templates TH3D(m_2B, B, pT), fit sliced in (B, pT),
// response 2D (B, pT)) is a large change to the chain and costs a full MC + data
// reprocessing. The binning has to be decided BEFORE that is spent, and it can be decided
// now: AggBHadronNtuple already stores recoPt1/recoPt2 and genPt1/genPt2, so both the reco
// and the gen momentum balance are already on disk for MC. Nothing is reprocessed here.
//
// It answers three questions:
//   1. what does B look like at particle level, and how much does the EEC weight reshape it;
//   2. how well is it measured -- the (B_reco - B_gen) resolution, and how it varies with B;
//   3. what binning to use -- purity and stability for several candidate binnings, including
//      equal-occupancy edges derived from the gen distribution itself.
//
// B is a RATIO of the two pTs, so the ordering convention does not matter (max/sum is
// symmetric) and the per-B pT scale largely cancels. That is the reason to expect it to be
// better measured than either pT alone -- see h_ptratio in the RMatrix files, where the
// reco/gen SV pT ratio is broad.
//
// SELECTIONS
// ----------
// Deliberately the SAME gates as the response-matrix branch of create_files_for_template_fit.cpp,
// so purity/stability here mean what they will mean in the real chain:
//
//   reco_pass = nRecoAgg == 2 && passRecoKin && passBtag && recoDr > 0.005
//   gen_pass  = nGenAgg >= 2  && passGenKin
//
// and the weight is w_reco = weight * recoEec, the response-matrix convention (gen-binned
// quantities are weighted with the RECO-side EEC weight there too -- see the comment on
// h_half0_eff_den). The unweighted versions are filled alongside, because a yield
// measurement in B would use those and the EEC weight reshapes B strongly.
//
// The ntuple is written for true 2b jets only (jtNbHad >= 2), which is exactly the
// population the response matrix is built from. Nothing here is a background study.
//
// RUN (on LLR, after source setup_roounfold_env.sh -- this file includes RooUnfold.h):
//
//   root -l -b -q -e '.L plot_purity_efficiency_response.cpp' \
//                 -e 'momentum_balance_mc_study("both","pythia")'
//
// Outputs go to
//   /data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/momentum_balance_study_<sample>_<generator>_upartv2/
// which is NOT on the sshfs mount -- use sync_plots_to_eos.sh or scp to look at the PDFs.
// ============================================================================

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

// MomBalance::value lives in observables.h, shared with the chain -- this macro must never
// carry its own copy of the definition, or the binning it recommends stops describing what
// create_files_for_template_fit.cpp actually fills. observables.h includes the binning
// header (guarded), so this is the only include needed for both.
#include "observables.h"

namespace {

// The study's own range, deliberately independent of B_binsVector: this macro exists to
// CHOOSE that binning, so it must be able to look outside it. [0.5, 1) is the observable's
// natural range, not a binning choice.
inline double kMin() { return 0.5; }
inline double kMax() { return 1.0; }

// Fold the open upper edge into the last bin rather than losing it to the overflow.
inline double Bfill(double B, double Bmax) { return (B >= Bmax) ? Bmax - 1e-6 : B; }

} // anonymous namespace

// ---- Input file discovery ------------------------------------------------------------
// Same conventions as aggChunkFiles() in apply_unfolding_2d.C: take what is on disk rather
// than a fixed block range, because the block count differs per sample/generator
// (10 qcd/pythia, 9 bjet/pythia, 8 qcd/herwig, 9 bjet/herwig).
namespace {

std::vector<TString> samplesIn_mb(const TString &sample)
{
    if (sample == "both") return {"qcd", "bjet"};
    return {sample};
}

TString sampleSubdir_mb(const TString &sample, const TString &generator)
{
    const bool herwig = (generator == "herwig");
    if (sample == "qcd")  return herwig ? "QCDHerwig"  : "QCD";
    if (sample == "bjet") return herwig ? "bJetHerwig" : "bJet";
    std::cerr << "ERROR: unknown sample '" << sample << "' (use qcd | bjet | both)" << std::endl;
    return "";
}

const char *kBase     = "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024";
const char *kBtagTag  = "btagWP0712";
const char *kOutSuffix = "_upartv2";   // the nominal production: no _trkdrop030, no _noeecw

int addChunks(TChain &chain, const TString &sample, const TString &generator,
              int max_blocks = 50)
{
    int total = 0;
    for (const TString &s : samplesIn_mb(sample)) {
        const TString subdir = sampleSubdir_mb(s, generator);
        if (subdir.Length() == 0) return -1;
        int found = 0;
        for (int b = 0; b < max_blocks; ++b) {
            TString f = TString::Format("%s/%s/agg_ntuple_chunks/block_%04d/"
                                        "AggBHadronNtuple_Run3_%s_template_for_fit_histos_3D_%s_f%s.root",
                                        kBase, subdir.Data(), b, kBtagTag, s.Data(), kOutSuffix);
            if (gSystem->AccessPathName(f) == kFALSE) { chain.Add(f); ++found; }
        }
        std::cout << "   " << subdir << ": " << found << " block files" << std::endl;
        if (found == 0)
            std::cerr << "ERROR: no AggBHadronNtuple block files under " << kBase << "/"
                      << subdir << "/agg_ntuple_chunks/ -- run run_agg_ntuple_chunks.sh"
                      << " for this sample" << std::endl;
        total += found;
    }
    return total;
}

// ---- Binning evaluation --------------------------------------------------------------
// purity_i    = N(reco in bin i AND gen in bin i) / N(reco in bin i)
// stability_i = N(reco in bin i AND gen in bin i) / N(gen in bin i)
//
// These are the diagonal fractions of the migration matrix -- the quantities that say
// whether a binning is coarse enough for the detector resolution. They are NOT the same as
// the "purity"/"efficiency" histograms the analysis unfolds with (those measure how many
// jets pass reco AND gen selection, not how many stay in their bin), which is why they are
// computed here rather than read from the RMatrix files.
struct BinQuality {
    std::vector<double> edges;
    std::vector<double> purity, stability, frac;   // per bin; frac = share of gen entries
    double min_purity = 1., min_stability = 1.;
};

BinQuality evaluate(const std::vector<double> &edges,
                    const std::vector<double> &B_reco,
                    const std::vector<double> &B_gen,
                    const std::vector<double> &w)
{
    const int nb = (int)edges.size() - 1;
    BinQuality q;
    q.edges = edges;
    q.purity.assign(nb, 0.);
    q.stability.assign(nb, 0.);
    q.frac.assign(nb, 0.);

    std::vector<double> n_reco(nb, 0.), n_gen(nb, 0.), n_diag(nb, 0.);

    auto bin_of = [&](double B) -> int {
        if (B < edges.front() || B >= edges.back()) return -1;
        for (int i = 0; i < nb; ++i)
            if (B >= edges[i] && B < edges[i + 1]) return i;
        return -1;
    };

    for (size_t k = 0; k < B_reco.size(); ++k) {
        const int ir = bin_of(B_reco[k]);
        const int ig = bin_of(B_gen[k]);
        if (ir >= 0) n_reco[ir] += w[k];
        if (ig >= 0) n_gen[ig]  += w[k];
        if (ir >= 0 && ir == ig) n_diag[ir] += w[k];
    }

    double gen_tot = 0.;
    for (int i = 0; i < nb; ++i) gen_tot += n_gen[i];

    for (int i = 0; i < nb; ++i) {
        q.purity[i]    = (n_reco[i] > 0.) ? n_diag[i] / n_reco[i] : 0.;
        q.stability[i] = (n_gen[i]  > 0.) ? n_diag[i] / n_gen[i]  : 0.;
        q.frac[i]      = (gen_tot   > 0.) ? n_gen[i]  / gen_tot   : 0.;
        q.min_purity    = std::min(q.min_purity,    q.purity[i]);
        q.min_stability = std::min(q.min_stability, q.stability[i]);
    }
    return q;
}

void printQuality(const char *name, const BinQuality &q)
{
    const int nb = (int)q.purity.size();
    std::cout << "\n  " << name << "  (" << nb << " bins)" << std::endl;
    std::cout << "    edges:";
    std::cout << std::fixed << std::setprecision(4);
    for (double e : q.edges) std::cout << " " << e;
    std::cout << std::endl;
    std::cout << "    bin |   range        |  gen frac |  purity | stability" << std::endl;
    std::cout << "    ----+----------------+-----------+---------+----------" << std::endl;
    for (int i = 0; i < nb; ++i) {
        std::cout << "    " << std::setw(3) << i + 1 << " | "
                  << std::setw(6) << q.edges[i] << " - " << std::setw(6) << q.edges[i + 1]
                  << " | " << std::setw(9) << q.frac[i]
                  << " | " << std::setw(7) << q.purity[i]
                  << " | " << std::setw(8) << q.stability[i] << std::endl;
    }
    std::cout << "    worst: purity " << q.min_purity
              << ", stability " << q.min_stability << std::endl;
}

// Equal-occupancy edges from the gen distribution: nb bins each holding ~1/nb of the
// weighted gen entries. First and last edges are pinned to the physical range.
std::vector<double> quantileEdges(TH1D *h_gen, int nb)
{
    std::vector<double> edges(nb + 1);
    edges.front() = kMin();
    edges.back()  = kMax();

    const double total = h_gen->Integral();
    if (total <= 0.) return edges;

    double running = 0.;
    int next = 1;
    for (int b = 1; b <= h_gen->GetNbinsX() && next < nb; ++b) {
        running += h_gen->GetBinContent(b);
        if (running / total >= (double)next / nb) {
            // Round to 3 decimals: an edge printed as 0.5734 is not more meaningful than
            // 0.573, and a tidy edge is what would actually go into binning_histos_small.h.
            edges[next] = std::round(h_gen->GetXaxis()->GetBinUpEdge(b) * 1000.) / 1000.;
            ++next;
        }
    }
    // Guard against ties collapsing two edges onto each other.
    for (int i = 1; i <= nb; ++i)
        if (edges[i] <= edges[i - 1]) edges[i] = edges[i - 1] + 0.001;
    return edges;
}

std::vector<double> uniformEdges(int nb)
{
    std::vector<double> edges(nb + 1);
    for (int i = 0; i <= nb; ++i)
        edges[i] = kMin()
                 + (kMax() - kMin()) * i / (double)nb;
    return edges;
}

} // anonymous namespace


void momentum_balance_mc_study(const TString &sample    = "qcd",
                               const TString &generator = "pythia",
                               Long64_t max_entries     = -1)
{
    if (sample != "qcd" && sample != "bjet" && sample != "both") {
        std::cerr << "ERROR: sample must be qcd | bjet | both (got '" << sample << "')" << std::endl;
        return;
    }
    if (generator != "pythia" && generator != "herwig") {
        std::cerr << "ERROR: generator must be pythia | herwig (got '" << generator << "')" << std::endl;
        return;
    }

    std::cout << "=== Momentum-balance MC study: B = pT_lead/(pT1+pT2) ===" << std::endl;
    std::cout << "sample " << sample << " " << generator << std::endl;

    TChain chain("aggBHadronKinematics");
    if (addChunks(chain, sample, generator) <= 0) return;

    // ---- Branches ----
    Float_t weight = 1.f, recoPt1 = 0.f, recoPt2 = 0.f, genPt1 = 0.f, genPt2 = 0.f;
    Float_t recoDr = -1.f, genDr = -1.f, recoEec = 0.f, genEec = 0.f, jtpt = 0.f, refpt = 0.f;
    Int_t   passRecoKin = 0, passGenKin = 0, passBtag = 0, nRecoAgg = 0, nGenAgg = 0;

    chain.SetBranchStatus("*", 0);
    auto on = [&](const char *n, void *addr) {
        chain.SetBranchStatus(n, 1);
        chain.SetBranchAddress(n, addr);
    };
    on("weight", &weight);
    on("recoPt1", &recoPt1);   on("recoPt2", &recoPt2);
    on("genPt1", &genPt1);     on("genPt2", &genPt2);
    on("recoDr", &recoDr);     on("genDr", &genDr);
    on("recoEec", &recoEec);   on("genEec", &genEec);
    on("jtpt", &jtpt);         on("refpt", &refpt);
    on("passRecoKin", &passRecoKin);
    on("passGenKin", &passGenKin);
    on("passBtag", &passBtag);
    on("nRecoAgg", &nRecoAgg);
    on("nGenAgg", &nGenAgg);

    // ---- Histograms ----
    const double Blo = kMin(), Bhi = kMax();
    const int    nfine = 100;

    // Particle level, with and without the EEC weight. The EEC weight is (pt1*pt2)^n, which
    // is LARGEST for a balanced pair, so it pulls B towards 0.5 -- the two curves are not
    // the same shape and the difference is part of what this study is for.
    TH1D *h_Bgen_w   = new TH1D("h_Bgen_w",   "gen B, EEC weighted;B_{gen};EEC weighted entries", nfine, Blo, Bhi);
    TH1D *h_Bgen_raw = new TH1D("h_Bgen_raw", "gen B, unweighted;B_{gen};jets",                   nfine, Blo, Bhi);
    TH1D *h_Breco_w  = new TH1D("h_Breco_w",  "reco B, EEC weighted;B_{reco};EEC weighted entries", nfine, Blo, Bhi);
    TH1D *h_Breco_raw= new TH1D("h_Breco_raw","reco B, unweighted;B_{reco};jets",                 nfine, Blo, Bhi);

    // Resolution. Absolute (B is already dimensionless and bounded, so a relative
    // resolution would only re-introduce a scale).
    TH1D *h_res      = new TH1D("h_res", "B resolution;B_{reco} - B_{gen};EEC weighted entries", 120, -0.3, 0.3);
    TH2D *h_res_vs_B = new TH2D("h_res_vs_B", "B resolution vs B_{gen};B_{gen};B_{reco} - B_{gen}",
                                25, Blo, Bhi, 120, -0.3, 0.3);
    TH2D *h_mig      = new TH2D("h_mig", "migration;B_{gen};B_{reco}", 50, Blo, Bhi, 50, Blo, Bhi);

    // For context: the same resolution for dR, the observable already measured. Puts the B
    // numbers on a scale -- if B is much worse than dR, a B measurement needs coarser bins.
    TH1D *h_Bres_by_dr = new TH1D("h_Bres_by_dr", "B resolution in dR bins;#DeltaR_{gen};RMS(B_{reco}-B_{gen})",
                                  bins_dr, dr_binsVector);

    for (TH1 *h : {(TH1 *)h_Bgen_w, (TH1 *)h_Bgen_raw, (TH1 *)h_Breco_w, (TH1 *)h_Breco_raw,
                   (TH1 *)h_res, (TH1 *)h_res_vs_B, (TH1 *)h_mig})
        h->Sumw2();

    // Kept in memory for the binning scan: one entry per jet passing BOTH gates.
    std::vector<double> v_Breco, v_Bgen, v_w, v_w_raw;
    // The same, for dR -- the observable already measured. Purity and stability are only
    // interpretable against a benchmark: "0.35" is disqualifying if dR scores 0.85 in its
    // own analysis binning and unremarkable if dR scores 0.40. Same events, same gates,
    // same weight, same metric -- the only difference is the variable.
    std::vector<double> v_drreco, v_drgen;

    // ---- Loop ----
    const Long64_t n = (max_entries > 0 && max_entries < chain.GetEntries())
                     ? max_entries : chain.GetEntries();
    std::cout << "Processing " << n << " ntuple rows (true 2b jets)" << std::endl;

    long n_reco_pass = 0, n_gen_pass = 0, n_both = 0, n_bad_B = 0;

    for (Long64_t i = 0; i < n; ++i) {
        if (i % 500000 == 0)
            std::cout << "\r  " << 100. * i / n << " %" << std::flush;
        chain.GetEntry(i);

        // The response-matrix gates, verbatim.
        const bool reco_pass = (nRecoAgg == 2) && passRecoKin && passBtag && (recoDr > 0.005);
        const bool gen_pass  = (nGenAgg >= 2) && passGenKin;

        const double B_reco = reco_pass ? MomBalance::value(recoPt1, recoPt2) : -1.;
        const double B_gen  = gen_pass  ? MomBalance::value(genPt1,  genPt2)  : -1.;

        // A pair that passed the gates but has an unusable pT is a real anomaly, not a
        // selection effect -- count it rather than dropping it silently.
        if ((reco_pass && B_reco < 0.) || (gen_pass && B_gen < 0.)) ++n_bad_B;

        const double w_reco = weight * recoEec;   // response-matrix convention
        const double w_gen  = weight * genEec;

        if (gen_pass && B_gen >= 0.) {
            ++n_gen_pass;
            h_Bgen_w  ->Fill(Bfill(B_gen, Bhi), w_gen);
            h_Bgen_raw->Fill(Bfill(B_gen, Bhi), weight);
        }
        if (reco_pass && B_reco >= 0.) {
            ++n_reco_pass;
            h_Breco_w  ->Fill(Bfill(B_reco, Bhi), w_reco);
            h_Breco_raw->Fill(Bfill(B_reco, Bhi), weight);
        }
        if (reco_pass && gen_pass && B_reco >= 0. && B_gen >= 0.) {
            ++n_both;
            h_res     ->Fill(B_reco - B_gen, w_reco);
            h_res_vs_B->Fill(Bfill(B_gen, Bhi), B_reco - B_gen, w_reco);
            h_mig     ->Fill(Bfill(B_gen, Bhi), Bfill(B_reco, Bhi), w_reco);

            v_Breco.push_back(Bfill(B_reco, Bhi));
            v_Bgen .push_back(Bfill(B_gen,  Bhi));
            v_w    .push_back(w_reco);
            v_w_raw.push_back(weight);

            // dR, with the same overflow folding the chain applies (dr_max_fill).
            v_drreco.push_back(recoDr >= dr_max ? dr_max_fill : recoDr);
            v_drgen .push_back(genDr  >= dr_max ? dr_max_fill : (double)genDr);
        }
    }
    std::cout << "\r  100 %" << std::endl;

    std::cout << "\n--- Jet statistics ---" << std::endl;
    std::cout << "  rows (true 2b jets):        " << n           << std::endl;
    std::cout << "  reco_pass:                  " << n_reco_pass << std::endl;
    std::cout << "  gen_pass:                   " << n_gen_pass  << std::endl;
    std::cout << "  both (enters the response): " << n_both      << std::endl;
    if (n_bad_B)
        std::cout << "  WARNING: " << n_bad_B << " pairs passed the gates with an unusable pT"
                  << " (non-positive or sentinel) and were dropped" << std::endl;
    if (n_both == 0) {
        std::cerr << "ERROR: no jet passes both gates -- nothing to study" << std::endl;
        return;
    }

    // ---- Resolution ----
    std::cout << "\n--- Resolution of B = pT_lead/(pT1+pT2) ---" << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  mean (B_reco - B_gen): " << h_res->GetMean()
              << "   RMS: " << h_res->GetRMS() << std::endl;
    std::cout << "  B_gen  mean " << h_Bgen_w->GetMean()  << "  RMS " << h_Bgen_w->GetRMS()
              << "   (EEC weighted)" << std::endl;
    std::cout << "  B_reco mean " << h_Breco_w->GetMean() << "  RMS " << h_Breco_w->GetRMS()
              << "   (EEC weighted)" << std::endl;
    std::cout << "  B_gen  mean " << h_Bgen_raw->GetMean() << "  RMS " << h_Bgen_raw->GetRMS()
              << "   (unweighted)" << std::endl;

    // Resolution as a function of B_gen: a bin narrower than the local resolution cannot be
    // unfolded stably, so this is what sets the minimum bin width.
    std::cout << "\n  resolution vs B_gen (slice RMS):" << std::endl;
    std::cout << "    B_gen range     |  entries  |   mean   |   RMS" << std::endl;
    std::cout << "    ----------------+-----------+----------+--------" << std::endl;
    for (int ix = 1; ix <= h_res_vs_B->GetNbinsX(); ++ix) {
        TH1D *slice = h_res_vs_B->ProjectionY(Form("_py%d", ix), ix, ix);
        if (slice->GetEntries() < 50) { delete slice; continue; }
        std::cout << "    " << std::setw(6) << h_res_vs_B->GetXaxis()->GetBinLowEdge(ix)
                  << " - " << std::setw(6) << h_res_vs_B->GetXaxis()->GetBinUpEdge(ix)
                  << " | " << std::setw(9) << (long)slice->GetEntries()
                  << " | " << std::setw(8) << slice->GetMean()
                  << " | " << std::setw(7) << slice->GetRMS() << std::endl;
        delete slice;
    }

    // ---- Binning scan ----
    std::cout << "\n--- Candidate binnings ---" << std::endl;
    std::cout << "  purity    = N(reco in bin AND gen in bin) / N(reco in bin)" << std::endl;
    std::cout << "  stability = N(reco in bin AND gen in bin) / N(gen  in bin)" << std::endl;
    std::cout << "  (EEC weighted, the weight the response matrix would use)" << std::endl;

    // The benchmark FIRST, so every B number below is read against it rather than against
    // an intuition about what purity "should" be.
    std::vector<double> dr_edges(dr_binsVector, dr_binsVector + dr_binsVectorSize);
    BinQuality qdr = evaluate(dr_edges, v_drreco, v_drgen, v_w);
    printQuality("BENCHMARK: dR, the analysis's own 9 bins", qdr);
    std::cout << "    ^ this is the observable already being measured and unfolded."
              << " Read the B numbers against it." << std::endl;

    std::vector<BinQuality> results;
    std::vector<TString>    names;

    for (int nb : {3, 4, 5, 6, 7, 8, 9}) {
        BinQuality q = evaluate(uniformEdges(nb), v_Breco, v_Bgen, v_w);
        printQuality(Form("uniform %d", nb), q);
        results.push_back(q);
        names.push_back(Form("uniform_%d", nb));
    }
    // Equal-occupancy edges naturally give the high-B tail one wide bin, which is where the
    // migration is worst -- so these tend to beat a uniform binning with the SAME number of
    // bins, not just a finer one.
    for (int nb : {3, 4, 5, 6, 7, 8, 9}) {
        BinQuality q = evaluate(quantileEdges(h_Bgen_w, nb), v_Breco, v_Bgen, v_w);
        printQuality(Form("equal-occupancy %d (EEC weighted gen)", nb), q);
        results.push_back(q);
        names.push_back(Form("quantile_%d", nb));
    }

    // ---- Recommendation -------------------------------------------------------------
    // The bar is the dR measurement itself, not a round number: a B binning whose worst bin
    // is at least as pure and as stable as dR's worst bin is, by construction, no harder to
    // unfold than what this analysis already does successfully. Among those, take the one
    // with the MOST bins -- the finest binning that does not cost resolution.
    //
    // Printed, not applied: the binning that goes into binning_histos_small.h is a physics
    // choice (it also has to make sense against the theory curves being compared to), and
    // this is the input to that choice, not a substitute for it.
    const double bar_p = qdr.min_purity;
    const double bar_s = qdr.min_stability;
    std::cout << "\n--- Recommendation ---" << std::endl;
    std::cout << "  bar (dR's worst bin): purity " << bar_p
              << ", stability " << bar_s << std::endl;

    int best = -1;
    for (size_t i = 0; i < results.size(); ++i) {
        if (results[i].min_purity < bar_p || results[i].min_stability < bar_s) continue;
        if (best < 0 || results[i].purity.size() > results[best].purity.size()) best = (int)i;
    }
    if (best < 0) {
        std::cout << "  NO candidate B binning matches dR's worst bin on both metrics."
                  << std::endl;
        int b2 = 0;
        for (size_t i = 1; i < results.size(); ++i)
            if (std::min(results[i].min_purity, results[i].min_stability) >
                std::min(results[b2].min_purity, results[b2].min_stability)) b2 = (int)i;
        std::cout << "  Closest: " << names[b2] << std::endl;
        printQuality(names[b2], results[b2]);
        std::cout << "  => B migrates more than dR does. Either accept coarser bins than the"
                  << " candidates scanned here, or expect a harder unfold than the dR one."
                  << std::endl;
    } else {
        std::cout << "  finest B binning that matches or beats it: " << names[best]
                  << std::endl;
        std::cout << "  B_binsVector for binning_histos_small.h:" << std::endl;
        std::cout << "    " << std::setprecision(3);
        for (size_t i = 0; i < results[best].edges.size(); ++i)
            std::cout << results[best].edges[i]
                      << (i + 1 < results[best].edges.size() ? ", " : "\n");
    }

    // ---- Output ----
    TString outdir = TString::Format("%s/results/momentum_balance_study_%s_%s_upartv2",
                                     kBase, sample.Data(), generator.Data());
    gSystem->mkdir(outdir, true);

    TFile fout(outdir + "/momentum_balance_study.root", "RECREATE");
    h_Bgen_w->Write();  h_Bgen_raw->Write();
    h_Breco_w->Write(); h_Breco_raw->Write();
    h_res->Write();     h_res_vs_B->Write();  h_mig->Write();
    h_Bres_by_dr->Write();
    fout.Close();

    // Plots. Palette convention of this analysis: red = measurement, blue = particle-level
    // MC. Nothing here is a measurement, so gen is blue (it IS particle level) and reco
    // takes green rather than borrowing the measurement's red.
    gStyle->SetOptStat(0);

    TCanvas c1("c_B", "", 900, 700);
    h_Bgen_w->SetLineColor(kBlue + 1);   h_Bgen_w->SetLineWidth(2);
    h_Breco_w->SetLineColor(kGreen + 2); h_Breco_w->SetLineWidth(2);
    h_Bgen_raw->SetLineColor(kBlue + 1); h_Bgen_raw->SetLineStyle(2); h_Bgen_raw->SetLineWidth(2);
    // Unit area: the shapes are the point, the normalisations are not comparable.
    for (TH1D *h : {h_Bgen_w, h_Breco_w, h_Bgen_raw})
        if (h->Integral() > 0) h->Scale(1. / h->Integral());
    // Same axis label as the chain uses, from obsB().axis, so the study and the results
    // plots cannot end up labelling the same quantity two different ways.
    h_Bgen_w->SetTitle(Form(";%s;normalised", obsB().axis.Data()));
    h_Bgen_w->SetMaximum(1.4 * std::max({h_Bgen_w->GetMaximum(), h_Breco_w->GetMaximum(),
                                         h_Bgen_raw->GetMaximum()}));
    h_Bgen_w->Draw("hist");
    h_Breco_w->Draw("hist same");
    h_Bgen_raw->Draw("hist same");
    TLegend leg(0.45, 0.68, 0.88, 0.88);
    leg.SetBorderSize(0); leg.SetFillStyle(0);
    leg.AddEntry(h_Bgen_w,   "gen, EEC weighted",  "l");
    leg.AddEntry(h_Breco_w,  "reco, EEC weighted", "l");
    leg.AddEntry(h_Bgen_raw, "gen, unweighted",    "l");
    leg.Draw();
    c1.SaveAs(outdir + "/B_distributions.pdf");
    c1.SaveAs(outdir + "/B_distributions.png");

    TCanvas c2("c_res", "", 900, 700);
    h_res->SetLineColor(kGreen + 2); h_res->SetLineWidth(2);
    h_res->Draw("hist");
    TLatex lat; lat.SetNDC(); lat.SetTextSize(0.035);
    lat.DrawLatex(0.16, 0.85, Form("mean %.4f", h_res->GetMean()));
    lat.DrawLatex(0.16, 0.80, Form("RMS %.4f",  h_res->GetRMS()));
    c2.SaveAs(outdir + "/B_resolution.pdf");
    c2.SaveAs(outdir + "/B_resolution.png");

    TCanvas c3("c_mig", "", 900, 800);
    c3.SetRightMargin(0.14);
    c3.SetLogz();
    h_mig->Draw("colz");
    c3.SaveAs(outdir + "/B_migration.pdf");
    c3.SaveAs(outdir + "/B_migration.png");

    std::cout << "\nWritten to " << outdir << std::endl;
    std::cout << "(/data_CMS is not on the sshfs mount -- use sync_plots_to_eos.sh or scp)"
              << std::endl;
}
