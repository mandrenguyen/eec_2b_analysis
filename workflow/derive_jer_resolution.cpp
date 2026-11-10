#include "../Run3HeavyIonJetCalibrations/include/JetCorrector.h"

#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

const std::vector<double> kAbsEtaEdges = {0.0, 0.5, 1.0, 1.5, 2.0};
const std::vector<double> kGenPtEdges = {
    30, 40, 50, 60, 80, 100, 120, 150, 200, 250, 300, 400, 500, 650, 800, 1000};

int findBin(double x, const std::vector<double> &edges) {
  const auto it = std::upper_bound(edges.begin(), edges.end(), x);
  if (it == edges.begin() || it == edges.end()) return -1;
  return static_cast<int>(it - edges.begin()) - 1;
}

struct WidthResult {
  double mean = 0;
  double sigma = 0;
  double sigmaError = 0;
  int entries = 0;
  bool valid = false;
};

WidthResult fitCore(TH1D &h) {
  WidthResult out;
  out.entries = static_cast<int>(h.GetEntries());
  if (out.entries < 200) return out;

  double probability[3] = {0.16, 0.50, 0.84};
  double quantile[3] = {0, 0, 0};
  h.GetQuantiles(3, quantile, probability);
  double mean = quantile[1];
  double sigma = 0.5 * (quantile[2] - quantile[0]);
  if (!(sigma > 0) || !std::isfinite(sigma)) return out;

  TF1 gaussian("core_gaussian", "gaus", mean - 1.5 * sigma, mean + 1.5 * sigma);
  gaussian.SetParameters(h.GetMaximum(), mean, sigma);
  for (int iteration = 0; iteration < 3; ++iteration) {
    gaussian.SetRange(mean - 1.5 * sigma, mean + 1.5 * sigma);
    const int status = h.Fit(&gaussian, "QNR");
    if (status != 0) return out;
    mean = gaussian.GetParameter(1);
    sigma = std::abs(gaussian.GetParameter(2));
    if (!(sigma > 0) || !std::isfinite(sigma)) return out;
  }

  out.mean = mean;
  out.sigma = sigma;
  out.sigmaError = gaussian.GetParError(2);
  out.valid = std::isfinite(out.sigmaError) && out.sigmaError > 0;
  return out;
}

std::vector<std::string> jecFiles() {
  const TString directory = TString(gSystem->DirName(__FILE__))
                          + "/../Run3HeavyIonJetCalibrations/txt/2024ppRef/";
  return {
      std::string(directory.Data()) + "Prompt24HIpp_V1_MC_L1FastJet_AK4PF.txt",
      std::string(directory.Data()) + "Prompt24HIpp_V2_MC_L2Relative_AK4PF.txt"};
}

}  // namespace

void derive_jer_resolution(
    const char *inputPattern =
        "/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/"
        "Pythia8_UParTV2_chunks/merged_block_0000_Pythia8_UParTV2_rho.root",
    const char *outputRoot = "calibrations/jer_resolution_diagnostics.root",
    const char *outputText =
        "calibrations/Prompt24HIpp_Derived_MC_PtResolution_AK4PF.txt",
    Long64_t maxEvents = -1) {
  TChain jets("ak4PFJetAnalyzer/t");
  TChain events("hiEvtAnalyzer/HiTree");
  if (jets.Add(inputPattern) <= 0 || events.Add(inputPattern) <= 0)
    throw std::runtime_error("No input files matched the requested pattern");
  if (jets.GetEntries() != events.GetEntries())
    throw std::runtime_error("Jet and event trees have different entry counts");

  Int_t nref = 0;
  Float_t rawpt[500] = {}, jteta[500] = {}, jtphi[500] = {}, jtarea[500] = {};
  Float_t refpt[500] = {}, refdrjt[500] = {};
  Float_t rho = 0;
  jets.SetBranchStatus("*", 0);
  for (const char *name : {"nref", "rawpt", "jteta", "jtphi", "jtarea", "refpt", "refdrjt"})
    jets.SetBranchStatus(name, 1);
  jets.SetBranchAddress("nref", &nref);
  jets.SetBranchAddress("rawpt", rawpt);
  jets.SetBranchAddress("jteta", jteta);
  jets.SetBranchAddress("jtphi", jtphi);
  jets.SetBranchAddress("jtarea", jtarea);
  jets.SetBranchAddress("refpt", refpt);
  jets.SetBranchAddress("refdrjt", refdrjt);
  events.SetBranchStatus("*", 0);
  events.SetBranchStatus("rho", 1);
  events.SetBranchAddress("rho", &rho);

  JetCorrector jec;
  const std::vector<std::string> files = jecFiles();
  for (const std::string &file : files)
    if (gSystem->AccessPathName(file.c_str()))
      throw std::runtime_error("Missing JEC file: " + file);
  jec.Initialize(files);

  const int nEta = static_cast<int>(kAbsEtaEdges.size()) - 1;
  const int nPt = static_cast<int>(kGenPtEdges.size()) - 1;
  std::vector<std::vector<std::unique_ptr<TH1D>>> response(nEta);
  std::vector<std::vector<std::unique_ptr<TH1D>>> genPt(nEta);
  for (int ie = 0; ie < nEta; ++ie) {
    for (int ip = 0; ip < nPt; ++ip) {
      response[ie].emplace_back(new TH1D(
          Form("response_eta%d_pt%d", ie, ip),
          Form("%.1f<|#eta|<%.1f, %.0f<p_{T}^{gen}<%.0f;"
               "(p_{T}^{reco}-p_{T}^{gen})/p_{T}^{gen};jets",
               kAbsEtaEdges[ie], kAbsEtaEdges[ie + 1], kGenPtEdges[ip],
               kGenPtEdges[ip + 1]),
          240, -0.8, 0.8));
      genPt[ie].emplace_back(new TH1D(Form("genpt_eta%d_pt%d", ie, ip), "", 100,
                                     kGenPtEdges[ip], kGenPtEdges[ip + 1]));
    }
  }

  Long64_t entries = jets.GetEntries();
  if (maxEvents >= 0) entries = std::min(entries, maxEvents);
  Long64_t accepted = 0;
  for (Long64_t entry = 0; entry < entries; ++entry) {
    jets.GetEntry(entry);
    events.GetEntry(entry);
    if (!std::isfinite(rho)) continue;
    for (int jet = 0; jet < nref; ++jet) {
      if (!(rawpt[jet] > 0) || !(refpt[jet] > 0) || refdrjt[jet] >= 0.2 ||
          !std::isfinite(jteta[jet]) || std::abs(jteta[jet]) >= 2.0) continue;
      const int ie = findBin(std::abs(jteta[jet]), kAbsEtaEdges);
      const int ip = findBin(refpt[jet], kGenPtEdges);
      if (ie < 0 || ip < 0) continue;
      jec.SetJetPT(rawpt[jet]);
      jec.SetJetEta(jteta[jet]);
      jec.SetJetPhi(jtphi[jet]);
      jec.SetJetArea(jtarea[jet]);
      jec.SetRho(rho);
      const double correctedPt = jec.GetCorrectedPT();
      if (!(correctedPt > 0) || !std::isfinite(correctedPt)) continue;
      const double r = (correctedPt - refpt[jet]) / refpt[jet];
      if (!std::isfinite(r) || std::abs(r) >= 0.8) continue;
      response[ie][ip]->Fill(r);
      genPt[ie][ip]->Fill(refpt[jet]);
      ++accepted;
    }
    if (entry && entry % 500000 == 0)
      std::cout << "Processed " << entry << "/" << entries << " events" << std::endl;
  }

  TFile output(outputRoot, "RECREATE");
  if (output.IsZombie()) throw std::runtime_error("Cannot create diagnostics ROOT file");
  std::vector<std::unique_ptr<TGraphErrors>> graphs;
  std::vector<std::unique_ptr<TF1>> fits;
  for (int ie = 0; ie < nEta; ++ie) {
    graphs.emplace_back(new TGraphErrors());
    TGraphErrors &graph = *graphs.back();
    graph.SetName(Form("resolution_eta%d", ie));
    graph.SetTitle(Form("%.1f<|#eta|<%.1f;p_{T}^{gen} [GeV];#sigma(p_{T})/p_{T}",
                        kAbsEtaEdges[ie], kAbsEtaEdges[ie + 1]));
    for (int ip = 0; ip < nPt; ++ip) {
      const WidthResult width = fitCore(*response[ie][ip]);
      response[ie][ip]->Write();
      genPt[ie][ip]->Write();
      if (!width.valid) continue;
      const int point = graph.GetN();
      const double x = genPt[ie][ip]->GetMean();
      graph.SetPoint(point, x, width.sigma);
      graph.SetPointError(point, genPt[ie][ip]->GetMeanError(), width.sigmaError);
    }
    fits.emplace_back(new TF1(
        Form("resolution_fit_eta%d", ie),
        "sqrt([0]*abs([0])/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2])", 30, 1000));
    TF1 &fit = *fits.back();
    fit.SetParameters(1.0, 0.5, 0.03, -0.5);
    fit.SetParLimits(3, -2.0, -0.05);
    graph.Fit(&fit, "QNR");
    graph.Write();
    fit.Write();
  }
  output.Close();

  std::ofstream payload(outputText);
  if (!payload) throw std::runtime_error("Cannot create resolution payload");
  payload << "{1 JetAbsEta 1 JetPt sqrt([0]*abs([0])/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2]) Resolution}\n";
  payload << std::setprecision(10);
  for (int ie = 0; ie < nEta; ++ie) {
    const TF1 &fit = *fits[ie];
    payload << kAbsEtaEdges[ie] << " " << kAbsEtaEdges[ie + 1]
            << " 6 30 1000 " << fit.GetParameter(0) << " "
            << fit.GetParameter(1) << " " << fit.GetParameter(2) << " "
            << fit.GetParameter(3) << "\n";
  }
  payload.close();
  std::cout << "Accepted matched jets: " << accepted << std::endl;
  std::cout << "Wrote " << outputRoot << " and " << outputText << std::endl;
}
