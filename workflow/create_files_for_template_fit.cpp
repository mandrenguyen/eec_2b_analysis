//For both
#include "binning_histos_small.h"
#include "tTree.h"
#include <iostream>
#include <Math/Vector4D.h>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <random>
#include "Math/VectorUtil.h"
#include "TString.h"
#include "TStyle.h"

//To create templates
#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <string>

// To create the response matrix
#include <cstdlib>
#include <cmath>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include <algorithm> // std::find
#include <functional> // std::not_equal_to<int>.
#include "TMatrixD.h"
#include "central_selections.h"
#include "tracking_efficiency_syst.h"
#include "sv_sf.h"         // SV b scale factor vs aggregate pT, applied per SV (MC only)
#include "track_bdt_sf.h"  // effective track-BDT WPs for signal/non-signal tracks (MC only)
#include "observables.h"   // the measured observables (dr, B) and the axis abstraction
#include "agg_nn.h"        // aggregation NN (AGG_NN=on): the two B candidates for the observables


namespace EecWeight {

inline bool &enabled()
{
    static bool on = true;
    return on;
}

inline void enable(bool on) { enabled() = on; }

// The weight for one jet, from the pTs of its two aggregated B hadrons; n is cfg.n.
inline double value(double pt1, double pt2, double n)
{
    return enabled() ? std::pow(pt1 * pt2, n) : 1.0;
}

// "" for the EEC analysis, so its filenames are unchanged; "_noeecw" for the yields.
// Must stay in step with eecWeightTag() in result_paths.h, which is what the readers use.
inline TString tag()
{
    return enabled() ? TString("") : TString("_noeecw");
}

inline void print()
{
    if (enabled()) {
        std::cout << "EEC weight: ON (nominal -- measuring the EEC)" << std::endl;
        return;
    }
    std::cout << "EEC weight: OFF -- every jet enters with weight 1, so this job measures "
              << "YIELDS (dN/dr), not the EEC. Output tagged '" << tag() << "'." << std::endl;
}

} 
#include "../Run3HeavyIonJetCalibrations/include/JetCorrector.h"
#include "../Run3HeavyIonJetCalibrations/include/JetSmearer.h"
#include "TLeaf.h"     // the event rho is read through its leaf
#include <stdexcept>
namespace Jec {

inline bool &enabled()
{
    static bool on = false;
    return on;
}

inline void enable(bool on) { enabled() = on; }

inline TString dir()
{
    return TString(gSystem->DirName(__FILE__)) + "/../Run3HeavyIonJetCalibrations/txt/2024ppRef/";
}

// The correction levels for this sample, in the order they are applied.
inline std::vector<std::string> files(bool isMC)
{
    std::vector<std::string> f = {
        (dir() + "Prompt24HIpp_V1_MC_L1FastJet_AK4PF.txt").Data(),
        (dir() + "Prompt24HIpp_V2_MC_L2Relative_AK4PF.txt").Data() };
    if (!isMC) f.push_back((dir() + "Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt").Data());
    return f;
}

inline void print(bool isMC)
{
    if (!enabled()) {
        std::cout << "JEC: OFF -- jtpt as stored in the forest" << std::endl;
        return;
    }
    std::cout << "JEC: ON -- jtpt recomputed from rawpt with:" << std::endl;
    for (const std::string &s : files(isMC)) std::cout << "   " << s << std::endl;
}

} 
namespace JesSyst {

enum Mode { Nominal = 0, Up = 1, Down = -1 };

inline Mode &mode() { static Mode value = Nominal; return value; }

inline TString uncertaintyFile()
{
    return TString(gSystem->DirName(__FILE__))
         + "/calibrations/Prompt24HIpp_Placeholder1Percent_JESUncertainty_AK4PF.txt";
}

inline void configure(const char *name, bool applyJec, bool isMC)
{
    const TString value = name ? TString(name) : TString("nominal");
    if (value == "nominal") mode() = Nominal;
    else if (value == "up") mode() = Up;
    else if (value == "down") mode() = Down;
    else throw std::runtime_error("JES variation must be nominal, up, or down; got: "
                                  + std::string(value.Data()));
    if (mode() != Nominal && !isMC)
        throw std::runtime_error("JES up/down is allowed only for MC; data must use nominal");
    if (mode() != Nominal && !applyJec)
        throw std::runtime_error("JES up/down requires apply_jec=true");
}

inline bool varied() { return mode() != Nominal; }
inline TString tag() { return mode() == Up ? "_jesUp" : mode() == Down ? "_jesDown" : ""; }

inline double apply(double nominalPt, const std::pair<double, double> &uncertainty)
{
    if (mode() == Up) return nominalPt * (1.0 + uncertainty.second);
    if (mode() == Down) return nominalPt * (1.0 - uncertainty.first);
    return nominalPt;
}

inline void print()
{
    if (mode() == Nominal) {
        std::cout << "JES variation: nominal" << std::endl;
        return;
    }
    std::cout << "JES variation: " << (mode() == Up ? "UP" : "DOWN")
              << " using " << uncertaintyFile() << std::endl;
}

}
namespace JerSyst {

enum Mode { Off = 0, Nominal, Up, Down };

inline Mode &mode() { static Mode value = Off; return value; }

inline TString resolutionFile()
{
    return TString(gSystem->DirName(__FILE__))
         + "/calibrations/Prompt24HIpp_Derived_MC_PtResolution_AK4PF.txt";
}

inline TString scaleFactorFile()
{
    return TString(gSystem->DirName(__FILE__))
         + "/calibrations/Prompt24HIpp_Estimated_MC_SF_AK4PF.txt";
}

inline void configure(const char *name, bool applyJec, bool isMC)
{
    const TString value = name ? TString(name) : TString("off");
    if (value == "off") mode() = Off;
    else if (value == "nominal") mode() = Nominal;
    else if (value == "up") mode() = Up;
    else if (value == "down") mode() = Down;
    else throw std::runtime_error("JER variation must be off, nominal, up, or down; got: "
                                  + std::string(value.Data()));
    if (mode() != Off && !isMC)
        throw std::runtime_error("JER smearing is allowed only for MC; data must use off");
    if (mode() != Off && !applyJec)
        throw std::runtime_error("JER smearing requires apply_jec=true");
    if ((mode() == Up || mode() == Down) && JesSyst::varied())
        throw std::runtime_error("JES and JER uncertainty shifts must not be varied simultaneously");
}

inline bool enabled() { return mode() != Off; }
inline TString tag()
{
    if (mode() == Nominal) return "_jerNominal";
    if (mode() == Up) return "_jerUp";
    if (mode() == Down) return "_jerDown";
    return "";
}

inline Variation variation()
{
    if (mode() == Up) return Variation::UP;
    if (mode() == Down) return Variation::DOWN;
    return Variation::NOMINAL;
}

inline void print()
{
    if (mode() == Off) {
        std::cout << "JER smearing: OFF" << std::endl;
        return;
    }
    const char *label = mode() == Up ? "UP" : mode() == Down ? "DOWN" : "nominal";
    std::cout << "JER smearing: " << label << " (hybrid method)" << std::endl
              << "   resolution: " << resolutionFile() << std::endl
              << "   scale factor: " << scaleFactorFile() << std::endl;
}

}
namespace JetVeto {

inline bool &enabled() { static bool on = false; return on; }
inline TH2D *&map() { static TH2D *m = nullptr; return m; }

inline TString file()
{
    return TString(gSystem->DirName(__FILE__)) + "/../Winter24Prompt24_2024BCDEFGHI.root";
}

inline void enable(bool on)
{
    enabled() = on;
    if (!on || map()) return;
    TFile *f = TFile::Open(file());
    if (!f || f->IsZombie()) { delete f; throw std::runtime_error(("Could not open jet veto map: " + file()).Data()); }
    TH2D *source = dynamic_cast<TH2D *>(f->Get("jetvetomap_all"));
    if (!source) { f->Close(); delete f; throw std::runtime_error(("Histogram jetvetomap_all is missing from: " + file()).Data()); }
    map() = (TH2D *) source->Clone("jetVetoMap_2024");
    map()->SetDirectory(nullptr);
    f->Close(); delete f;
}

// true = the jet is kept. Always true when off.
inline bool pass(double eta, double phi)
{
    if (!enabled()) return true;
    if (!std::isfinite(eta) || !std::isfinite(phi)) return false;
    return map()->GetBinContent(map()->FindBin(eta, phi)) == 0.0;
}

inline void print()
{
    if (enabled()) std::cout << "Jet veto map: ON -- jetvetomap_all from " << file() << std::endl;
    else           std::cout << "Jet veto map: OFF" << std::endl;
}

} // namespace JetVeto

// RooUnfold comes from the private build; setup_roounfold_env.sh puts it on the paths.

#if __has_include("RooUnfold.h") && __has_include("RooUnfoldResponse.h")
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#define HAS_ROOUNFOLD 1
#else
#define HAS_ROOUNFOLD 0
// Named distinctly from the real class so ROOT autoload can never clash with it.
class RooUnfoldResponseFallback {
 public:
  template <typename... Args>
  RooUnfoldResponseFallback(Args...) {}
  template <typename... Args>
  void Fill(Args...) {}
  void Write() {}
  TMatrixD Mresponse() const { return TMatrixD(); }
};
#define RooUnfoldResponse RooUnfoldResponseFallback
#endif
// No longer part of the output file names; only the disabled make_templates /
// create_response_templatefit paths below still call it.
TString kinLabel(const AnalysisConfig& cfg) {
  return Form("_%g_%g_%d", cfg.kin.ptLow, jtpt_max, jtpt_binsVectorSize - 1);
}

struct Vertex {
    ROOT::Math::PtEtaPhiMVector p4;
    std::vector<ROOT::Math::PtEtaPhiMVector> tracks;
    std::vector<int> trkMatchSta, trkSvtxId;
};

//given two vectors, computes deltaR squared
auto deltaR2 = [](const ROOT::Math::PtEtaPhiMVector &a,
                  const ROOT::Math::PtEtaPhiMVector &b) {
    double dEta = a.Eta() - b.Eta();
    double dPhi = std::remainder(a.Phi() - b.Phi(), 2*M_PI);
    return dEta*dEta + dPhi*dPhi;
};

// Mass of a charged particle from its |PDG id|; -1 for an id not in the list.
inline double pdgMass(int pid) {
  switch (pid) {
    case 211:  return 0.139570; case 13:   return 0.105658; case 11:   return 0.000510;
    case 2212: return 0.938272; case 321:  return 0.493677; case 3112: return 1.19744;
    case 3222: return 1.18937;  case 3312: return 1.32171;  case 3334: return 1.67245;
    default:   return -1.;
  }
}

// Four-vector of reco track itrk; an id not in pdgMass() gets the pion mass (safety).
inline ROOT::Math::PtEtaPhiMVector recoTrackP4(tTree& t, Int_t itrk) {
  ROOT::Math::PtEtaPhiMVector v1;
  v1.SetEta(t.trkEta[itrk]); v1.SetPt(t.trkPt[itrk]); v1.SetPhi(t.trkPhi[itrk]);
  const double m = pdgMass(std::abs(t.trkPdgId[itrk]));
  v1.SetM(m >= 0. ? m : 0.139570);
  return v1;
}

// ---- disabled (kept for reference): printvtx, SV-merging debug printout ----
#if 0
//printing function for a Vertex object
void printvtx (const Vertex& vertex,
              const TString& txt,
	       int vtxnb,  Int_t& ijet, Long64_t& ient) {
  for (size_t itrk = 0; itrk < vertex.trkSvtxId.size(); itrk++) {
    std::cout<<"ient : "<<ient <<"  ijet : "<<ijet << txt << "vertex number: "<<vtxnb<<"  track original vertex : "<< vertex.trkSvtxId[itrk] <<"    and track original status : "<<vertex.trkMatchSta[itrk] << " pT = " << vertex.tracks[itrk].Pt() <<  std::endl;
  }
}
#endif // ---- end disabled: printvtx ----

bool skipMC_event(double jtpt, double pthat) { // -- Since it is applied all the time, at event level, keep it general please.
  return (pthat < 0.40 * jtpt);}

bool skipMC_gen(double refpt) {// Indeed, this is just a saftey check like refpt > 0
  return !(refpt > 0);}

void PartialBsAggregation(std::vector<ROOT::Math::PtEtaPhiMVector>& hadrons_4vec, std::vector<Int_t>& hadrons_stat, tTree& t, Int_t ijet){
  // Only using gen level info
  // Add tacks of status >= 100 to make Bs.
  hadrons_4vec.clear();
  hadrons_stat.clear();
  for (Int_t itrk = 0; itrk < t.nrefTrk; itrk++) {
    // Track must belong to this jet
    if (t.refTrkJetId[itrk] != ijet) continue;
    // pT cut
    if (t.refTrkPt[itrk] < 1) continue;
    // Assign mass
    const int pid = std::abs(t.refTrkPdgId[itrk]);
    double mass = pdgMass(pid);
    if (mass < 0.) { std::cout << "PDG:" << pid << std::endl; mass = 0.0; }
    // 4 vector
    ROOT::Math::PtEtaPhiMVector v(t.refTrkPt [itrk], t.refTrkEta[itrk], t.refTrkPhi[itrk], mass);

    Int_t status = t.refTrkSta[itrk];
    if (status < 100) continue;
    else {
      auto it = std::find(hadrons_stat.begin(), hadrons_stat.end(), status);
      if (it == hadrons_stat.end()) {
        hadrons_stat.push_back(status);
	      hadrons_4vec.push_back(v);
      }
      else {
        size_t index = std::distance(hadrons_stat.begin(), it);
        hadrons_4vec[index] += v;}
    }
  }
}
double genNonBTrackPt(tTree& t, Int_t ijet){
  double sum = 0.;
  for (Int_t itrk = 0; itrk < t.nrefTrk; itrk++) {
    if (t.refTrkJetId[itrk] != ijet) continue;
    if (t.refTrkPt[itrk] < 1) continue;
    if (t.refTrkSta[itrk] >= 100) continue;   // from a B decay: that is the b side
    sum += t.refTrkPt[itrk];
  }
  return sum;
}

// ---- MatchingTracksAggregation: RECO tracks grouped by their TRUE B (trkMatchSta >= 100) ----
// perfect assignment, i.e. the best any aggregation could do with these tracks. Only feeds the
// orPt*/orEta*/orPhi* branches of the aggregation ntuple
// no histogram uses it
// #if 0
void MatchingTracksAggregation(std::vector<ROOT::Math::PtEtaPhiMVector>& hadrons_4vec, std::vector<Int_t>& hadrons_stat, tTree& t, Int_t ijet, Int_t ient){
  hadrons_4vec.clear();
  hadrons_stat.clear();

  for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {       //loop over trks
    // Track must belong to this jet
    if (t.trkJetId[itrk] != ijet) continue;
    // pT cut
    if (t.trkPt[itrk] < 1) continue;
    if (t.trkMatchSta[itrk] <100 ) continue;
    //if (t.trkSvtxId[itrk] <0 ) continue ;      //trk must be in SV
// Assign mass
    double mass = 0.0;
    int pid = std::abs(t.trkMatchPdgId[itrk]);
    if      (pid == 211)  mass = 0.139570;
    else if (pid == 13)   mass = 0.105658;
    else if (pid == 11)   mass = 0.000510;
    else if (pid == 2212) mass = 0.938272;
    else if (pid == 321)  mass = 0.493677;
    else if (pid == 3112) mass = 1.19744;
    else if (pid == 3222) mass = 1.18937;
    else if (pid == 3312) mass = 1.32171;
    else if (pid == 3334) mass = 1.67245;
    else if (pid == 1)   continue;
        //skip not matched
    // else std::cout << "PDG:" << pid << std::endl;
    

    // 4 vector
    ROOT::Math::PtEtaPhiMVector v(t.trkPt [itrk], t.trkEta[itrk], t.trkPhi[itrk], mass);

    Int_t status = t.trkMatchSta[itrk];

    if (status < 100) {   continue;
      // I don't aggregate, just store each particle separately
      //hadrons_stat.push_back(status);
      //hadrons_4vec.push_back(v);
    }
    else {
      // I aggregate only if status >= 100
      auto it = std::find(hadrons_stat.begin(), hadrons_stat.end(), status);      //check if status alreday exists in hadrons_stat
      if (it == hadrons_stat.end()) {                     //if not, create new hadrons
        hadrons_stat.push_back(status);
        hadrons_4vec.push_back(v);
      }
      else {
        size_t index = std::distance(hadrons_stat.begin(), it);
        hadrons_4vec[index] += v;}                                                     //if yes, add to existing
    }
  }

}
// #endif 

void groupVertexes1(std::vector<Vertex>& vertices)
{
    if (vertices.size() < 2) return;

    while (vertices.size() > 2) {

        double min_distance = std::numeric_limits<double>::infinity();
        size_t index1 = 0;
        size_t index2 = 1;

        // find closest pair in (eta, phi)
        for (size_t i = 0; i < vertices.size(); ++i) {
            for (size_t j = i + 1; j < vertices.size(); ++j) {

                double dEta = vertices[i].p4.Eta() - vertices[j].p4.Eta();
                double dPhi = std::remainder(vertices[i].p4.Phi() - vertices[j].p4.Phi(), 2*M_PI);

                double dist = dEta*dEta + dPhi*dPhi;

                if (dist < min_distance) {
                    min_distance = dist;
                    index1 = i;
                    index2 = j;
                }
            }
        }

	// ---- merge vertex index2 into index1 ----
	vertices[index1].p4 += vertices[index2].p4;

	// merge track lists
	auto &t1 = vertices[index1].tracks;
	auto &t2 = vertices[index2].tracks;
	t1.insert(t1.end(), t2.begin(), t2.end());

	// merge match status lists
	auto &s1 = vertices[index1].trkMatchSta;
	auto &s2 = vertices[index2].trkMatchSta;
	s1.insert(s1.end(), s2.begin(), s2.end());

	//merge vertex id
	vertices[index1].trkSvtxId.insert(
	vertices[index1].trkSvtxId.end(),
        vertices[index2].trkSvtxId.begin(),
        vertices[index2].trkSvtxId.end());

	// remove merged vertex
	vertices.erase(vertices.begin() + index2);

    }
}

// ---- Shared by makeSvtxs and makeSvtxs_withBDT ----
typedef std::unordered_map<Int_t, std::vector<ROOT::Math::PtEtaPhiMVector>> SvTrackMap;
typedef std::unordered_map<Int_t, std::vector<int>> SvStatusMap;

// One Vertex per original SV id; tracks[i] <-> trkMatchSta[i] <-> trkSvtxId[i] correspond.
// An SV whose tracks do not all come from the same B counts in sv_fail.
std::vector<Vertex> buildVertices(SvTrackMap& secVtxs, SvStatusMap& secVtxsMatchSta, double& sv_fail) {
    std::vector<Vertex> vertices;
    for (auto &it : secVtxs) {           //loop over the map
      int originalSvId = it.first;       //saves it.first = vertex ID (map key)
      auto& v = secVtxsMatchSta[originalSvId];
      if (std::adjacent_find(v.begin(), v.end(),std::not_equal_to<int>()) != v.end()) {      // if not all tracks in the sv come from the same B
        sv_fail += 1.0;
        //return empty;
      }
      Vertex vtx;
      vtx.p4 = ROOT::Math::PtEtaPhiMVector();
      vtx.tracks.reserve(it.second.size()); vtx.trkMatchSta.reserve(it.second.size()); vtx.trkSvtxId.reserve(it.second.size());
      for (size_t i = 0; i < it.second.size(); ++i) {      //loop over tracks in this vertex (it.second)
        vtx.p4 += it.second[i];
        vtx.tracks.push_back(it.second[i]);
        vtx.trkMatchSta.push_back(v[i]);  // reuse already-looked-up reference
        vtx.trkSvtxId.push_back(originalSvId);
      }
      vertices.push_back(std::move(vtx));
    }
    return vertices;
}

// The most frequent match status among a vertex's SV tracks (the lowest on a tie) and, in
// maxCount, how many tracks carry it. maxCount < number of tracks = a merging failure.
int dominantStatus(const Vertex& vtx, int& maxCount) {
    std::map<int,int> countMap;   //count nb of tracks of the same MatchSta within a reconstructed vertex
    for (size_t i = 0; i < vtx.tracks.size(); ++i) countMap[vtx.trkMatchSta[i]]++;
    int domSta = -1;
    maxCount = 0;
    for (const auto &kv : countMap)
      if (kv.second > maxCount) { maxCount = kv.second; domSta = kv.first; }
    return domSta;
}

// Add the particles that are not part of an SV to the nearer vertex (a tie adds to neither).
// agg_fail counts those whose match status differs from the vertex's first SV track.
void absorbNoSvTracks(const std::vector<ROOT::Math::PtEtaPhiMVector>& no_sv_list, const std::vector<int>& no_sv_sta_list,
                      Vertex& vertex0, Vertex& vertex1, double& agg_fail) {
    int i = 0;
    for (const auto& v1 : no_sv_list) {
      double d0 = deltaR2(v1, vertex0.p4), d1 = deltaR2(v1, vertex1.p4);
      if (d0 < d1) {
        vertex0.p4 += v1;
        if (no_sv_sta_list[i] != vertex0.trkMatchSta[0]) agg_fail += 1;
      }
      else if (d0 > d1) {
        vertex1.p4 += v1;
        if (no_sv_sta_list[i] != vertex1.trkMatchSta[0]) agg_fail += 1;
      }
      i += 1;
    }
}

// Track aggregation: makeSvtxs uses the status, makeSvtxs_withBDT the BDT score.

vector<ROOT::Math::PtEtaPhiMVector> makeSvtxs(
  tTree& t,     //event
  Int_t& ijet,  //jet nb
  Long64_t& ient,
  double& agg_fail,     //dominant trkMatchSta fraction for vertex 0
  double& nb_sv ,    //   ..      ..           ..     ..    ..   1
  double& sv_fail,
  double& merge_fail){

  SvTrackMap secVtxs;
  SvStatusMap secVtxsMatchSta;
  std::vector<ROOT::Math::PtEtaPhiMVector> no_sv_list;
  std::vector<int> no_sv_sta_list;
  no_sv_list.reserve(t.ntrk); no_sv_sta_list.reserve(t.ntrk);

    // Loop over tracks
    for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {
      //cuts for tracks
        if (t.trkJetId[itrk] != ijet) continue;
        if (t.trkPt[itrk] < 1) continue;
        if (t.trkMatchSta[itrk] <100 ) continue;
        // Tracking-efficiency systematic, after the track cuts -- same convention as
        // makeSvtxs_withBDT. No-op in the nominal. See tracking_efficiency_syst.h.
        if (TrkEffSyst::dropTrack(ient, itrk)) continue;

        //build v1, the 4 vec of the track
        const ROOT::Math::PtEtaPhiMVector v1 = recoTrackP4(t, itrk);
        // no_sv_tracks, else map vertex id to list of tracks
        if (t.trkSvtxId[itrk] < 0) { no_sv_list.push_back(v1); no_sv_sta_list.push_back(t.trkMatchSta[itrk]); }
        else { secVtxs[t.trkSvtxId[itrk]].push_back(v1); secVtxsMatchSta[t.trkSvtxId[itrk]].push_back(t.trkMatchSta[itrk]); }
    }

    if (secVtxs.size() <  2) {
      nb_sv += 1 ;
      if (nb_sv < 15 ) {}
      return {};    }

    // Build Vertex objects
    std::vector<Vertex> vertices = buildVertices(secVtxs, secVtxsMatchSta, sv_fail);

    // Merge until only 2 vertices remain
    if (vertices.size() > 2) groupVertexes1(vertices);
    if (vertices.size() != 2) {
       std::cout<<" merging pb"<<std::endl;
       //return empty;
    }

    //vertices names
    auto &vertex0 = vertices[0], &vertex1 = vertices[1];
    int maxCount0, maxCount1;
    dominantStatus(vertex0, maxCount0); dominantStatus(vertex1, maxCount1);
    if (double(maxCount0) / double(vertex0.tracks.size()) != 1) merge_fail += 1;
    if (double(maxCount1) / double(vertex1.tracks.size()) != 1) {
      merge_fail += 1;
      //printvtx( vertex1 , " wrong merging ", 1);
      //return empty ;
    }

    absorbNoSvTracks(no_sv_list, no_sv_sta_list, vertex0, vertex1, agg_fail);

    // Return summed 4-vectors of the 2 vertices
    return {vertices[0].p4, vertices[1].p4};
}

// Status is used for quality plots and BDT cut for the merging.
vector<ROOT::Math::PtEtaPhiMVector> makeSvtxs_withBDT(
  tTree& t,     //event
  Int_t& ijet,  //jet nb
  Long64_t& ient,
  double& agg_fail,     //dominant trkMatchSta fraction for vertex 0
  double& nb_sv ,    //   ..      ..           ..     ..    ..   1
  double& sv_fail,
  double& merge_fail,
  TH1D* h_score_bkg, TH1D* h_score_sg,
  double* pt_nonb = nullptr,
  std::vector<int>* sv_sta = nullptr,
  bool allowSvSeed = true   // false: the plain BDT selection even with SV_SEED=on (nominalVertexMass)
){
  if (pt_nonb) *pt_nonb = 0.;
  if (sv_sta) sv_sta->clear();

  // Failure counters: nb_sv (fewer than 2 SVs), sv_fail / merge_fail / agg_fail (match status).
  SvTrackMap secVtxs;
  SvStatusMap secVtxsMatchSta;
  std::vector<ROOT::Math::PtEtaPhiMVector> no_sv_list;
  std::vector<int> no_sv_sta_list;
  no_sv_list.reserve(t.ntrk); no_sv_sta_list.reserve(t.ntrk);

  // SV_SEED=on (2026-11-05, Zoe): the vertices with >= 1 track passing the BDT keep ALL their tracks below.
  // No vertex appears that the nominal selection did not have. A track dropped by the tracking systematic
  // cannot seed one. Off: seedSvs stays empty, the loop below is the nominal one.
  static const bool kSvSeed = [] { const char *e = std::getenv("SV_SEED"); return e && TString(e) == "on"; }();   // svSeedOn() in result_paths.h
  std::set<Int_t> seedSvs;
  if (kSvSeed && allowSvSeed)
    for (Int_t itrk = 0; itrk < t.ntrk; itrk++)
      if (t.trkJetId[itrk] == ijet && t.trkPt[itrk] >= 1 && t.trkSvtxId[itrk] >= 0
          && TrackBdtSF::passes(t.trkBdtScore[itrk], t.trkMatchSta[itrk], t.trkSvtxId[itrk]) && !TrkEffSyst::dropTrack(ient, itrk))
        seedSvs.insert(t.trkSvtxId[itrk]);

    // Loop over tracks
    for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {
      //cuts for tracks
        if (t.trkJetId[itrk] != ijet) continue;
        if (t.trkPt[itrk] < 1) continue;
    //if (t.trkMatchSta[itrk] <100 ) continue;
    //if (t.trkMatchSta[itrk] <= 1) continue; //activate or not
    //control plots for the BDT score distribution in 0b and 1b jets
    //if (t.trkMatchSta[itrk] == 1)
      //{h_score_bkg->Fill(t.trkBdtScore[itrk], weight_tree);}
    //else if (t.trkMatchSta[itrk] >= 100)
     // {h_score_sg->Fill(t.trkBdtScore[itrk], weight_tree);}

    // -- Add trkbdt score
    // if (!TrackBdtSF::passes(t.trkBdtScore[itrk], t.trkMatchSta[itrk], t.trkSvtxId[itrk])) {   // until 2026-11-05
    if (!TrackBdtSF::passes(t.trkBdtScore[itrk], t.trkMatchSta[itrk], t.trkSvtxId[itrk]) && !seedSvs.count(t.trkSvtxId[itrk])) {
       if (pt_nonb && !TrkEffSyst::dropTrack(ient, itrk)) *pt_nonb += t.trkPt[itrk];
      continue;
    }

        // Tracking-efficiency systematic: throw the track away before anything is built
        // from it, so the SVs, m_2B, dR and the EEC weight all see a jet genuinely
        // missing 3% of the tracks that would have entered the B reconstruction.
        // No-op in the nominal. See tracking_efficiency_syst.h.
  
        if (TrkEffSyst::dropTrack(ient, itrk)) continue;

        //build v1, the 4 vec of the track
        const ROOT::Math::PtEtaPhiMVector v1 = recoTrackP4(t, itrk);
        // no_sv_tracks, else map vertex id to list of tracks
        if (t.trkSvtxId[itrk] < 0) { no_sv_list.push_back(v1); no_sv_sta_list.push_back(t.trkMatchSta[itrk]); }
        else { secVtxs[t.trkSvtxId[itrk]].push_back(v1); secVtxsMatchSta[t.trkSvtxId[itrk]].push_back(t.trkMatchSta[itrk]); }
    } // tracks loop

    if (secVtxs.size() <  2) {
      nb_sv += 1 ;
      // test: if (nb_sv < 15 ) {} // not needed!
      return {};
    }

    // Build Vertex objects
    std::vector<Vertex> vertices = buildVertices(secVtxs, secVtxsMatchSta, sv_fail);

    // Merge until only 2 vertices remain
    if (vertices.size() > 2) groupVertexes1(vertices); // does it give 2 vertices?
    if (vertices.size() != 2) {
       std::cout<<" merging pb"<<std::endl;
       //return empty;
    }

    //vertices names; dominant status of each, for the SV scale factor
    auto &vertex0 = vertices[0], &vertex1 = vertices[1];
    int maxCount0, maxCount1;
    const int domSta0 = dominantStatus(vertex0, maxCount0), domSta1 = dominantStatus(vertex1, maxCount1);
    if (double(maxCount0) / double(vertex0.tracks.size()) != 1) merge_fail += 1;
    if (sv_sta) { sv_sta->push_back(domSta0); sv_sta->push_back(domSta1); }
    if (double(maxCount1) / double(vertex1.tracks.size()) != 1) {
      merge_fail += 1;
      //printvtx( vertex1 , " wrong merging ", 1);
      //return empty ;
    }

    absorbNoSvTracks(no_sv_list, no_sv_sta_list, vertex0, vertex1, agg_fail);

    // Return summed 4-vectors of the 2 vertices
    return {vertices[0].p4, vertices[1].p4};
}

// SV_SEED_M2B=nominal (2026-11-06, with SV_SEED=on): the fit variable m2B from the vertices of the plain BDT
// selection, for the EEC and A alike; pT, dR and the EEC weight keep the SV-seed vertices. Same jets: the
// rule creates no vertex, so both passes find the same number of them.
inline bool svSeedNominalM2B() { static const bool on = [] { const char *e = std::getenv("SV_SEED_M2B"); return e && TString(e) == "nominal"; }(); return on; }
inline double nominalVertexMass(tTree& t, Int_t& ijet, Long64_t& ient, double fallback) {
  double a = 0., b = 0., c = 0., d = 0.;   // this pass does not touch the failure counters
  const vector<ROOT::Math::PtEtaPhiMVector> v = makeSvtxs_withBDT(t, ijet, ient, a, b, c, d, nullptr, nullptr, nullptr, nullptr, false);
  return v.size() == 2 ? v[0].M() + v[1].M() : fallback;
}

// SV_SEED_DIAG=on (2026-11-06, with SV_SEED=on): the data/MC check of the SV seed rule. Same jets and same
// weight as the templates, in a fine jet pT axis so the reader can apply the fit's jet-pT weights. Per
// flavour (data | 0b | b | bb), names h_diag_<var>_<flavour>:
//   per jet      what the response sees: seeded vs plain-BDT vertices (A, dR, pT sum, pT1*pT2, m2B)
//   rescued trk  tracks kept only because of the rule. NOT the same set in data and MC when the track-BDT SF
//                is on: MC's effective WPs differ from data's 0.365, so read them with the next group
//   sv trk       ALL kept tracks of the seeded vertices: independent of that WP difference once a vertex is seeded
struct SvSeedDiag {
    bool on = false;
    std::map<TString, TH2D *> h;
    static bool requested() { const char *e = std::getenv("SV_SEED_DIAG"); return e && TString(e) == "on"; }
    void book(bool isMC) {
        on = requested();
        if (!on) return;
        struct V { const char *n; int nb; double lo, hi; };
        const std::vector<V> vars = {
            {"nresc", 10, 0, 10}, {"ptresc", 40, 0, 40}, {"share", 20, 0, 1}, {"ptratio", 42, 0.9, 3.0}, {"dA", 50, -0.5, 0.5},
            {"ddR", 40, -0.2, 0.2}, {"eecratio", 42, 0.8, 5.0}, {"m2Bnom", 10, 0, 10}, {"m2Bseed", 10, 0, 10},
            {"Anom", 5, 0.5, 1.0}, {"Aseed", 5, 0.5, 1.0}, {"drnom", 9, 0, 0.45}, {"drseed", 9, 0, 0.45},
            {"nsvtrk", 10, 1, 11}, {"svshare", 20, 0, 1}, {"svmass", 24, 0, 6},
            {"resc_pt", 40, 1, 21}, {"resc_bdt", 40, -1, 1}, {"resc_ip", 60, -20, 100}, {"resc_dr", 40, 0, 0.4}, {"resc_z", 30, 0, 0.3},
            {"svtrk_pt", 40, 1, 21}, {"svtrk_ip", 60, -20, 100}, {"svtrk_dr", 40, 0, 0.4}, {"svtrk_bdt", 40, -1, 1} };
        const int nfine = 64; double e[nfine + 1];
        for (int i = 0; i <= nfine; ++i) e[i] = 80. + 5. * i;
        std::vector<TString> fl; if (isMC) fl = {"0b", "b", "bb"}; else fl = {"data"};
        for (const V &v : vars) for (const TString &f : fl) {
            std::vector<double> x(v.nb + 1); for (int i = 0; i <= v.nb; ++i) x[i] = v.lo + (v.hi - v.lo) * i / v.nb;
            TH2D *o = new TH2D(TString::Format("h_diag_%s_%s", v.n, f.Data()), TString::Format(";%s;p_{T} [GeV]", v.n), v.nb, x.data(), nfine, e);
            o->Sumw2(); h[TString::Format("%s_%s", v.n, f.Data())] = o;
        }
    }
    void put(const char *var, const TString &f, double x, double pt, double w) {
        TH2D *o = h.at(TString::Format("%s_%s", var, f.Data()));
        const TAxis *a = o->GetXaxis();
        o->Fill(std::max(a->GetXmin() + 1e-9, std::min(a->GetXmax() - 1e-9, x)), pt, w);   // over/underflow in the edge bins
    }
    void fill(tTree &t, Int_t &ijet, Long64_t &ient, bool isMC, Int_t nb, const ROOT::Math::PtEtaPhiMVector &s0,
              const ROOT::Math::PtEtaPhiMVector &s1, double pt, double w) {
        if (!on || pt >= 400.) return;
        const TString f = !isMC ? "data" : nb == 0 ? "0b" : nb == 1 ? "b" : "bb";
        double a = 0., b = 0., c = 0., d = 0.;
        const vector<ROOT::Math::PtEtaPhiMVector> n = makeSvtxs_withBDT(t, ijet, ient, a, b, c, d, nullptr, nullptr, nullptr, nullptr, false);
        if (n.size() != 2) return;
        std::set<Int_t> seeds;
        for (Int_t k = 0; k < t.ntrk; ++k)
            if (t.trkJetId[k] == ijet && t.trkPt[k] >= 1 && t.trkSvtxId[k] >= 0 && TrackBdtSF::passes(t.trkBdtScore[k], t.trkMatchSta[k], t.trkSvtxId[k]))
                seeds.insert(t.trkSvtxId[k]);
        std::map<Int_t, ROOT::Math::PtEtaPhiMVector> axis; std::map<Int_t, int> nsv; std::map<Int_t, double> ptsv, rescsv;
        double ptkept = 0., ptresc = 0.; int nresc = 0;
        for (Int_t k = 0; k < t.ntrk; ++k) {
            if (t.trkJetId[k] != ijet || t.trkPt[k] < 1) continue;
            const Int_t sv = t.trkSvtxId[k];
            const bool pass = TrackBdtSF::passes(t.trkBdtScore[k], t.trkMatchSta[k], sv), seeded = sv >= 0 && seeds.count(sv);
            if (pass || seeded) ptkept += t.trkPt[k];
            if (seeded) { axis[sv] += recoTrackP4(t, k); ++nsv[sv]; ptsv[sv] += t.trkPt[k]; }
            if (seeded && !pass) { ptresc += t.trkPt[k]; ++nresc; rescsv[sv] += t.trkPt[k]; }
        }
        for (Int_t k = 0; k < t.ntrk; ++k) {
            if (t.trkJetId[k] != ijet || t.trkPt[k] < 1) continue;
            const Int_t sv = t.trkSvtxId[k];
            if (sv < 0 || !seeds.count(sv)) continue;
            const ROOT::Math::PtEtaPhiMVector tr = recoTrackP4(t, k);
            const double ip = std::isfinite(t.trkIp3dSig[k]) ? t.trkIp3dSig[k] : 0.;
            const double dra = t.calc_dr(tr.Eta(), tr.Phi(), axis[sv].Eta(), axis[sv].Phi());
            put("svtrk_pt", f, t.trkPt[k], pt, w); put("svtrk_ip", f, ip, pt, w); put("svtrk_dr", f, dra, pt, w); put("svtrk_bdt", f, t.trkBdtScore[k], pt, w);
            if (!TrackBdtSF::passes(t.trkBdtScore[k], t.trkMatchSta[k], sv)) {
                put("resc_pt", f, t.trkPt[k], pt, w); put("resc_bdt", f, t.trkBdtScore[k], pt, w); put("resc_ip", f, ip, pt, w);
                put("resc_dr", f, dra, pt, w); put("resc_z", f, t.trkPt[k] / pt, pt, w);
            }
        }
        for (const auto &kv : nsv) {
            put("nsvtrk", f, kv.second, pt, w); put("svshare", f, rescsv[kv.first] / ptsv[kv.first], pt, w); put("svmass", f, axis[kv.first].M(), pt, w);
        }
        const double ns = s0.Pt() + s1.Pt(), nn = n[0].Pt() + n[1].Pt();
        const double As = std::max(s0.Pt(), s1.Pt()) / ns, An = std::max(n[0].Pt(), n[1].Pt()) / nn;
        const double ds = t.calc_dr(s0.Eta(), s0.Phi(), s1.Eta(), s1.Phi()), dn = t.calc_dr(n[0].Eta(), n[0].Phi(), n[1].Eta(), n[1].Phi());
        put("nresc", f, nresc, pt, w); put("ptresc", f, ptresc, pt, w); put("share", f, ptkept > 0. ? ptresc / ptkept : 0., pt, w);
        put("ptratio", f, ns / nn, pt, w); put("dA", f, As - An, pt, w); put("ddR", f, ds - dn, pt, w);
        put("eecratio", f, (s0.Pt() * s1.Pt()) / (n[0].Pt() * n[1].Pt()), pt, w);
        put("m2Bnom", f, n[0].M() + n[1].M(), pt, w); put("m2Bseed", f, s0.M() + s1.M(), pt, w);
        put("Anom", f, An, pt, w); put("Aseed", f, As, pt, w); put("drnom", f, dn, pt, w); put("drseed", f, ds, pt, w);
    }
    void write() const { for (const auto &kv : h) kv.second->Write(); }
};

//fill_jk_resampling_response to be added
// Jackknife resampling: fill all but one of the 10 objects; num in [0,1) picks the one left
// out. Works for TH1D/TH2D/TH3D (x..., w) and RooUnfoldResponse (reco..., gen..., w) alike --
// it replaces the six fill_jk_resampling{,_response}_{1D,2D,3D} functions. Not called at the
// moment.
template <class H, class... X>
void fill_jk_resampling(const std::vector<H *> &objs, double num, X... args) {
    int numRescaled = (int) (num*10); // from 0 to 9
    for (int i=0; i<10; i++) {
        if (i==numRescaled) continue; // fill all but 1
        objs[i]->Fill(args...);
    }
}

struct AggBHadronNtupleRow {
  Long64_t entry;
  // Same widths as the input branches (see tTree.h): evt does not fit in 32 bits.
  UInt_t run, lumi;
  ULong64_t evt;
  Int_t jetIndex;
  Float_t weight, jtpt, jteta, refpt, refeta;
  Int_t jtNbHad;
  Float_t btagScore;
  Int_t passRecoKin, passGenKin, passBtag, passJetVeto, nRecoAgg, nGenAgg, genStatus1, genStatus2;
  Float_t recoPt1, recoEta1, recoPhi1, recoM1, recoPt2, recoEta2, recoPhi2, recoM2, recoDr, recoMB, recoEec, genPt1, genEta1, genPhi1, genM1, genPt2, genEta2, genPhi2, genM2, genDr, genMB, genEec;
  Float_t tsPt1, tsEta1, tsPhi1, tsPt2, tsEta2, tsPhi2, orPt1, orEta1, orPhi1, orPt2, orEta2, orPhi2;
  Float_t jtphi;
  Float_t cxPt1, cxEta1, cxPhi1, cxPt2, cxEta2, cxPhi2; Int_t cxStatus;
  std::vector<float> trk_pt, trk_eta, trk_phi, trk_ip3d, trk_ip3dsig, trk_dist, trk_bdt,
                     trk_sv_m, trk_sv_mcorr, trk_sv_pt, trk_sv_dl, trk_sv_dls, trk_sv_dl2d, trk_sv_dls2d, trk_sv_chi2;
  std::vector<int> trk_pdg, trk_svidx, trk_sv_ntrk, trk_nomsel, truth_sta, truth_label;

  void reset() {
    entry = -1;
    run = lumi = 0;   // unsigned: 0 is the "unset" marker, not -1
    evt = 0;
    jetIndex = -1;
    weight = 1.0;
    jtpt = jteta = refpt = refeta = -999.0;
    jtNbHad = -1;
    btagScore = -999.0;
    passRecoKin = passGenKin = passBtag = passJetVeto = 0;
    nRecoAgg = nGenAgg = 0;
    genStatus1 = genStatus2 = 0;
    recoPt1 = recoEta1 = recoPhi1 = recoM1 = -999.0;
    recoPt2 = recoEta2 = recoPhi2 = recoM2 = -999.0;
    recoDr = recoMB = recoEec = -999.0;
    genPt1 = genEta1 = genPhi1 = genM1 = -999.0;
    genPt2 = genEta2 = genPhi2 = genM2 = -999.0;
    genDr = genMB = genEec = -999.0;
    tsPt1 = tsEta1 = tsPhi1 = tsPt2 = tsEta2 = tsPhi2 = -999.0;
    orPt1 = orEta1 = orPhi1 = orPt2 = orEta2 = orPhi2 = -999.0;
    jtphi = -999.0;
    cxPt1 = cxEta1 = cxPhi1 = cxPt2 = cxEta2 = cxPhi2 = -999.0; cxStatus = -1;
    for (auto *v : {&trk_pt, &trk_eta, &trk_phi, &trk_ip3d, &trk_ip3dsig, &trk_dist, &trk_bdt, &trk_sv_m, &trk_sv_mcorr,
                    &trk_sv_pt, &trk_sv_dl, &trk_sv_dls, &trk_sv_dl2d, &trk_sv_dls2d, &trk_sv_chi2}) v->clear();
    for (auto *v : {&trk_pdg, &trk_svidx, &trk_sv_ntrk, &trk_nomsel, &truth_sta, &truth_label}) v->clear();
  }

  // Kinematics of the two aggregated B hadrons: the reco SVs, or the chosen gen pair.
  void setReco(const ROOT::Math::PtEtaPhiMVector &a, const ROOT::Math::PtEtaPhiMVector &b) {
    recoPt1 = a.Pt(); recoEta1 = a.Eta(); recoPhi1 = a.Phi(); recoM1 = a.M();
    recoPt2 = b.Pt(); recoEta2 = b.Eta(); recoPhi2 = b.Phi(); recoM2 = b.M();
  }
  void setGen(const ROOT::Math::PtEtaPhiMVector &a, const ROOT::Math::PtEtaPhiMVector &b) {
    genPt1 = a.Pt(); genEta1 = a.Eta(); genPhi1 = a.Phi(); genM1 = a.M();
    genPt2 = b.Pt(); genEta2 = b.Eta(); genPhi2 = b.Phi(); genM2 = b.M();
  }
};

void makeAggBHadronBranches(TTree* tree, AggBHadronNtupleRow& row) {
  // Leaf types: /L Long64_t, /i UInt_t, /l ULong64_t, /I Int_t, /F Float_t. Branch order = file order.
  auto br = [tree](const char *name, void *addr, const char *type) { tree->Branch(name, addr, Form("%s/%s", name, type)); };
  br("entry", &row.entry, "L"); br("run", &row.run, "i"); br("lumi", &row.lumi, "i"); br("evt", &row.evt, "l");
  br("jetIndex", &row.jetIndex, "I");
  br("weight", &row.weight, "F"); br("jtpt", &row.jtpt, "F"); br("jteta", &row.jteta, "F");
  br("refpt", &row.refpt, "F"); br("refeta", &row.refeta, "F");
  br("jtNbHad", &row.jtNbHad, "I"); br("btagScore", &row.btagScore, "F");
  br("passRecoKin", &row.passRecoKin, "I"); br("passGenKin", &row.passGenKin, "I"); br("passBtag", &row.passBtag, "I");
  br("passJetVeto", &row.passJetVeto, "I");
  br("nRecoAgg", &row.nRecoAgg, "I"); br("nGenAgg", &row.nGenAgg, "I");
  br("genStatus1", &row.genStatus1, "I"); br("genStatus2", &row.genStatus2, "I");
  br("recoPt1", &row.recoPt1, "F"); br("recoEta1", &row.recoEta1, "F"); br("recoPhi1", &row.recoPhi1, "F"); br("recoM1", &row.recoM1, "F");
  br("recoPt2", &row.recoPt2, "F"); br("recoEta2", &row.recoEta2, "F"); br("recoPhi2", &row.recoPhi2, "F"); br("recoM2", &row.recoM2, "F");
  br("recoDr", &row.recoDr, "F"); br("recoMB", &row.recoMB, "F"); br("recoEec", &row.recoEec, "F");
  br("genPt1", &row.genPt1, "F"); br("genEta1", &row.genEta1, "F"); br("genPhi1", &row.genPhi1, "F"); br("genM1", &row.genM1, "F");
  br("genPt2", &row.genPt2, "F"); br("genEta2", &row.genEta2, "F"); br("genPhi2", &row.genPhi2, "F"); br("genM2", &row.genM2, "F");
  br("genDr", &row.genDr, "F"); br("genMB", &row.genMB, "F"); br("genEec", &row.genEec, "F");
  br("tsPt1", &row.tsPt1, "F"); br("tsEta1", &row.tsEta1, "F"); br("tsPhi1", &row.tsPhi1, "F");
  br("tsPt2", &row.tsPt2, "F"); br("tsEta2", &row.tsEta2, "F"); br("tsPhi2", &row.tsPhi2, "F");
  br("orPt1", &row.orPt1, "F"); br("orEta1", &row.orEta1, "F"); br("orPhi1", &row.orPhi1, "F");
  br("orPt2", &row.orPt2, "F"); br("orEta2", &row.orEta2, "F"); br("orPhi2", &row.orPhi2, "F");
  br("jtphi", &row.jtphi, "F");
  br("cxPt1", &row.cxPt1, "F"); br("cxEta1", &row.cxEta1, "F"); br("cxPhi1", &row.cxPhi1, "F");
  br("cxPt2", &row.cxPt2, "F"); br("cxEta2", &row.cxEta2, "F"); br("cxPhi2", &row.cxPhi2, "F");
  br("cxStatus", &row.cxStatus, "I");
  // Per-track arrays (std::vector branches; uproot reads them as jagged arrays).
  tree->Branch("trk_pt", &row.trk_pt); tree->Branch("trk_eta", &row.trk_eta); tree->Branch("trk_phi", &row.trk_phi);
  tree->Branch("trk_ip3d", &row.trk_ip3d); tree->Branch("trk_ip3dsig", &row.trk_ip3dsig); tree->Branch("trk_dist", &row.trk_dist);
  tree->Branch("trk_bdt", &row.trk_bdt); tree->Branch("trk_pdg", &row.trk_pdg); tree->Branch("trk_nomsel", &row.trk_nomsel);
  tree->Branch("trk_svidx", &row.trk_svidx); tree->Branch("trk_sv_ntrk", &row.trk_sv_ntrk);
  tree->Branch("trk_sv_m", &row.trk_sv_m); tree->Branch("trk_sv_mcorr", &row.trk_sv_mcorr); tree->Branch("trk_sv_pt", &row.trk_sv_pt);
  tree->Branch("trk_sv_dl", &row.trk_sv_dl); tree->Branch("trk_sv_dls", &row.trk_sv_dls);
  tree->Branch("trk_sv_dl2d", &row.trk_sv_dl2d); tree->Branch("trk_sv_dls2d", &row.trk_sv_dls2d); tree->Branch("trk_sv_chi2", &row.trk_sv_chi2);
  tree->Branch("truth_sta", &row.truth_sta); tree->Branch("truth_label", &row.truth_label);
}

#if 0
// Build templates for the template fit.
// MC:   fills h3D_b (jtNbHad==1) and h3D_bb (jtNbHad==2) using reco-level cuts and reco SV reconstruction.
// Data: fills h3D_data with the same reco logic — no truth classification.
void make_templates(const AnalysisConfig& cfg, Long64_t ev_first = 0, Long64_t ev_last = -1, Int_t job_idx = -1) {

  bool isMC = cfg.dataset.isMC;

  tTree t;
  t.Init(cfg.dataset.filename, cfg.dataset.isMC, cfg.dataset.RunN);
  t.SetBranchStatus("*", 0);
  auto active_branches = getActiveBranches(cfg);
  t.SetBranchStatus(active_branches, 1);

  double prescale = cfg.dataset.data_prescale;

  //-- Variables
  double agg_fail = 0, nb_sv = 0, sv_fail = 0, merge_fail = 0;

  // MC: separate 0b, b and bb templates
  TH3D *h3D_0b = new TH3D("h3D_0b", ";m_{2B} [GeV];#DeltaR;p_{T} [GeV]", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h3D_b  = new TH3D("h3D_b",  ";m_{2B} [GeV];#DeltaR;p_{T} [GeV]", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h3D_bb = new TH3D("h3D_bb", ";m_{2B} [GeV];#DeltaR;p_{T} [GeV]", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  // Data: single distribution to be fit
  TH3D *h3D_data = new TH3D("h3D_data", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);

  h3D_0b->Sumw2();   h3D_0b->SetCanExtend(TH1::kNoAxis);
  h3D_b->Sumw2();    h3D_b->SetCanExtend(TH1::kNoAxis);
  h3D_bb->Sumw2();   h3D_bb->SetCanExtend(TH1::kNoAxis);
  h3D_data->Sumw2(); h3D_data->SetCanExtend(TH1::kNoAxis);

  // Jet counts (no EEC weight): 3D (mB, dr, jtpt) — same axes as the EEC histograms
  TH3D *h_count_0b   = new TH3D("h_count_0b",   "jet counts 0b;m_{2B} [GeV];#DeltaR;p_{T} [GeV]",   bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h_count_b    = new TH3D("h_count_b",    "jet counts 1b;m_{2B} [GeV];#DeltaR;p_{T} [GeV]",   bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h_count_bb   = new TH3D("h_count_bb",   "jet counts 2b;m_{2B} [GeV];#DeltaR;p_{T} [GeV]",   bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  TH3D *h_count_data = new TH3D("h_count_data", "jet counts data;m_{2B} [GeV];#DeltaR;p_{T} [GeV]", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  h_count_0b->Sumw2();   h_count_0b->SetCanExtend(TH1::kNoAxis);
  h_count_b->Sumw2();    h_count_b->SetCanExtend(TH1::kNoAxis);
  h_count_bb->Sumw2();   h_count_bb->SetCanExtend(TH1::kNoAxis);
  h_count_data->Sumw2(); h_count_data->SetCanExtend(TH1::kNoAxis);

  Long64_t n_events = t.GetEntries();
  if (ev_last < 0 || ev_last > n_events) ev_last = n_events;
  if (ev_first < 0) ev_first = 0;
  std::cout << "Processing events [" << ev_first << ", " << ev_last << ") of " << n_events << std::endl;

  for (Long64_t ient = ev_first; ient < ev_last; ient++) {
    if (ient % 50000 == 0)
      std::cout << "\rProcessing: " << 100.0 * ient / n_events << " %" << std::flush;
      //std::cout << "\rProcessing: " << ient  <<std::endl;

      t.GetEntry(ient);
      TrkEffSyst::beginEvent(t.run, t.lumi, t.evt);
      // -- test tree branches are read:
      //cout << "jtpt = " << t.jtpt[0] << endl;

      double weight_tree = cfg.dataset.isMC ? t.weight : 1.0;

    // -- Trigger selections
    if (! passEventSelection(t, cfg)) continue;
    // cout << "event pass ok" << endl;

    for (Int_t ijet = 0; ijet < t.nref; ijet++) {

      // -- jet knematiocs + skip event of large weight based on jet pt
      if(! passRecoJetKinematics(t, ijet, cfg)) continue;
          // cout << "Jet kin pass ok" << endl;

	  	// batgging
      if (! passBtag(t, ijet, cfg)) continue;
          // cout << "Jet btag pass ok" << endl;

      // reco SV reconstruction — same for data and MC
      vector<ROOT::Math::PtEtaPhiMVector> reco_sv = makeSvtxs_withBDT(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail, nullptr, nullptr);
      if (reco_sv.size() != 2) continue; // #sv must = 2 (default of makeSvtxs_withBDT())

      double dr   = t.calc_dr(reco_sv[0].Eta(), reco_sv[0].Phi(), reco_sv[1].Eta(), reco_sv[1].Phi());
      double pt1  = reco_sv[0].Pt();
      double pt2  = reco_sv[1].Pt();
      double eec  = EecWeight::value(pt1, pt2, cfg.n);
      double jtpt = jtpt_fill(t.jtpt[ijet]);
      double mB   = reco_sv[0].M() + reco_sv[1].M();

      //Fix the under/overflow
      //if(dr < dr_min) dr = dr_min_fill;
      if(dr >= dr_max) dr = dr_max_fill;
      if(mB >= mb_max) mB = mb_max_fill;

      //std::cout << "weight: " << weight_tree << std::endl;
      //std::cout << "eec: " << eec << std::endl;

      if (cfg.dataset.RunN == 2 && !cfg.dataset.isMC && t.HLT_HIAK4PFJet40_v1 && !(t.HLT_HIAK4PFJet60_v1 || t.HLT_HIAK4PFJet80_v1 || t.HLT_HIAK4PFJet100_v1))
      {eec *= prescale;}

      if (cfg.dataset.RunN == 3 && !cfg.dataset.isMC && t.HLT_AK4PFJet60_v8 && !(t.HLT_AK4PFJet80_v8 || t.HLT_AK4PFJet100_v8 || t.HLT_AK4PFJet120_v8) )
      {eec *= prescale;}

      if (cfg.dataset.isMC && dr > 0.005) {
        // use truth to classify: fill separate 0b, b and bb templates
        if      (t.jtNbHad[ijet] == 0) { h3D_0b->Fill(mB, dr, jtpt, eec * weight_tree); h_count_0b->Fill(mB, dr, jtpt, weight_tree); }
        else if (t.jtNbHad[ijet] == 1) { h3D_b ->Fill(mB, dr, jtpt, eec * weight_tree); h_count_b ->Fill(mB, dr, jtpt, weight_tree); }
        else if (t.jtNbHad[ijet] >= 2) { h3D_bb->Fill(mB, dr, jtpt, eec * weight_tree); h_count_bb->Fill(mB, dr, jtpt, weight_tree); }
      } else {
        h3D_data->Fill(mB, dr, jtpt, eec * weight_tree);
        h_count_data->Fill(mB, dr, jtpt, weight_tree);
      }
    }
  }
  std::cout << std::endl;

  TString label = cfg.physics.useBtag ? "_btag" : "_nobtag";
  TString job_suffix = (job_idx >= 0) ? Form("_job%d", job_idx) : "";
  TFile outFile((cfg.dataset.output_folder + cfg.dataset.output_hist + label + kinLabel(cfg) + job_suffix + cfg.dataset.domain).Data(), "RECREATE");
  if (isMC) {
    h3D_0b->Write();
    h3D_b->Write();
    h3D_bb->Write();
    h_count_0b->Write();
    h_count_b->Write();
    h_count_bb->Write();
  } else {
    h3D_data->Write();
    h_count_data->Write();
  }
  outFile.Close();
}
#endif 
#if 0
// Creates a 3D response matrix in (mB, dr_SV, jtpt) for unfolding the bb template
// fit distributions. One entry per jet (not per track pair).
//
// Reco level : makeSvtxs_withBDT -> two reconstructed SVs
// Gen level  : PartialBsAggregation -> all gen b hadrons from refTrk* branches;
//              the pair with the LARGEST EEC weight (pt_i * pt_j)^n is chosen.
//
// Purity   = (reco_pass && gen_pass) / reco_pass  [binned in reco observables]
// Efficiency = (reco_pass && gen_pass) / gen_pass [binned in gen observables]
// reco_pass uses reco jet kinematics (jtpt, jteta, btag, reco SVs, reco mB and dr).
// gen_pass  uses gen  jet kinematics (refpt, refeta) and gen observable range.
//
// No fakes: the jet always has a real gen b-hadron pair — only migration corrections needed.
void create_response_templatefit(
    const AnalysisConfig& cfg,
    TString  output_hist_response, // is it needed to be different from that of template fits?
    Long64_t ev_first = 0,
    Long64_t ev_last  = -1)
{

  bool isMC = cfg.dataset.isMC;

  // -- Only for MC------------
  if(!isMC) { // For response matrix
    std::cerr<< "ERROR: Input sample is not MC!, check sample or MC tag."<<endl;
    return;
  }

    TString fout_name = cfg.dataset.output_folder + output_hist_response + (cfg.physics.useBtag ? "_btag" : "_nobtag") + kinLabel(cfg) + ".root"; // -- can be changed to be in the same file as templates?(but with Response keyword if they are produced seperately)

    // -------------------------------------------------------------------
    // -- Change this to be used from the default header ?
    // Re-derive bin counts from const sizes (Cling init-order workaround for non-const globals)
    Int_t n_mb = mb_binsVectorSize - 1;   // same as mb_bins in binning_histos_small.h
    Int_t n_dr = dr_binsVectorSize - 1;   // same as dr_bins
    Int_t n_pt = jtpt_binsVectorSize - 1; // same as jtpt_bins
    Double_t mb_min = mb_binsVector[0];
    Double_t mb_max = mb_binsVector[n_mb];
    Double_t mb_max_fill = 9.9;
    Double_t dr_min = dr_binsVector[0];
    Double_t dr_max = dr_binsVector[n_dr];
    std::cout << "Histogram binning: mB=" << n_mb << " dr=" << n_dr << " jtpt=" << n_pt << std::endl;
    std::cout << "mb range: [" << mb_min << ", " << mb_max << "], fill cap=" << mb_max_fill << std::endl;
    std::cout << "dr range: [" << dr_min << ", " << dr_max << "]" << std::endl;
    // -------------------------------------------------------------------

    tTree t;
    t.Init(cfg.dataset.filename, cfg.dataset.isMC, cfg.dataset.RunN);
      t.SetBranchStatus("*", 0);
      auto active_branches = getActiveBranches(cfg);
      t.SetBranchStatus(active_branches, 1);

    std::mt19937 generator(12345u);
    std::uniform_real_distribution<double> distr(0., 1.);

    // -- needed parameters for reco SV func()
    double agg_fail = 0, nb_sv = 0, sv_fail = 0, merge_fail = 0;

    //-- debugging variables
    long n_bb_jets = 0, n_reco_pass = 0, n_gen_pass = 0, n_both_pass = 0;
    long n_gen_bh_ok = 0;  // jets passing gen_bh.size() >= 2

    // per-condition failure counters for reco_pass
    long n_fail_reco_sv   = 0, n_fail_reco_pt  = 0, n_fail_reco_eta = 0;
    long n_fail_reco_btag = 0, n_fail_reco_mb  = 0, n_fail_reco_dr  = 0;
    // per-condition failure counters for gen_pass
    long n_fail_gen_pt  = 0, n_fail_gen_eta = 0;
    long n_fail_gen_mb  = 0, n_fail_gen_dr  = 0;
    int  n_debug_printed = 0;

    //purity
    TH2D *h_half0_purity_num = new TH2D("h_half0_purity_numerator_tf",   "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half0_purity_den = new TH2D("h_half0_purity_denominator_tf", "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half1_purity_num = new TH2D("h_half1_purity_numerator_tf",   "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half1_purity_den = new TH2D("h_half1_purity_denominator_tf", "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    //efficiency
    TH2D *h_half0_eff_num = new TH2D("h_half0_efficiency_numerator_tf",   "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half0_eff_den = new TH2D("h_half0_efficiency_denominator_tf", "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half1_eff_num = new TH2D("h_half1_efficiency_numerator_tf",   "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half1_eff_den = new TH2D("h_half1_efficiency_denominator_tf", "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    //btag efficiency
    TH2D *h_half0_btag_eff_num = new TH2D("h_half0_btag_efficiency_numerator_tf",   "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half0_btag_eff_den = new TH2D("h_half0_btag_efficiency_denominator_tf", "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half1_btag_eff_num = new TH2D("h_half1_btag_efficiency_numerator_tf",   "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);
    TH2D *h_half1_btag_eff_den = new TH2D("h_half1_btag_efficiency_denominator_tf", "x=dr_SV, y=jtpt", n_dr, dr_binsVector, n_pt, jtpt_binsVector);

    RooUnfoldResponse *response_half0 = new RooUnfoldResponse(h_half0_purity_den, h_half0_eff_den, "response_tf_half0", "tf response half0");
    RooUnfoldResponse *response_half1 = new RooUnfoldResponse(h_half1_purity_den, h_half1_eff_den, "response_tf_half1", "tf response half1");
    RooUnfoldResponse *response_full  = new RooUnfoldResponse(h_half0_purity_den, h_half0_eff_den, "response_tf_full",  "tf response full");

    Long64_t n_events = t.GetEntries();
    if (ev_last < 0 || ev_last > n_events) ev_last = n_events;
    std::cout << "Processing events [" << ev_first << ", " << ev_last << ") of " << n_events << std::endl;

    for (Long64_t ient = ev_first; ient < ev_last; ient++) {
        if (ient % 50000 == 0) std::cout << "\rProcessing: " << 100.0 * ient / n_events << " %" << std::flush;
        t.GetEntry(ient);

        double weight_tree = cfg.dataset.isMC ? t.weight : 1.0;

        // MC trigger selection (run-dependent)
        if (! passEventSelection(t, cfg)) continue;

        for (Int_t ijet = 0; ijet < t.nref; ijet++) {

            if (t.jtNbHad[ijet] < 2) continue; // select jets of 2b (truth)
            n_bb_jets++;

            // step1: for Response matrix ---- Gen b hadrons ----
            std::vector<ROOT::Math::PtEtaPhiMVector> gen_bh;
            std::vector<Int_t> gen_bh_sta;
            PartialBsAggregation(gen_bh, gen_bh_sta, t, ijet);
            if (gen_bh.size() < 2) continue;
            n_gen_bh_ok++;

            // From aggregated gen B: Pick gen pair with largest EEC weight (pt_i * pt_j)^n
            int best_i = 0, best_j = 1;
            double best_pt_prod = -1;
            for (size_t gi = 0; gi < gen_bh.size(); gi++)
                for (size_t gj = gi+1; gj < gen_bh.size(); gj++) {
                    double pp = gen_bh[gi].Pt() * gen_bh[gj].Pt();
                    if (pp > best_pt_prod) { best_pt_prod = pp; best_i = gi; best_j = gj; }
            }

            double eec_gen = EecWeight::value(gen_bh[best_i].Pt(), gen_bh[best_j].Pt(), cfg.n);
            double mB_gen  = gen_bh[best_i].M() + gen_bh[best_j].M();
            double dr_gen  = t.calc_dr(gen_bh[best_i].Eta(), gen_bh[best_i].Phi(),
                                       gen_bh[best_j].Eta(), gen_bh[best_j].Phi());

            // -- common variables repeatdly used in fill histograms
            // jtpt_fill() folds jets above the last bin edge into that bin, so they are
            // not lost to the overflow bin (see binning_histos_small.h).
            double jpt_reco = jtpt_fill(reco_jet_pt(t, ijet));
            // double jpt_gen = jtpt_fill(gen_jet_pt(t, ijet));   // until 2026-10-29
            double jpt_gen = jtptGenFill(gen_jet_pt(t, ijet));   // = jtpt_fill unless JET_PT_MAX (gen overflow)
            double jeta_reco = reco_jet_eta(t, ijet);
            double jeta_gen = gen_jet_eta(t, ijet);
              // -- for cout only
            double btagVar =  (cfg.dataset.RunN == 2) ? (t.discr_particleNet_BvsAll[ijet]) :
                              ( (cfg.dataset.RunN == 3) ?
                                  (t.discr_unifiedParticleTransformer_probb[ijet] +
                                  t.discr_unifiedParticleTransformer_problepb[ijet] +
                                  t.discr_unifiedParticleTransformer_probbb[ijet]) :
                                  -1);

            // step2: for Response matrix ---- Reco SVs ----
            vector<ROOT::Math::PtEtaPhiMVector> reco_sv =
                makeSvtxs_withBDT(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail, nullptr, nullptr);

            double mB_reco = -1, dr_reco = -1, eec_reco = -1;
            bool reco_sv_ok = (reco_sv.size() == 2);// must equal 2
             if (reco_sv_ok) {
                mB_reco  = reco_sv[0].M() + reco_sv[1].M();
                dr_reco  = t.calc_dr(reco_sv[0].Eta(), reco_sv[0].Phi(),
                                     reco_sv[1].Eta(), reco_sv[1].Phi());
                eec_reco = std::pow(reco_sv[0].Pt() * reco_sv[1].Pt(), cfg.n);
            }

            // Overflow protection
            double mB_reco_fill = mB_reco, dr_reco_fill = dr_reco;
            if (reco_sv_ok) {
                if (mB_reco_fill >= mb_max) mB_reco_fill = mb_max_fill;
                if (dr_reco_fill >= dr_max) dr_reco_fill = dr_max_fill;
                //if (dr_reco_fill  < dr_min) dr_reco_fill = dr_min_fill; // underflow use or not ?
            }
            double mB_gen_fill = mB_gen, dr_gen_fill = dr_gen;
              if (mB_gen_fill >= mb_max) mB_gen_fill = mb_max_fill;
              if (dr_gen_fill >= dr_max)  dr_gen_fill = dr_max_fill;
              //if (dr_gen_fill  < dr_min)  dr_gen_fill = dr_min_fill;  // underflow use or not ?

            // reco_pass: full detector-level selection
            //bool reco_pass = reco_sv_ok &&
                             //passRecoJetKinematics(t, ijet, cfg) &&
                             //passBtag(t, ijet, cfg);
                             //  && (mB_reco_fill >= mb_min && mB_reco_fill < mb_max) && (dr_reco_fill < dr_max); // Not needed since reco_sv_ok is already required.

            bool reco_pass = reco_sv_ok && passRecoJetKinematics(t, ijet, cfg) && dr_reco_fill > 0.005;

            bool reco_pass_btag = reco_sv_ok && passBtag(t, ijet, cfg); //check

            /*
            // -- Debugging paragraph-------
            //std::cout << "reco_sv_ok: " << reco_sv_ok << std::endl;
            //if (reco_sv_ok) {std::cout << "jpt_reco: " << jpt_reco << std::endl;}
            //if (reco_sv_ok && (jpt_reco >= pT_low && jpt_reco < pT_high)) {
                         //std::cout << "jteta_reco: " << t.jteta[ijet] << std::endl; }
            //if (reco_sv_ok && (jpt_reco >= pT_low && jpt_reco < pT_high) && (std::abs(t.jteta[ijet]) < 1.6)) {
             //std::cout << "btag: " << (t.discr_particleNet_BvsAll[ijet] ) << std::endl;}

            //if (reco_sv_ok && (jpt_reco >= pT_low && jpt_reco < pT_high) && (std::abs(t.jteta[ijet]) < 1.6) && (!btag || t.discr_particleNet_BvsAll[ijet] > 0.898)){
             //std::cout << "mB_reco: " << mb_min << " " << mB_reco_fill << " " << mb_max << std::endl;
            //}
            //if (reco_sv_ok && (jpt_reco >= pT_low && jpt_reco < pT_high) && (std::abs(t.jteta[ijet]) < 1.6) && (!btag || t.discr_particleNet_BvsAll[ijet] > 0.898 ) && (mB_reco_fill >= mb_min && mB_reco_fill < mb_max)){
                //std::cout << "dr_reco: " << dr_reco << "MAX" <<dr_max << std::endl;}
            */

            // gen_pass: particle-level jet kinematics + gen observable range
            bool gen_pass  = passGenJetKinematics(t, ijet, cfg); //  && (mB_gen_fill >= mb_min && mB_gen_fill < mb_max) &&(dr_gen_fill < dr_max); //dr_gen_fill >= dr_min &&

            bool gen_pass_btag = passBtag(t, ijet, cfg); //check
            // --- debugging paragraph ----
            //std::cout << "gen_pass: " << gen_pass << std::endl;
            //std::cout << "jpt_gen: " << jpt_gen << std::endl;
            //std::cout << "pT_low: " << pT_low << ", pT_high: " << pT_high << std::endl;
            //std::cout << "refeta_gen: " << t.refeta[ijet] << std::endl;
            //std::cout << "mB_gen: " << mB_gen << std::endl;
            //std::cout << "mB_low: " << mb_min << ", mB_high: " << mb_max << std::endl;
            //std::cout << "dr_gen: " << dr_gen << std::endl;
            //---------------------------

            // -- For debugging only (not needed to fill histograms) ---
            // --- per-condition failure tallies ---
            if (reco_pass) n_reco_pass++;
            if (gen_pass)  n_gen_pass++;
            if (reco_pass && gen_pass) n_both_pass++;
            if (!reco_sv_ok)                                                                              n_fail_reco_sv++;
            else if (!( jpt_reco >= cfg.kin.ptLow ))  n_fail_reco_pt++;
            else if (!(std::abs(jeta_reco) < cfg.kin.etaMax))                                 n_fail_reco_eta++;
            else if (cfg.physics.useBtag && !passBtag(t, ijet, cfg))                                      n_fail_reco_btag++;
            else if (!(mB_reco_fill >= mb_min && mB_reco_fill < mb_max))                                  n_fail_reco_mb++;
            else if (!(dr_reco_fill < dr_max))                                                            n_fail_reco_dr++;

            if (!(jpt_gen >= cfg.kin.ptLow))       n_fail_gen_pt++;
            else if (!(std::abs(jeta_gen) < cfg.kin.etaMax))                               n_fail_gen_eta++;
            else if (!(mB_gen_fill >= mb_min && mB_gen_fill < mb_max))                                 n_fail_gen_mb++;
            else if (!(dr_gen_fill < dr_max))                                                          n_fail_gen_dr++;

            // --- debug: print first 10 jets that reach cut evaluation ---
            if (n_debug_printed < 10) {
                std::cout << "[DBG jet " << n_debug_printed << "]"
                    << " jpt_reco=" << jpt_reco
                    << " jpt_gen="  <<  t.refpt[ijet]  // jpt_gen
                    << " reco jteta="    << jeta_reco
                    << " gen jteta="   << t.refeta[ijet] //  jeta_gen
                    << " discr="    << btagVar
                    << " mB_reco="  << mB_reco
                    << " dr_reco="  << dr_reco
                    << " mB_gen="   << mB_gen
                    << " dr_gen="   << dr_gen
                    << " reco_sv="  << reco_sv.size()
                    << " gen_bh="   << gen_bh.size()
                    << " reco_pass=" << reco_pass
                    << " gen_pass="  << gen_pass
                    << std::endl;
                n_debug_printed++;
            }
            //---------------------------

            // -------------------------------------------
            // ---- Fill Response matrix -----------------
            // -------------------------------------------
            double num    = distr(generator);
            double w_reco = weight_tree * eec_reco;
            //double w_gen  = weight_tree * eec_gen;

            if (reco_pass) {
              if (num < 0.5) h_half0_purity_den->Fill(dr_reco_fill, jpt_reco, w_reco);
              else           h_half1_purity_den->Fill(dr_reco_fill, jpt_reco, w_reco);

            }

            if (gen_pass) {
                // Intentional: the efficiency is binned at gen level but weighted with
                // the reco-side EEC weight, matching the response matrix convention.
                if (num < 0.5) h_half0_eff_den->Fill(dr_gen_fill, jpt_gen, w_reco);
                else           h_half1_eff_den->Fill(dr_gen_fill, jpt_gen, w_reco);
            }

            if (reco_pass && gen_pass) {
                if (num < 0.5) {
                    h_half0_purity_num->Fill(dr_reco_fill, jpt_reco, w_reco);
                    h_half0_eff_num   ->Fill(dr_gen_fill,  jpt_gen,  w_reco);
                    response_half0->Fill(dr_reco_fill, jpt_reco,
                                         dr_gen_fill,  jpt_gen,  w_reco);
                } else {
                    h_half1_purity_num->Fill( dr_reco_fill, jpt_reco, w_reco);
                    h_half1_eff_num   ->Fill(dr_gen_fill,  jpt_gen,  w_reco);
                    response_half1->Fill(dr_reco_fill, jpt_reco,
                                         dr_gen_fill,  jpt_gen,  w_reco);
                }
                response_full->Fill(dr_reco_fill, jpt_reco,
                                    dr_gen_fill,  jpt_gen,  w_reco);
          }

        } // jet loop
    } // event loop
    // -- Debug output ----------
    std::cout << std::endl;
    std::cout << "--- Jet statistics (bb jets, jtNbHad >= 2) ---" << std::endl;
    std::cout << "  selected bb jets (after triggers):      " << n_bb_jets    << std::endl;
    std::cout << "  Gen bh ok (>= 2 gen bh):     " << n_gen_bh_ok  << std::endl;
    std::cout << "  Passing reco cuts:           " << n_reco_pass  << std::endl;
    std::cout << "  Passing gen cuts:            " << n_gen_pass   << std::endl;
    std::cout << "  Passing both (numerator):    " << n_both_pass  << std::endl;
    std::cout << "  Reco SV failures (< 2 SVs): " << nb_sv        << std::endl;
    std::cout << "  SV purity failures:          " << sv_fail      << std::endl;
    std::cout << "  SV merging failures:         " << merge_fail   << std::endl;
    std::cout << "  No-SV track agg failures:    " << agg_fail     << std::endl;
    std::cout << "--- Reco-pass per-condition failures (first failing cond shown) ---" << std::endl;
    std::cout << "  fail reco_sv_ok:   " << n_fail_reco_sv   << std::endl;
    std::cout << "  fail reco jpt:     " << n_fail_reco_pt   << std::endl;
    std::cout << "  fail reco eta:     " << n_fail_reco_eta  << std::endl;
    std::cout << "  fail reco btag:    " << n_fail_reco_btag << std::endl;
    std::cout << "  fail reco mB:      " << n_fail_reco_mb   << std::endl;
    std::cout << "  fail reco dr:      " << n_fail_reco_dr   << std::endl;
    std::cout << "--- Gen-pass per-condition failures (first failing cond shown) ---" << std::endl;
    std::cout << "  fail gen jpt:      " << n_fail_gen_pt    << std::endl;
    std::cout << "  fail gen eta:      " << n_fail_gen_eta   << std::endl;
    std::cout << "  fail gen mB:       " << n_fail_gen_mb    << std::endl;
    std::cout << "  fail gen dr:       " << n_fail_gen_dr    << std::endl;
    //----------------------------

    // --- Compte Purity, Eff.
    auto divide = [](TH2D* num, TH2D* den, const char* name) -> TH2D* {
        TH2D *h = (TH2D*) num->Clone(name);
        h->Divide(num, den, 1., 1., "b"); // Bionmial error propogation
        return h;
    };
      // For halfs
    TH2D *h_half0_purity = divide(h_half0_purity_num, h_half0_purity_den, "h_half0_purity_tf");
    TH2D *h_half1_purity = divide(h_half1_purity_num, h_half1_purity_den, "h_half1_purity_tf");
    TH2D *h_half0_eff    = divide(h_half0_eff_num,    h_half0_eff_den,    "h_half0_efficiency_tf");
    TH2D *h_half1_eff    = divide(h_half1_eff_num,    h_half1_eff_den,    "h_half1_efficiency_tf");

      // For full
    TH2D *h_full_purity_num = (TH2D*) h_half0_purity_num->Clone("h_full_purity_numerator_tf");
    h_full_purity_num->Add(h_half1_purity_num);
    TH2D *h_full_purity_den = (TH2D*) h_half0_purity_den->Clone("h_full_purity_denominator_tf");
    h_full_purity_den->Add(h_half1_purity_den);
    TH2D *h_full_purity = divide(h_full_purity_num, h_full_purity_den, "h_full_purity_tf");
    TH2D *h_full_eff_num = (TH2D*) h_half0_eff_num->Clone("h_full_efficiency_numerator_tf");
    h_full_eff_num->Add(h_half1_eff_num);
    TH2D *h_full_eff_den = (TH2D*) h_half0_eff_den->Clone("h_full_efficiency_denominator_tf");
    h_full_eff_den->Add(h_half1_eff_den);
    TH2D *h_full_eff = divide(h_full_eff_num, h_full_eff_den, "h_full_efficiency_tf");

    std::cout << "Creating: " << fout_name << std::endl;
    TFile *fout = new TFile(fout_name, "recreate");
    h_half0_purity_num->Write(); h_half0_purity_den->Write(); h_half0_purity->Write();
    h_half0_eff_num->Write();    h_half0_eff_den->Write();    h_half0_eff->Write();
    response_half0->Write();
    h_half1_purity_num->Write(); h_half1_purity_den->Write(); h_half1_purity->Write();
    h_half1_eff_num->Write();    h_half1_eff_den->Write();    h_half1_eff->Write();
    response_half1->Write();
    h_full_purity_num->Write(); h_full_purity_den->Write(); h_full_purity->Write();
    h_full_eff_num->Write();    h_full_eff_den->Write();    h_full_eff->Write();
    response_full->Write();
    fout->Close();
    delete fout;
}
#endif // ---- end disabled: create_response_templatefit ----


inline std::vector<ObsDef> obsList()
{
    return { obsDr(), obsB() };

}

inline std::vector<double> obsValues(double dr, double B)
{
    return { dr, B };
 }

// ---- Per-observable histogram sets ------------------------------------------------------
// One TemplateHistos and one ResponseHistos per measured observable (dr, B). Every name
// comes from ObsDef::n(), whose suffix is "" for dr -- so the dr objects are booked, filled
// and written under exactly the names they have always had, and template_fit.cpp,
// apply_unfolding_2d.C and every result file already on disk are unaffected.
//
// The fill entry points take the observable VALUE as an argument (observables.h explains why
// the value is not computed inside ObsDef), so the caller computes dr and B once each and
// hands every set its own number. Binning, overflow folding and naming are the set's own.

struct TemplateHistos {
    ObsDef o;
    TH3D *h3D_0b = nullptr, *h3D_b = nullptr, *h3D_bb = nullptr;
    TH3D *h3D_pseudo_0b = nullptr, *h3D_pseudo_b = nullptr, *h3D_pseudo_bb = nullptr;
    TH3D *h3D_data = nullptr, *h3D_pseudodata = nullptr, *h3D_pseudodata_bb = nullptr;
    TH3D *h_count_0b = nullptr, *h_count_b = nullptr, *h_count_bb = nullptr, *h_count_data = nullptr;
    // Same templates in 5 GeV jet-pT bins, 80-300 GeV (above 300 folded into the last bin), for the
    // jet-pT reweighting of the MC in template_fit.cpp (2026-10-27, Zoe). Same fill weights.
    // 2026-10-28 (Zoe: "dont keep jets with more than 400 GeV, adding the correspondant bins"):
    // 80-400 GeV, 64 bins, no fold -- jets >= 400 GeV are not filled (kPtFineMax).
    TH3D *h3D_0b_ptfine = nullptr, *h3D_b_ptfine = nullptr, *h3D_bb_ptfine = nullptr, *h3D_data_ptfine = nullptr;
    static constexpr double kPtFineMax = 400.;

    TH3D *book3Dptfine(const char *base, const char *what) const {
        // const int nfine = 44; double e[nfine + 1];   // 80-300 GeV, until 2026-10-28
        const int nfine = 64; double e[nfine + 1];
        for (int i = 0; i <= nfine; ++i) e[i] = 80. + 5. * i;
        TH3D *h = new TH3D(o.n(base),
                           Form("%s;m_{2B} [GeV];%s;p_{T} [GeV]", what, o.axis.Data()),
                           bins_mb, mb_binsVector, o.nbins, o.bins, nfine, e);
        h->Sumw2();
        h->SetCanExtend(TH1::kNoAxis);
        return h;
    }

    TH3D *book3D(const char *base, const char *what) const {
        TH3D *h = new TH3D(o.n(base),
                           Form("%s;m_{2B} [GeV];%s;p_{T} [GeV]", what, o.axis.Data()),
                           bins_mb, mb_binsVector, o.nbins, o.bins, jtpt_bins, jtpt_binsVector);
        h->Sumw2();
        h->SetCanExtend(TH1::kNoAxis);
        return h;
    }

    void book(const ObsDef &def) {
        o = def;
        h3D_0b = book3D("h3D_0b", "EEC 0b");
        h3D_b  = book3D("h3D_b",  "EEC 1b");
        h3D_bb = book3D("h3D_bb", "EEC 2b");
        h3D_pseudo_0b = book3D("h3D_pseudo_0b", "EEC 0b, odd half");
        h3D_pseudo_b  = book3D("h3D_pseudo_b",  "EEC 1b, odd half");
        h3D_pseudo_bb = book3D("h3D_pseudo_bb", "EEC 2b, odd half");
        h3D_data          = book3D("h3D_data",          "EEC data");
        h3D_pseudodata    = book3D("h3D_pseudodata",    "EEC pseudodata, even half");
        h3D_pseudodata_bb = book3D("h3D_pseudodata_bb", "EEC pseudodata 2b, even half");
        h_count_0b   = book3D("h_count_0b",   "jet counts 0b");
        h_count_b    = book3D("h_count_b",    "jet counts 1b");
        h_count_bb   = book3D("h_count_bb",   "jet counts 2b");
        h_count_data = book3D("h_count_data", "jet counts data");
        h3D_0b_ptfine   = book3Dptfine("h3D_0b_ptfine",   "EEC 0b, fine jet pT");
        h3D_b_ptfine    = book3Dptfine("h3D_b_ptfine",    "EEC 1b, fine jet pT");
        h3D_bb_ptfine   = book3Dptfine("h3D_bb_ptfine",   "EEC 2b, fine jet pT");
        h3D_data_ptfine = book3Dptfine("h3D_data_ptfine", "EEC data, fine jet pT");
    }

    // x is this observable's RECO value, already past o.passReco(). eec carries the data
    // prescale the caller applied; w is the tree weight. jtpt is the folded template pT,
    // jtpt_raw the jet pT itself (for the _ptfine histograms).
    void fill(bool isMC, Int_t jtNbHad, Long64_t ient,
              double mB, double x, double jtpt, double eec, double w, double jtpt_raw) {
        const double xf = o.fill(x);
        // const double ptf = std::min(jtpt_raw, 299.9);   // fold above 300 GeV into the last fine bin (until 2026-10-28)
        const double ptf = jtpt_raw;
        const bool fine_in = jtpt_raw < kPtFineMax;    // jets >= 400 GeV: not in the _ptfine histograms
        if (isMC) {
            if      (jtNbHad == 0) { h3D_0b->Fill(mB, xf, jtpt, eec * w); h_count_0b->Fill(mB, xf, jtpt, w); }
            else if (jtNbHad == 1) { h3D_b ->Fill(mB, xf, jtpt, eec * w); h_count_b ->Fill(mB, xf, jtpt, w); }
            else if (jtNbHad >= 2) { h3D_bb->Fill(mB, xf, jtpt, eec * w); h_count_bb->Fill(mB, xf, jtpt, w); }
            if      (!fine_in)     ;
            else if (jtNbHad == 0) h3D_0b_ptfine->Fill(mB, xf, ptf, eec * w);
            else if (jtNbHad == 1) h3D_b_ptfine ->Fill(mB, xf, ptf, eec * w);
            else if (jtNbHad >= 2) h3D_bb_ptfine->Fill(mB, xf, ptf, eec * w);

            if (ient % 2 == 0) {
                h3D_pseudodata->Fill(mB, xf, jtpt, eec * w);
                // >=2 matches the response matrix gate, not the ==2 of the templates
                if (jtNbHad >= 2) h3D_pseudodata_bb->Fill(mB, xf, jtpt, eec * w);
            } else {
                if      (jtNbHad == 0) h3D_pseudo_0b->Fill(mB, xf, jtpt, eec * w);
                else if (jtNbHad == 1) h3D_pseudo_b ->Fill(mB, xf, jtpt, eec * w);
                else if (jtNbHad >= 2) h3D_pseudo_bb->Fill(mB, xf, jtpt, eec * w);
            }
        } else {
            h3D_data    ->Fill(mB, xf, jtpt, eec * w);
            h_count_data->Fill(mB, xf, jtpt, w);
            // h3D_data_ptfine->Fill(mB, xf, ptf, eec * w);   // until 2026-10-28
            if (fine_in) h3D_data_ptfine->Fill(mB, xf, ptf, eec * w);
        }
    }

    void write(bool isMC) const {
        if (isMC) {
            h3D_0b->Write(); h3D_b->Write(); h3D_bb->Write();
            h3D_pseudodata->Write(); h3D_pseudodata_bb->Write();
            // Odd-half templates, for fitting the even-half pseudodata in the split test
            h3D_pseudo_0b->Write(); h3D_pseudo_b->Write(); h3D_pseudo_bb->Write();
            h_count_0b->Write(); h_count_b->Write(); h_count_bb->Write();
            h3D_0b_ptfine->Write(); h3D_b_ptfine->Write(); h3D_bb_ptfine->Write();
        } else {
            h3D_data->Write(); h_count_data->Write();
            h3D_data_ptfine->Write();
        }
    }
};

struct ResponseHistos {
    ObsDef o;
    // [0] = half0, [1] = half1
    TH2D *purity_num[2] = {nullptr, nullptr}, *purity_den[2] = {nullptr, nullptr};
    TH2D *eff_num[2]    = {nullptr, nullptr}, *eff_den[2]    = {nullptr, nullptr};
    TH2D *ps_purity_num[2] = {nullptr, nullptr}, *ps_purity_den[2] = {nullptr, nullptr};
    TH2D *ps_eff_num[2]    = {nullptr, nullptr}, *ps_eff_den[2]    = {nullptr, nullptr};
    TH2D *pseudodata_truth = nullptr;
    RooUnfoldResponse *response[2] = {nullptr, nullptr}, *response_full = nullptr;
    RooUnfoldResponse *ps_response[2] = {nullptr, nullptr}, *ps_response_full = nullptr;
    TH2D *hgenjet_2b = nullptr, *hgenjet_2b_passbtag = nullptr;
    TH2D *hgenjet_2b_all = nullptr, *hgenjet_2b_reco_btag = nullptr;

    TH2D *book2D(const char *base, const char *title) const {
        return new TH2D(o.n(base), title, o.nbins, o.bins, jtpt_bins, jtpt_binsVector);
    }
    // gen-level histograms: the gen pT binning (= the reco one unless JET_PT_MAX, then + a gen overflow bin)
    TH2D *book2Dgen(const char *base, const char *title) const {
        return new TH2D(o.n(base), title, o.nbins, o.bins, jtptGenNbins(), jtptGenEdges());
    }

    void book(const ObsDef &def) {
        o = def;
        const TString xr = Form("x=%s_reco, y=jtpt", o.name.Data());
        const TString xg = Form("x=%s_gen, y=jtpt",  o.name.Data());
        for (int h = 0; h < 2; ++h) {
            purity_num[h] = book2D(Form("h_half%d_purity_numerator_tf",   h), xr);
            purity_den[h] = book2D(Form("h_half%d_purity_denominator_tf", h), xr);
            eff_num[h]    = book2Dgen(Form("h_half%d_efficiency_numerator_tf",   h), xg);   // book2D until 2026-10-29
            eff_den[h]    = book2Dgen(Form("h_half%d_efficiency_denominator_tf", h), xg);
            ps_purity_num[h] = book2D(Form("h_half%d_pseudo_purity_numerator_tf",   h), xr);
            ps_purity_den[h] = book2D(Form("h_half%d_pseudo_purity_denominator_tf", h), xr);
            ps_eff_num[h]    = book2Dgen(Form("h_half%d_pseudo_efficiency_numerator_tf",   h), xg);
            ps_eff_den[h]    = book2Dgen(Form("h_half%d_pseudo_efficiency_denominator_tf", h), xg);
        }
        // Closure target for the split test: the pseudodata half's gen-level distribution,
        // with the same gate (gen_pass) and weight (w_reco) as the purity/response chain.
        pseudodata_truth = book2Dgen("h_pseudodata_truth_tf", xg);
        pseudodata_truth->Sumw2();

        for (int h = 0; h < 2; ++h) {
            response[h] = new RooUnfoldResponse(purity_den[h], eff_den[h],
                                                o.n(Form("response_tf_half%d", h)),
                                                Form("tf response half%d", h));
            ps_response[h] = new RooUnfoldResponse(ps_purity_den[h], ps_eff_den[h],
                                                   o.n(Form("response_tf_pseudo_half%d", h)),
                                                   Form("tf response pseudo half%d", h));
        }
        response_full    = new RooUnfoldResponse(purity_den[0], eff_den[0],
                                                 o.n("response_tf_full"), "tf response full");
        ps_response_full = new RooUnfoldResponse(ps_purity_den[0], ps_eff_den[0],
                                                 o.n("response_tf_pseudo_full"),
                                                 "tf response pseudo full");

        // -- For b-tagging eff. correction after unfolding (at particle level)
        hgenjet_2b          = book2Dgen("hgenjet_2b",          Form("b-tagging eff. DENO;%s;p_{T} [GeV]", o.axis.Data()));
        hgenjet_2b_passbtag = book2Dgen("hgenjet_2b_passbtag", Form("b-tagging eff. NUM;%s;p_{T} [GeV]",  o.axis.Data()));
        // Combined SV-reco + b-tag efficiency, gen-EEC weighted (w_gen) so the denominator
        // is defined even for true 2b jets without 2 reco SVs.
        hgenjet_2b_all       = book2Dgen("hgenjet_2b_all",       Form("SV+btag eff. DENO;%s;p_{T} [GeV]", o.axis.Data()));
        hgenjet_2b_reco_btag = book2Dgen("hgenjet_2b_reco_btag", Form("SV+btag eff. NUM;%s;p_{T} [GeV]",  o.axis.Data()));
    }

  
    // Gen values are folded with o.fill(), like the response/efficiency fills, so the last
    // bin means the same thing here as in the unfolded result (dr_gen >= 0.45 -> last bin,
    // not ROOT's overflow). Changed 2026-10-26; before, these were filled unfolded.
    void fillGenAll(double x_gen_raw, double jpt_gen, double w_gen) const {
        // hgenjet_2b_all->Fill(x_gen_raw, jpt_gen, w_gen);
        hgenjet_2b_all->Fill(o.fill(x_gen_raw), jpt_gen, w_gen);
    }
    void fillGenPreBtag(double x_gen_raw, double jpt_gen, double w_reco) const {
        // hgenjet_2b->Fill(x_gen_raw, jpt_gen, w_reco);
        hgenjet_2b->Fill(o.fill(x_gen_raw), jpt_gen, w_reco);
    }
    void fillGenPostBtag(double x_gen_raw, double jpt_gen, double w_reco, double w_gen) const {
        // hgenjet_2b_passbtag  ->Fill(x_gen_raw, jpt_gen, w_reco);
        // hgenjet_2b_reco_btag ->Fill(x_gen_raw, jpt_gen, w_gen);
        hgenjet_2b_passbtag  ->Fill(o.fill(x_gen_raw), jpt_gen, w_reco);
        hgenjet_2b_reco_btag ->Fill(o.fill(x_gen_raw), jpt_gen, w_gen);
    }

    // The response/purity/efficiency fills. x_reco_fill and x_gen_fill are already folded.
    void fillResponse(bool reco_pass, bool gen_pass, double num, Long64_t ient,
                      double x_reco_fill, double jpt_reco,
                      double x_gen_fill,  double jpt_gen, double w_reco) const {
        const int h = (num < 0.5) ? 0 : 1;

        if (reco_pass) {
            purity_den[h]->Fill(x_reco_fill, jpt_reco, w_reco);
            // Odd half: same half as the pseudo numerators and response, so purity is a
            // self-consistent ratio. The even half is reserved for the pseudodata.
            if (ient % 2 == 1) ps_purity_den[h]->Fill(x_reco_fill, jpt_reco, w_reco);
        }

        if (gen_pass) {
            // Intentional: the efficiency is binned at gen level but weighted with the
            // reco-side EEC weight, matching the response matrix convention.
            eff_den[h]->Fill(x_gen_fill, jpt_gen, w_reco);
            if (ient % 2 == 1) ps_eff_den[h]->Fill(x_gen_fill, jpt_gen, w_reco);
            // Even half = the pseudodata half. This is what the corrected chain must reproduce.
            else pseudodata_truth->Fill(x_gen_fill, jpt_gen, w_reco);
        }

        if (reco_pass && gen_pass) {
            purity_num[h]->Fill(x_reco_fill, jpt_reco, w_reco);
            eff_num[h]   ->Fill(x_gen_fill,  jpt_gen,  w_reco);
            response[h]  ->Fill(x_reco_fill, jpt_reco, x_gen_fill, jpt_gen, w_reco);
            response_full->Fill(x_reco_fill, jpt_reco, x_gen_fill, jpt_gen, w_reco);

            if (ient % 2 == 1) {
                ps_purity_num[h]->Fill(x_reco_fill, jpt_reco, w_reco);
                ps_eff_num[h]   ->Fill(x_gen_fill,  jpt_gen,  w_reco);
                ps_response[h]  ->Fill(x_reco_fill, jpt_reco, x_gen_fill, jpt_gen, w_reco);
                ps_response_full->Fill(x_reco_fill, jpt_reco, x_gen_fill, jpt_gen, w_reco);
            }
        }
    }

    // Derive the ratios and write everything. Same operations, same names as before.
    void write() const {
        auto divide = [](TH2D *num, TH2D *den, const TString &name) -> TH2D * {
            TH2D *h = (TH2D *) num->Clone(name);
            h->Divide(num, den, 1., 1., "b");   // Binomial error propagation
            return h;
        };
        auto sum2 = [](TH2D *a, TH2D *b, const TString &name) -> TH2D * {
            TH2D *h = (TH2D *) a->Clone(name);
            h->Add(b);
            return h;
        };

        for (int h = 0; h < 2; ++h) {
            purity_num[h]->Write(); purity_den[h]->Write();
            divide(purity_num[h], purity_den[h], o.n(Form("h_half%d_purity_tf", h)))->Write();
            eff_num[h]->Write(); eff_den[h]->Write();
            divide(eff_num[h], eff_den[h], o.n(Form("h_half%d_efficiency_tf", h)))->Write();
            response[h]->Write();
        }
        TH2D *full_purity_num = sum2(purity_num[0], purity_num[1], o.n("h_full_purity_numerator_tf"));
        TH2D *full_purity_den = sum2(purity_den[0], purity_den[1], o.n("h_full_purity_denominator_tf"));
        TH2D *full_eff_num    = sum2(eff_num[0],    eff_num[1],    o.n("h_full_efficiency_numerator_tf"));
        TH2D *full_eff_den    = sum2(eff_den[0],    eff_den[1],    o.n("h_full_efficiency_denominator_tf"));
        full_purity_num->Write(); full_purity_den->Write();
        divide(full_purity_num, full_purity_den, o.n("h_full_purity_tf"))->Write();
        full_eff_num->Write(); full_eff_den->Write();
        divide(full_eff_num, full_eff_den, o.n("h_full_efficiency_tf"))->Write();
        response_full->Write();

        for (int h = 0; h < 2; ++h) {
            ps_purity_num[h]->Write(); ps_purity_den[h]->Write();
            divide(ps_purity_num[h], ps_purity_den[h], o.n(Form("h_half%d_pseudo_purity_tf", h)))->Write();
            ps_eff_num[h]->Write(); ps_eff_den[h]->Write();
            divide(ps_eff_num[h], ps_eff_den[h], o.n(Form("h_half%d_pseudo_efficiency_tf", h)))->Write();
            ps_response[h]->Write();
        }
        TH2D *full_ps_purity_num = sum2(ps_purity_num[0], ps_purity_num[1], o.n("h_full_pseudo_purity_numerator_tf"));
        TH2D *full_ps_purity_den = sum2(ps_purity_den[0], ps_purity_den[1], o.n("h_full_pseudo_purity_denominator_tf"));
        TH2D *full_ps_eff_num    = sum2(ps_eff_num[0],    ps_eff_num[1],    o.n("h_full_pseudo_efficiency_numerator_tf"));
        TH2D *full_ps_eff_den    = sum2(ps_eff_den[0],    ps_eff_den[1],    o.n("h_full_pseudo_efficiency_denominator_tf"));
        full_ps_purity_num->Write(); full_ps_purity_den->Write();
        divide(full_ps_purity_num, full_ps_purity_den, o.n("h_full_pseudo_purity_tf"))->Write();
        full_ps_eff_num->Write(); full_ps_eff_den->Write();
        divide(full_ps_eff_num, full_ps_eff_den, o.n("h_full_pseudo_efficiency_tf"))->Write();
        ps_response_full->Write();
        pseudodata_truth->Write();

        // -- b-tagging efficiency correction
        TH2D *btagEff = divide(hgenjet_2b_passbtag, hgenjet_2b, o.n("hbtagEff_correction_plevel"));
        // -- combined SV-reco + b-tag efficiency (num = 2 reco SV + btag, den = all true 2b)
        TH2D *svbtagEff = divide(hgenjet_2b_reco_btag, hgenjet_2b_all, o.n("hSVbtagEff_correction_plevel"));
        // EEC-weight correction (reco -> gen EEC weight): r_eec = sum(eec_gen)/sum(eec_reco).
        // Same jets, different weight, so NOT a subset ratio -> normal (not binomial) errors.
        // Applied by MULTIPLYING the unfolded result.
        TH2D *eecwEff = (TH2D *) hgenjet_2b_reco_btag->Clone(o.n("hEECweightEff_correction_plevel"));
        eecwEff->Divide(hgenjet_2b_reco_btag, hgenjet_2b_passbtag, 1., 1., "");

        // DEBUG: global scale of the EEC-weight correction (should be ~1 if reco pt ~ gen pt).
        const double I_reco_eec = hgenjet_2b_passbtag->Integral();   // reco-EEC weighted
        const double I_gen_eec  = hgenjet_2b_reco_btag->Integral();  // gen-EEC weighted
        std::cout << "[EEC-weight debug " << o.name << "] Integral reco-EEC (hgenjet_2b_passbtag)  = " << I_reco_eec << std::endl;
        std::cout << "[EEC-weight debug " << o.name << "] Integral gen-EEC  (hgenjet_2b_reco_btag) = " << I_gen_eec  << std::endl;
        std::cout << "[EEC-weight debug " << o.name << "] global eec_gen/eec_reco = "
                  << (I_reco_eec > 0 ? I_gen_eec / I_reco_eec : 0) << std::endl;

        hgenjet_2b->Write(); hgenjet_2b_passbtag->Write(); btagEff->Write();
        hgenjet_2b_all->Write(); hgenjet_2b_reco_btag->Write(); svbtagEff->Write(); eecwEff->Write();
    }
};

// ---- 3D response in (dR, A, jet pT), added 2026-10-27 ---------------------------------------
// For the 3D Bayesian unfolding. Reco and gen have DIFFERENT binnings (the gen-only side bins,
// see u3d_* in binning_histos_small.h), so the response is non-square. Axes: x = dR, y = A,
// z = jet pT -- the 2D convention (observable first, pT last) extended by one.
//
// Same fill pattern and weights as ResponseHistos, except:
//  - the gates are the 3D ones: gen pT >= u3d_gen_pt[0] (60, not ptLow) and the reco hard
//    edges (u3d_reco_in), so these purity/efficiency are NOT the 2D ones;
//  - full sample + split test only (odd half = response, even half = pseudodata reco and
//    truth), no jackknife halves;
//  - values come in RAW (no jtpt_fill / o.fill), the clamps are applied here;
//  - only num/den are written, the ratios are formed in the unfolding from the block sums.
// One weight per production, like everything else: EEC-weighted in the nominal production,
// unweighted in the _noeecw one.
struct Response3D {
    TH3D *purity_num = nullptr, *purity_den = nullptr, *eff_num = nullptr, *eff_den = nullptr;
    TH3D *ps_purity_num = nullptr, *ps_purity_den = nullptr, *ps_eff_num = nullptr, *ps_eff_den = nullptr;
    TH3D *pseudodata_reco = nullptr, *pseudodata_truth = nullptr;
    TH3D *hgenjet_2b = nullptr, *hgenjet_2b_passbtag = nullptr;
    TH3D *hgenjet_2b_all = nullptr, *hgenjet_2b_reco_btag = nullptr;
    RooUnfoldResponse *response_full = nullptr, *ps_response_full = nullptr;

    TH3D *bookReco(const char *name, const char *title) const {
        TH3D *h = new TH3D(name, title, u3d_reco_dr_size - 1, u3d_reco_dr,
                           B_binsVectorSize - 1, B_binsVector, u3d_reco_pt_size - 1, u3d_reco_pt);
        h->Sumw2();
        return h;
    }
    TH3D *bookGen(const char *name, const char *title) const {
        TH3D *h = new TH3D(name, title, u3d_gen_dr_size - 1, u3d_gen_dr,
                           B_binsVectorSize - 1, B_binsVector, u3d_gen_pt_size - 1, u3d_gen_pt);
        h->Sumw2();
        return h;
    }

    void book() {
        const char *xr = "x=dr_reco, y=A_reco, z=jtpt_reco";
        const char *xg = "x=dr_gen, y=A_gen, z=jtpt_gen";
        purity_num = bookReco("h_full_purity_numerator_3d",   xr);
        purity_den = bookReco("h_full_purity_denominator_3d", xr);
        eff_num    = bookGen("h_full_efficiency_numerator_3d",   xg);
        eff_den    = bookGen("h_full_efficiency_denominator_3d", xg);
        ps_purity_num = bookReco("h_full_pseudo_purity_numerator_3d",   xr);
        ps_purity_den = bookReco("h_full_pseudo_purity_denominator_3d", xr);
        ps_eff_num    = bookGen("h_full_pseudo_efficiency_numerator_3d",   xg);
        ps_eff_den    = bookGen("h_full_pseudo_efficiency_denominator_3d", xg);
        // Split test: the even half, what the odd-half chain must recover.
        pseudodata_reco  = bookReco("h_pseudodata_reco_3d",  xr);
        pseudodata_truth = bookGen("h_pseudodata_truth_3d", xg);
        response_full    = new RooUnfoldResponse(purity_den, eff_den, "response_full_3d",
                                                 "3D response full");
        ps_response_full = new RooUnfoldResponse(ps_purity_den, ps_eff_den, "response_pseudo_full_3d",
                                                 "3D response pseudo (odd half)");
        // Post-unfolding corrections, gen binning, same definitions as the 2D hgenjet_2b_*.
        hgenjet_2b           = bookGen("hgenjet_2b_3d",           "b-tagging eff. DENO");
        hgenjet_2b_passbtag  = bookGen("hgenjet_2b_passbtag_3d",  "b-tagging eff. NUM");
        hgenjet_2b_all       = bookGen("hgenjet_2b_all_3d",       "SV+btag eff. DENO");
        hgenjet_2b_reco_btag = bookGen("hgenjet_2b_reco_btag_3d", "SV+btag eff. NUM");
    }

    void fillGen(TH3D *h, double dr, double A, double pt, double w) const {
        h->Fill(u3d_gen_dr_fill(dr), u3d_A_fill(A), u3d_gen_pt_fill(pt), w);
    }
    void fillGenAll(double dr, double A, double pt, double w_gen) const {
        fillGen(hgenjet_2b_all, dr, A, pt, w_gen);
    }
    void fillGenPreBtag(double dr, double A, double pt, double w_reco) const {
        fillGen(hgenjet_2b, dr, A, pt, w_reco);
    }
    void fillGenPostBtag(double dr, double A, double pt, double w_reco, double w_gen) const {
        fillGen(hgenjet_2b_passbtag,  dr, A, pt, w_reco);
        fillGen(hgenjet_2b_reco_btag, dr, A, pt, w_gen);
    }

    // reco_pass must already include u3d_reco_in, gen_pass the 60 GeV gen threshold.
    void fillResponse(bool reco_pass, bool gen_pass, Long64_t ient,
                      double dr_reco, double A_reco, double pt_reco,
                      double dr_gen,  double A_gen,  double pt_gen, double w_reco) const {
        const bool odd = (ient % 2 == 1);
        const double ar = u3d_A_fill(A_reco);
        const double dg = u3d_gen_dr_fill(dr_gen), ag = u3d_A_fill(A_gen), pg = u3d_gen_pt_fill(pt_gen);
        if (reco_pass) {
            purity_den->Fill(dr_reco, ar, pt_reco, w_reco);
            if (odd) ps_purity_den->Fill(dr_reco, ar, pt_reco, w_reco);
            else     pseudodata_reco->Fill(dr_reco, ar, pt_reco, w_reco);
        }
        if (gen_pass) {
            eff_den->Fill(dg, ag, pg, w_reco);   // reco-side weight, as in the 2D efficiency
            if (odd) ps_eff_den->Fill(dg, ag, pg, w_reco);
            else     pseudodata_truth->Fill(dg, ag, pg, w_reco);
        }
        if (reco_pass && gen_pass) {
            purity_num->Fill(dr_reco, ar, pt_reco, w_reco);
            eff_num->Fill(dg, ag, pg, w_reco);
            response_full->Fill(dr_reco, ar, pt_reco, dg, ag, pg, w_reco);
            if (odd) {
                ps_purity_num->Fill(dr_reco, ar, pt_reco, w_reco);
                ps_eff_num->Fill(dg, ag, pg, w_reco);
                ps_response_full->Fill(dr_reco, ar, pt_reco, dg, ag, pg, w_reco);
            }
        }
    }

    void write() const {
        for (TH3D *h : {purity_num, purity_den, eff_num, eff_den,
                        ps_purity_num, ps_purity_den, ps_eff_num, ps_eff_den,
                        pseudodata_reco, pseudodata_truth,
                        hgenjet_2b, hgenjet_2b_passbtag, hgenjet_2b_all, hgenjet_2b_reco_btag})
            h->Write();
        response_full->Write();
        ps_response_full->Write();
    }
};

void Build_templates(const AnalysisConfig& cfg, bool isMakeTemplates = true, bool isCreateRmatrix = true,
                     Long64_t ev_first = 0, Long64_t ev_last = -1, Int_t job_idx = -1, bool makeAggNtuple = true) {
  // -- make templates of Data/MC, and Response matrix for MC

#if !HAS_ROOUNFOLD
  if (isCreateRmatrix) {
    std::cout << "Warning: RooUnfold headers are not available. "
              << "RooUnfoldResponse objects will be no-op placeholders; "
              << "histograms and aggBHadronKinematics will still be written."
              << std::endl;
  }
#endif

  // -- Output files name
  TString job_suffix = (job_idx >= 0) ? Form("_job%d", job_idx) : "";
 
  TString syst_suffix = TrkEffSyst::tag() + SvSF::tag() + TrackBdtSF::tag()
                       + JesSyst::tag() + JerSyst::tag() + EecWeight::tag();
  TString fout_name = cfg.dataset.output_folder + cfg.dataset.output_hist + job_suffix + "MCGEN" + syst_suffix + ".root"; // for reposnse matrix: has Prefix: Response
  TString ResponseMatrix_fout_name =  cfg.dataset.output_folder + "RMatrix_" + cfg.dataset.output_hist + job_suffix + syst_suffix + ".root";
  TString AggBHadronNtuple_fout_name = cfg.dataset.output_folder + "AggBHadronNtuple_" + cfg.dataset.output_hist + job_suffix + syst_suffix + ".root";
  gSystem->mkdir(cfg.dataset.output_folder, true);

  // ---- Histograms ----
  const std::vector<ObsDef> observables = obsList();
 
  if (obsValues(0., 0.).size() != observables.size()) {
    std::cerr << "FATAL: obsList() has " << observables.size() << " observables but "
              << "obsValues() returns " << obsValues(0., 0.).size() << " values. "
              << "Add the new observable to BOTH." << std::endl;
    return;
  }
  std::cout << "Measured observables (second histogram axis):" << std::endl;
  for (const ObsDef &d : observables)
    std::cout << "  " << d.name << ": " << d.nbins << " bins ["
              << d.bins[0] << ", " << d.bins[d.nbins] << "]"
              << (d.suffix.Length() ? Form(", name suffix '%s'", d.suffix.Data())
                                    : ", original names")
              << std::endl;

  std::vector<TemplateHistos> tmpl(observables.size());
  SvSeedDiag svDiag; svDiag.book(cfg.dataset.isMC);   // 2026-11-06, empty unless SV_SEED_DIAG=on
  for (size_t i = 0; i < observables.size(); ++i) tmpl[i].book(observables[i]);

  // Selected jets after btagging + >=2 SV, within the kinematic selection
  TH1D *hpt_selectedJets   = new TH1D("hpt_selectedJets", "For templates: Selected jets (btagged, NSvx >=2, within kienmatics(jet pt and eta + skipped events of large weights) of anlalysis);p_{T} [GeV]; weighted counts(event weight)",jtpt_bins, jtpt_binsVector);
  TH1D *hpt_selectedJets_noweight   = new TH1D("hpt_selectedJets_noweight", "For templates: Selected jets (btagged, NSvx >=2, within kienmatics(jet pt and eta + skipped events of large weights) of anlalysis);p_{T} [GeV]; counts",jtpt_bins, jtpt_binsVector);

  hpt_selectedJets->Sumw2(); hpt_selectedJets_noweight->Sumw2();

  // -- For response matrix

    // Re-derive bin counts from const sizes (Cling init-order workaround);
    // same as mb_bins, dr_bins, jtpt_bins in binning_histos_small.h
    Int_t n_mb = mb_binsVectorSize - 1, n_dr = dr_binsVectorSize - 1, n_pt = jtpt_binsVectorSize - 1;
    Double_t mb_min = mb_binsVector[0], mb_max = mb_binsVector[n_mb], mb_max_fill = 9.9;
    Double_t dr_min = dr_binsVector[0], dr_max = dr_binsVector[n_dr];
    std::cout << "Histogram binning: mB=" << n_mb << " dr=" << n_dr << " jtpt=" << n_pt << std::endl;
    std::cout << "mb range: [" << mb_min << ", " << mb_max << "], fill cap=" << mb_max_fill << std::endl;
    std::cout << "dr range: [" << dr_min << ", " << dr_max << "]" << std::endl;
    // -------------------------------------------------------------------

    // One response set per observable, same ObsDef list as the templates above. dr's names
    // are unchanged; B's carry "_B".
    std::vector<ResponseHistos> resp(observables.size());
    for (size_t i = 0; i < observables.size(); ++i) resp[i].book(observables[i]);
    // The 3D (dR, A, jet pT) response, filled next to the 2D ones (2026-10-27).
    Response3D resp3d;
    resp3d.book();
    // Its gen pT threshold is the bottom of the 60-80 underflow bin, not ptLow.
    AnalysisConfig cfg3d = cfg;
    cfg3d.kin.ptLow = u3d_gen_pt[0];
    cfg3d.kin.ptHigh = 0.;   // the 3D keeps its own edges (gen overflow to 1000, reco hard edge at 200)

	// DEBUG: reco SV pt vs gen partial-B pt (EEC weight = pt_i * pt_j).
	TH2D* h_svpt_vs_bpt = new TH2D("h_svpt_vs_bpt", "reco SV p_{T} vs gen B p_{T};gen partial-B p_{T} [GeV];reco SV p_{T} [GeV]", 60, 0, 120, 60, 0, 120);
	TH1D* h_ptratio     = new TH1D("h_ptratio", "reco/gen SV p_{T};p_{T}^{reco}/p_{T}^{gen};entries", 60, 0., 2.);
	h_svpt_vs_bpt->Sumw2(); h_ptratio->Sumw2();
	const bool doAggNtuple = makeAggNtuple && cfg.dataset.isMC;
	TFile* fout_agg = nullptr;
	TTree* aggBHadronTree = nullptr;
	AggBHadronNtupleRow aggRow;
	if (doAggNtuple) {
	  std::cout << "Creating aggregation diagnostic ntuple: "
	            << AggBHadronNtuple_fout_name << std::endl;
	  fout_agg = new TFile(AggBHadronNtuple_fout_name, "RECREATE");
	  if (!fout_agg || fout_agg->IsZombie()) {
	    std::cerr << "ERROR: could not create aggregation ntuple file: "
	              << AggBHadronNtuple_fout_name << std::endl;
	    delete fout_agg;
	    fout_agg = nullptr;
	  } else {
	    aggBHadronTree = new TTree("aggBHadronKinematics",
	                               "Reco and charged-truth aggregated b-hadron kinematics");
	    makeAggBHadronBranches(aggBHadronTree, aggRow);
	  }
	}

    // Data / reco MC variables
    // -- needed parameters for reco SV func()
      double agg_fail = 0, nb_sv = 0, sv_fail = 0, merge_fail = 0; // for templates conditions
      double agg_fail_rm = 0, nb_sv_rm = 0, sv_fail_rm = 0, merge_fail_rm = 0; // for Response matrix conditions

  // Response-matrix variables
    
    //-- debugging variables
    long n_bb_jets = 0, n_reco_pass = 0, n_gen_pass = 0, n_both_pass = 0;
    long n_gen_bh_ok = 0;  // jets passing gen_bh.size() >= 2
    // per-condition failure counters for reco_pass, then gen_pass
    long n_fail_reco_sv = 0, n_fail_reco_pt = 0, n_fail_reco_eta = 0, n_fail_reco_btag = 0, n_fail_reco_mb = 0, n_fail_reco_dr = 0;
    long n_fail_gen_pt = 0, n_fail_gen_eta = 0, n_fail_gen_mb = 0, n_fail_gen_dr = 0;
    int  n_debug_printed = 10, n_weight_debug_printed = 10;
    long n_weight_checked = 0, n_zero_weight_tree = 0, n_zero_eec_gen = 0, n_zero_eec_reco = 0, n_zero_w_reco = 0;
    //long n_zero_w_gen = 0;
    double sum_weight_tree_checked = 0.0, sum_w_reco_checked = 0.0;
    //double sum_w_gen_checked = 0.0;

  // ---- Event loop ----
  tTree t;
  t.Init(cfg.dataset.filename, cfg.dataset.isMC, cfg.dataset.RunN);
  t.SetBranchStatus("*", 0);
  auto active_branches = getActiveBranches(cfg);
  t.SetBranchStatus(active_branches, 1);

  // JEC (see the Jec namespace at the top): MC-truth levels on everything, residual on data.
  JetCorrector jec;
  if (Jec::enabled()) {
    const std::vector<std::string> jec_files = Jec::files(cfg.dataset.isMC);
    for (const std::string &f : jec_files)
      if (gSystem->AccessPathName(f.c_str()))
        throw std::runtime_error("Missing JEC file: " + f);
    jec.Initialize(jec_files);
  }
  JetUncertainty jesUncertainty;
  if (JesSyst::varied()) {
    const TString uncertaintyFile = JesSyst::uncertaintyFile();
    if (gSystem->AccessPathName(uncertaintyFile))
      throw std::runtime_error("Missing JES uncertainty file: "
                               + std::string(uncertaintyFile.Data()));
    jesUncertainty.Initialize(uncertaintyFile.Data());
  }
  std::unique_ptr<JetSmearer> jerSmearer;
  if (JerSyst::enabled()) {
    const TString resolutionFile = JerSyst::resolutionFile();
    const TString scaleFactorFile = JerSyst::scaleFactorFile();
    if (gSystem->AccessPathName(resolutionFile))
      throw std::runtime_error("Missing JER resolution file: "
                               + std::string(resolutionFile.Data()));
    if (gSystem->AccessPathName(scaleFactorFile))
      throw std::runtime_error("Missing JER scale-factor file: "
                               + std::string(scaleFactorFile.Data()));
    jerSmearer.reset(new JetSmearer(resolutionFile.Data(), scaleFactorFile.Data(),
                                    JetSmearing::Method::Hybrid));
  }
  TLeaf* hiEvtRhoLeaf = nullptr;

  TFile* eventInfoFile = nullptr;
  TTree *hiEvtTree = nullptr, *hltTree = nullptr, *skimTree = nullptr, *weightTree = nullptr;
  Float_t hiEvtVz = 999.0, hiEvtWeight = 1.0;
  Int_t skimPprimaryVertexFilter = 1;
  Int_t hltAK4PFJet40 = 0, hltAK4PFJet60 = 0, hltAK4PFJet80 = 0, hltAK4PFJet100 = 0, hltAK4PFJet120 = 0;

  eventInfoFile = TFile::Open(cfg.dataset.filename);
  if (eventInfoFile && !eventInfoFile->IsZombie()) {
    hiEvtTree = (TTree*)eventInfoFile->Get("hiEvtAnalyzer/HiTree");
    if (hiEvtTree) {
      hiEvtTree->SetBranchStatus("*", 0);
      if (Jec::enabled()) {
        if (!hiEvtTree->GetBranch("rho"))
          throw std::runtime_error("JEC requires hiEvtAnalyzer/HiTree::rho in the input forest");
        hiEvtTree->SetBranchStatus("rho", 1);
        hiEvtRhoLeaf = hiEvtTree->GetLeaf("rho");
        if (!hiEvtRhoLeaf)
          throw std::runtime_error("Could not read hiEvtAnalyzer/HiTree::rho leaf");
      }
      if (hiEvtTree->GetBranch("vz")) {
        hiEvtTree->SetBranchStatus("vz", 1);
        hiEvtTree->SetBranchAddress("vz", &hiEvtVz);
        std::cout << "Reading event vz from hiEvtAnalyzer/HiTree::vz" << std::endl;
      }
      if (cfg.dataset.isMC) {
        weightTree = hiEvtTree;
        if (weightTree && weightTree->GetBranch("weight")) {
          weightTree->SetBranchStatus("weight", 1);
          weightTree->SetBranchAddress("weight", &hiEvtWeight);
          std::cout << "Reading MC event weights from hiEvtAnalyzer/HiTree::weight" << std::endl;
        } else {
          std::cout << "WARNING: could not find hiEvtAnalyzer/HiTree::weight; falling back to t.weight" << std::endl;
          weightTree = nullptr;
        }
      }
    }
    hltTree = (TTree*)eventInfoFile->Get("hltanalysis/HltTree");
    if (hltTree && cfg.dataset.RunN == 3) {
      hltTree->SetBranchStatus("*", 0);
      const std::vector<std::pair<const char *, Int_t *>> hltBits = {
        {"HLT_AK4PFJet40_v8", &hltAK4PFJet40}, {"HLT_AK4PFJet60_v8", &hltAK4PFJet60}, {"HLT_AK4PFJet80_v8", &hltAK4PFJet80},
        {"HLT_AK4PFJet100_v8", &hltAK4PFJet100}, {"HLT_AK4PFJet120_v8", &hltAK4PFJet120} };
      for (const auto &b : hltBits)
        if (hltTree->GetBranch(b.first)) { hltTree->SetBranchStatus(b.first, 1); hltTree->SetBranchAddress(b.first, b.second); }
      std::cout << "Reading Run3 HLT bits from hltanalysis/HltTree" << std::endl;
    }
    skimTree = (TTree*)eventInfoFile->Get("skimanalysis/HltTree");
    if (skimTree && cfg.dataset.RunN == 3 && !cfg.dataset.isMC) {
      skimTree->SetBranchStatus("*", 0);
      if (skimTree->GetBranch("pprimaryVertexFilter")) {
        skimTree->SetBranchStatus("pprimaryVertexFilter", 1);
        skimTree->SetBranchAddress("pprimaryVertexFilter", &skimPprimaryVertexFilter);
        std::cout << "Reading pprimaryVertexFilter from skimanalysis/HltTree" << std::endl;
      }
    }
  } else {
    std::cout << "WARNING: could not reopen input file for event-level friend trees; falling back to tTree friends" << std::endl;
    eventInfoFile = nullptr;
  }
  if (!cfg.dataset.isMC) weightTree = nullptr;

  // Tree variables
  double prescale = cfg.dataset.data_prescale;

  Long64_t n_events = t.GetEntries();
  if (ev_last < 0 || ev_last > n_events) ev_last = n_events;
  if (ev_first < 0) ev_first = 0;
  std::cout << "Processing events [" << ev_first << ", " << ev_last << ") of " << n_events << std::endl;

  for (Long64_t ient = ev_first; ient < ev_last; ient++) {
    // Keep batch logs compact; final counters are printed after the loop.

      t.GetEntry(ient);
      TrkEffSyst::beginEvent(t.run, t.lumi, t.evt);
      if (hiEvtTree) {
        hiEvtTree->GetEntry(ient);
        t.vz = hiEvtVz;
      }
      // JEC: replace jtpt by the corrected pT BEFORE any selection reads it. A jet whose
      // inputs are unusable gets -999, which every pT cut downstream rejects.
      if (Jec::enabled()) {
        const double eventRho = hiEvtRhoLeaf->GetValue();
        if (!std::isfinite(eventRho))
          throw std::runtime_error("Non-finite event rho encountered while applying JEC");
        for (Int_t ijet = 0; ijet < t.nref; ++ijet) {
          if (!std::isfinite(t.rawpt[ijet]) || t.rawpt[ijet] <= 0 ||
              !std::isfinite(t.jteta[ijet]) || !std::isfinite(t.jtarea[ijet])) {
            t.jtpt[ijet] = -999.0f;
            continue;
          }
          jec.SetJetPT(t.rawpt[ijet]);
          jec.SetJetEta(t.jteta[ijet]);
          jec.SetJetPhi(t.jtphi[ijet]);
          jec.SetJetArea(t.jtarea[ijet]);
          jec.SetRho(eventRho);
          const double correctedPt = jec.GetCorrectedPT();
          if (!std::isfinite(correctedPt) || correctedPt <= 0) {
            t.jtpt[ijet] = -999.0f;
            continue;
          }
          double variedPt = correctedPt;
          if (JesSyst::varied()) {
            jesUncertainty.SetJetPT(correctedPt);
            jesUncertainty.SetJetEta(t.jteta[ijet]);
            jesUncertainty.SetJetPhi(t.jtphi[ijet]);
            jesUncertainty.SetJetArea(t.jtarea[ijet]);
            jesUncertainty.SetRho(eventRho);
            const std::pair<double, double> uncertainty = jesUncertainty.GetUncertainty();
            if (!std::isfinite(uncertainty.first) || !std::isfinite(uncertainty.second)
                || uncertainty.first < 0 || uncertainty.second < 0)
              throw std::runtime_error("Invalid JES uncertainty for a reconstructed jet");
            variedPt = JesSyst::apply(correctedPt, uncertainty);
          }
          t.jtpt[ijet] = std::isfinite(variedPt) && variedPt > 0 ? variedPt : -999.0f;
        }
        if (JerSyst::enabled()) {
          const std::vector<int> genMatch = jerSmearer->Match(
              t.nref, t.jtpt, t.jteta, t.jtphi,
              t.ngen, t.genpt, t.geneta, t.genphi,
              0.4, eventRho, JetSmearing::Mode::JME);
          for (Int_t ijet = 0; ijet < t.nref; ++ijet) {
            if (!(t.jtpt[ijet] > 0) || !std::isfinite(t.jtpt[ijet])) continue;
            const int igen = genMatch[ijet];
            const double genPt = igen >= 0 ? t.genpt[igen] : -1.0;
            JetSmearing::Result smear = jerSmearer->Smear(
                t.jtpt[ijet], t.jteta[ijet], eventRho, genPt,
                static_cast<std::int64_t>(t.evt), JerSyst::variation());
            // Correlate the stochastic JER draw across JES variations.  The
            // correctionlib-compatible default hashes reco pT into the random
            // seed, so shifting JES would otherwise give the same jet a new
            // Gaussian draw.  Run/lumi/event and reco-jet index are stable under
            // these variations and distinguish jets within an event.
            if (!smear.matched) {
              const std::uint64_t eventKey = TrkEffSyst::eventKey();
              const std::uint64_t jetKey =
                  (static_cast<std::uint64_t>(ijet) + 1ULL)
                  * 0x9e3779b97f4a7c15ULL;
              const std::int64_t stableSeed =
                  static_cast<std::int64_t>(eventKey ^ jetKey
                                            ^ 0x4a45525f32303234ULL);
              const double gaussian =
                  JetSmearing::HashNormal(0.0, 0.0, 0.0, stableSeed);
              const double variance =
                  std::max(smear.scaleFactor * smear.scaleFactor - 1.0, 0.0);
              smear.smearFactor =
                  1.0 + smear.resolution * std::sqrt(variance) * gaussian;
            }
            const double smearedPt =
                JetSmearing::SmearedPt(t.jtpt[ijet], smear.smearFactor);
            t.jtpt[ijet] = std::isfinite(smearedPt) && smearedPt > 0
                               ? smearedPt : -999.0f;
          }
        }
      }
      if (hltTree && cfg.dataset.RunN == 3) {
        hltTree->GetEntry(ient);
        t.HLT_AK4PFJet40_v8 = hltAK4PFJet40;
        t.HLT_AK4PFJet60_v8 = hltAK4PFJet60;
        t.HLT_AK4PFJet80_v8 = hltAK4PFJet80;
        t.HLT_AK4PFJet100_v8 = hltAK4PFJet100;
        t.HLT_AK4PFJet120_v8 = hltAK4PFJet120;
      }
      if (skimTree && cfg.dataset.RunN == 3 && !cfg.dataset.isMC) {
        skimTree->GetEntry(ient);
        t.pprimaryVertexFilter = skimPprimaryVertexFilter;
      }

      // -- test tree branches are read:
        // cout << "jtpt = " << t.jtpt[0] << endl;

      double weight_tree = cfg.dataset.isMC ? (weightTree ? hiEvtWeight : t.weight) : 1.0;

    
      if (! passPVQuality_EventSelection(t, cfg)) continue;

        // cout << "PV quality pass ok: t.vz = "<< t.vz  << endl;
          // if (cfg.dataset.RunN == 3 &&  !isMC ) cout << " and t.pprimaryVertexFilter = "<< t.pprimaryVertexFilter << endl;

      // -- Trigger selections for all
      if (! passEventSelection(t, cfg)) continue;
      // cout << "event pass ok" << endl;

    for (Int_t ijet = 0; ijet < t.nref; ijet++) { // Jet loop

      // 2024 jet veto map (JetVeto at the top): true when off.
      const bool passJetVeto = JetVeto::pass(t.jteta[ijet], t.jtphi[ijet]);

      // Fill templates from DATA or RECO MC: jet kinematics + btagging

      if(isMakeTemplates && passJetVeto && passRecoJetKinematics(t, ijet, cfg) &&  passBtag(t, ijet, cfg)){

              // -- Fill here Selected jets histogram: selected jets after btagging + >=2svx (cut)
              if (t.jtNsvtx[ijet] >= 2){
                hpt_selectedJets ->Fill(t.jtpt[ijet], weight_tree);
                hpt_selectedJets_noweight ->Fill(t.jtpt[ijet]);
              }

              // reco SV reconstruction — same for data and MC
              // pt_nonb_reco: the jet's non-b track pT, filled by the same pass that builds
              // the SVs (it is the only place that knows which tracks were rejected).
              double pt_nonb_reco = 0.;
              std::vector<int> reco_sv_sta;
              vector<ROOT::Math::PtEtaPhiMVector> reco_sv = makeSvtxs_withBDT(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail, nullptr, nullptr, &pt_nonb_reco, &reco_sv_sta);

              if (reco_sv.size() == 2){ // #sv must = 2 (default of makeSvtxs_withBDT())

                // SV b scale factor, one per b-matched SV (sv_sf.h). 1 for data and when off.
                // Always from the NOMINAL vertices: they make the selection the SF corrects.
                double w_svsf = SvSF::jetWeight(reco_sv[0].Pt(), reco_sv_sta[0], reco_sv[1].Pt(), reco_sv_sta[1]);
                double jtpt = jtpt_fill(t.jtpt[ijet]);
                // The two B candidates of the observables: the nominal vertices, or with AGG_NN=on
                // the aggregation NN's two clusters (agg_nn.h). Selection unchanged either way.
                ROOT::Math::PtEtaPhiMVector cand0 = reco_sv[0], cand1 = reco_sv[1];
                if (AggNN::enabled()) {
                  const AggNN::Result nn = AggNN::run(t, ijet, ient, jtpt);
                  if (nn.status >= 0) {
                    cand0 = cand1 = ROOT::Math::PtEtaPhiMVector();
                    for (int k : nn.g1) cand0 += recoTrackP4(t, k);
                    for (int k : nn.g2) cand1 += recoTrackP4(t, k);
                  }
                }
                double dr   = t.calc_dr(cand0.Eta(), cand0.Phi(), cand1.Eta(), cand1.Phi());
                double pt1  = cand0.Pt();
                double pt2  = cand1.Pt();
                double eec  = EecWeight::value(pt1, pt2, cfg.n);
                double mB   = cand0.M() + cand1.M();
                if (svSeedNominalM2B()) mB = nominalVertexMass(t, ijet, ient, mB);   // 2026-11-06
               
                double Bbal = MomBalance::value(pt1, pt2);
                 (void)pt_nonb_reco;

                //Fix the under/overflow
                //if(dr < dr_min) dr = dr_min_fill;
                if(dr >= dr_max) dr = dr_max_fill;
                if(mB >= mb_max) mB = mb_max_fill;

                // Match the response-matrix reco definition: drop the dr<=0.005 spike so
                // templates and response cover the same reco domain. Data and all MC.
          
                if (obsDr().passReco(dr)) {

                // Prescale factor for data
                if (cfg.dataset.RunN == 2 && !cfg.dataset.isMC && t.HLT_HIAK4PFJet40_v1 && !(t.HLT_HIAK4PFJet60_v1 || t.HLT_HIAK4PFJet80_v1 || t.HLT_HIAK4PFJet100_v1))
                {eec *= prescale;}

                if (cfg.dataset.RunN == 3 && !cfg.dataset.isMC && t.HLT_AK4PFJet60_v8 && !(t.HLT_AK4PFJet80_v8 || t.HLT_AK4PFJet100_v8 || t.HLT_AK4PFJet120_v8) )
                {eec *= prescale;}

                // Same jet, same weight, one fill per observable. The values are in the
                // SAME ORDER as the `observables` vector -- see obsValues() for the contract.
                const std::vector<double> xreco = obsValues(dr, Bbal);
                for (size_t io = 0; io < tmpl.size(); ++io) {
                  if (xreco[io] < 0.) continue;   // unusable pair for this observable
                  tmpl[io].fill(cfg.dataset.isMC, t.jtNbHad[ijet], ient,
                                mB, xreco[io], jtpt, eec, weight_tree * w_svsf, t.jtpt[ijet]);
                }
                svDiag.fill(t, ijet, ient, cfg.dataset.isMC, t.jtNbHad[ijet], cand0, cand1, t.jtpt[ijet], eec * weight_tree * w_svsf);

                } // end if (dr > 0.005)
                }// Fill Templates
              } // end if 2 SV

    /////////---- To Prepare Response matrix (of true >=2B) ---- ONLY for MC (both RECO, GEN) ----
	          if((isCreateRmatrix || doAggNtuple) && cfg.dataset.isMC && t.jtNbHad[ijet] >= 2)  // -- select jets of 2b (truth)
          {

            // -- common variables repeatdly used in fill histograms
            // jtpt_fill() folds jets above the last bin edge into that bin, so they are
            // not lost to the overflow bin (see binning_histos_small.h).
            double jpt_reco = jtpt_fill(reco_jet_pt(t, ijet));
            // double jpt_gen = jtpt_fill(gen_jet_pt(t, ijet));   // until 2026-10-29
            double jpt_gen = jtptGenFill(gen_jet_pt(t, ijet));   // = jtpt_fill unless JET_PT_MAX (gen overflow)
            double jeta_reco = reco_jet_eta(t, ijet);
            double jeta_gen = gen_jet_eta(t, ijet);
            // -- for cout only
            double btagVar =  (cfg.dataset.RunN == 2) ? (t.discr_particleNet_BvsAll[ijet]) :
                              ( (cfg.dataset.RunN == 3) ?
                                  (t.discr_unifiedParticleTransformer_probb[ijet] +
                                  t.discr_unifiedParticleTransformer_problepb[ijet] +
                                  t.discr_unifiedParticleTransformer_probbb[ijet]) :
	                                  -1);

	            // -- counts stats
	            n_bb_jets++; // true >=2b jets

	            aggRow.reset();
	            aggRow.entry = ient; aggRow.run = t.run; aggRow.lumi = t.lumi; aggRow.evt = t.evt; aggRow.jetIndex = ijet;
	            aggRow.weight = weight_tree; aggRow.jtNbHad = t.jtNbHad[ijet]; aggRow.btagScore = btagVar;
	            aggRow.jtpt = jpt_reco; aggRow.jteta = jeta_reco; aggRow.refpt = jpt_gen; aggRow.refeta = jeta_gen;
	            aggRow.passRecoKin = passRecoJetKinematics(t, ijet, cfg);
	            aggRow.passGenKin = passGenJetKinematics(t, ijet, cfg);
	            aggRow.passBtag = passBtag(t, ijet, cfg);
	            aggRow.passJetVeto = passJetVeto;

	                // step1: for Response matrix ---- Gen b hadrons ----
	                std::vector<ROOT::Math::PtEtaPhiMVector> gen_bh;
	                std::vector<Int_t> gen_bh_sta;
	                PartialBsAggregation(gen_bh, gen_bh_sta, t, ijet);
	                aggRow.nGenAgg = gen_bh.size();
	                // MEASUREMENT DEFINITION: a reco jet with < 2 aggregatable gen B-hadrons
	                // is a fake and must stay in the purity DENOMINATOR, so do NOT drop it.
	                // gen_ok only guards the gen-pair math and the gen-level fills.
	                const bool gen_ok = (gen_bh.size() >= 2);
	                if (gen_ok) n_gen_bh_ok++;

	                // gen-pair observables: assigned only when gen_ok
	                double eec_gen = -1., mB_gen = -1., dr_gen = -1.;
	                double mB_gen_fill = -1., dr_gen_fill = -1.;
	                // The second observable at gen level, from the same best-EEC pair.
	                double B_gen = -1.;
	                // DISABLED 2026-09-22 -- the gen b side is that same pair, the gen non-b
	                // side every other gen track of the jet above 1 GeV:
	                //   double fnb_gen = -1., fb_gen = -1., lnfb_gen = -1.;
	                if (gen_ok) {

                // From aggregated gen B: Pick gen pair with largest EEC weight (pt_i * pt_j)^n
                int best_i = 0, best_j = 1;
                double best_pt_prod = -1;
                for (size_t gi = 0; gi < gen_bh.size(); gi++)
                    for (size_t gj = gi+1; gj < gen_bh.size(); gj++) {
                        double pp = gen_bh[gi].Pt() * gen_bh[gj].Pt();
                        if (pp > best_pt_prod) { best_pt_prod = pp; best_i = gi; best_j = gj; }
                }
					// gen pair info
	                eec_gen = EecWeight::value(gen_bh[best_i].Pt(), gen_bh[best_j].Pt(), cfg.n);
	                mB_gen  = gen_bh[best_i].M() + gen_bh[best_j].M();
	                dr_gen  = t.calc_dr(gen_bh[best_i].Eta(), gen_bh[best_i].Phi(),
	                                           gen_bh[best_j].Eta(), gen_bh[best_j].Phi());
	                B_gen   = MomBalance::value(gen_bh[best_i].Pt(), gen_bh[best_j].Pt());
	                // DISABLED 2026-09-22 with the fraction observables. genNonBTrackPt() is
	                // the only caller of that gen-track loop, so nothing walks the gen tracks
	                // while these three are commented out.
	                // {
	                //   const double pt_b_gen    = gen_bh[best_i].Pt() + gen_bh[best_j].Pt();
	                //   const double pt_nonb_gen = genNonBTrackPt(t, ijet);
	                //   fnb_gen  = NonBFraction::value(pt_b_gen, pt_nonb_gen);
	                //   fb_gen   = BFraction::value(pt_b_gen, pt_nonb_gen);
	                //   lnfb_gen = LogBFraction::value(pt_b_gen, pt_nonb_gen);
	                // }
	                aggRow.genStatus1 = gen_bh_sta[best_i]; aggRow.genStatus2 = gen_bh_sta[best_j];
	                aggRow.setGen(gen_bh[best_i], gen_bh[best_j]);
	                aggRow.genDr = dr_gen; aggRow.genMB = mB_gen; aggRow.genEec = eec_gen;
					// overflow treatement at gen level
				    mB_gen_fill = mB_gen; dr_gen_fill = dr_gen;
                if (mB_gen_fill >= mb_max) mB_gen_fill = mb_max_fill;
                if (dr_gen_fill >= dr_max)  dr_gen_fill = dr_max_fill;
                } // end if (gen_ok): gen-pair observables ready

			  // -- Prepare combined b-tagger:
        // ---- Reco SVs ----
        // Same call as the template block above, and the non-b pT comes back the same way.
        double pt_nonb_reco_rm = 0.;
        std::vector<int> reco_sv_sta_rm;
        vector<ROOT::Math::PtEtaPhiMVector> reco_sv_rm = makeSvtxs_withBDT(t, ijet, ient, agg_fail_rm, nb_sv_rm, sv_fail_rm, merge_fail_rm, nullptr, nullptr, &pt_nonb_reco_rm, &reco_sv_sta_rm);

	      bool reco_sv_ok = (reco_sv_rm.size() == 2);
	      aggRow.nRecoAgg = reco_sv_rm.size();
	      if (reco_sv_ok) {
	        const double mB_reco_diag = reco_sv_rm[0].M() + reco_sv_rm[1].M();
	        const double dr_reco_diag = t.calc_dr(reco_sv_rm[0].Eta(), reco_sv_rm[0].Phi(),
	                                                        reco_sv_rm[1].Eta(), reco_sv_rm[1].Phi());
	        const double eec_reco_diag =
	            EecWeight::value(reco_sv_rm[0].Pt(), reco_sv_rm[1].Pt(), cfg.n);
	        aggRow.setReco(reco_sv_rm[0], reco_sv_rm[1]);
	        aggRow.recoDr = dr_reco_diag; aggRow.recoMB = mB_reco_diag; aggRow.recoEec = eec_reco_diag;
	      }
	      // Aggregation NN (AGG_NN=on, agg_nn.h): the two candidates the observables use below, for
	      // the jets the nominal aggregation selects. Kept in the ntuple (cx* branches) as well.
	      ROOT::Math::PtEtaPhiMVector cand0_rm, cand1_rm;
	      if (reco_sv_ok) { cand0_rm = reco_sv_rm[0]; cand1_rm = reco_sv_rm[1]; }
	      if (reco_sv_ok && AggNN::enabled()) {
	        const AggNN::Result nn = AggNN::run(t, ijet, ient, jpt_reco);
	        aggRow.cxStatus = nn.status;
	        if (nn.status >= 0) {
	          cand0_rm = cand1_rm = ROOT::Math::PtEtaPhiMVector();
	          for (int k : nn.g1) cand0_rm += recoTrackP4(t, k);
	          for (int k : nn.g2) cand1_rm += recoTrackP4(t, k);
	          aggRow.cxPt1 = cand0_rm.Pt(); aggRow.cxEta1 = cand0_rm.Eta(); aggRow.cxPhi1 = cand0_rm.Phi();
	          aggRow.cxPt2 = cand1_rm.Pt(); aggRow.cxEta2 = cand1_rm.Eta(); aggRow.cxPhi2 = cand1_rm.Phi();
	        }
	      }
	    // Reconstruction-ceiling study (2026-10-21): the same jet's reco tracks aggregated two
	    // other ways (see the ts*/or* branches). Ntuple only -- no histogram reads these.
	    if (aggBHadronTree) {
	      double ts_af = 0., ts_nb = 0., ts_sf = 0., ts_mf = 0.;
	      const std::vector<ROOT::Math::PtEtaPhiMVector> ts = makeSvtxs(t, ijet, ient, ts_af, ts_nb, ts_sf, ts_mf);
	      if (ts.size() == 2) {
	        aggRow.tsPt1 = ts[0].Pt(); aggRow.tsEta1 = ts[0].Eta(); aggRow.tsPhi1 = ts[0].Phi();
	        aggRow.tsPt2 = ts[1].Pt(); aggRow.tsEta2 = ts[1].Eta(); aggRow.tsPhi2 = ts[1].Phi();
	      }
	      if (gen_ok) {
	        std::vector<ROOT::Math::PtEtaPhiMVector> orc; std::vector<Int_t> orc_sta;
	        MatchingTracksAggregation(orc, orc_sta, t, ijet, ient);
	        aggRow.orPt1 = aggRow.orPt2 = 0.;
	        for (size_t k = 0; k < orc.size(); ++k) {
	          if (orc_sta[k] == aggRow.genStatus1) { aggRow.orPt1 = orc[k].Pt(); aggRow.orEta1 = orc[k].Eta(); aggRow.orPhi1 = orc[k].Phi(); }
	          if (orc_sta[k] == aggRow.genStatus2) { aggRow.orPt2 = orc[k].Pt(); aggRow.orEta2 = orc[k].Eta(); aggRow.orPhi2 = orc[k].Phi(); }
	        }
	      }
	      // Aggregation-NN training set: the per-track arrays (see AggBHadronNtupleRow).
	      aggRow.jtphi = t.jtphi[ijet];
	      for (Int_t itrk = 0; itrk < t.ntrk; ++itrk) {
	        if (t.trkJetId[itrk] != ijet || t.trkPt[itrk] < 1) continue;
	        const Int_t sv = t.trkSvtxId[itrk];
	        const bool in_sv = (sv >= 0 && sv < t.nsvtx);
	        aggRow.trk_pt.push_back(t.trkPt[itrk]); aggRow.trk_eta.push_back(t.trkEta[itrk]); aggRow.trk_phi.push_back(t.trkPhi[itrk]);
	        aggRow.trk_ip3d.push_back(t.trkIp3d[itrk]); aggRow.trk_ip3dsig.push_back(t.trkIp3dSig[itrk]); aggRow.trk_dist.push_back(t.trkDistToAxis[itrk]);
	        aggRow.trk_bdt.push_back(t.trkBdtScore[itrk]); aggRow.trk_pdg.push_back(t.trkPdgId[itrk]);
	        aggRow.trk_nomsel.push_back(TrackBdtSF::passes(t.trkBdtScore[itrk], t.trkMatchSta[itrk], sv) ? 1 : 0);
	        aggRow.trk_svidx.push_back(in_sv ? sv : -1);
	        aggRow.trk_sv_ntrk.push_back(in_sv ? t.svtxNtrk[sv] : -1);
	        aggRow.trk_sv_m.push_back(in_sv ? t.svtxm[sv] : -999.f);       aggRow.trk_sv_mcorr.push_back(in_sv ? t.svtxmcorr[sv] : -999.f);
	        aggRow.trk_sv_pt.push_back(in_sv ? t.svtxpt[sv] : -999.f);     aggRow.trk_sv_dl.push_back(in_sv ? t.svtxdl[sv] : -999.f);
	        aggRow.trk_sv_dls.push_back(in_sv ? t.svtxdls[sv] : -999.f);   aggRow.trk_sv_dl2d.push_back(in_sv ? t.svtxdl2d[sv] : -999.f);
	        aggRow.trk_sv_dls2d.push_back(in_sv ? t.svtxdls2d[sv] : -999.f); aggRow.trk_sv_chi2.push_back(in_sv ? t.svtxnormchi2[sv] : -999.f);
	        const Int_t sta = t.trkMatchSta[itrk];
	        aggRow.truth_sta.push_back(sta);
	        aggRow.truth_label.push_back(sta < 100 ? 0 : (gen_ok && sta == aggRow.genStatus1) ? 1 : (gen_ok && sta == aggRow.genStatus2) ? 2 : 3);
	      }
	    }
	    if (aggBHadronTree) aggBHadronTree->Fill();
	    if (!isCreateRmatrix) continue;

	    // Combined SV-reco + b-tag efficiency: DENOMINATOR = all true 2b jets passing the gen
	    // cuts (gen_kin), before the 2-SV cut. Gen-EEC weighted (w_gen), reused for the numerator below.
	    double w_gen = gen_ok ? weight_tree * eec_gen : 0.0;
	    // Gen values in observable order. Reused by every hgenjet_* fill below.
	    const std::vector<double> xgen_raw = obsValues(dr_gen, B_gen);
	    // The hgenjet_* fills (SV+b-tag, b-tag and EEC-weight corrections) use the same gen
	    // cuts as gen_pass below (pT_gen, |eta_gen|, pthat cleaning): they are applied after
	    // h_full_efficiency, whose result is already restricted to that phase space.
	    // Changed 2026-10-26; before, these three fills were gated on gen_ok only.
	    // const bool gen_kin = gen_ok && passGenJetKinematics(t, ijet, cfg);   // until 2026-10-30
    const bool gen_kin = gen_ok && passGenJetKinematics(t, ijet, cfg) && passDrMinGen(dr_gen);   // + DR_MIN, as gen_pass
	    // if (gen_ok)
	    if (gen_kin)
	      for (size_t io = 0; io < resp.size(); ++io)
	        if (xgen_raw[io] >= 0.) resp[io].fillGenAll(xgen_raw[io], jpt_gen, w_gen);
	    // 3D: raw pT (no jtpt_fill), gen threshold 60, and the jet needs BOTH dR and A.
	    const double pt_gen3d = gen_jet_pt(t, ijet), pt_reco3d = reco_jet_pt(t, ijet);
	    const bool gen_kin3d = gen_ok && passGenJetKinematics(t, ijet, cfg3d) && dr_gen >= 0. && B_gen >= 0.;
	    if (gen_kin3d) resp3d.fillGenAll(dr_gen, B_gen, pt_gen3d, w_gen);

	    // Guard the reco pair access: makeSvtxs_withBDT returns an empty vector when the
	    // jet has < 2 secondary vertices, and reco_sv_rm[0] below would then segfault.
	    if (!reco_sv_ok) continue;

      // -- Prepare to Fill Rmatrix --
		  // -- Reco pair info
      // The candidates: the nominal vertices, or the NN's two clusters with AGG_NN=on (above).
      // double mB_reco  = reco_sv_rm[0].M() + reco_sv_rm[1].M();
      // double dr_reco  = t.calc_dr(reco_sv_rm[0].Eta(), reco_sv_rm[0].Phi(), reco_sv_rm[1].Eta(), reco_sv_rm[1].Phi());
      // double eec_reco = EecWeight::value(reco_sv_rm[0].Pt(), reco_sv_rm[1].Pt(), cfg.n);
      // double B_reco = MomBalance::value(reco_sv_rm[0].Pt(), reco_sv_rm[1].Pt());
      double mB_reco  = cand0_rm.M() + cand1_rm.M();
      if (svSeedNominalM2B()) mB_reco = nominalVertexMass(t, ijet, ient, mB_reco);   // 2026-11-06, as the templates
      double dr_reco  = t.calc_dr(cand0_rm.Eta(), cand0_rm.Phi(), cand1_rm.Eta(), cand1_rm.Phi());
      double eec_reco = EecWeight::value(cand0_rm.Pt(), cand1_rm.Pt(), cfg.n);
      double B_reco = MomBalance::value(cand0_rm.Pt(), cand1_rm.Pt());
      // DISABLED 2026-09-22 with the fraction observables:
      //   double fnb_reco  = NonBFraction::value(reco_sv_rm[0].Pt() + reco_sv_rm[1].Pt(),
      //                                          pt_nonb_reco_rm);
      //   double fb_reco   = BFraction::value(reco_sv_rm[0].Pt() + reco_sv_rm[1].Pt(),
      //                                       pt_nonb_reco_rm);
      //   double lnfb_reco = LogBFraction::value(reco_sv_rm[0].Pt() + reco_sv_rm[1].Pt(),
      //                                          pt_nonb_reco_rm);
      (void)pt_nonb_reco_rm;
      // SV b scale factor (sv_sf.h): every fill from here on is made from these two reco SVs,
      // so it carries the SF -- including w_gen in the SV+b-tag efficiency NUMERATOR below.
      // The all-true-2b denominator (fillGenAll, above) is before the SVs and stays unweighted.
      double w_svsf = SvSF::jetWeight(reco_sv_rm[0].Pt(), reco_sv_sta_rm[0],
                                      reco_sv_rm[1].Pt(), reco_sv_sta_rm[1]);
		  double w_reco = weight_tree * eec_reco * w_svsf;

      // Fill total number of True 2b Jets (Deno of b-tagging eff. correction) before b-tagger: simialr axis to what we unfold to. - used after unfolding
			// if (gen_ok)
			if (gen_kin)
			  for (size_t io = 0; io < resp.size(); ++io)
			    if (xgen_raw[io] >= 0.) resp[io].fillGenPreBtag(xgen_raw[io], jpt_gen, w_reco);
			if (gen_kin3d) resp3d.fillGenPreBtag(dr_gen, B_gen, pt_gen3d, w_reco);

			// -- Use combined b-tagger: Upart tagger + 2SV of aggreagted Bs., for Rmatrix unfolding, purity, reconstruction eff.
			if (!passBtag(t, ijet, cfg)) continue; // select btagged jets only (reco_sv_ok already required above)

			// Fill total number of True 2b Jets that survive after btagging condition (Num of b tagging eff. correction - used after unfolding)
			// and the combined SV-reco + b-tag efficiency NUMERATOR = 2b jets reconstructed
			// with 2 SV AND tagged (this point is past both cuts). Same gen-EEC weight as
			// hgenjet_2b_all.
			// if (gen_ok)
			if (gen_kin)
			  for (size_t io = 0; io < resp.size(); ++io)
			    if (xgen_raw[io] >= 0.) resp[io].fillGenPostBtag(xgen_raw[io], jpt_gen, w_reco, w_gen * w_svsf);
			if (gen_kin3d) resp3d.fillGenPostBtag(dr_gen, B_gen, pt_gen3d, w_reco, w_gen * w_svsf);

			// DEBUG: reco SV pt vs gen partial-B pt (pt-rank matched), to explain the EEC-weight factor.
			if (gen_ok) {
			  double g_hi = std::max(aggRow.genPt1,  aggRow.genPt2),  g_lo = std::min(aggRow.genPt1,  aggRow.genPt2);
			  double r_hi = std::max(aggRow.recoPt1, aggRow.recoPt2), r_lo = std::min(aggRow.recoPt1, aggRow.recoPt2);
			  h_svpt_vs_bpt->Fill(g_hi, r_hi, weight_tree);
			  h_svpt_vs_bpt->Fill(g_lo, r_lo, weight_tree);
			  if (g_hi > 0) h_ptratio->Fill(r_hi / g_hi, weight_tree);
			  if (g_lo > 0) h_ptratio->Fill(r_lo / g_lo, weight_tree);
			}

	      // Overflow protection at reco level
	              double mB_reco_fill = mB_reco, dr_reco_fill = dr_reco;
                  if (mB_reco_fill >= mb_max) mB_reco_fill = mb_max_fill;
                  if (dr_reco_fill >= dr_max) dr_reco_fill = dr_max_fill;

            // Define reco_pass: full detector-level selection
            // bool reco_pass = passJetVeto && passRecoJetKinematics(t, ijet, cfg) && dr_reco_fill > 0.005;   // until 2026-10-30
            // obsDr().passReco: dr > 0.005, or dr > DR_MIN when set (the same gate as the templates)
            bool reco_pass = passJetVeto && passRecoJetKinematics(t, ijet, cfg) && obsDr().passReco(dr_reco_fill);

            // Define gen_pass: particle-level jet kinematics + gen observable range.
            // Requires gen_ok, so a jet with < 2 gen B-hadrons stays out of num/eff/response.
            // DR_MIN (2026-10-30): the dR cut at gen too (Zoe: cut at reco and gen), for every observable.
            bool gen_pass  = gen_ok && passGenJetKinematics(t, ijet, cfg) && passDrMinGen(dr_gen);

            ++n_weight_checked;
            sum_weight_tree_checked += weight_tree;
            //sum_w_gen_checked += w_gen;
            sum_w_reco_checked += w_reco;
            if (weight_tree == 0.0) ++n_zero_weight_tree;
            if (eec_gen == 0.0) ++n_zero_eec_gen;
            if (eec_reco == 0.0) ++n_zero_eec_reco;
            //if (w_gen == 0.0) ++n_zero_w_gen;
            if (w_reco == 0.0) ++n_zero_w_reco;
            if (n_weight_debug_printed < 10) {
                std::cout << "[WEIGHT DBG " << n_weight_debug_printed << "]"
                    << " entry=" << ient
                    << " jet=" << ijet
                    << " cfg.isMC=" << cfg.dataset.isMC
                    << " sample=" << cfg.dataset.sample
                    << " raw t.weight=" << t.weight
                    << " weight_tree=" << weight_tree
                    << " eec_gen=" << eec_gen
                    << " eec_reco=" << eec_reco
                    //<< " w_gen=" << w_gen
                    << " w_reco=" << w_reco
                    << " reco_pass=" << reco_pass
                    << " gen_pass=" << gen_pass
                    << std::endl;
                ++n_weight_debug_printed;
            }

            // --- debugging paragraph ----
            //std::cout << "gen_pass: " << gen_pass << std::endl;
            //std::cout << "jpt_gen: " << jpt_gen << std::endl;
            //std::cout << "pT_low: " << pT_low << ", pT_high: " << pT_high << std::endl;
            //std::cout << "refeta_gen: " << t.refeta[ijet] << std::endl;
            //std::cout << "mB_gen: " << mB_gen << std::endl;
            //std::cout << "mB_low: " << mb_min << ", mB_high: " << mb_max << std::endl;
            //std::cout << "dr_gen: " << dr_gen << std::endl;
            // --- per-condition failure tallies (debug only) ---
            if (reco_pass) n_reco_pass++;
            if (gen_pass)  n_gen_pass++;
            if (reco_pass && gen_pass) n_both_pass++;
            if (!reco_sv_ok)                                                                              n_fail_reco_sv++;
            else if (!( jpt_reco >= cfg.kin.ptLow ))                         n_fail_reco_pt++;
            else if (!(std::abs(jeta_reco) < cfg.kin.etaMax))                                             n_fail_reco_eta++;
            else if (cfg.physics.useBtag && !passBtag(t, ijet, cfg))                                      n_fail_reco_btag++;
            else if (!(mB_reco_fill >= mb_min && mB_reco_fill < mb_max))                                  n_fail_reco_mb++;
            else if (!(dr_reco_fill < dr_max))                                                            n_fail_reco_dr++;
            if (!(jpt_gen >= cfg.kin.ptLow))                                   n_fail_gen_pt++;
            else if (!(std::abs(jeta_gen) < cfg.kin.etaMax))                                               n_fail_gen_eta++;
            else if (!(mB_gen_fill >= mb_min && mB_gen_fill < mb_max))                                     n_fail_gen_mb++;
            else if (!(dr_gen_fill < dr_max))                                                              n_fail_gen_dr++;
            // --- debug: print first 10 jets that reach cut evaluation ---
            if (n_debug_printed < 10) {
                std::cout << "[DBG jet " << n_debug_printed << "]"
                    << " jpt_reco=" << jpt_reco
                    << " jpt_gen="  <<  t.refpt[ijet]  // jpt_gen
                    << " reco jteta="    << jeta_reco
                    << " gen jteta="   << t.refeta[ijet] //  jeta_gen
                    << " discr="    << btagVar
                    << " mB_reco="  << mB_reco
                    << " dr_reco="  << dr_reco
                    << " mB_gen="   << mB_gen
                    << " dr_gen="   << dr_gen
                    << " reco_sv="  << reco_sv_rm.size()
                    << " gen_bh="   << gen_bh.size()
                    << " reco_pass=" << reco_pass
                    << " gen_pass="  << gen_pass
                    << std::endl;
                n_debug_printed++;
            }//debug if

            // ---- Fill response matrix ----
            double num    = TrkEffSyst::jetSplitRandom(ient, ijet);

            // One response set per observable, same jet, same weight, same half.
            // reco_pass/gen_pass are the shared jet-level gates computed above; the only
            // per-observable addition is that the jet must HAVE a usable value for this
            // observable (B is undefined if a B hadron has non-positive pT).
            const std::vector<double> xreco_raw = obsValues(dr_reco, B_reco);
            for (size_t io = 0; io < resp.size(); ++io) {
              const ObsDef &o = resp[io].o;
              const bool has_reco = (xreco_raw[io] >= 0.);
              const bool has_gen  = (xgen_raw[io]  >= 0.);
              resp[io].fillResponse(reco_pass && has_reco, gen_pass && has_gen, num, ient,
                                    has_reco ? o.fill(xreco_raw[io]) : 0., jpt_reco,
                                    has_gen  ? o.fill(xgen_raw[io])  : 0., jpt_gen, w_reco);
            }
            // 3D: same jet, weight and halves; reco must also sit inside the reco hard edges.
            const bool reco_pass3d = reco_pass && B_reco >= 0. && u3d_reco_in(dr_reco, pt_reco3d);
            resp3d.fillResponse(reco_pass3d, gen_kin3d, ient,
                                dr_reco, B_reco, pt_reco3d, dr_gen, B_gen, pt_gen3d, w_reco);

          } // END if MC()

      }// JET LOOP

	  } // EVENT LOOP
	  std::cout << std::endl;

	  if (aggBHadronTree && fout_agg) { fout_agg->cd(); aggBHadronTree->Write(); fout_agg->Close(); delete fout_agg; }
	  if (eventInfoFile) { eventInfoFile->Close(); delete eventInfoFile; }

	  // Response-matrix efficiency and purity, then write the output file
  if(isCreateRmatrix && cfg.dataset.isMC)
  {
  // -- Debug output ----------
    std::cout << std::endl;
    std::cout << "--- Jet statistics (bb jets, jtNbHad >= 2) ---" << std::endl;
    std::cout << "  selected bb jets (after triggers):      " << n_bb_jets    << std::endl;
    std::cout << "  Gen bh ok (>= 2 gen bh):     " << n_gen_bh_ok  << std::endl;
    std::cout << "  Passing reco cuts:           " << n_reco_pass  << std::endl;
    std::cout << "  Passing gen cuts:            " << n_gen_pass   << std::endl;
    std::cout << "  Passing both (numerator):    " << n_both_pass  << std::endl;
    std::cout << "  Reco SV failures (< 2 SVs): " << nb_sv_rm        << std::endl;
    std::cout << "  SV purity failures:          " << sv_fail_rm      << std::endl;
    std::cout << "  SV merging failures:         " << merge_fail_rm   << std::endl;
    std::cout << "  No-SV track agg failures:    " << agg_fail_rm     << std::endl;
    std::cout << "--- Reco-pass per-condition failures (first failing cond shown) ---" << std::endl;
    std::cout << "  fail reco_sv_ok:   " << n_fail_reco_sv   << std::endl;
    std::cout << "  fail reco jpt:     " << n_fail_reco_pt   << std::endl;
    std::cout << "  fail reco eta:     " << n_fail_reco_eta  << std::endl;
    std::cout << "  fail reco btag:    " << n_fail_reco_btag << std::endl;
    std::cout << "  fail reco mB:      " << n_fail_reco_mb   << std::endl;
    std::cout << "  fail reco dr:      " << n_fail_reco_dr   << std::endl;
    std::cout << "--- Gen-pass per-condition failures (first failing cond shown) ---" << std::endl;
    std::cout << "  fail gen jpt:      " << n_fail_gen_pt    << std::endl;
    std::cout << "  fail gen eta:      " << n_fail_gen_eta   << std::endl;
    std::cout << "  fail gen mB:       " << n_fail_gen_mb    << std::endl;
    std::cout << "  fail gen dr:       " << n_fail_gen_dr    << std::endl;
    std::cout << "--- Response-matrix weight diagnostics ---" << std::endl;
    std::cout << "  checked selected reco-aggregate jets: " << n_weight_checked << std::endl;
    std::cout << "  zero weight_tree: " << n_zero_weight_tree << std::endl;
    std::cout << "  zero eec_gen:     " << n_zero_eec_gen << std::endl;
    std::cout << "  zero eec_reco:    " << n_zero_eec_reco << std::endl;
    //std::cout << "  zero w_gen:       " << n_zero_w_gen << std::endl;
    std::cout << "  zero w_reco:      " << n_zero_w_reco << std::endl;
    std::cout << "  sum weight_tree:  " << sum_weight_tree_checked << std::endl;
    //std::cout << "  sum w_gen:        " << sum_w_gen_checked << std::endl;
    std::cout << "  sum w_reco:       " << sum_w_reco_checked << std::endl;

    // --- Compute purity and efficiency, then write, once per observable ---
    // ResponseHistos::write() derives the half/full purity and efficiency ratios, the
    // b-tag and SV+b-tag efficiency corrections and the EEC-weight correction, and writes
    // them under this observable's names. dr's suffix is "", so the dr objects land in the
    // file with exactly the names they have always had.
    std::cout << "Creating: " << ResponseMatrix_fout_name << std::endl;
    TFile *fout_rm = new TFile(ResponseMatrix_fout_name, "recreate");
    for (const ResponseHistos &r : resp) r.write();
    resp3d.write();

    // Observable-independent debug histograms, written once.
    std::cout << "[EEC-weight debug] mean reco/gen SV pt = " << h_ptratio->GetMean() << std::endl;
    h_svpt_vs_bpt->Write(); h_ptratio->Write();
    fout_rm->Close(); delete fout_rm;
  } // end if MC () after event loop  -- For Response matrix related hists.

// Write templates to the output root file
if(isMakeTemplates)
{
  std::cout << "Creating: " << fout_name << std::endl;
  TFile* outFile = new TFile(fout_name, "RECREATE");
  hpt_selectedJets->Write(); hpt_selectedJets_noweight->Write();

  // One template set per observable; dr's names are unchanged, B's carry "_B".
  for (const TemplateHistos &h : tmpl) h.write(cfg.dataset.isMC);
  svDiag.write();

  outFile->Close(); delete outFile;
}

} // end Function

#if 0
void filter_b_bb_as_data_and_mc(const AnalysisConfig& cfg) {
  TString output_hist =  Form("SPLIT_TEST_templatefit_n%d_bjet_Run%d", 1, cfg.dataset.RunN);

  double prescale = cfg.dataset.data_prescale;
  TString fout_name = cfg.dataset.output_folder + cfg.dataset.output_hist + "MCGEN.root"; // for reposnse matrix: has Prefix: Response

  // true level information aggregated to partial Bs
  tTree t;
  t.Init(cfg.dataset.filename, cfg.dataset.isMC, cfg.dataset.RunN);
  t.SetBranchStatus("*", 0);
  auto active_branches = getActiveBranches(cfg);
  t.SetBranchStatus(active_branches, 1);

  double agg_fail = 0;
  double nb_sv=0 ;
  double lt2sv = 0 ;
  double totjt = 0;
  double beforeskip = 0 ;
  double vertskip = 0 ;
  double matchskip = 0;
  double merge_fail = 0 ;
  double sv_fail = 0 ;
  double global_count = 0 ;

  // Plots or histograms
  TH3D *h3D_bb_bkg = new TH3D("h3D_bb_bkg", "#DeltaR;EEC",bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);
  TH3D *h3D_b_bkg = new TH3D("h3D_b_bkg", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);
  TH3D *h3D_bb_signal = new TH3D("h3D_bb_signal", "#DeltaR;EEC",bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);
  TH3D *h3D_b_signal = new TH3D("h3D_b_signal", "#DeltaR;EEC", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);

  TH2F* global_histo = new TH2F("global_histo", "distribution of wrong reco with pt; pT; reco proportion", 40, 0, 180, 40, 0, 1.2);

  TH1D *h_0b_score = new TH1D("h_0b_score", "Distribution of Score", 50, -1.5, 1.5);
  TH1D *h_1b_score = new TH1D("h_1b_score", "Distribution of Score", 50, -1.5, 1.5);

  h3D_bb_bkg->Sumw2();
  h3D_b_bkg->Sumw2();
  h3D_bb_signal->Sumw2();
  h3D_b_signal->Sumw2();

  TRandom3 rnd(12345);
  Long64_t n_events = t.GetEntries();
  for (Long64_t ient = 0; ient < n_events; ient++) {
    //Progress
    if (ient % 50000 == 0) {
      float percent = 100.0 * ient / n_events;
      std::cout << "\rProcessing: "  << percent << " %" << std::flush;
    }
    t.GetEntry(ient);
    double weight_tree = cfg.dataset.isMC ? t.weight : 1.0;

    bool isTemplate = (rnd.Uniform() < 0.5);

    // Loop over jets — identical logic for data and MC
    for (Int_t ijet = 0; ijet < t.nref; ijet++) {

      if (std::abs(t.jteta[ijet]) > cfg.kin.etaMax) continue;
      if (cfg.dataset.isMC && skipMC_event(t.jtpt[ijet], t.pthat)) continue;
      if (t.jtpt[ijet] < cfg.kin.ptLow) continue;
      /* ---- disabled (kept for reference): upper pT cut ----
      if (!cfg.kin.foldPtOverflow && t.jtpt[ijet] > cfg.kin.ptHigh) continue;
      ---- */
      if (cfg.physics.useBtag && !passBtag(t, ijet, cfg)) continue;

      // classify by number of reconstructed secondary vertices (same for data and MC)
      int nSV = t.jtNsvtx[ijet];
      if (nSV < 1) continue;

      beforeskip += 1;
      vector<ROOT::Math::PtEtaPhiMVector> reco_hadrons_4vec = makeSvtxs(t, ijet, ient, agg_fail, nb_sv, sv_fail, merge_fail);
      totjt += 1;
      if (reco_hadrons_4vec.empty()) lt2sv += 1;
      if (reco_hadrons_4vec.size() != 2) continue;

      double dr   = t.calc_dr(reco_hadrons_4vec[0].Eta(), reco_hadrons_4vec[0].Phi(), reco_hadrons_4vec[1].Eta(), reco_hadrons_4vec[1].Phi());
      double pt1  = reco_hadrons_4vec[0].Pt();
      double pt2  = reco_hadrons_4vec[1].Pt();
      double eec  = EecWeight::value(pt1, pt2, cfg.n);
      double jtpt = jtpt_fill(t.jtpt[ijet]);
      double mB   = reco_hadrons_4vec[0].M() + reco_hadrons_4vec[1].M();

      if (nSV == 1) {
        if (!isTemplate) h3D_b_signal->Fill(mB, dr, jtpt, eec * weight_tree);
        else             h3D_b_bkg->Fill(mB, dr, jtpt, eec * weight_tree);
      } else {  // nSV >= 2
        if (!isTemplate) h3D_bb_signal->Fill(mB, dr, jtpt, eec * weight_tree);
        else             h3D_bb_bkg->Fill(mB, dr, jtpt, eec * weight_tree);
      }

  } //close jet loop
  } //close event loop

  // -------- Save histograms -------
  TString label = "_btag";
  if(!cfg.physics.useBtag) label = "_nobtag";

  TFile outFile( (cfg.dataset.output_folder + output_hist + label + cfg.dataset.domain).Data(), "RECREATE");
  h3D_b_signal->Write();
  h3D_b_bkg->Write();
  h3D_bb_signal->Write();
  h3D_bb_bkg->Write();

  h_0b_score->Write();
  h_1b_score->Write();

  outFile.Close();

}
#endif 

void create_files_for_template_fit(Int_t RunN = 3, Float_t pT_low = 80, Float_t etaCut = 2, Int_t n = 1,
                                   bool btag = true, bool isMC = true, Double_t btagWP = 0.712, bool makeTemplates = true, bool createRmatrix = true,
                                   bool makeAggNtuple = true, Long64_t ev_first = 0, Long64_t ev_last = -1, const char* inputFileOverride = "", const char* outputFolderOverride = "",
                                   const char* sampleTag = "",
                                   bool track_eff_unc = false, ULong64_t trkSeed = 20260908,
                                   bool eec_weight_off = false,
                                   bool apply_jec = false,
                                   const char* sv_sf = "off",
                                   bool apply_jetveto = false,
                                   const char* track_bdt_sf = "nominal",
                                   const char* jes_variation = "nominal",
                                   const char* jer_variation = "off"){
 // load at prompt: gSystem->Load("libGenVector");
 // -- test use of central configuration
  AnalysisConfig cfg = buildConfig(RunN, pT_low, etaCut, n, btag, isMC, btagWP, sampleTag);
  // JET_PT_MAX (env, 2026-10-29): hard upper jet-pT cut at reco and gen (KinematicConfig::ptHigh).
  // Unset = no cut. The OUT_TAG must then say "ptmax<value>" (checked by the run scripts).
  if (const char *pm = std::getenv("JET_PT_MAX")) if (std::atof(pm) > 0.) cfg.kin.ptHigh = std::atof(pm);
  std::cout << "Jet pT upper cut (JET_PT_MAX): " << (cfg.kin.ptHigh > 0. ? Form("%g GeV, reco and gen", cfg.kin.ptHigh) : "none") << std::endl;
  // DR_MIN (env, 2026-10-30): lower dR cut at reco and gen, a jet selection for every observable
  // (binning_histos_small.h). It must sit on a dR bin edge, or the new first bin would be cut in two.
  if (drMinCut() > 0. && std::abs(dr_binsVector[drFirstEdge()] - drMinCut()) > 1e-9) {
    std::cerr << "FATAL: DR_MIN=" << drMinCut() << " is not a dR bin edge -- refusing" << std::endl;
    return;
  }
  std::cout << "Jet dR lower cut (DR_MIN): " << (drFirstEdge() > 0 ? Form("dR > %g, reco and gen; dR axis %d bins from %g", dr_binsVector[drFirstEdge()], obsDr().nbins, obsDr().bins[0]) : "none (reco dR > 0.005)") << std::endl;
  { const char *ss = std::getenv("SV_SEED");   // 2026-11-05, see makeSvtxs_withBDT
    std::cout << "SV seed rule (SV_SEED): " << ((ss && TString(ss) == "on") ? "ON, an SV with a track passing the BDT keeps all its tracks" : "off") << std::endl;
    std::cout << "m2B from the nominal vertices (SV_SEED_M2B): " << (svSeedNominalM2B() ? "ON" : "off") << std::endl; }

	  if (inputFileOverride && std::string(inputFileOverride).size() > 0) {
	    cfg.dataset.filename = inputFileOverride;
	    std::cout << "Overriding input file: " << cfg.dataset.filename << std::endl;
	  }
	  if (outputFolderOverride && std::string(outputFolderOverride).size() > 0) {
	    cfg.dataset.output_folder = outputFolderOverride;
	    if (!cfg.dataset.output_folder.EndsWith("/")) cfg.dataset.output_folder += "/";
	    std::cout << "Overriding output folder: " << cfg.dataset.output_folder << std::endl;
	  }

	  // ---- Tracking-efficiency systematic ----

	  TrkEffSyst::enable(track_eff_unc, cfg.dataset.isMC, trkSeed);
	  TrkEffSyst::print();

	  // ---- EEC weight on/off ----
	  EecWeight::enable(!eec_weight_off);
	  EecWeight::print();

	  // ---- Jet energy corrections on/off (see the Jec namespace at the top) ----
	  Jec::enable(apply_jec);
	  Jec::print(cfg.dataset.isMC);

	  // ---- JES uncertainty: nominal, up, or down on reconstructed MC only ----
	  JesSyst::configure(jes_variation, apply_jec, cfg.dataset.isMC);
	  JesSyst::print();

	  // ---- JER smearing: hybrid method on reconstructed MC only ----
	  JerSyst::configure(jer_variation, apply_jec, cfg.dataset.isMC);
	  JerSyst::print();

	  // ---- SV b scale factor (see sv_sf.h) ----
	  SvSF::enable(sv_sf, cfg.dataset.isMC);
	  SvSF::print();

	  // ---- Track-BDT efficiency/mistag scale factors ----
	  TrackBdtSF::enable(track_bdt_sf, cfg.dataset.isMC);
	  TrackBdtSF::print();

	  // ---- 2024 jet veto map on/off (see the JetVeto namespace at the top) ----
	  JetVeto::enable(apply_jetveto && cfg.dataset.RunN == 3);
	  JetVeto::print();

      // make_templates(cfg, 1, 1e+04);
      // TString output_hist_response =  Form("response_templatefit_n1_bjet_Run%d", cfg.dataset.RunN);
      // create_response_templatefit(cfg, output_hist_response, 0, 1e+04);

    cout << "Build_templates " << endl;
    Build_templates(cfg, makeTemplates, createRmatrix, ev_first, ev_last, -1, makeAggNtuple); // switches exposed for diagnostic/template-only runs

   std::cout << "finished :D " << std::endl;
  
  }
