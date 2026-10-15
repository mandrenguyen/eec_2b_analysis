//For both
#include "binning_histos_small.h"
#include "tTree.h"
#include <iostream>
#include <Math/Vector4D.h>
#include <vector>
#include <map>
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

// ============================================================================
// The EEC weight, on or off
// ============================================================================
// The analysis measures the energy-energy correlator: every selected 2b jet enters its
// histograms weighted by (pt1*pt2)^n, the pTs being those of its two aggregated B hadrons.
// Turning that weight off turns the SAME chain into a plain YIELD measurement -- the final
// unfolded distribution becomes dN/dr, the number of 2b jets per dr bin, with nothing else
// about the selection, the corrections or the unfolding changed.
//
// This is one switch rather than a second analysis because the weight is computed in exactly
// six places below, all the same expression, and everything downstream is built from those
// six values:
//
//   h3D_0b / h3D_b / h3D_bb / h3D_data   the template-fit and unfolding inputs
//   w_reco = weight_tree * eec_reco      the response matrix, the purity and efficiency
//   w_gen  = weight_tree * eec_gen       corrections, and the MC truth
//   hgenjet_2b_*                         the b-tag and SV efficiency corrections
//
// So value() returning 1 makes all of them unweighted at once, consistently, with no second
// code path to keep in step. Two consequences worth knowing:
//
//   * the unweighted h3D_bb is then identical to h_count_bb, which is filled with
//     weight_tree alone and already exists in every production -- a free cross-check that
//     the switch does what it claims;
//   * the EEC-weight correction (r_eec = sum(eec_gen)/sum(eec_reco), applied by multiplying
//     the unfolded result) collapses to exactly 1 bin by bin, because its numerator and
//     denominator become the same unweighted histogram. It needs no special case: it turns
//     itself off.
//
// NOT a systematic variation -- a different observable measured with the same machinery, so
// its outputs carry their own "_noeecw" tag (eecWeightTag() in result_paths.h) and never mix
// with the EEC ones. Default ON, so a caller that never mentions this gets what it always
// got.
// ============================================================================
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

} // namespace EecWeight

// ============================================================================
// Jet energy corrections (added 2026-10-01)
// ============================================================================
// The 2024 pp-reference JEC, exactly as in Matt's JEC_Prompt24HIpp productions (the top-level
// eec_2b_analysis/create_files_for_template_fit.cpp): rho-dependent L1FastJet + MC-truth
// L2Relative on MC and data, plus the dijet L2Residual on data only, applied to rawpt. The
// corrected pT REPLACES t.jtpt right after each event is read, so the jet-pT cut, the pT
// bins, the response and every histogram see corrected jets; gen jets are untouched.
//
// Default OFF, so a caller that never mentions it reproduces the uncorrected productions
// bit for bit. Switched on by the last argument of create_files_for_template_fit(); the run
// scripts insist the OUT_TAG then says "jec" (eraTag() in result_paths.h), so corrected and
// uncorrected files can never share a name.
//
// Text files: the Run3HeavyIonJetCalibrations git submodule at the repo top, pinned commit
// c64293a (fetched 2026-10-01 with `git submodule update --init`).
#include "../Run3HeavyIonJetCalibrations/include/JetCorrector.h"
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
    return "/home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/Run3HeavyIonJetCalibrations/txt/2024ppRef/";
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

} // namespace Jec

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

// ⚠️ UNCALLED since 2026-09-22: its only caller was the fnb/fb/lnfb gen block, which is
// commented out with the rest of the fraction observables. Left compiled rather than
// disabled so that re-enabling them is an uncomment at the fill site alone.
//
// The gen half of the fnb observable (observables.h): the scalar pT sum of this jet's gen
// tracks that are NOT from a B decay. Exactly the complement of PartialBsAggregation just
// above -- same jet, same pT > 1 cut, status < 100 instead of >= 100 -- so pT_b and pT_nonb
// between them account for every gen track of the jet above 1 GeV, and neither double counts.
//
// No tracking-efficiency variation here, deliberately: the systematic varies RECO only
// (TrkEffSyst is never applied to refTrk*), because a variation that moved truth as well
// would not change the response matrix, which is the thing it exists to probe.
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

// ---- disabled (kept for reference): MatchingTracksAggregation, unused ----
#if 0
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
    else std::cout << "PDG:" << pid << std::endl;

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
#endif // ---- end disabled: MatchingTracksAggregation ----

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
// pt_nonb (optional, may be null): the scalar pT sum of this jet's tracks that do NOT end up
// in either of the two B hadrons -- the reco half of the fnb observable (observables.h).
// ⚠️ Every caller passes it since 2026-09-22 but nothing CONSUMES it: fnb, fb and lnfb are
// commented out. It stays filled because it is an out-parameter of a pass being made anyway.
// It is computed HERE and nowhere else because this is the only place that knows which
// tracks the B reconstruction kept: a track fails the BDT cut below and is gone, and every
// track that passes but has no SV is absorbed into the nearer vertex further down. So the
// complement is exactly "passed the jet and pT cuts, failed the BDT cut", which cannot be
// reconstructed from the two returned four-vectors afterwards.
// sv_sta (optional, may be null): the dominant trkMatchSta of the SV tracks of each returned
// vertex, same order as the four-vectors. >= 100 means the SV is from a B hadron -- that is
// what decides whether the SV b scale factor applies to it (sv_sf.h). MC only; meaningless
// in data, where nothing reads it.
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
  std::vector<int>* sv_sta = nullptr
){
  if (pt_nonb) *pt_nonb = 0.;
  if (sv_sta) sv_sta->clear();

  // Failure counters: nb_sv (fewer than 2 SVs), sv_fail / merge_fail / agg_fail (match status).
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
    //if (t.trkMatchSta[itrk] <100 ) continue;
    //if (t.trkMatchSta[itrk] <= 1) continue; //activate or not
    //control plots for the BDT score distribution in 0b and 1b jets
    //if (t.trkMatchSta[itrk] == 1)
      //{h_score_bkg->Fill(t.trkBdtScore[itrk], weight_tree);}
    //else if (t.trkMatchSta[itrk] >= 100)
     // {h_score_sg->Fill(t.trkBdtScore[itrk], weight_tree);}

    // -- Add trkbdt score
    if (!TrackBdtSF::passes(t.trkBdtScore[itrk], t.trkMatchSta[itrk], t.trkSvtxId[itrk])) {
      // Rejected as not b-like: this track is the NON-B part of the jet. Sum its pT for the
      // fnb observable before dropping it. The tracking-efficiency variation is applied here
      // too, so a 3% track drop thins the non-b side exactly as it thins the Bs -- the
      // detector does not lose tracks according to a BDT score. Safe to call dropTrack() on
      // a track the nominal path never reaches: it hashes (seed, ient, itrk) and reads no
      // track property, so asking about this track changes no other track's verdict.
      if (pt_nonb && !TrkEffSyst::dropTrack(ient, itrk)) *pt_nonb += t.trkPt[itrk];
      continue;
    }

        // Tracking-efficiency systematic: throw the track away before anything is built
        // from it, so the SVs, m_2B, dR and the EEC weight all see a jet genuinely
        // missing 3% of the tracks that would have entered the B reconstruction.
        // No-op in the nominal. See tracking_efficiency_syst.h.
        //
        // Placed AFTER the pT and BDT cuts on purpose: "3% of the tracks used to build
        // the Bs". The position is in fact irrelevant -- dropTrack() hashes
        // (seed, ient, itrk) and reads no track property, so the surviving set is the
        // same intersection either way and the two orderings give bit-identical output
        // (verified). That is only true because the number is hashed rather than drawn
        // from a stream: a stream would assign different numbers to different tracks
        // depending on how many were consumed before the cuts.
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
  Int_t passRecoKin, passGenKin, passBtag, nRecoAgg, nGenAgg, genStatus1, genStatus2;
  Float_t recoPt1, recoEta1, recoPhi1, recoM1, recoPt2, recoEta2, recoPhi2, recoM2, recoDr, recoMB, recoEec, genPt1, genEta1, genPhi1, genM1, genPt2, genEta2, genPhi2, genM2, genDr, genMB, genEec;

  void reset() {
    entry = -1;
    run = lumi = 0;   // unsigned: 0 is the "unset" marker, not -1
    evt = 0;
    jetIndex = -1;
    weight = 1.0;
    jtpt = jteta = refpt = refeta = -999.0;
    jtNbHad = -1;
    btagScore = -999.0;
    passRecoKin = passGenKin = passBtag = 0;
    nRecoAgg = nGenAgg = 0;
    genStatus1 = genStatus2 = 0;
    recoPt1 = recoEta1 = recoPhi1 = recoM1 = -999.0;
    recoPt2 = recoEta2 = recoPhi2 = recoM2 = -999.0;
    recoDr = recoMB = recoEec = -999.0;
    genPt1 = genEta1 = genPhi1 = genM1 = -999.0;
    genPt2 = genEta2 = genPhi2 = genM2 = -999.0;
    genDr = genMB = genEec = -999.0;
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
  br("nRecoAgg", &row.nRecoAgg, "I"); br("nGenAgg", &row.nGenAgg, "I");
  br("genStatus1", &row.genStatus1, "I"); br("genStatus2", &row.genStatus2, "I");
  br("recoPt1", &row.recoPt1, "F"); br("recoEta1", &row.recoEta1, "F"); br("recoPhi1", &row.recoPhi1, "F"); br("recoM1", &row.recoM1, "F");
  br("recoPt2", &row.recoPt2, "F"); br("recoEta2", &row.recoEta2, "F"); br("recoPhi2", &row.recoPhi2, "F"); br("recoM2", &row.recoM2, "F");
  br("recoDr", &row.recoDr, "F"); br("recoMB", &row.recoMB, "F"); br("recoEec", &row.recoEec, "F");
  br("genPt1", &row.genPt1, "F"); br("genEta1", &row.genEta1, "F"); br("genPhi1", &row.genPhi1, "F"); br("genM1", &row.genM1, "F");
  br("genPt2", &row.genPt2, "F"); br("genEta2", &row.genEta2, "F"); br("genPhi2", &row.genPhi2, "F"); br("genM2", &row.genM2, "F");
  br("genDr", &row.genDr, "F"); br("genMB", &row.genMB, "F"); br("genEec", &row.genEec, "F");
}

// ---- disabled (kept for reference): make_templates, superseded by Build_templates ----
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
#endif // ---- end disabled: make_templates ----

// ---- disabled (kept for reference): create_response_templatefit, now in Build_templates ----
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

    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
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
            double jpt_gen = jtpt_fill(gen_jet_pt(t, ijet));
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

// ---- The observable list, and the one place that fixes its order -----------------------
// obsList() and obsValues() are a matched pair: the Nth entry of obsValues() must be the
// value of the Nth observable in obsList(). Everything downstream indexes the two together,
// so this is the single place where "which observable is index 1" is decided. Adding a third
// observable means one more entry in EACH -- and the size check in Build_templates fails
// loudly if only one of them is edited.
//
// ⚠️ 2026-09-22: this pass fills TWO observables, dr and the momentum balance B. The three
// b/non-b fractions are commented out here and in observables.h -- they are one quantity
// three ways (fb = 1 - fnb, lnfb = ln(1 + fb)) and none of them is being measured. The
// non-b pT sums they need are still computed at the fill sites, also commented out; see
// the DISABLED markers there.
//
// ⚠️ B was called z until 2026-09-22. Its histograms now carry "_B", so this pass must be
// re-run into a fresh production (OUT_TAG=_upartv2_B) before anything downstream can read
// the balance axis -- obsProdTag() in result_paths.h.
inline std::vector<ObsDef> obsList()
{
    return { obsDr(), obsB() };
    // DISABLED 2026-09-22: { obsDr(), obsB(), obsFnb(), obsFb(), obsLnFb() }
    // fnb, fb and lnfb were kept together on purpose: each costs one extra set of histograms
    // in a pass that is already being made, and having all three lets the binnings be judged
    // against each other. fb and fnb on matching uniform edges are exact mirrors, which is a
    // free closure check on the whole chain. Re-enable all three together, or the mirror
    // check goes away with them.
}

inline std::vector<double> obsValues(double dr, double B)
{
    return { dr, B };
    // DISABLED 2026-09-22, matching obsList() above -- the signature was
    //   obsValues(double dr, double B, double fnb, double fb, double lnfb)
    // returning { dr, B, fnb, fb, lnfb }. Restoring it means restoring the extra arguments
    // at all four call sites, which the size check in Build_templates will catch.
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
    }

    // x is this observable's RECO value, already past o.passReco(). eec carries the data
    // prescale the caller applied; w is the tree weight.
    void fill(bool isMC, Int_t jtNbHad, Long64_t ient,
              double mB, double x, double jtpt, double eec, double w) {
        const double xf = o.fill(x);
        if (isMC) {
            if      (jtNbHad == 0) { h3D_0b->Fill(mB, xf, jtpt, eec * w); h_count_0b->Fill(mB, xf, jtpt, w); }
            else if (jtNbHad == 1) { h3D_b ->Fill(mB, xf, jtpt, eec * w); h_count_b ->Fill(mB, xf, jtpt, w); }
            else if (jtNbHad >= 2) { h3D_bb->Fill(mB, xf, jtpt, eec * w); h_count_bb->Fill(mB, xf, jtpt, w); }

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
        }
    }

    void write(bool isMC) const {
        if (isMC) {
            h3D_0b->Write(); h3D_b->Write(); h3D_bb->Write();
            h3D_pseudodata->Write(); h3D_pseudodata_bb->Write();
            // Odd-half templates, for fitting the even-half pseudodata in the split test
            h3D_pseudo_0b->Write(); h3D_pseudo_b->Write(); h3D_pseudo_bb->Write();
            h_count_0b->Write(); h_count_b->Write(); h_count_bb->Write();
        } else {
            h3D_data->Write(); h_count_data->Write();
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

    void book(const ObsDef &def) {
        o = def;
        const TString xr = Form("x=%s_reco, y=jtpt", o.name.Data());
        const TString xg = Form("x=%s_gen, y=jtpt",  o.name.Data());
        for (int h = 0; h < 2; ++h) {
            purity_num[h] = book2D(Form("h_half%d_purity_numerator_tf",   h), xr);
            purity_den[h] = book2D(Form("h_half%d_purity_denominator_tf", h), xr);
            eff_num[h]    = book2D(Form("h_half%d_efficiency_numerator_tf",   h), xg);
            eff_den[h]    = book2D(Form("h_half%d_efficiency_denominator_tf", h), xg);
            ps_purity_num[h] = book2D(Form("h_half%d_pseudo_purity_numerator_tf",   h), xr);
            ps_purity_den[h] = book2D(Form("h_half%d_pseudo_purity_denominator_tf", h), xr);
            ps_eff_num[h]    = book2D(Form("h_half%d_pseudo_efficiency_numerator_tf",   h), xg);
            ps_eff_den[h]    = book2D(Form("h_half%d_pseudo_efficiency_denominator_tf", h), xg);
        }
        // Closure target for the split test: the pseudodata half's gen-level distribution,
        // with the same gate (gen_pass) and weight (w_reco) as the purity/response chain.
        pseudodata_truth = book2D("h_pseudodata_truth_tf", xg);
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
        hgenjet_2b          = book2D("hgenjet_2b",          Form("b-tagging eff. DENO;%s;p_{T} [GeV]", o.axis.Data()));
        hgenjet_2b_passbtag = book2D("hgenjet_2b_passbtag", Form("b-tagging eff. NUM;%s;p_{T} [GeV]",  o.axis.Data()));
        // Combined SV-reco + b-tag efficiency, gen-EEC weighted (w_gen) so the denominator
        // is defined even for true 2b jets without 2 reco SVs.
        hgenjet_2b_all       = book2D("hgenjet_2b_all",       Form("SV+btag eff. DENO;%s;p_{T} [GeV]", o.axis.Data()));
        hgenjet_2b_reco_btag = book2D("hgenjet_2b_reco_btag", Form("SV+btag eff. NUM;%s;p_{T} [GeV]",  o.axis.Data()));
    }

    // ⚠️ The four hgenjet_* histograms are filled with the RAW gen value, not the folded one,
    // while the response and purity/efficiency use the folded value. That asymmetry is
    // inherited from the original dr code: a gen value at or above the top edge lands in
    // ROOT's overflow here but in the last bin there. Preserved deliberately -- "fixing" it
    // would silently change every existing dr correction.
    void fillGenAll(double x_gen_raw, double jpt_gen, double w_gen) const {
        hgenjet_2b_all->Fill(x_gen_raw, jpt_gen, w_gen);
    }
    void fillGenPreBtag(double x_gen_raw, double jpt_gen, double w_reco) const {
        hgenjet_2b->Fill(x_gen_raw, jpt_gen, w_reco);
    }
    void fillGenPostBtag(double x_gen_raw, double jpt_gen, double w_reco, double w_gen) const {
        hgenjet_2b_passbtag  ->Fill(x_gen_raw, jpt_gen, w_reco);
        hgenjet_2b_reco_btag ->Fill(x_gen_raw, jpt_gen, w_gen);
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
  // No kinematic label in the names: the pT range is fixed by jtpt_binsVector in
  // binning_histos_small.h, so spelling it into every filename only invited drift.
  // Systematic-variation tag: "" for the nominal, so nominal names are unchanged;
  // "_trkdrop030" for the 3% tracking-efficiency variation. It goes last, right before
  // ".root", which is where apply_unfolding_2d.C's MC_TAG expects it (the run scripts
  // then append their own OUT_TAG after it).
  // EecWeight::tag() is "_noeecw" for the yield run: not a systematic, a different
  // observable, so it must never land in the same file as the EEC histograms.
  // SvSF::tag() is "_svsf" / "_svsfup" / "_svsfdn" when the SV b scale factor is applied
  // (MC only), "" when off -- see sv_sf.h.
  TString syst_suffix = TrkEffSyst::tag() + SvSF::tag() + TrackBdtSF::tag() + EecWeight::tag();
  TString fout_name = cfg.dataset.output_folder + cfg.dataset.output_hist + job_suffix + "MCGEN" + syst_suffix + ".root"; // for reposnse matrix: has Prefix: Response
  TString ResponseMatrix_fout_name =  cfg.dataset.output_folder + "RMatrix_" + cfg.dataset.output_hist + job_suffix + syst_suffix + ".root";
  TString AggBHadronNtuple_fout_name = cfg.dataset.output_folder + "AggBHadronNtuple_" + cfg.dataset.output_hist + job_suffix + syst_suffix + ".root";
  gSystem->mkdir(cfg.dataset.output_folder, true);

  // ---- Histograms ----
  // One set per measured observable. dr's ObsDef has an empty suffix, so its histograms are
  // booked under exactly the names they have always had; B's carry "_B". Both sets are
  // filled from the SAME jets in the SAME pass -- adding an observable costs no extra run.
  const std::vector<ObsDef> observables = obsList();
  // obsList() and obsValues() must stay the same length, or the per-jet value of one
  // observable would be filled into another's histograms. Caught here, at startup, rather
  // than as a plausible-looking distribution three hours into a production.
  // The zero-argument list below has to be widened when an observable is added, so the
  // COMPILER catches a half-done addition before this runtime check ever runs. That is the
  // stronger guard of the two; keep both.
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
    // -- To split any sample into two parts randomly.
    // Was std::mt19937 seeded from std::random_device, which made two runs of the same
    // job disagree in the half0/half1 histograms and made the split depend on how the
    // sample was chunked. TrkEffSyst::jetSplitRandom() hashes (ient, ijet) instead, so
    // the half a jet lands in is a property of the jet: reproducible across reruns, and
    // identical between the nominal and a systematic variation, which is what keeps the
    // split noise out of the systematic. See tracking_efficiency_syst.h.

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
          t.jtpt[ijet] = std::isfinite(correctedPt) && correctedPt > 0 ? correctedPt : -999.0f;
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

      // -- NEW:
      if (! passPVQuality_EventSelection(t, cfg)) continue;

        // cout << "PV quality pass ok: t.vz = "<< t.vz  << endl;
          // if (cfg.dataset.RunN == 3 &&  !isMC ) cout << " and t.pprimaryVertexFilter = "<< t.pprimaryVertexFilter << endl;

      // -- Trigger selections for all
      if (! passEventSelection(t, cfg)) continue;
      // cout << "event pass ok" << endl;

    for (Int_t ijet = 0; ijet < t.nref; ijet++) { // Jet loop

      // Fill templates from DATA or RECO MC: jet kinematics + btagging

      if(isMakeTemplates && passRecoJetKinematics(t, ijet, cfg) &&  passBtag(t, ijet, cfg)){

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

                double dr   = t.calc_dr(reco_sv[0].Eta(), reco_sv[0].Phi(), reco_sv[1].Eta(), reco_sv[1].Phi());
                double pt1  = reco_sv[0].Pt();
                double pt2  = reco_sv[1].Pt();
                // SV b scale factor, one per b-matched SV (sv_sf.h). 1 for data and when off.
                double w_svsf = SvSF::jetWeight(pt1, reco_sv_sta[0], pt2, reco_sv_sta[1]);
                double eec  = EecWeight::value(pt1, pt2, cfg.n);
                double jtpt = jtpt_fill(t.jtpt[ijet]);
                double mB   = reco_sv[0].M() + reco_sv[1].M();
                // The second observable, from the same two B hadrons. Symmetric in pt1/pt2,
                // so the SV ordering does not matter.
                double Bbal = MomBalance::value(pt1, pt2);
                // DISABLED 2026-09-22 -- fnb/fb/lnfb, how much of the jet is NOT in those two
                // B hadrons, the same split the other way up, and the log of 1 + that. No
                // overflow fold here; ObsDef::fill() does it inside the histogram set, which
                // is where every observable's fold belongs (the dr line above predates that).
                //   double fnb  = NonBFraction::value(pt1 + pt2, pt_nonb_reco);
                //   double fb   = BFraction::value(pt1 + pt2, pt_nonb_reco);
                //   double lnfb = LogBFraction::value(pt1 + pt2, pt_nonb_reco);
                // pt_nonb_reco above is still computed by makeSvtxs_withBDT -- it is an
                // out-parameter of the SV pass, not a separate loop, so leaving it costs
                // nothing and keeps the three lines above a straight uncomment.
                (void)pt_nonb_reco;

                //Fix the under/overflow
                //if(dr < dr_min) dr = dr_min_fill;
                if(dr >= dr_max) dr = dr_max_fill;
                if(mB >= mb_max) mB = mb_max_fill;

                // Match the response-matrix reco definition: drop the dr<=0.005 spike so
                // templates and response cover the same reco domain. Data and all MC.
                //
                // ⚠️ This gate is dr's, and it applies to EVERY observable: the B templates
                // are filled from the same jets as the dr ones, so a jet dropped here is
                // dropped from B too. That is deliberate -- the two measurements must share
                // a reco domain, or the signal fraction fitted on one does not apply to the
                // other. It is ObsDef::passReco() for dr, spelled out here because it gates
                // the whole block rather than one fill.
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
                                mB, xreco[io], jtpt, eec, weight_tree * w_svsf);
                }

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
            double jpt_gen = jtpt_fill(gen_jet_pt(t, ijet));
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
	    if (aggBHadronTree) aggBHadronTree->Fill();
	    if (!isCreateRmatrix) continue;

	    // Combined SV-reco + b-tag efficiency: DENOMINATOR = all true 2b jets, before the
	    // 2-SV cut. Gen-EEC weighted (w_gen), reused for the numerator below.
	    double w_gen = gen_ok ? weight_tree * eec_gen : 0.0;
	    // Gen values in observable order. Reused by every hgenjet_* fill below.
	    const std::vector<double> xgen_raw = obsValues(dr_gen, B_gen);
	    if (gen_ok)
	      for (size_t io = 0; io < resp.size(); ++io)
	        if (xgen_raw[io] >= 0.) resp[io].fillGenAll(xgen_raw[io], jpt_gen, w_gen);

	    // Guard the reco pair access: makeSvtxs_withBDT returns an empty vector when the
	    // jet has < 2 secondary vertices, and reco_sv_rm[0] below would then segfault.
	    if (!reco_sv_ok) continue;

      // -- Prepare to Fill Rmatrix --
		  // -- Reco pair info
      double mB_reco  = reco_sv_rm[0].M() + reco_sv_rm[1].M();
      double dr_reco  = t.calc_dr(reco_sv_rm[0].Eta(), reco_sv_rm[0].Phi(), reco_sv_rm[1].Eta(), reco_sv_rm[1].Phi());
      double eec_reco = EecWeight::value(reco_sv_rm[0].Pt(), reco_sv_rm[1].Pt(), cfg.n);
      double B_reco = MomBalance::value(reco_sv_rm[0].Pt(), reco_sv_rm[1].Pt());
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
			if (gen_ok)
			  for (size_t io = 0; io < resp.size(); ++io)
			    if (xgen_raw[io] >= 0.) resp[io].fillGenPreBtag(xgen_raw[io], jpt_gen, w_reco);

			// -- Use combined b-tagger: Upart tagger + 2SV of aggreagted Bs., for Rmatrix unfolding, purity, reconstruction eff.
			if (!passBtag(t, ijet, cfg)) continue; // select btagged jets only (reco_sv_ok already required above)

			// Fill total number of True 2b Jets that survive after btagging condition (Num of b tagging eff. correction - used after unfolding)
			// and the combined SV-reco + b-tag efficiency NUMERATOR = 2b jets reconstructed
			// with 2 SV AND tagged (this point is past both cuts). Same gen-EEC weight as
			// hgenjet_2b_all.
			if (gen_ok)
			  for (size_t io = 0; io < resp.size(); ++io)
			    if (xgen_raw[io] >= 0.) resp[io].fillGenPostBtag(xgen_raw[io], jpt_gen, w_reco, w_gen * w_svsf);

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
            bool reco_pass = passRecoJetKinematics(t, ijet, cfg) && dr_reco_fill > 0.005;

            // Define gen_pass: particle-level jet kinematics + gen observable range.
            // Requires gen_ok, so a jet with < 2 gen B-hadrons stays out of num/eff/response.
            bool gen_pass  = gen_ok && passGenJetKinematics(t, ijet, cfg);

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

  outFile->Close(); delete outFile;
}

} // end Function

// ---- disabled (kept for reference): filter_b_bb_as_data_and_mc, split-sample variant ----
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
#endif // ---- end disabled: filter_b_bb_as_data_and_mc ----

// NOTE: pT_high (4th) and foldPtOverflow (last) were removed; the upper pT edge lives in
// jtpt_binsVector in binning_histos_small.h. Old callers must drop the 4th argument.
// NOTE: dataType (2nd) was removed too — isMC says data vs MC, and sampleTag (last, set by
// the run scripts) says which sample the output files are named for.
void create_files_for_template_fit(Int_t RunN = 3, Float_t pT_low = 80, Float_t etaCut = 2, Int_t n = 1,
                                   bool btag = true, bool isMC = true, Double_t btagWP = 0.712, bool makeTemplates = true, bool createRmatrix = true,
                                   bool makeAggNtuple = true, Long64_t ev_first = 0, Long64_t ev_last = -1, const char* inputFileOverride = "", const char* outputFolderOverride = "",
                                   const char* sampleTag = "",
                                   bool track_eff_unc = false, ULong64_t trkSeed = 20260908,
                                   bool eec_weight_off = false,
                                   bool apply_jec = false,
                                   const char* sv_sf = "off",
                                   const char* track_bdt_sf = "nominal"){
 // load at prompt: gSystem->Load("libGenVector");
 // -- test use of central configuration
  AnalysisConfig cfg = buildConfig(RunN, pT_low, etaCut, n, btag, isMC, btagWP, sampleTag);

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
	  // track_eff_unc = false is the nominal and changes nothing. true throws away 3%
	  // (TrkEffSyst::kDropFraction) of the reconstructed tracks in MC during the B
	  // reconstruction, chosen per track from a hash of (seed, entry, track index) so
	  // the same job rerun -- or the same sample re-chunked, or run in parallel instead
	  // of serially -- drops exactly the same tracks.
	  // Data is never varied, hence the isMC guard. See tracking_efficiency_syst.h.
	  TrkEffSyst::enable(track_eff_unc, cfg.dataset.isMC, trkSeed);
	  TrkEffSyst::print();

	  // ---- EEC weight on/off ----
	  // eec_weight_off = false is the nominal EEC analysis and changes nothing. true drops
	  // the (pt1*pt2)^n weight everywhere, turning the same chain into a yield measurement
	  // whose unfolded result is dN/dr. Applies to MC and data alike -- it is the observable,
	  // not a correction -- so there is no isMC guard here. See EecWeight at the top.
	  EecWeight::enable(!eec_weight_off);
	  EecWeight::print();

	  // ---- Jet energy corrections on/off (see the Jec namespace at the top) ----
	  Jec::enable(apply_jec);
	  Jec::print(cfg.dataset.isMC);

	  // ---- SV b scale factor (see sv_sf.h) ----
	  // "off" (default) changes nothing. "nominal" weights every b-matched reco SV in MC by
	  // SF(aggregate pT); "up"/"down" shift that SF by its stat (+) syst uncertainty. Data is
	  // never weighted, hence the isMC guard inside.
	  SvSF::enable(sv_sf, cfg.dataset.isMC);
	  SvSF::print();

	  // ---- Track-BDT efficiency/mistag scale factors ----
	  // MC uses category-dependent effective working points derived from the measured SFb
	  // and SFl. "boff" and "loff" independently remove the corresponding correction;
	  // "off" is the uncorrected cross-check. Data always keeps the nominal 0.365 cut.
	  TrackBdtSF::enable(track_bdt_sf, cfg.dataset.isMC);
	  TrackBdtSF::print();

      // make_templates(cfg, 1, 1e+04);
      // TString output_hist_response =  Form("response_templatefit_n1_bjet_Run%d", cfg.dataset.RunN);
      // create_response_templatefit(cfg, output_hist_response, 0, 1e+04);

    cout << "Build_templates " << endl;
    Build_templates(cfg, makeTemplates, createRmatrix, ev_first, ev_last, -1, makeAggNtuple); // switches exposed for diagnostic/template-only runs

   std::cout << "finished :) :D " << std::endl;
  //filter_b_bb_as_data_and_mc(cfg);
  }
