#ifndef SV_SF_H
#define SV_SF_H

// ============================================================================
// Secondary-vertex b scale factor (added 2026-10-08)
// ============================================================================
//
// Data/MC scale factor for reconstructing a b-hadron secondary vertex, measured on
// single-b jets (JP, muon pT > 5, UParT WP 0.712) as a function of the SV's aggregate pT
// under the single-b hypothesis:
//
//     sv_sfb_singleb_aggpt_jp_wp712_mu5_systematics.root   (workflow folder)
//     key h_SFb_singleb_sv_aggpt_central      SF, bin error = statistical
//     key h_SFb_singleb_sv_aggpt_syst_total   systematic, = no_muon (+) bb_model (+)
//                                             jpcalib_hf (+) qqrate in quadrature (checked
//                                             2026-10-08, all 7 bins)
//     7 bins in aggregate pT: 0-20-40-60-80-100-140-200 GeV, SF 0.97 -> 1.04.
//
// ⚠️ As with the UParT SF file, KEYS and internal histogram names differ (key
// "h_SFb_singleb_sv_aggpt_central" holds "h_SFb_singleb_sv_aggpt_jp_mu5_central").
// Get() takes the key.
//
// ---------------------------------------------------------------------------
// WHERE IT IS APPLIED
// ---------------------------------------------------------------------------
// At the SV reconstruction, per vertex: makeSvtxs_withBDT() returns the two aggregated
// B candidates, and each one that is TRUTH-MATCHED to a B hadron (most of its SV tracks
// have trkMatchSta >= 100) gets SF(its pT). An SV that is not from a B (light / fake,
// dominant status < 100) gets 1 -- a b-efficiency SF says nothing about it, and the light
// side is the mistag systematic's job. The jet weight is the product over its two SVs.
//
// MC only (data never). The weight multiplies every reco-level MC fill made from those
// SVs: the 0b/1b/2b templates, the response, purity and efficiency, and the numerators of
// the b-tag and SV+b-tag efficiency corrections. The all-true-2b DENOMINATOR
// (hgenjet_2b_all) is filled before any SV is reconstructed and is NOT weighted, so the
// SV+b-tag efficiency correction becomes the data-like SF * eff_MC.
//
// ---------------------------------------------------------------------------
// VARIATIONS
// ---------------------------------------------------------------------------
//   off      no SF (default -- output bit-identical to before this header existed)
//   nominal  central SF                                         file tag "_svsf"
//   up       central + sqrt(stat^2 + syst_total^2), every bin   file tag "_svsfup"
//   down     central - sqrt(stat^2 + syst_total^2), every bin   file tag "_svsfdn"
// The shift is COHERENT in pT (all bins up, or all down), the conservative choice for a
// single-number uncertainty.
//
// PER-CONTRIBUTION VARIATIONS (2026-10-13, Zoe: the systematic is each contribution varied on
// its own and added in quadrature). Each is coherent in pT, file tag "_svsf<name>":
//   statup / statdn       central +/- its bin error (statistical)
//   nomuonup / nomuondn   central +/- h_..._syst_no_muon
//   bbup / bbdn           central +/- h_..._syst_bb_model
//   jpup / jpdn           central +/- h_..._syst_jpcalib_hf  (= |JP HF swap - central|, checked)
//   qqup / qqdn           the h_..._qqrate_up / _qqrate_down curves themselves
// jp* and qq* are varied TOGETHER with the UParT SF's jpsyst_up/dn and qqrate_up/down in the
// same chain (the two SFs share those calibrations) -- that pairing is done by the driver
// (run_svsf_components.sh) and the band, not here. The tag is written by the macro into its own filenames, the
// same way as TrkEffSyst::tag(), so no variation can overwrite another or the nominal.
//
// An SV above 200 GeV takes the last bin, below 0 cannot happen.
// ============================================================================

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"

namespace SvSF {

inline TString file()
{
    return "/home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/workflow/"
           "sv_sfb_singleb_aggpt_jp_wp712_mu5_systematics.root";
}

enum class Mode { Off, Nominal, Up, Down,
                  StatUp, StatDn, NoMuonUp, NoMuonDn, BbUp, BbDn, JpUp, JpDn, QqUp, QqDn };

// variation name <-> mode <-> file tag, in one table.
struct ModeName { const char *name; Mode mode; const char *tag; };
inline const std::vector<ModeName> &modeNames()
{
    static const std::vector<ModeName> t = {
        {"off", Mode::Off, ""}, {"nominal", Mode::Nominal, "_svsf"},
        {"up", Mode::Up, "_svsfup"}, {"down", Mode::Down, "_svsfdn"},
        {"statup", Mode::StatUp, "_svsfstatup"},     {"statdn", Mode::StatDn, "_svsfstatdn"},
        {"nomuonup", Mode::NoMuonUp, "_svsfnomuonup"}, {"nomuondn", Mode::NoMuonDn, "_svsfnomuondn"},
        {"bbup", Mode::BbUp, "_svsfbbup"},           {"bbdn", Mode::BbDn, "_svsfbbdn"},
        {"jpup", Mode::JpUp, "_svsfjpup"},           {"jpdn", Mode::JpDn, "_svsfjpdn"},
        {"qqup", Mode::QqUp, "_svsfqqup"},           {"qqdn", Mode::QqDn, "_svsfqqdn"} };
    return t;
}

struct Config {
    Mode  mode = Mode::Off;
    bool  isMC = false;
    TH1D *central = nullptr;   // owned, detached from the file
    TH1D *syst    = nullptr;
    TH1D *shift   = nullptr;   // per-contribution variations: the +/- histogram, or the curve
};

inline Config &config()
{
    static Config c;
    return c;
}

inline bool active()
{
    const Config &c = config();
    return c.mode != Mode::Off && c.isMC;
}

// variation: "off" | "nominal" | "up" | "down". Anything else throws -- an unknown value
// silently becoming "off" (or nominal) is exactly the aliasing result_paths.h guards against.
inline void enable(const TString &variation, bool isMC)
{
    Config &c = config();
    bool known = false;
    for (const ModeName &m : modeNames())
        if (variation == m.name || (variation == "" && m.mode == Mode::Off)) { c.mode = m.mode; known = true; }
    if (!known)
        throw std::runtime_error(("unknown SV SF variation '" + std::string(variation.Data())
                                   + "' (use off | nominal | up | down | stat|nomuon|bb|jp|qq + up|dn)").c_str());
    c.isMC = isMC;
    if (c.mode == Mode::Off) return;

    TFile *f = TFile::Open(file());
    if (!f || f->IsZombie())
        throw std::runtime_error(("cannot open SV SF file " + std::string(file().Data())).c_str());
    TH1D *hc = dynamic_cast<TH1D *>(f->Get("h_SFb_singleb_sv_aggpt_central"));
    TH1D *hs = dynamic_cast<TH1D *>(f->Get("h_SFb_singleb_sv_aggpt_syst_total"));
    if (!hc || !hs)
        throw std::runtime_error("SV SF file is missing h_SFb_singleb_sv_aggpt_central or _syst_total");
    c.central = (TH1D *) hc->Clone("svsf_central");
    c.syst    = (TH1D *) hs->Clone("svsf_syst_total");
    c.central->SetDirectory(nullptr);
    c.syst->SetDirectory(nullptr);
    // The histogram a per-contribution variation reads.
    const char *key = nullptr;
    switch (c.mode) {
        case Mode::NoMuonUp: case Mode::NoMuonDn: key = "h_SFb_singleb_sv_aggpt_syst_no_muon";    break;
        case Mode::BbUp:     case Mode::BbDn:     key = "h_SFb_singleb_sv_aggpt_syst_bb_model";   break;
        case Mode::JpUp:     case Mode::JpDn:     key = "h_SFb_singleb_sv_aggpt_syst_jpcalib_hf"; break;
        case Mode::QqUp:                          key = "h_SFb_singleb_sv_aggpt_qqrate_up";       break;
        case Mode::QqDn:                          key = "h_SFb_singleb_sv_aggpt_qqrate_down";     break;
        default: break;
    }
    if (key) {
        TH1D *h = dynamic_cast<TH1D *>(f->Get(key));
        if (!h) throw std::runtime_error(std::string("SV SF file is missing ") + key);
        c.shift = (TH1D *) h->Clone("svsf_shift");
        c.shift->SetDirectory(nullptr);
    }
    f->Close();
    delete f;
}

// SF for one SV of aggregate pT `pt`, in the configured variation. 1 when off or data.
inline double value(double pt)
{
    if (!active()) return 1.0;
    const Config &c = config();
    const int nb = c.central->GetNbinsX();
    int b = c.central->FindFixBin(pt);
    if (b < 1)  b = 1;
    if (b > nb) b = nb;   // above 200 GeV: last bin
    const double sf = c.central->GetBinContent(b);
    switch (c.mode) {
        case Mode::Nominal: return sf;
        case Mode::Up:   return sf + std::hypot(c.central->GetBinError(b), c.syst->GetBinContent(b));
        case Mode::Down: return sf - std::hypot(c.central->GetBinError(b), c.syst->GetBinContent(b));
        case Mode::StatUp: return sf + c.central->GetBinError(b);
        case Mode::StatDn: return sf - c.central->GetBinError(b);
        case Mode::NoMuonUp: case Mode::BbUp: case Mode::JpUp: return sf + c.shift->GetBinContent(b);
        case Mode::NoMuonDn: case Mode::BbDn: case Mode::JpDn: return sf - c.shift->GetBinContent(b);
        case Mode::QqUp: case Mode::QqDn: return c.shift->GetBinContent(b);   // the curve itself
        default: return 1.0;
    }
}

// The jet weight: product over the two reconstructed SVs, b-matched ones only.
// sta = dominant trkMatchSta of each SV (makeSvtxs_withBDT's sv_sta out-parameter).
inline double jetWeight(double pt1, int sta1, double pt2, int sta2)
{
    if (!active()) return 1.0;
    return (sta1 >= 100 ? value(pt1) : 1.0) * (sta2 >= 100 ? value(pt2) : 1.0);
}

// "" when off, so every existing filename is unchanged.
inline TString tag()
{
    for (const ModeName &m : modeNames()) if (m.mode == config().mode) return m.tag;
    return "";
}

inline void print()
{
    const Config &c = config();
    if (c.mode == Mode::Off) {
        std::cout << "SV b scale factor: OFF" << std::endl;
        return;
    }
    std::cout << "SV b scale factor: variation '" << tag() << "'"
              << " per b-matched SV vs aggregate pT, from " << file()
              << " (file tag '" << tag() << "')" << std::endl;
    for (int b = 1; b <= c.central->GetNbinsX(); ++b)
        std::cout << Form("   [%5.0f,%5.0f) GeV  SF = %.4f  -> %.4f",
                          c.central->GetXaxis()->GetBinLowEdge(b),
                          c.central->GetXaxis()->GetBinUpEdge(b),
                          c.central->GetBinContent(b),
                          value(c.central->GetBinCenter(b)))
                  << std::endl;
    if (!c.isMC)
        std::cout << "   -> this job is DATA, so no SF is applied" << std::endl;
}

} // namespace SvSF

#endif // SV_SF_H
