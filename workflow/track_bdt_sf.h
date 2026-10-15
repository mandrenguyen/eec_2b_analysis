#ifndef TRACK_BDT_SF_H
#define TRACK_BDT_SF_H

#include <iostream>
#include <stdexcept>
#include <string>

#include "Rtypes.h"
#include "TString.h"

namespace TrackBdtSF {

constexpr double kNominalWP = 0.365;

// Effective working points derived in selected double-b jets
// (100 < jet pT < 200 GeV, |eta| < 2, UParTV2 > 0.712, at least two SVs).
constexpr double kSignalInSvWP     = 0.5645440855;
constexpr double kSignalOutSvWP    = 0.7255153474;
constexpr double kNonSignalInSvWP  = 0.3601433651;
constexpr double kNonSignalOutSvWP = 0.3246904670;

enum class Mode { Off, Nominal, BOff, LOff };

struct Config {
    Mode mode = Mode::Off;
    bool isMC = false;
};

inline Config &config()
{
    static Config c;
    return c;
}

inline const char *modeName(Mode mode)
{
    switch (mode) {
        case Mode::Off:     return "off";
        case Mode::Nominal: return "nominal";
        case Mode::BOff:    return "boff";
        case Mode::LOff:    return "loff";
    }
    return "invalid";
}

inline Mode parseMode(const char *value)
{
    const std::string mode = value ? value : "";
    if (mode == "off")     return Mode::Off;
    if (mode == "nominal") return Mode::Nominal;
    if (mode == "boff")    return Mode::BOff;
    if (mode == "loff")    return Mode::LOff;
    throw std::runtime_error("TRACK_BDT_SF must be off, nominal, boff, or loff; got '" + mode + "'");
}

inline void enable(const char *mode, bool isMC)
{
    Config &c = config();
    c.mode = parseMode(mode);
    c.isMC = isMC;
}

inline bool applySignal()
{
    const Config &c = config();
    return c.isMC && (c.mode == Mode::Nominal || c.mode == Mode::LOff);
}

inline bool applyNonSignal()
{
    const Config &c = config();
    return c.isMC && (c.mode == Mode::Nominal || c.mode == Mode::BOff);
}

inline double effectiveWP(Int_t matchStatus, Int_t svtxId)
{
    const bool signal = matchStatus >= 100;
    const bool inSv = svtxId >= 0;
    if (signal && applySignal())
        return inSv ? kSignalInSvWP : kSignalOutSvWP;
    if (!signal && applyNonSignal())
        return inSv ? kNonSignalInSvWP : kNonSignalOutSvWP;
    return kNominalWP;
}

inline bool passes(double score, Int_t matchStatus, Int_t svtxId)
{
    return score > effectiveWP(matchStatus, svtxId);
}

// Data always has an empty tag because truth categories are unavailable and no shift is
// applied. For MC, "off" retains the historical filenames; corrected productions are
// explicitly tagged so they cannot overwrite the uncorrected cross-check.
inline TString tag()
{
    const Config &c = config();
    if (!c.isMC || c.mode == Mode::Off) return "";
    if (c.mode == Mode::Nominal) return "_trkbdtsf";
    if (c.mode == Mode::BOff)    return "_trkbdtsf_boff";
    if (c.mode == Mode::LOff)    return "_trkbdtsf_loff";
    return "";
}

inline void print()
{
    const Config &c = config();
    std::cout << "Track-BDT scale-factor mode: " << modeName(c.mode)
              << " (isMC " << c.isMC << ", tag '" << tag() << "')" << std::endl;
    if (!c.isMC) {
        std::cout << "   -> data scores use the nominal 0.365 working point" << std::endl;
        return;
    }
    std::cout << "   signal WP in/out SV: "
              << (applySignal() ? kSignalInSvWP : kNominalWP) << " / "
              << (applySignal() ? kSignalOutSvWP : kNominalWP) << std::endl;
    std::cout << "   non-signal WP in/out SV: "
              << (applyNonSignal() ? kNonSignalInSvWP : kNominalWP) << " / "
              << (applyNonSignal() ? kNonSignalOutSvWP : kNominalWP) << std::endl;
}

} // namespace TrackBdtSF

#endif
