#ifndef TRACKING_EFFICIENCY_SYST_H
#define TRACKING_EFFICIENCY_SYST_H

// ============================================================================
// Tracking-efficiency systematic
// ============================================================================
//
// Track reconstruction in the core of a jet is prone to mismodelling in
// simulation. To propagate that through the unfolding we randomly throw away a
// fraction (nominally 3%) of the RECONSTRUCTED tracks in MC during the
// substructure extraction, redo the whole chain, and take the symmetrised
// difference of the unfolded distribution against the nominal one as the
// uncertainty.
//
// Applies to MC only, and only to the reco-level tracks (trk*). The gen-level
// b hadrons (refTrk*) are the truth and are never touched -- the variation is
// a detector effect, so it must move reco while gen stays put, which is exactly
// what makes the response matrix and the corrections change.
//
// ---------------------------------------------------------------------------
// WHY THE RANDOM NUMBER IS HASHED RATHER THAN DRAWN FROM A STREAM
// ---------------------------------------------------------------------------
// "Repeatable" has to mean more than "seeded". A seeded std::mt19937 consumed
// as a stream gives the same answer only if the tracks are consumed in exactly
// the same order, by exactly one process, over exactly the same event range.
// This analysis violates all three: run_agg_ntuple_chunks.sh runs one job per
// block in parallel, the number of blocks changes per sample (10 qcd/pythia,
// 9 bjet/pythia, 8 qcd/herwig, 9 bjet/herwig), and the Condor data path splits
// differently again. With a stream, re-chunking the sample silently changes
// which tracks are dropped, and the systematic moves with it.
//
// So instead of drawing the n-th number of a stream, the random number for a
// track is COMPUTED from the track's identity:
//
//     u = hash(seed, run/lumi/event, track index)  ->  uniform in [0,1)
//     drop the track if u >= 1 - dropFraction    (3% -> drop if u >= 0.97)
//
// The number attached to a given track is then a property of that track, not of
// the processing history. Consequences, all of them wanted:
//
//   * running the same job twice gives bit-identical output;
//   * running with a different number of blocks, in a different order, or in
//     parallel vs serial gives the same output;
//   * resuming or re-running one failed block reproduces exactly what the
//     first attempt would have written;
//   * the nominal and the varied production see the *same* events, so their
//     difference is the track dropping and nothing else.
//
// The hash is splitmix64: fixed integer arithmetic, so it is also identical
// across compilers and machines, which std::mt19937 seeding conventions and
// std::hash are not.
//
// (ient, itrk) identifies a track within one input file. The same pair in two
// different blocks gets the same u, which is harmless: those are unrelated
// events and the decision is uncorrelated with anything physical. Deliberately
// NOT salted with the file name, so re-running over moved or renamed inputs
// still reproduces the published numbers.
//
// ---------------------------------------------------------------------------
// HOW TO USE
// ---------------------------------------------------------------------------
// One boolean, everywhere: track_eff_unc = false is the nominal (every function
// below is a no-op, output bit-identical to before this header existed),
// track_eff_unc = true is the systematic. Set it with
//
//     TrkEffSyst::enable(track_eff_unc, isMC);
//
// The 3% itself is kDropFraction below -- one constant, not a number spelled
// into each caller, so there is nothing to keep in step by hand.
//
// create_files_for_template_fit() takes both as arguments; run_agg_ntuple_chunks.sh
// and make_hardprobes_condor_scripts.sh set them from one knob at the top and
// tag the output files so a variation can never overwrite the nominal.
// ============================================================================

#include <cstdint>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "Rtypes.h"
#include "TString.h"

namespace TrkEffSyst {

// ---- The prescription ------------------------------------------------------
// 3% covers the residual data/MC difference for tracking in the jet core seen in
// the jet energy scale determination. Change it here and nowhere else -- but note
// that the output-file tag below is derived from it, so a different value writes
// differently-named files (_trkdrop050 for 5%) and the readers in
// apply_unfolding_2d.C / apply_weights_and_systematics.C need the new tag too.
constexpr double kDropFraction = 0.03;

// ---- Configuration, set once per job before the event loop -----------------
struct Config {
    double             dropFraction = 0.0;         // 0 = nominal (no dropping)
    unsigned long long seed         = 20260908ULL; // fixed: change it only to
                                                   // produce an independent throw
    bool               isMC         = false;       // data is never varied
    unsigned long long splitSeed    = 12345ULL;    // half0/half1 jet split
    unsigned long long eventKey     = 0ULL;         // hash(run,lumi,event)
    bool               hasEventKey  = false;
};

inline Config &config()
{
    static Config c;
    return c;
}

// The single switch. track_eff_unc = false leaves everything nominal.
inline void enable(bool track_eff_unc, bool isMC, unsigned long long seed = 20260908ULL)
{
    Config &c = config();
    c.dropFraction = track_eff_unc ? kDropFraction : 0.0;
    c.seed         = seed;
    c.isMC         = isMC;
}

// ---- splitmix64: 64 bits in, 64 well-mixed bits out ------------------------
inline uint64_t mix64(uint64_t x)
{
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

// Uniform in [0,1) from three integers. Every bit of every input is mixed into
// the result, so neighbouring (ient, itrk) do not give neighbouring u.
inline double uniformFrom(uint64_t a, uint64_t b, uint64_t c)
{
    uint64_t h = mix64(a + 0x9E3779B97F4A7C15ULL);
    h = mix64(h ^ mix64(b + 0x632BE59BD9B4E019ULL));
    h = mix64(h ^ mix64(c + 0xD6E8FEB86659FD93ULL));
    // Top 53 bits -> a double with full mantissa resolution, in [0,1).
    return static_cast<double>(h >> 11) * (1.0 / 9007199254740992.0);
}

// Set once after reading each event.  Unlike a local TTree entry number, this
// identity is unchanged by chunking, merging, or processing only a sub-range.
inline void beginEvent(UInt_t run, UInt_t lumi, ULong64_t event)
{
    Config &c = config();
    c.eventKey = mix64(static_cast<uint64_t>(run));
    c.eventKey = mix64(c.eventKey ^ mix64(static_cast<uint64_t>(lumi)));
    c.eventKey = mix64(c.eventKey ^ mix64(static_cast<uint64_t>(event)));
    c.hasEventKey = true;
}

inline uint64_t eventKey()
{
    if (!config().hasEventKey)
        throw std::logic_error("TrkEffSyst::beginEvent must be called before drawing random numbers");
    return config().eventKey;
}

// ---- The decision ----------------------------------------------------------
// true  = this reconstructed track is thrown away by the variation.
// beginEvent() supplies the stable event identity; itrk indexes the trk* arrays.
// Deliberately independent of the jet index: a track keeps its verdict no
// matter which jet's loop happens to be looking at it.
inline bool dropTrack(Long64_t ient, Int_t itrk)
{
    const Config &c = config();
    if (c.dropFraction <= 0.) return false;   // nominal: no-op
    if (!c.isMC)              return false;   // data is never varied
    const double u = uniformFrom(c.seed, eventKey(),
                                 static_cast<uint64_t>(static_cast<uint32_t>(itrk)));
    return u >= (1.0 - c.dropFraction);       // 3% -> dropped when u >= 0.97
}

// ---- Deterministic jet half-split ------------------------------------------
// Replaces a std::random_device-seeded mt19937 that made two runs of the SAME
// job disagree in the half0/half1 histograms. Hashing (event identity, ijet) makes the
// half a jet lands in a property of the jet, so it is stable across reruns and
// across re-chunking, and identical between the nominal and the varied
// production -- otherwise the split noise would leak into the systematic.
inline double jetSplitRandom(Long64_t ient, Int_t ijet)
{
    const Config &c = config();
    return uniformFrom(c.splitSeed, eventKey(),
                       0x5FF1CEULL ^ static_cast<uint64_t>(static_cast<uint32_t>(ijet)));
}

// ---- Output-file tag -------------------------------------------------------
// "" for the nominal, so nominal filenames are unchanged. Otherwise the drop
// fraction in PER MILLE: 3% -> "_trkdrop030".
inline TString tag()
{
    const Config &c = config();
    if (c.dropFraction <= 0.) return TString("");
    return TString(Form("_trkdrop%03d", static_cast<int>(std::lround(1000. * c.dropFraction))));
}

inline void print()
{
    const Config &c = config();
    if (c.dropFraction <= 0.) {
        std::cout << "Tracking-efficiency systematic: OFF (nominal)" << std::endl;
        return;
    }
    std::cout << "Tracking-efficiency systematic: dropping "
              << 100. * c.dropFraction << "% of reco tracks"
              << " (seed " << c.seed << ", isMC " << c.isMC
              << ", tag '" << tag() << "')" << std::endl;
    if (!c.isMC)
        std::cout << "   -> this job is DATA, so nothing will actually be dropped"
                  << std::endl;
}

} // namespace TrkEffSyst

#endif // TRACKING_EFFICIENCY_SYST_H
