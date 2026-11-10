#pragma once
// Aggregation NN in step 1 (2026-10-22): the C++ forward pass of the network trained by
// aggregation_nn/agg_nn.py, so data and MC productions can build the two B candidates with it.
//
// Switch: environment variable AGG_NN = on | off (default off: nothing changes). With it on, the
// OBSERVABLES of a selected jet (dR, A, the EEC weight pT1 pT2, mB for the template fit) come from
// the NN's two candidates; the jet SELECTION (2 vertices from makeSvtxs_withBDT, b-tag) and the SV
// b SF weight stay those of the nominal aggregation, exactly as in the NN study. Productions made
// with it on must carry the "vetonn" era tag (aggnnInputs() in result_paths.h; the run scripts check).
//
// Mirrors agg_nn.py exactly: build() (inputs, transforms, standardisation, the 40 highest-pT
// tracks with pT >= 1 GeV), AggNet.forward() (post-norm transformer encoder, node and symmetric edge
// heads) and cluster() (node threshold, average-linkage agglomeration to two groups). Weights from
// cmd_export (agg_nn_weights.txt). Validated against the python friend files (see the note).
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include "TString.h"
// Needs tTree (tTree.h) and TrkEffSyst (tracking_efficiency_syst.h) declared BEFORE this header:
// it is included from create_files_for_template_fit.cpp after both. tTree.h has no include guard,
// so including it here a second time breaks the build.

namespace AggNN {

// DISABLED 2026-10-22 (Zoe): back to the nominal aggregation (BDT + 2SV) while the template fit
// and unfolding are reviewed. Always off; AGG_NN=on now refuses instead of silently doing nothing.
inline bool enabled()
{
    const char *e = std::getenv("AGG_NN");
    if (e && std::string(e) == "on") throw std::runtime_error("AGG_NN=on: the aggregation NN is disabled (agg_nn.h)");
    return false;
}
/* ---- disabled 2026-10-22: the env switch ----
inline bool enabled()
{
    static const bool on = [] {
        const char *e = std::getenv("AGG_NN");
        const std::string v = e ? e : "off";
        if (v != "on" && v != "off") throw std::runtime_error("AGG_NN must be on or off, got '" + v + "'");
        return v == "on";
    }();
    return on;
}
---- end disabled ---- */

inline std::string weightsFile()
{
    const char *e = std::getenv("AGG_NN_WEIGHTS");
    return e ? e : "/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/agg_nn/agg_nn_weights.txt";
}

struct Model {
    int nmax = 40; double thr = 0.5;
    std::vector<double> mu, sd;
    std::map<std::string, std::vector<double>> W;
    bool loaded = false;
    const std::vector<double> &w(const std::string &k) const {
        auto it = W.find(k);
        if (it == W.end()) throw std::runtime_error("AggNN: tensor " + k + " missing from " + weightsFile());
        return it->second;
    }
};

inline Model &model()
{
    static Model m;
    if (m.loaded) return m;
    std::ifstream in(weightsFile());
    if (!in) throw std::runtime_error("AggNN: cannot open " + weightsFile());
    std::string line;
    auto nums = [](const std::string &s) { std::vector<double> v; std::istringstream is(s); double x; while (is >> x) v.push_back(x); return v; };
    while (std::getline(in, line)) {
        std::istringstream is(line); std::string key; is >> key;
        std::string rest; std::getline(is, rest);
        if (key == "nmax") m.nmax = std::stoi(rest);
        else if (key == "thr") m.thr = std::stod(rest);
        else if (key == "mu") m.mu = nums(rest);
        else if (key == "sd") m.sd = nums(rest);
        else { std::string vals; std::getline(in, vals); m.W[key] = nums(vals); }
    }
    if (m.mu.size() != 24 || m.sd.size() != 24) throw std::runtime_error("AggNN: expected 24 input features");
    m.loaded = true;
    std::cout << "AggNN: loaded " << m.W.size() << " tensors from " << weightsFile() << " (node threshold " << m.thr << ")" << std::endl;
    return m;
}

constexpr int D = 64, H = 4, DH = 16, FF = 128, NF = 24, NL = 3;

// y = W x + b, W row-major (nout x nin)
inline void linear(const std::vector<double> &Wm, const std::vector<double> &b, int nout, int nin, const double *x, double *y)
{
    for (int o = 0; o < nout; ++o) { double s = b[o]; const double *r = &Wm[o * nin]; for (int i = 0; i < nin; ++i) s += r[i] * x[i]; y[o] = s; }
}

inline void layerNorm(double *x, const std::vector<double> &g, const std::vector<double> &b)
{
    double m = 0., v = 0.;
    for (int i = 0; i < D; ++i) m += x[i];
    m /= D;
    for (int i = 0; i < D; ++i) v += (x[i] - m) * (x[i] - m);
    v /= D;
    const double inv = 1. / std::sqrt(v + 1e-5);
    for (int i = 0; i < D; ++i) x[i] = (x[i] - m) * inv * g[i] + b[i];
}

inline double wrapPy(double x)          // numpy (x + pi) % (2 pi) - pi
{
    const double tp = 2. * M_PI; double y = std::fmod(x + M_PI, tp); if (y < 0.) y += tp; return y - M_PI;
}

inline double sigmoid(double x) { return 1. / (1. + std::exp(-x)); }

struct Result { std::vector<int> g1, g2; int status = -1; };   // track indices in the event's arrays

// The two B candidates of jet ijet. jtpt = jtpt_fill(t.jtpt[ijet]) (the folded value, as in training).
inline Result run(tTree &t, Int_t ijet, Long64_t ient, double jtpt)
{
    const Model &M = model();
    Result res;
    // tracks: pT >= 1 GeV in the jet, the tracking-efficiency variation applied (a no-op nominally),
    // sorted by pT (descending), the 40 highest kept
    std::vector<int> trk;
    for (Int_t i = 0; i < t.ntrk; ++i) {
        if (t.trkJetId[i] != ijet || t.trkPt[i] < 1) continue;
        if (TrkEffSyst::dropTrack(ient, i)) continue;
        trk.push_back(i);
    }
    std::stable_sort(trk.begin(), trk.end(), [&](int a, int b) { return t.trkPt[a] > t.trkPt[b]; });
    if ((int) trk.size() > M.nmax) trk.resize(M.nmax);
    const int n = trk.size();
    if (n < 2) return res;

    // ---- inputs (agg_nn.py build(): same order, same transforms; NaN/inf -> 0 before scaling)
    const double jeta = t.jteta[ijet], jphi = t.jtphi[ijet];
    std::vector<double> X(n * NF), deta(n), dphi(n);
    std::vector<int> svid(n);
    for (int k = 0; k < n; ++k) {
        const int i = trk[k];
        const double pt = t.trkPt[i];
        const int sv = t.trkSvtxId[i];
        const bool insv = (sv >= 0 && sv < t.nsvtx);
        svid[k] = insv ? sv : -1;
        deta[k] = t.trkEta[i] - jeta; dphi[k] = wrapPy(t.trkPhi[i] - jphi);
        const int pdg = std::abs(t.trkPdgId[i]);
        const double bdt = t.trkBdtScore[i];
        const double f[NF] = {
            std::log(pt), pt / jtpt, deta[k], dphi[k], std::hypot(deta[k], dphi[k]),
            std::asinh(t.trkIp3d[i] * 100.), std::asinh(t.trkIp3dSig[i]), std::asinh(t.trkDistToAxis[i] * 100.),
            bdt > -900 ? bdt : 0., bdt > -900 ? 1. : 0.,
            pdg == 11 ? 1. : 0., pdg == 13 ? 1. : 0., pdg == 211 ? 1. : 0., insv ? 1. : 0.,
            insv ? t.svtxm[sv] : 0., insv ? t.svtxmcorr[sv] : 0., insv ? std::asinh(t.svtxpt[sv]) : 0.,
            insv ? std::asinh(t.svtxdl[sv] * 10.) : 0., insv ? std::asinh(t.svtxdls[sv]) : 0.,
            insv ? std::asinh(t.svtxdl2d[sv] * 10.) : 0., insv ? std::asinh(t.svtxdls2d[sv]) : 0.,
            insv ? double(t.svtxNtrk[sv]) : 0., insv ? std::asinh(t.svtxnormchi2[sv]) : 0., std::log(jtpt)};
        for (int j = 0; j < NF; ++j) {
            const double v = std::isfinite(f[j]) ? f[j] : 0.;
            // float32 inputs, as in the python training (pad(...).astype(float32))
            X[k * NF + j] = (double(float(v)) - M.mu[j]) / M.sd[j];
        }
    }

    // ---- embedding: Linear(24, 64) -> ReLU -> Linear(64, 64)
    std::vector<double> Hm(n * D), tmp(D);
    for (int k = 0; k < n; ++k) {
        linear(M.w("inp.0.weight"), M.w("inp.0.bias"), D, NF, &X[k * NF], tmp.data());
        for (auto &v : tmp) v = std::max(v, 0.);
        linear(M.w("inp.2.weight"), M.w("inp.2.bias"), D, D, tmp.data(), &Hm[k * D]);
    }
    // ---- encoder: 3 post-norm layers, x = norm1(x + SA(x)); x = norm2(x + FF(x))
    std::vector<double> qkv(n * 3 * D), att(n * D), sa(n * D), ff1(FF), ff2(D), sc(n);
    for (int L = 0; L < NL; ++L) {
        const std::string p = "enc.layers." + std::to_string(L) + ".";
        for (int k = 0; k < n; ++k) linear(M.w(p + "self_attn.in_proj_weight"), M.w(p + "self_attn.in_proj_bias"), 3 * D, D, &Hm[k * D], &qkv[k * 3 * D]);
        for (int h = 0; h < H; ++h)
            for (int i = 0; i < n; ++i) {
                double mx = -std::numeric_limits<double>::infinity();
                for (int j = 0; j < n; ++j) {
                    double s = 0.;
                    for (int c = 0; c < DH; ++c) s += qkv[i * 3 * D + h * DH + c] * qkv[j * 3 * D + D + h * DH + c];
                    sc[j] = s / std::sqrt(double(DH)); mx = std::max(mx, sc[j]);
                }
                double z = 0.;
                for (int j = 0; j < n; ++j) { sc[j] = std::exp(sc[j] - mx); z += sc[j]; }
                for (int c = 0; c < DH; ++c) {
                    double s = 0.;
                    for (int j = 0; j < n; ++j) s += sc[j] * qkv[j * 3 * D + 2 * D + h * DH + c];
                    att[i * D + h * DH + c] = s / z;
                }
            }
        for (int k = 0; k < n; ++k) {
            linear(M.w(p + "self_attn.out_proj.weight"), M.w(p + "self_attn.out_proj.bias"), D, D, &att[k * D], &sa[k * D]);
            double *x = &Hm[k * D];
            for (int c = 0; c < D; ++c) x[c] += sa[k * D + c];
            layerNorm(x, M.w(p + "norm1.weight"), M.w(p + "norm1.bias"));
            linear(M.w(p + "linear1.weight"), M.w(p + "linear1.bias"), FF, D, x, ff1.data());
            for (auto &v : ff1) v = std::max(v, 0.);
            linear(M.w(p + "linear2.weight"), M.w(p + "linear2.bias"), D, FF, ff1.data(), ff2.data());
            for (int c = 0; c < D; ++c) x[c] += ff2[c];
            layerNorm(x, M.w(p + "norm2.weight"), M.w(p + "norm2.bias"));
        }
    }
    // ---- heads
    std::vector<double> pn(n), h1(64), e(3 * D + 2);
    for (int k = 0; k < n; ++k) {
        linear(M.w("node.0.weight"), M.w("node.0.bias"), 64, D, &Hm[k * D], h1.data());
        for (auto &v : h1) v = std::max(v, 0.);
        double o; linear(M.w("node.2.weight"), M.w("node.2.bias"), 1, 64, h1.data(), &o);
        pn[k] = sigmoid(o);
    }
    // track selection; fewer than two -> the two highest-scoring tracks
    std::vector<int> sel;
    for (int k = 0; k < n; ++k) if (pn[k] > M.thr) sel.push_back(k);
    res.status = 0;
    if (sel.size() < 2) {
        std::vector<int> o(n); for (int k = 0; k < n; ++k) o[k] = k;
        std::stable_sort(o.begin(), o.end(), [&](int a, int b) { return pn[a] > pn[b]; });
        sel = {o[0], o[1]}; res.status = 1;
    }
    const int m = sel.size();
    std::vector<double> S(m * m, -std::numeric_limits<double>::infinity());
    for (int a = 0; a < m; ++a)
        for (int b = a + 1; b < m; ++b) {
            const int i = sel[a], j = sel[b];
            for (int c = 0; c < D; ++c) {
                const double hi = Hm[i * D + c], hj = Hm[j * D + c];
                e[c] = hi * hj; e[D + c] = hi + hj; e[2 * D + c] = std::fabs(hi - hj);
            }
            e[3 * D] = (svid[i] >= 0 && svid[i] == svid[j]) ? 1. : 0.;
            double dd = dphi[i] - dphi[j] + M_PI; dd -= 2. * M_PI * std::floor(dd / (2. * M_PI)); dd -= M_PI;   // torch.remainder
            e[3 * D + 1] = std::hypot(deta[i] - deta[j], dd);
            linear(M.w("edge.0.weight"), M.w("edge.0.bias"), 64, 3 * D + 2, e.data(), h1.data());
            for (auto &v : h1) v = std::max(v, 0.);
            double o; linear(M.w("edge.2.weight"), M.w("edge.2.bias"), 1, 64, h1.data(), &o);
            S[a * m + b] = S[b * m + a] = sigmoid(o);
        }
    // ---- average-linkage agglomeration to two clusters (agg_nn.py cluster())
    std::vector<std::vector<int>> cl(m);
    std::vector<double> size(m, 1.);
    std::vector<bool> alive(m, true);
    for (int a = 0; a < m; ++a) cl[a] = {a};
    int nal = m;
    while (nal > 2) {
        double best = -std::numeric_limits<double>::infinity(); int ba = -1, bb = -1;
        for (int a = 0; a < m; ++a) { if (!alive[a]) continue;
            for (int b = 0; b < m; ++b) { if (!alive[b] || a == b) continue;
                if (S[a * m + b] > best) { best = S[a * m + b]; ba = a; bb = b; } } }
        if (ba < 0) { ba = -1; for (int a = 0; a < m && bb < 0; ++a) if (alive[a]) { if (ba < 0) ba = a; else bb = a; } }
        for (int c = 0; c < m; ++c) S[ba * m + c] = (size[ba] * S[ba * m + c] + size[bb] * S[bb * m + c]) / (size[ba] + size[bb]);
        for (int c = 0; c < m; ++c) S[c * m + ba] = S[ba * m + c];
        S[ba * m + ba] = -std::numeric_limits<double>::infinity();
        cl[ba].insert(cl[ba].end(), cl[bb].begin(), cl[bb].end()); size[ba] += size[bb]; alive[bb] = false; --nal;
    }
    std::vector<std::vector<int>> groups;
    for (int a = 0; a < m; ++a) if (alive[a]) { std::vector<int> g; for (int k : cl[a]) g.push_back(trk[sel[k]]); groups.push_back(g); }
    res.g1 = groups[0]; res.g2 = groups[1];
    return res;
}

} // namespace AggNN
