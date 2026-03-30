#pragma once
#include "dspar.hpp"
#include "rng.hpp"
#include <vector>

// Returns keep[i] = true if edge i (canonical index) is retained.
// Retention probability: p(i) = alpha * scores.score[i] / scores.mean
// Clamped to [0, 1] — edges where p(i) >= 1 are always kept.
std::vector<bool> sample_edges(
    const DsparScores& scores,
    double alpha,
    Xoshiro256pp& rng);

// In-place variant — writes into a pre-allocated buffer (avoids allocation per trial)
void sample_edges_inplace(
    const DsparScores& scores,
    double alpha,
    Xoshiro256pp& rng,
    std::vector<bool>& keep_buf);

// ── Opt B: Bernoulli sampler with pre-computed thresholds ─────────────────────
// threshold[i] = min(1.0, alpha * score[i] / mean) — computed once per alpha.
// Eliminates the per-edge division inside the hot loop.
void sample_edges_bernoulli(
    const std::vector<double>& threshold,
    Xoshiro256pp& rng,
    std::vector<bool>& keep_buf);

// ── Opt D: prefix-sum (PPS) sampler ──────────────────────────────────────────
// Retains exactly k edges sampled proportional to their scores.
// Faster than m Bernoulli draws when alpha is very small (most edges rejected).
//   score_cdf : cumulative sum of scores, size m+1, cdf[0]=0, cdf[m]=total_score
//   total_score: sum of all edge scores
//   k          : number of edges to retain = round(alpha * m)
// Draw k uniform values in [0, total_score], binary-search in CDF → edge index.
// Duplicates are collapsed (an edge is either kept or not).
void sample_edges_prefix_sum(
    const std::vector<double>& score_cdf,
    double total_score,
    int k,
    Xoshiro256pp& rng,
    std::vector<bool>& keep_buf);
