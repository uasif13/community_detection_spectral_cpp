#include "sampler.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>

void sample_edges_inplace(
    const DsparScores& scores,
    double alpha,
    Xoshiro256pp& rng,
    std::vector<bool>& keep_buf)
{
    const int m = static_cast<int>(scores.score.size());
    keep_buf.resize(m);
    for (int i = 0; i < m; i++) {
        double p = alpha * scores.score[i] / scores.mean;
        assert(p <= 1.0 + 1e-12);   // clamping guard (debug builds)
        if (p >= 1.0) {
            keep_buf[i] = true;
        } else {
            keep_buf[i] = rng.next_double() < p;
        }
    }
}

std::vector<bool> sample_edges(
    const DsparScores& scores,
    double alpha,
    Xoshiro256pp& rng)
{
    std::vector<bool> keep;
    sample_edges_inplace(scores, alpha, rng, keep);
    return keep;
}

// ── Opt B ─────────────────────────────────────────────────────────────────────

void sample_edges_bernoulli(
    const std::vector<double>& threshold,
    Xoshiro256pp& rng,
    std::vector<bool>& keep_buf)
{
    const int m = static_cast<int>(threshold.size());
    keep_buf.resize(m);
    for (int i = 0; i < m; i++)
        keep_buf[i] = rng.next_double() < threshold[i];
}

// ── Opt D ─────────────────────────────────────────────────────────────────────

void sample_edges_prefix_sum(
    const std::vector<double>& score_cdf,
    double total_score,
    int k,
    Xoshiro256pp& rng,
    std::vector<bool>& keep_buf)
{
    const int m = static_cast<int>(score_cdf.size()) - 1;
    keep_buf.assign(m, false);
    if (k <= 0 || total_score <= 0.0) return;

    // Draw k uniform values, binary-search in CDF to find edge index.
    // Duplicates are collapsed naturally (keep_buf is bool).
    for (int draw = 0; draw < k; draw++) {
        double u = rng.next_double() * total_score;
        // lower_bound finds first cdf entry > u → edge index = that position - 1
        int idx = static_cast<int>(
            std::upper_bound(score_cdf.begin(), score_cdf.end(), u)
            - score_cdf.begin()) - 1;
        if (idx < 0)  idx = 0;
        if (idx >= m) idx = m - 1;
        keep_buf[idx] = true;
    }
}
