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
