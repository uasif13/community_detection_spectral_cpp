#include "sampler.hpp"
#include <cassert>

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
