#include "dspar.hpp"
#include <stdexcept>

DsparScores compute_scores(const Graph& g) {
    if (g.m == 0)
        throw std::runtime_error("compute_scores: graph has no edges");

    DsparScores s;
    s.score.resize(g.m);

    // Compute per-edge scores and accumulate mean with Kahan summation
    double kahan_sum  = 0.0;
    double kahan_comp = 0.0;
    for (int i = 0; i < g.m; i++) {
        auto [u, v] = g.edges[i];
        s.score[i] = 1.0 / g.deg[u] + 1.0 / g.deg[v];

        double y = s.score[i] - kahan_comp;
        double t = kahan_sum + y;
        kahan_comp = (t - kahan_sum) - y;
        kahan_sum  = t;
    }
    s.mean = kahan_sum / g.m;
    return s;
}

void compute_structural_metrics(
    DsparScores& scores,
    const Graph& g,
    const std::vector<int>& membership)
{
    double sum_intra  = 0.0, comp_intra = 0.0;
    double sum_inter  = 0.0, comp_inter = 0.0;
    int n1 = 0, n2 = 0;

    for (int i = 0; i < g.m; i++) {
        auto [u, v] = g.edges[i];
        double sc = scores.score[i];
        if (membership[u] == membership[v]) {
            // Kahan accumulation for intra
            double y = sc - comp_intra;
            double t = sum_intra + y;
            comp_intra = (t - sum_intra) - y;
            sum_intra  = t;
            n1++;
        } else {
            // Kahan accumulation for inter
            double y = sc - comp_inter;
            double t = sum_inter + y;
            comp_inter = (t - sum_inter) - y;
            sum_inter  = t;
            n2++;
        }
    }

    scores.n1       = n1;
    scores.n2       = n2;
    scores.mu_intra = (n1 > 0) ? sum_intra / n1 : 0.0;
    scores.mu_inter = (n2 > 0) ? sum_inter / n2 : 0.0;
    scores.delta    = scores.mu_intra - scores.mu_inter;
}
