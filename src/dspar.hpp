#pragma once
#include "graph.hpp"
#include <vector>

struct DsparScores {
    std::vector<double> score;  // score[i] = 1/deg(u) + 1/deg(v) for edges[i]
    double mean;                // s̄ = mean of all scores (Kahan summation)
    // Structural metrics w.r.t. a partition
    double mu_intra = 0.0;     // mean score of intra-community edges
    double mu_inter = 0.0;     // mean score of inter-community edges
    double delta    = 0.0;     // mu_intra - mu_inter
    int n1 = 0;                // intra-community edge count
    int n2 = 0;                // inter-community edge count
};

// Compute per-edge scores and mean from graph degrees (no partition needed)
DsparScores compute_scores(const Graph& g);

// Fill structural metrics (mu_intra, mu_inter, delta, n1, n2) given a partition
void compute_structural_metrics(
    DsparScores& scores,
    const Graph& g,
    const std::vector<int>& membership);
