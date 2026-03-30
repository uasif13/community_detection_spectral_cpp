#pragma once
#include "graph.hpp"
#include <vector>

struct ModResult {
    double F;             // intra-edge fraction
    double G;             // null-model penalty
    double Q;             // F - G
    int m_sparse;         // edge count in sparse graph
    int preserved_intra;  // intra edges kept
    int preserved_inter;  // inter edges kept
    double intra_rate;    // preserved_intra / n1  (0 if n1==0)
    double inter_rate;    // preserved_inter / n2  (0 if n2==0)
    double ratio_obs;     // inter_rate / intra_rate  (0 if intra_rate==0)
};

// Compute modularity with fixed membership on a (possibly sparse) graph.
//   keep == nullptr  → treat all edges as present
//   keep[i] == true  → edge i is present
//   n1, n2           → intra/inter edge counts in the ORIGINAL graph (from DsparScores)
ModResult compute_modularity(
    const Graph& g,
    const std::vector<bool>* keep,
    const std::vector<int>& member,
    int n1,
    int n2);
