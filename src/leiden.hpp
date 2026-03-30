#pragma once
#include "graph.hpp"
#include <vector>

struct LeidenResult {
    std::vector<int> membership;  // membership[v] = community ID, 0-indexed contiguous
    double modularity;
    int n_communities;
};

// Build an igraph from the Graph with an optional edge mask, then run Leiden.
//   keep == nullptr  → use all edges
//   keep[i] == true  → include edge i
// Community IDs are renumbered to a contiguous 0..n_communities-1 range.
LeidenResult run_leiden(
    const Graph& g,
    const std::vector<bool>* keep = nullptr);
