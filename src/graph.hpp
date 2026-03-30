#pragma once
#include <string>
#include <vector>
#include <utility>

struct Graph {
    int n;                                   // number of nodes
    int m;                                   // number of undirected edges
    std::vector<int> deg;                    // deg[v] = degree of node v
    std::vector<int> adj_off;               // CSR offsets, size n+1
    std::vector<int> adj;                   // CSR neighbor list, size 2*m
    std::vector<std::pair<int,int>> edges;  // canonical edge list, u < v, size m
};

// Load a SNAP-format edge list file.
// - Skips lines starting with '#'
// - Treats edges as undirected: stores (min(s,d), max(s,d))
// - Remaps node IDs to contiguous 0-indexed range
// - Removes self-loops
Graph load_snap(const std::string& path);
