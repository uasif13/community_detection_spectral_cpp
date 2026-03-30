#include "graph.hpp"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <set>

Graph load_snap(const std::string& path) {
    std::ifstream fin(path);
    if (!fin.is_open())
        throw std::runtime_error("Cannot open file: " + path);

    // Collect unique undirected edges (no self-loops) and all node IDs
    std::set<std::pair<int,int>> edge_set;
    std::set<int> node_set;

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        int s, d;
        if (!(iss >> s >> d)) continue;
        if (s == d) continue;  // self-loop
        node_set.insert(s);
        node_set.insert(d);
        edge_set.insert({std::min(s, d), std::max(s, d)});
    }

    // Build contiguous 0-indexed remapping
    std::unordered_map<int,int> remap;
    remap.reserve(node_set.size());
    int idx = 0;
    for (int v : node_set)
        remap[v] = idx++;

    Graph g;
    g.n = static_cast<int>(node_set.size());
    g.m = static_cast<int>(edge_set.size());

    // Build canonical edge list with remapped IDs
    g.edges.reserve(g.m);
    for (auto& [u, v] : edge_set)
        g.edges.push_back({remap[u], remap[v]});

    // Sort edges for deterministic order (already sorted because edge_set is ordered,
    // but remap may reorder — sort again)
    std::sort(g.edges.begin(), g.edges.end());

    // Build degree array
    g.deg.assign(g.n, 0);
    for (auto& [u, v] : g.edges) {
        g.deg[u]++;
        g.deg[v]++;
    }

    // Build CSR
    g.adj_off.resize(g.n + 1, 0);
    for (int v = 0; v < g.n; v++)
        g.adj_off[v + 1] = g.adj_off[v] + g.deg[v];

    g.adj.resize(2 * g.m);
    std::vector<int> pos(g.n, 0);  // current insertion position per node
    for (auto& [u, v] : g.edges) {
        g.adj[g.adj_off[u] + pos[u]++] = v;
        g.adj[g.adj_off[v] + pos[v]++] = u;
    }

    return g;
}
