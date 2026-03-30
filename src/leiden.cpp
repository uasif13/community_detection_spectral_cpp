#include "leiden.hpp"
#include "modularity.hpp"
#include <igraph/igraph.h>
#include <map>
#include <stdexcept>

LeidenResult run_leiden(const Graph& g, const std::vector<bool>* keep) {
    // ── 1. Build flat edge list ───────────────────────────────────────────────
    std::vector<igraph_integer_t> edge_list;
    edge_list.reserve(2 * g.m);
    for (int i = 0; i < g.m; i++) {
        if (keep && !(*keep)[i]) continue;
        edge_list.push_back(static_cast<igraph_integer_t>(g.edges[i].first));
        edge_list.push_back(static_cast<igraph_integer_t>(g.edges[i].second));
    }

    // ── 2. Construct igraph ───────────────────────────────────────────────────
    igraph_t ig;
    igraph_vector_int_t edges_vec;
    igraph_vector_int_init_array(
        &edges_vec,
        edge_list.data(),
        static_cast<igraph_integer_t>(edge_list.size()));
    igraph_error_t rc = igraph_create(
        &ig, &edges_vec,
        static_cast<igraph_integer_t>(g.n),
        IGRAPH_UNDIRECTED);
    igraph_vector_int_destroy(&edges_vec);
    if (rc != IGRAPH_SUCCESS)
        throw std::runtime_error("igraph_create failed");

    // ── 3. Build degree-based node weights for modularity optimisation ────────
    // igraph's Leiden quality: Q = Σ_ij (A_ij - γ·n_i·n_j)·δ(σ_i,σ_j)
    // Setting n_i = deg(i) and γ = 1/(2m) recovers standard modularity × 2m.
    // This makes Leiden maximise standard modularity rather than CPM.
    int m_edges = static_cast<int>(edge_list.size() / 2);
    igraph_vector_t node_weights;
    igraph_vector_init(&node_weights, g.n);
    igraph_vector_fill(&node_weights, 0.0);
    for (int i = 0; i < g.m; i++) {
        if (keep && !(*keep)[i]) continue;
        VECTOR(node_weights)[g.edges[i].first]  += 1.0;
        VECTOR(node_weights)[g.edges[i].second] += 1.0;
    }
    double resolution = (m_edges > 0) ? 1.0 / (2.0 * m_edges) : 1.0;

    // ── 4. Run Leiden ─────────────────────────────────────────────────────────
    igraph_vector_int_t membership_vec;
    igraph_vector_int_init(&membership_vec, 0);
    igraph_integer_t nb_clusters = 0;
    igraph_real_t quality_val = 0.0;

    rc = igraph_community_leiden(
        &ig,
        /*edge_weights=*/nullptr,
        &node_weights,
        /*resolution_parameter=*/resolution,
        /*beta=*/0.01,
        /*start=*/false,
        /*n_iterations=*/2,
        &membership_vec,
        &nb_clusters,
        &quality_val);

    igraph_vector_destroy(&node_weights);
    igraph_destroy(&ig);

    if (rc != IGRAPH_SUCCESS) {
        igraph_vector_int_destroy(&membership_vec);
        throw std::runtime_error("igraph_community_leiden failed");
    }

    // ── 5. Extract and renumber membership ───────────────────────────────────
    LeidenResult result;
    result.membership.resize(g.n);
    for (int v = 0; v < g.n; v++)
        result.membership[v] = static_cast<int>(VECTOR(membership_vec)[v]);
    igraph_vector_int_destroy(&membership_vec);

    // Renumber to contiguous 0..k-1 (igraph 0.10 should already do this,
    // but we guarantee it here to be safe)
    std::map<int,int> remap;
    int next = 0;
    for (auto& c : result.membership)
        if (!remap.count(c)) remap[c] = next++;
    for (auto& c : result.membership) c = remap[c];

    result.n_communities = next;

    // Compute standard modularity from our formula (authoritative, not igraph's
    // quality which is in different units when node_weights != uniform)
    int n1 = 0, n2 = 0;
    for (int i = 0; i < g.m; i++) {
        if (keep && !(*keep)[i]) continue;
        auto [u, v] = g.edges[i];
        if (result.membership[u] == result.membership[v]) n1++; else n2++;
    }
    ModResult mod = compute_modularity(g, keep, result.membership, n1, n2);
    result.modularity = mod.Q;
    return result;
}
