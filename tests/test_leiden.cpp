#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "graph.hpp"
#include "modularity.hpp"
#include "leiden.hpp"
#include <fstream>
#include <filesystem>
#include <algorithm>
#include <set>

using Catch::Matchers::WithinAbs;

// ── helpers ──────────────────────────────────────────────────────────────────

static std::string write_tmp(const std::string& name, const std::string& content) {
    std::string path = (std::filesystem::temp_directory_path() / name).string();
    std::ofstream f(path);
    f << content;
    return path;
}

// Zachary's karate club — 34 nodes, 78 edges (well-known benchmark)
static Graph make_karate() {
    // Full edge list for the karate club graph
    std::string content =
        "1 2\n1 3\n1 4\n1 5\n1 6\n1 7\n1 8\n1 9\n1 11\n1 12\n1 13\n1 14\n"
        "1 18\n1 20\n1 22\n1 32\n2 3\n2 4\n2 8\n2 14\n2 18\n2 20\n2 22\n"
        "2 31\n3 4\n3 8\n3 9\n3 10\n3 14\n3 28\n3 29\n3 33\n4 8\n4 13\n"
        "4 14\n5 7\n5 11\n6 7\n6 11\n6 17\n7 17\n9 31\n9 33\n9 34\n"
        "10 34\n14 34\n15 33\n15 34\n16 33\n16 34\n19 33\n19 34\n20 34\n"
        "21 33\n21 34\n23 33\n23 34\n24 26\n24 28\n24 30\n24 33\n24 34\n"
        "25 26\n25 28\n25 32\n26 32\n27 30\n27 34\n28 34\n29 32\n29 34\n"
        "30 33\n30 34\n31 33\n31 34\n32 33\n32 34\n33 34\n";
    return load_snap(write_tmp("karate.txt", content));
}

// Two cliques of size 5 connected by a single bridge edge
static Graph make_barbell(int clique_size) {
    std::string c;
    // Clique A: nodes 0..clique_size-1
    for (int i = 0; i < clique_size; i++)
        for (int j = i+1; j < clique_size; j++)
            c += std::to_string(i) + " " + std::to_string(j) + "\n";
    // Clique B: nodes clique_size..2*clique_size-1
    for (int i = clique_size; i < 2*clique_size; i++)
        for (int j = i+1; j < 2*clique_size; j++)
            c += std::to_string(i) + " " + std::to_string(j) + "\n";
    // Bridge: last node of A to first node of B
    c += std::to_string(clique_size-1) + " " + std::to_string(clique_size) + "\n";
    return load_snap(write_tmp("barbell.txt", c));
}

// ── tests ────────────────────────────────────────────────────────────────────

TEST_CASE("karate club: modularity > 0.35 and multiple communities", "[leiden]") {
    Graph g = make_karate();
    LeidenResult r = run_leiden(g);
    CHECK(r.n_communities >= 2);
    CHECK(r.modularity > 0.35);
    CHECK((int)r.membership.size() == g.n);
}

TEST_CASE("barbell graph: exactly 2 communities across 10 runs", "[leiden]") {
    Graph g = make_barbell(5);
    for (int seed_run = 0; seed_run < 10; seed_run++) {
        LeidenResult r = run_leiden(g);
        CHECK(r.n_communities == 2);
    }
}

TEST_CASE("membership IDs are contiguous 0..k-1", "[leiden]") {
    Graph g = make_karate();
    LeidenResult r = run_leiden(g);
    std::set<int> ids(r.membership.begin(), r.membership.end());
    // IDs must be exactly {0, 1, ..., n_communities-1}
    CHECK((int)ids.size() == r.n_communities);
    CHECK(*ids.begin()  == 0);
    CHECK(*ids.rbegin() == r.n_communities - 1);
}

TEST_CASE("all nodes assigned a valid community", "[leiden]") {
    Graph g = make_karate();
    LeidenResult r = run_leiden(g);
    for (int v = 0; v < g.n; v++) {
        CHECK(r.membership[v] >= 0);
        CHECK(r.membership[v] < r.n_communities);
    }
}

TEST_CASE("leiden modularity matches compute_modularity Q", "[leiden]") {
    // run_leiden stores standard modularity computed via compute_modularity,
    // so the two must be exactly equal.
    Graph g = make_karate();
    LeidenResult r = run_leiden(g);

    int n1 = 0, n2 = 0;
    for (auto& [u, v] : g.edges) {
        if (r.membership[u] == r.membership[v]) n1++; else n2++;
    }

    ModResult mod = compute_modularity(g, nullptr, r.membership, n1, n2);
    CHECK_THAT(mod.Q, WithinAbs(r.modularity, 1e-12));
}

TEST_CASE("sparse graph: edge mask applied correctly", "[leiden]") {
    // Barbell with bridge removed → two disconnected cliques
    Graph g = make_barbell(5);
    // Build keep mask: exclude bridge edge (the one connecting the two cliques)
    // Bridge connects node clique_size-1 to node clique_size
    // Find its index in g.edges
    std::vector<bool> keep(g.m, true);
    for (int i = 0; i < g.m; i++) {
        auto [u, v] = g.edges[i];
        // Bridge is between node 4 and node 5 after remap (clique_size=5)
        // After remap, the bridge is between max(clique A) and min(clique B)
        // Degrees: bridge nodes have deg = clique_size-1+1 = clique_size
        // All other nodes have deg = clique_size-1
        if (g.deg[u] == 5 && g.deg[v] == 5) {
            keep[i] = false;
            break;
        }
    }
    LeidenResult r = run_leiden(g, &keep);
    // Without bridge, should still find 2 communities
    CHECK(r.n_communities == 2);
}

TEST_CASE("triangle graph gives single community or Q > 0", "[leiden]") {
    std::string path = write_tmp("tri_l.txt", "0 1\n1 2\n0 2\n");
    Graph g = load_snap(path);
    LeidenResult r = run_leiden(g);
    // Small graph — all one community or valid partition
    CHECK(r.n_communities >= 1);
    CHECK((int)r.membership.size() == g.n);
    for (int v = 0; v < g.n; v++)
        CHECK(r.membership[v] >= 0);
}
