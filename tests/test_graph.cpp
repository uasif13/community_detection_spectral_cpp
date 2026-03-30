#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "graph.hpp"
#include <fstream>
#include <numeric>
#include <filesystem>
#include <cstdlib>

static std::string write_tmp(const std::string& content) {
    std::string path = std::filesystem::temp_directory_path() / "test_graph_tmp.txt";
    std::ofstream f(path);
    f << content;
    return path;
}

TEST_CASE("triangle graph", "[graph]") {
    std::string path = write_tmp(
        "# comment\n"
        "0 1\n"
        "1 2\n"
        "0 2\n"
    );
    Graph g = load_snap(path);
    CHECK(g.n == 3);
    CHECK(g.m == 3);
    CHECK(g.deg[0] == 2);
    CHECK(g.deg[1] == 2);
    CHECK(g.deg[2] == 2);
    // adj is symmetric: each node appears in others' lists
    CHECK(g.adj.size() == 6u);
    // edges: all u < v
    for (auto& [u, v] : g.edges)
        CHECK(u < v);
    CHECK(g.edges.size() == 3u);
}

TEST_CASE("star K1_10", "[graph]") {
    std::string content = "# star\n";
    for (int i = 1; i <= 10; i++)
        content += "0 " + std::to_string(i) + "\n";
    std::string path = write_tmp(content);
    Graph g = load_snap(path);
    CHECK(g.n == 11);
    CHECK(g.m == 10);
    // Hub has degree 10 — but hub may be remapped; find which node has deg 10
    int hub = -1;
    for (int v = 0; v < g.n; v++)
        if (g.deg[v] == 10) { hub = v; break; }
    CHECK(hub != -1);
    for (int v = 0; v < g.n; v++)
        if (v != hub) CHECK(g.deg[v] == 1);
    // All edges canonical
    for (auto& [u, v] : g.edges)
        CHECK(u < v);
    CHECK(g.edges.size() == 10u);
}

TEST_CASE("self-loop rejection", "[graph]") {
    std::string path = write_tmp(
        "0 1\n"
        "1 2\n"
        "3 3\n"  // self-loop — must be discarded
        "0 2\n"
    );
    Graph g = load_snap(path);
    // node 3 only appears in a self-loop → excluded (no valid edges)
    // nodes from valid edges: 0,1,2 → n=3
    CHECK(g.n == 3);
    CHECK(g.m == 3);
    for (auto& [u, v] : g.edges)
        CHECK(u != v);
}

TEST_CASE("degree sum invariant", "[graph]") {
    std::string path = write_tmp(
        "10 20\n"
        "20 30\n"
        "10 30\n"
        "10 40\n"
    );
    Graph g = load_snap(path);
    int deg_sum = std::accumulate(g.deg.begin(), g.deg.end(), 0);
    CHECK(deg_sum == 2 * g.m);
}

TEST_CASE("duplicate edge deduplication", "[graph]") {
    // Same edge listed twice → should appear once
    std::string path = write_tmp(
        "0 1\n"
        "1 0\n"  // duplicate (reversed)
        "0 2\n"
    );
    Graph g = load_snap(path);
    CHECK(g.m == 2);
}

TEST_CASE("CSR adj_off consistency", "[graph]") {
    std::string path = write_tmp(
        "0 1\n"
        "1 2\n"
        "2 3\n"
        "3 0\n"
    );
    Graph g = load_snap(path);
    // adj_off[v+1] - adj_off[v] == deg[v] for all v
    for (int v = 0; v < g.n; v++)
        CHECK(g.adj_off[v+1] - g.adj_off[v] == g.deg[v]);
    CHECK(g.adj_off[g.n] == 2 * g.m);
}

TEST_CASE("SNAP file ca-GrQc", "[graph][snap]") {
    const char* data_path = std::getenv("SNAP_DATA_DIR");
    std::string path = data_path
        ? std::string(data_path) + "/ca-GrQc.txt"
        : "data/ca-GrQc.txt";
    if (!std::filesystem::exists(path)) {
        SKIP("ca-GrQc.txt not found — set SNAP_DATA_DIR or place in data/");
    }
    Graph g = load_snap(path);
    CHECK(g.n == 5242);
    CHECK(g.m == 14496);
    int deg_sum = std::accumulate(g.deg.begin(), g.deg.end(), 0);
    CHECK(deg_sum == 2 * g.m);
}
