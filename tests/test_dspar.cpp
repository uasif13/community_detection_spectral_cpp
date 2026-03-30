#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "graph.hpp"
#include "dspar.hpp"
#include <fstream>
#include <filesystem>
#include <numeric>
#include <cstdlib>

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

// ── helpers ──────────────────────────────────────────────────────────────────

static std::string write_tmp(const std::string& name, const std::string& content) {
    std::string path = (std::filesystem::temp_directory_path() / name).string();
    std::ofstream f(path);
    f << content;
    return path;
}

// Build a star graph K1_N: hub=0, leaves=1..N
static Graph make_star(int N) {
    std::string content = "# star\n";
    for (int i = 1; i <= N; i++)
        content += "0 " + std::to_string(i) + "\n";
    return load_snap(write_tmp("star.txt", content));
}

// Build complete graph K_N
static Graph make_complete(int N) {
    std::string content = "# K" + std::to_string(N) + "\n";
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++)
            content += std::to_string(i) + " " + std::to_string(j) + "\n";
    return load_snap(write_tmp("complete.txt", content));
}

// ── tests ────────────────────────────────────────────────────────────────────

TEST_CASE("star graph edge score", "[dspar]") {
    // K1_10: hub deg=10, each leaf deg=1
    // Score of hub-leaf edge = 1/10 + 1/1 = 1.1 exactly
    Graph g = make_star(10);
    DsparScores s = compute_scores(g);
    REQUIRE(s.score.size() == 10u);
    for (double sc : s.score)
        CHECK_THAT(sc, WithinAbs(1.1, 1e-15));
}

TEST_CASE("complete graph K5 uniform scores", "[dspar]") {
    // All nodes deg=4 → all edge scores = 1/4 + 1/4 = 0.5
    Graph g = make_complete(5);
    DsparScores s = compute_scores(g);
    REQUIRE((int)s.score.size() == g.m);
    for (double sc : s.score)
        CHECK_THAT(sc, WithinAbs(0.5, 1e-15));
    CHECK_THAT(s.mean, WithinAbs(0.5, 1e-15));
}

TEST_CASE("mean invariant", "[dspar]") {
    // sum(scores) / m == mean, to 1e-14
    Graph g = make_star(20);
    DsparScores s = compute_scores(g);
    double direct_mean = 0.0;
    for (double sc : s.score) direct_mean += sc;
    direct_mean /= g.m;
    CHECK_THAT(s.mean, WithinAbs(direct_mean, 1e-14));
}

TEST_CASE("K5 mean invariant", "[dspar]") {
    Graph g = make_complete(5);
    DsparScores s = compute_scores(g);
    double direct_mean = 0.0;
    for (double sc : s.score) direct_mean += sc;
    direct_mean /= g.m;
    CHECK_THAT(s.mean, WithinAbs(direct_mean, 1e-14));
}

TEST_CASE("two-community structural metrics", "[dspar]") {
    // Two triangles: nodes 0-1-2 (community 0) and 3-4-5 (community 1),
    // plus one inter-community edge 2-3
    std::string path = write_tmp("two_tri.txt",
        "0 1\n1 2\n0 2\n"   // triangle 1
        "3 4\n4 5\n3 5\n"   // triangle 2
        "2 3\n"             // bridge
    );
    Graph g = load_snap(path);
    DsparScores s = compute_scores(g);

    // membership: nodes 0,1,2 → 0; nodes 3,4,5 → 1
    // After remap the IDs may be different — find them by degree pattern
    // All nodes have deg 2 except bridge endpoints (deg 3)
    std::vector<int> membership(g.n);
    // Identify the two bridge nodes (deg 3) and their communities
    // Actually all nodes in each triangle have deg 2 except bridge endpoints deg 3
    // Build membership based on which component each node is closer to:
    // We rely on the graph being loaded deterministically (node IDs 0..5)
    // and the remap being stable (node 0 maps to 0, etc.)
    for (int v = 0; v < g.n; v++)
        membership[v] = (v < 3) ? 0 : 1;  // assumes stable remap of 0..5

    compute_structural_metrics(s, g, membership);
    CHECK(s.n1 + s.n2 == g.m);
    // 6 intra edges, 1 inter edge
    CHECK(s.n1 == 6);
    CHECK(s.n2 == 1);
    CHECK(s.delta == s.mu_intra - s.mu_inter);
}

TEST_CASE("single community gives mu_inter = 0", "[dspar]") {
    // All edges intra → n2=0, mu_inter=0, delta=mu_intra
    Graph g = make_complete(4);
    DsparScores s = compute_scores(g);
    std::vector<int> membership(g.n, 0);  // everyone in community 0
    compute_structural_metrics(s, g, membership);
    CHECK(s.n2 == 0);
    CHECK_THAT(s.mu_inter, WithinAbs(0.0, 1e-15));
    CHECK_THAT(s.delta, WithinAbs(s.mu_intra, 1e-15));
    CHECK(s.n1 == g.m);
}

TEST_CASE("SNAP ca-GrQc reference mean", "[dspar][snap]") {
    const char* dp = std::getenv("SNAP_DATA_DIR");
    std::string path = dp ? std::string(dp) + "/ca-GrQc.txt" : "data/ca-GrQc.txt";
    if (!std::filesystem::exists(path))
        SKIP("ca-GrQc.txt not found");

    Graph g = load_snap(path);
    DsparScores s = compute_scores(g);
    // Known mean score from Python reference (~0.07845)
    CHECK_THAT(s.mean, WithinAbs(0.07845, 0.001));
    // Mean must equal sum/m to within 1e-10
    double direct = 0.0;
    for (double sc : s.score) direct += sc;
    direct /= g.m;
    CHECK_THAT(s.mean, WithinAbs(direct, 1e-10));
}
