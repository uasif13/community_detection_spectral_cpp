#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "graph.hpp"
#include "dspar.hpp"
#include "sampler.hpp"
#include "modularity.hpp"
#include "rng.hpp"
#include <fstream>
#include <filesystem>
#include <cmath>

using Catch::Matchers::WithinAbs;

// ── helpers ──────────────────────────────────────────────────────────────────

static std::string write_tmp(const std::string& name, const std::string& content) {
    std::string path = (std::filesystem::temp_directory_path() / name).string();
    std::ofstream f(path);
    f << content;
    return path;
}

static Graph make_complete(int N) {
    std::string c;
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++)
            c += std::to_string(i) + " " + std::to_string(j) + "\n";
    return load_snap(write_tmp("comp_m.txt", c));
}

// ── tests ────────────────────────────────────────────────────────────────────

TEST_CASE("two perfect cliques — no inter edges", "[modularity]") {
    // K3 ∪ K3: 6 nodes, 6 edges, 2 communities
    // community 0: nodes 0,1,2  community 1: nodes 3,4,5
    std::string path = write_tmp("two_k3.txt",
        "0 1\n0 2\n1 2\n"    // clique 0
        "3 4\n3 5\n4 5\n");  // clique 1
    Graph g = load_snap(path);
    REQUIRE(g.n == 6);
    REQUIRE(g.m == 6);

    std::vector<int> member(6);
    for (int v = 0; v < 6; v++) member[v] = (v < 3) ? 0 : 1;

    ModResult r = compute_modularity(g, nullptr, member, 6, 0);
    CHECK_THAT(r.F, WithinAbs(1.0, 1e-15));      // all intra
    // vol[0] = vol[1] = 2+2+2 = 6 (each node deg 2, 3 nodes per community)
    // G = (36+36)/(4*36) = 72/144 = 0.5
    CHECK_THAT(r.G, WithinAbs(0.5, 1e-15));
    CHECK_THAT(r.Q, WithinAbs(0.5, 1e-15));
    CHECK(r.m_sparse == 6);
    CHECK(r.preserved_intra == 6);
    CHECK(r.preserved_inter == 0);
}

TEST_CASE("all inter-community edges", "[modularity]") {
    // Complete bipartite K2_2: 4 nodes, 4 edges, all inter
    // member: 0,1 → community 0;  2,3 → community 1
    std::string path = write_tmp("k22.txt",
        "0 2\n0 3\n1 2\n1 3\n");
    Graph g = load_snap(path);
    std::vector<int> member = {0, 0, 1, 1};
    ModResult r = compute_modularity(g, nullptr, member, 0, 4);
    CHECK_THAT(r.F, WithinAbs(0.0, 1e-15));   // no intra edges
    // vol[0] = deg(0)+deg(1) = 2+2 = 4; vol[1] = 4
    // G = (16+16)/(4*16) = 32/64 = 0.5
    CHECK_THAT(r.G, WithinAbs(0.5, 1e-15));
    CHECK_THAT(r.Q, WithinAbs(-0.5, 1e-15));
}

TEST_CASE("single community — Q == 1 - 1/n^2 * something", "[modularity]") {
    // K4, all one community: F=1, G = (4*3/2 * 2 / n)^2... compute directly
    // Each node deg=3, vol = 4*3=12, m=6
    // G = 12^2 / (4*36) = 144/144 = 1.0
    // Q = 1 - 1 = 0
    Graph g = make_complete(4);
    std::vector<int> member(4, 0);
    ModResult r = compute_modularity(g, nullptr, member, 6, 0);
    CHECK_THAT(r.F, WithinAbs(1.0, 1e-15));
    CHECK_THAT(r.G, WithinAbs(1.0, 1e-15));
    CHECK_THAT(r.Q, WithinAbs(0.0, 1e-15));
}

TEST_CASE("keep=nullptr matches keep=all-true", "[modularity]") {
    std::string path = write_tmp("tri_m.txt", "0 1\n1 2\n0 2\n");
    Graph g = load_snap(path);
    std::vector<int> member = {0, 0, 1};
    std::vector<bool> all_true(g.m, true);
    ModResult r1 = compute_modularity(g, nullptr,    member, 2, 1);
    ModResult r2 = compute_modularity(g, &all_true,  member, 2, 1);
    CHECK_THAT(r1.F, WithinAbs(r2.F, 1e-15));
    CHECK_THAT(r1.G, WithinAbs(r2.G, 1e-15));
    CHECK_THAT(r1.Q, WithinAbs(r2.Q, 1e-15));
    CHECK(r1.m_sparse == r2.m_sparse);
}

TEST_CASE("sparse graph uses sparse degrees for G", "[modularity]") {
    // Triangle 0-1-2 with member={0,0,1}
    // Remove edge (1,2) — inter edge
    // Remaining: (0,1) intra, (0,2) inter — m_sparse=2
    // Sparse deg: 0→2, 1→1, 2→1
    // vol[0] = 2+1 = 3, vol[1] = 1
    // G = (9+1)/(4*4) = 10/16 = 0.625
    // F = 1/2 = 0.5
    std::string path = write_tmp("tri2.txt", "0 1\n1 2\n0 2\n");
    Graph g = load_snap(path);
    // Identify edges: sort gives (0,1),(0,2),(1,2)
    // member: 0→0, 1→0, 2→1
    std::vector<int> member = {0, 0, 1};
    // Keep only (0,1) and (0,2): indices 0 and 1
    std::vector<bool> keep = {true, true, false};
    ModResult r = compute_modularity(g, &keep, member, 2, 1);
    CHECK(r.m_sparse == 2);
    CHECK_THAT(r.F, WithinAbs(0.5, 1e-15));
    CHECK_THAT(r.G, WithinAbs(10.0/16.0, 1e-12));
    CHECK_THAT(r.Q, WithinAbs(0.5 - 10.0/16.0, 1e-12));
}

TEST_CASE("identity ΔQ = ΔF - ΔG holds for 100 random masks", "[modularity]") {
    // Karate-like: use K10 with two communities
    Graph g = make_complete(10);
    // Communities: 0..4 → 0, 5..9 → 1
    std::vector<int> member(10);
    for (int v = 0; v < 10; v++) member[v] = (v < 5) ? 0 : 1;
    // n1 = C(5,2)*2 = 20, n2 = 5*5 = 25
    ModResult orig = compute_modularity(g, nullptr, member, 20, 25);

    Xoshiro256pp rng(2024);
    DsparScores scores = compute_scores(g);

    for (int t = 0; t < 100; t++) {
        auto keep = sample_edges(scores, 0.6, rng);
        ModResult sp = compute_modularity(g, &keep, member, 20, 25);

        double dQ = sp.Q - orig.Q;
        double dF = sp.F - orig.F;
        double dG = sp.G - orig.G;
        double err = std::abs(dQ - (dF - dG));
        CHECK(err < 1e-12);
    }
}

TEST_CASE("empty sparse graph (all edges removed)", "[modularity]") {
    std::string path = write_tmp("tri3.txt", "0 1\n1 2\n0 2\n");
    Graph g = load_snap(path);
    std::vector<int> member = {0, 1, 1};
    std::vector<bool> keep(g.m, false);
    ModResult r = compute_modularity(g, &keep, member, 1, 2);
    CHECK(r.m_sparse == 0);
    CHECK_THAT(r.F, WithinAbs(0.0, 1e-15));
    CHECK_THAT(r.G, WithinAbs(0.0, 1e-15));
    CHECK_THAT(r.Q, WithinAbs(0.0, 1e-15));
}

TEST_CASE("preserved counts and rates", "[modularity]") {
    // Triangle: (0,1) intra, (0,2) inter, (1,2) inter  (member={0,0,1})
    std::string path = write_tmp("tri4.txt", "0 1\n1 2\n0 2\n");
    Graph g = load_snap(path);
    std::vector<int> member = {0, 0, 1};
    // keep: (0,1)=true intra, (0,2)=true inter, (1,2)=false inter
    std::vector<bool> keep = {true, true, false};
    ModResult r = compute_modularity(g, &keep, member, 1, 2);
    CHECK(r.preserved_intra == 1);
    CHECK(r.preserved_inter == 1);
    CHECK_THAT(r.intra_rate, WithinAbs(1.0, 1e-15));   // 1/1
    CHECK_THAT(r.inter_rate, WithinAbs(0.5, 1e-15));   // 1/2
    CHECK_THAT(r.ratio_obs,  WithinAbs(0.5, 1e-15));   // 0.5/1.0
}
