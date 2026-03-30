#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "graph.hpp"
#include "dspar.hpp"
#include "sampler.hpp"
#include "rng.hpp"
#include <fstream>
#include <filesystem>
#include <numeric>
#include <cmath>
#include <algorithm>

using Catch::Matchers::WithinAbs;

// ── helpers ──────────────────────────────────────────────────────────────────

static std::string write_tmp(const std::string& name, const std::string& content) {
    std::string path = (std::filesystem::temp_directory_path() / name).string();
    std::ofstream f(path);
    f << content;
    return path;
}

static Graph make_complete(int N) {
    std::string content;
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++)
            content += std::to_string(i) + " " + std::to_string(j) + "\n";
    return load_snap(write_tmp("complete_s.txt", content));
}

static Graph make_star(int N) {
    std::string content;
    for (int i = 1; i <= N; i++)
        content += "0 " + std::to_string(i) + "\n";
    return load_snap(write_tmp("star_s.txt", content));
}

// ── tests ────────────────────────────────────────────────────────────────────

TEST_CASE("alpha=1.0 keeps all edges on uniform graph", "[sampler]") {
    // K5: all scores equal → all p(e) = 1.0 → all kept
    Graph g = make_complete(5);
    DsparScores s = compute_scores(g);
    Xoshiro256pp rng(12345);
    auto keep = sample_edges(s, 1.0, rng);
    REQUIRE((int)keep.size() == g.m);
    for (int i = 0; i < g.m; i++)
        CHECK(keep[i] == true);
}

TEST_CASE("alpha=0.0 keeps no edges", "[sampler]") {
    Graph g = make_complete(5);
    DsparScores s = compute_scores(g);
    Xoshiro256pp rng(99);
    auto keep = sample_edges(s, 0.0, rng);
    int kept = 0;
    for (bool k : keep) if (k) kept++;
    CHECK(kept == 0);
}

TEST_CASE("expected retention rate matches alpha", "[sampler]") {
    // 100-edge path graph → run many trials, mean kept/100 ≈ alpha
    // Use K15 (105 edges, close to 100) — all scores equal → exact p = alpha
    Graph g = make_complete(15);
    DsparScores s = compute_scores(g);
    const double alpha = 0.6;
    const int trials = 5000;
    double total = 0.0;
    Xoshiro256pp rng(777);
    for (int t = 0; t < trials; t++) {
        auto keep = sample_edges(s, alpha, rng);
        int cnt = 0;
        for (bool k : keep) if (k) cnt++;
        total += cnt;
    }
    double mean_rate = total / (trials * g.m);
    // Within 3σ: σ = sqrt(alpha*(1-alpha)/m) ≈ sqrt(0.24/105) ≈ 0.048 per trial
    // Over 5000 trials: 3σ/sqrt(trials) ≈ 3*0.048/sqrt(5000) ≈ 0.002
    CHECK_THAT(mean_rate, WithinAbs(alpha, 0.01));
}

TEST_CASE("seed reproducibility", "[sampler]") {
    Graph g = make_star(20);
    DsparScores s = compute_scores(g);
    Xoshiro256pp rng1(42);
    Xoshiro256pp rng2(42);
    auto keep1 = sample_edges(s, 0.5, rng1);
    auto keep2 = sample_edges(s, 0.5, rng2);
    REQUIRE(keep1.size() == keep2.size());
    for (int i = 0; i < (int)keep1.size(); i++)
        CHECK(keep1[i] == keep2[i]);
}

TEST_CASE("different seeds give different results", "[sampler]") {
    Graph g = make_complete(15);
    DsparScores s = compute_scores(g);
    Xoshiro256pp rng1(1);
    Xoshiro256pp rng2(2);
    auto keep1 = sample_edges(s, 0.5, rng1);
    auto keep2 = sample_edges(s, 0.5, rng2);
    bool any_diff = false;
    for (int i = 0; i < (int)keep1.size(); i++)
        if (keep1[i] != keep2[i]) { any_diff = true; break; }
    CHECK(any_diff);
}

TEST_CASE("inplace variant matches allocating variant", "[sampler]") {
    Graph g = make_star(15);
    DsparScores s = compute_scores(g);
    Xoshiro256pp rng1(55), rng2(55);
    auto keep1 = sample_edges(s, 0.7, rng1);
    std::vector<bool> keep2;
    sample_edges_inplace(s, 0.7, rng2, keep2);
    REQUIRE(keep1.size() == keep2.size());
    for (int i = 0; i < (int)keep1.size(); i++)
        CHECK(keep1[i] == keep2[i]);
}

TEST_CASE("star mean retention rate matches alpha", "[sampler]") {
    // Star K1_10: all edges have same score (1/10+1/1=1.1) and same mean,
    // so p(e) = alpha for all edges
    Graph g = make_star(10);
    DsparScores s = compute_scores(g);
    const double alpha = 0.5;
    const int trials = 5000;
    double total = 0.0;
    Xoshiro256pp rng(314159);
    for (int t = 0; t < trials; t++) {
        auto keep = sample_edges(s, alpha, rng);
        int cnt = 0;
        for (bool k : keep) if (k) cnt++;
        total += cnt;
    }
    double mean_rate = total / (trials * g.m);
    CHECK_THAT(mean_rate, WithinAbs(alpha, 0.02));
}

TEST_CASE("p(e) <= 1 clamping — edge with score >> mean", "[sampler]") {
    // Star K1_3: hub has deg=3, leaves deg=1
    // Hub-leaf score = 1/3 + 1 = 1.333; mean = 1.333 (only 3 such edges)
    // At alpha=1.0: p = 1.0 exactly → no assertion fires
    // At alpha=0.9 with score/mean > 1 scenario: test at alpha slightly above
    // where some edges would exceed 1 if not clamped.
    // Use a 2-node path (m=1, score=1+1=2, mean=2) at alpha=1.0 → p=1 exactly
    std::string path = write_tmp("two.txt", "0 1\n");
    Graph g = load_snap(path);
    DsparScores s = compute_scores(g);
    // score[0]=2, mean=2, so p=alpha*1=alpha; no clamping needed at alpha<=1
    Xoshiro256pp rng(0);
    // alpha=1.0 → keep everything, no assertion
    auto keep = sample_edges(s, 1.0, rng);
    CHECK(keep[0] == true);

    // Construct a graph where one edge has score > mean: star K1_2
    // Hub deg=2, leaves deg=1; edge scores both = 1/2+1/1 = 1.5; mean=1.5
    // p = alpha * 1.5 / 1.5 = alpha → no clamping needed
    // But a heterogeneous case: path 0-1-2 (m=2)
    // deg[1]=2, deg[0]=deg[2]=1
    // edge(0,1) score = 1+0.5=1.5, edge(1,2) score=0.5+1=1.5, mean=1.5
    // p = alpha * 1.5/1.5 = alpha → still no > 1 case
    // True > 1 case requires a node with very low degree relative to mean
    // Use triangle + pendant: edges 0-1, 1-2, 0-2, 0-3
    // deg: 0→3, 1→2, 2→2, 3→1
    // scores: (0,1)=1/3+1/2=5/6, (0,2)=5/6, (1,2)=1/2+1/2=1, (0,3)=1/3+1=4/3
    // mean = (5/6+5/6+1+4/3)/4 = (5/6+5/6+6/6+8/6)/4 = 24/6/4 = 4/4=1
    // p(0,3) at alpha=0.9: 0.9*4/3/1 = 1.2 > 1 → must be clamped
    std::string path2 = write_tmp("clamp.txt", "0 1\n1 2\n0 2\n0 3\n");
    Graph g2 = load_snap(path2);
    DsparScores s2 = compute_scores(g2);
    Xoshiro256pp rng2(7);
    // alpha=0.9 → p for pendant edge = 0.9*1.333/1.0 = 1.2 → clamped → no assert
    auto keep2 = sample_edges(s2, 0.9, rng2);
    CHECK((int)keep2.size() == g2.m);
    // The pendant edge (highest score) should be always retained at alpha=0.9
    // since p was clamped to 1.0
    // Run many times to confirm pendant always kept
    int always_kept = 0;
    // Find pendant edge index (connects node with deg=1)
    int pendant_idx = -1;
    for (int i = 0; i < g2.m; i++) {
        auto [u, v] = g2.edges[i];
        if (g2.deg[u] == 1 || g2.deg[v] == 1) { pendant_idx = i; break; }
    }
    REQUIRE(pendant_idx != -1);
    for (int t = 0; t < 200; t++) {
        Xoshiro256pp r(t);
        auto k = sample_edges(s2, 0.9, r);
        if (k[pendant_idx]) always_kept++;
    }
    CHECK(always_kept == 200);  // always kept because p was clamped to 1.0
}

// ── Opt B: pre-computed thresholds ───────────────────────────────────────────

TEST_CASE("bernoulli threshold matches inplace for same RNG state", "[sampler][optB]") {
    // sample_edges_bernoulli with pre-computed thresholds must produce
    // the same result as sample_edges_inplace for the same RNG seed.
    Graph g = make_complete(10);
    DsparScores s = compute_scores(g);
    const double alpha = 0.6;

    // Build thresholds the same way run_experiment does
    std::vector<double> threshold(g.m);
    for (int i = 0; i < g.m; i++)
        threshold[i] = std::min(1.0, alpha * s.score[i] / s.mean);

    Xoshiro256pp rng1(42), rng2(42);
    std::vector<bool> keep1, keep2;
    sample_edges_inplace(s, alpha, rng1, keep1);
    sample_edges_bernoulli(threshold, rng2, keep2);

    REQUIRE(keep1.size() == keep2.size());
    for (int i = 0; i < (int)keep1.size(); i++)
        CHECK(keep1[i] == keep2[i]);
}

TEST_CASE("bernoulli threshold mean retention matches alpha", "[sampler][optB]") {
    Graph g = make_complete(15);  // uniform scores → p(e) = alpha for all e
    DsparScores s = compute_scores(g);
    const double alpha = 0.7;

    std::vector<double> threshold(g.m);
    for (int i = 0; i < g.m; i++)
        threshold[i] = std::min(1.0, alpha * s.score[i] / s.mean);

    const int trials = 3000;
    double total = 0.0;
    Xoshiro256pp rng(99);
    std::vector<bool> keep;
    for (int t = 0; t < trials; t++) {
        sample_edges_bernoulli(threshold, rng, keep);
        for (bool k : keep) if (k) total++;
    }
    double mean_rate = total / (trials * g.m);
    CHECK_THAT(mean_rate, WithinAbs(alpha, 0.02));
}

// ── Opt D: prefix-sum sampler ─────────────────────────────────────────────────

TEST_CASE("prefix-sum never exceeds k unique edges", "[sampler][optD]") {
    // Duplicate draws are collapsed into a single keep; result is always <= k.
    Graph g = make_complete(20);  // 190 edges
    DsparScores s = compute_scores(g);

    std::vector<double> cdf(g.m + 1, 0.0);
    for (int i = 0; i < g.m; i++)
        cdf[i + 1] = cdf[i] + s.score[i];
    double total = cdf[g.m];

    const int k = 50;
    Xoshiro256pp rng(7);
    std::vector<bool> keep;
    for (int t = 0; t < 20; t++) {
        sample_edges_prefix_sum(cdf, total, k, rng, keep);
        int kept = std::count(keep.begin(), keep.end(), true);
        CHECK(kept <= k);
        CHECK(kept > 0);
    }
}

TEST_CASE("prefix-sum k=0 keeps nothing", "[sampler][optD]") {
    Graph g = make_complete(10);
    DsparScores s = compute_scores(g);
    std::vector<double> cdf(g.m + 1, 0.0);
    for (int i = 0; i < g.m; i++) cdf[i+1] = cdf[i] + s.score[i];
    Xoshiro256pp rng(1);
    std::vector<bool> keep;
    sample_edges_prefix_sum(cdf, cdf[g.m], 0, rng, keep);
    int kept = std::count(keep.begin(), keep.end(), true);
    CHECK(kept == 0);
}

TEST_CASE("prefix-sum mean unique edges matches birthday formula", "[sampler][optD]") {
    // When scores are uniform the expected number of unique edges retained is
    //   E[U] = m * (1 - (1 - 1/m)^k)  ≈  m * (1 - exp(-k/m))
    // Build a path of 200 nodes (199 edges) — non-uniform scores exercise CDF code.
    std::string content;
    for (int i = 0; i < 199; i++)
        content += std::to_string(i) + " " + std::to_string(i+1) + "\n";
    Graph g = load_snap(write_tmp("path200.txt", content));
    DsparScores s = compute_scores(g);

    std::vector<double> cdf(g.m + 1, 0.0);
    for (int i = 0; i < g.m; i++) cdf[i+1] = cdf[i] + s.score[i];
    double total = cdf[g.m];

    const double alpha = 0.15;
    const int k = static_cast<int>(std::round(alpha * g.m));
    // Expected unique = m * (1 - exp(-k/m))
    double expected_unique = g.m * (1.0 - std::exp(-static_cast<double>(k) / g.m));

    const int trials = 2000;
    double total_kept = 0.0;
    Xoshiro256pp rng(42);
    std::vector<bool> keep;
    for (int t = 0; t < trials; t++) {
        sample_edges_prefix_sum(cdf, total, k, rng, keep);
        total_kept += std::count(keep.begin(), keep.end(), true);
    }
    double mean_kept = total_kept / trials;
    // Tolerate ±10% of expected_unique (sampling variance + non-uniform scores)
    CHECK_THAT(mean_kept, WithinAbs(expected_unique, expected_unique * 0.10 + 2.0));
}
