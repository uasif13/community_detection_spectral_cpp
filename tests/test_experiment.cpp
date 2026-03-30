#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "graph.hpp"
#include "experiment.hpp"
#include <fstream>
#include <filesystem>
#include <sstream>
#include <algorithm>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

using Catch::Matchers::WithinAbs;

// ── helpers ──────────────────────────────────────────────────────────────────

static std::string write_tmp(const std::string& name, const std::string& content) {
    std::string path = (std::filesystem::temp_directory_path() / name).string();
    std::ofstream f(path);
    f << content;
    return path;
}

// Zachary's karate club
static Graph make_karate() {
    std::string content =
        "1 2\n1 3\n1 4\n1 5\n1 6\n1 7\n1 8\n1 9\n1 11\n1 12\n1 13\n1 14\n"
        "1 18\n1 20\n1 22\n1 32\n2 3\n2 4\n2 8\n2 14\n2 18\n2 20\n2 22\n"
        "2 31\n3 4\n3 8\n3 9\n3 10\n3 14\n3 28\n3 29\n3 33\n4 8\n4 13\n"
        "4 14\n5 7\n5 11\n6 7\n6 11\n6 17\n7 17\n9 31\n9 33\n9 34\n"
        "10 34\n14 34\n15 33\n15 34\n16 33\n16 34\n19 33\n19 34\n20 34\n"
        "21 33\n21 34\n23 33\n23 34\n24 26\n24 28\n24 30\n24 33\n24 34\n"
        "25 26\n25 28\n25 32\n26 32\n27 30\n27 34\n28 34\n29 32\n29 34\n"
        "30 33\n30 34\n31 33\n31 34\n32 33\n32 34\n33 34\n";
    return load_snap(write_tmp("karate_exp.txt", content));
}

static std::vector<std::string> csv_header(const std::string& path) {
    std::ifstream f(path);
    std::string line;
    std::getline(f, line);
    std::vector<std::string> cols;
    std::stringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, ',')) cols.push_back(tok);
    return cols;
}

// ── tests ────────────────────────────────────────────────────────────────────

TEST_CASE("identity holds for all karate trials", "[experiment]") {
    Graph g = make_karate();
    ExperimentConfig cfg;
    cfg.alphas   = {0.5, 0.7, 0.9};
    cfg.n_trials = 5;
    cfg.base_seed = 42;
    cfg.run_leiden_resample = false;  // faster

    auto results = run_experiment(g, "karate", cfg);
    REQUIRE((int)results.size() == 3 * 5);
    for (auto& r : results)
        CHECK(r.identity_error < 1e-12);
}

TEST_CASE("result count matches alphas * n_trials", "[experiment]") {
    Graph g = make_karate();
    ExperimentConfig cfg;
    cfg.alphas   = {0.5, 0.6, 0.7, 0.8, 0.9};
    cfg.n_trials = 4;
    cfg.run_leiden_resample = false;

    auto results = run_experiment(g, "karate", cfg);
    CHECK((int)results.size() == 5 * 4);
}

TEST_CASE("alpha=1.0 keeps all edges on regular graph", "[experiment]") {
    // K6: every node has degree 5, so all edge scores equal mean → p(e)=1 for all e
    std::string content;
    for (int i = 0; i < 6; i++)
        for (int j = i+1; j < 6; j++)
            content += std::to_string(i) + " " + std::to_string(j) + "\n";
    Graph g = load_snap(write_tmp("k6.txt", content));

    ExperimentConfig cfg;
    cfg.alphas   = {1.0};
    cfg.n_trials = 3;
    cfg.run_leiden_resample = false;

    auto results = run_experiment(g, "k6", cfg);
    for (auto& r : results) {
        CHECK(r.m_sparse == g.m);
        CHECK_THAT(r.modularity_fixed_change, WithinAbs(0.0, 1e-12));
        CHECK_THAT(r.dQ_reconstructed,        WithinAbs(0.0, 1e-12));
    }
}

TEST_CASE("seed reproducibility", "[experiment]") {
    Graph g = make_karate();
    ExperimentConfig cfg;
    cfg.alphas   = {0.6, 0.8};
    cfg.n_trials = 3;
    cfg.run_leiden_resample = false;

    auto r1 = run_experiment(g, "karate", cfg);
    auto r2 = run_experiment(g, "karate", cfg);
    REQUIRE(r1.size() == r2.size());
    for (int i = 0; i < (int)r1.size(); i++) {
        CHECK(r1[i].seed    == r2[i].seed);
        CHECK(r1[i].m_sparse == r2[i].m_sparse);
        CHECK_THAT(r1[i].modularity_fixed_change,
                   WithinAbs(r2[i].modularity_fixed_change, 1e-15));
    }
}

TEST_CASE("raw CSV columns match spec exactly", "[experiment]") {
    Graph g = make_karate();
    ExperimentConfig cfg;
    cfg.alphas   = {0.5};
    cfg.n_trials = 2;
    cfg.run_leiden_resample = false;

    auto results = run_experiment(g, "karate", cfg);
    std::string path = (std::filesystem::temp_directory_path() / "karate_raw.csv").string();
    write_raw_csv(results, path);

    auto cols = csv_header(path);
    std::vector<std::string> expected = {
        "retention", "seed", "m", "n1", "n2",
        "mu_intra", "mu_inter", "delta",
        "F_original", "G_original", "modularity_fixed_original",
        "m_sparse", "F_observed", "G_sparse_observed",
        "dG_observed", "F_improvement_observed",
        "modularity_fixed_sparse", "modularity_fixed_change",
        "modularity_leiden_sparse", "modularity_leiden_change",
        "dQ_reconstructed",
        "preserved_intra", "preserved_inter",
        "intra_rate", "inter_rate", "ratio_observed"
    };
    CHECK(cols == expected);
}

TEST_CASE("summary CSV has correct columns", "[experiment]") {
    Graph g = make_karate();
    ExperimentConfig cfg;
    cfg.alphas   = {0.5, 0.7};
    cfg.n_trials = 3;
    cfg.run_leiden_resample = false;

    auto results = run_experiment(g, "karate", cfg);
    std::string path = (std::filesystem::temp_directory_path() / "karate_summary.csv").string();
    write_summary_csv(results, path);

    auto cols = csv_header(path);
    std::vector<std::string> expected = {
        "retention",
        "modularity_fixed_change_mean", "modularity_fixed_change_std",
        "modularity_leiden_change_mean", "modularity_leiden_change_std",
        "F_improvement_observed_mean", "F_improvement_observed_std",
        "dG_observed_mean", "dG_observed_std",
        "delta_mean"
    };
    CHECK(cols == expected);

    // One summary row per alpha
    std::ifstream f(path);
    int line_count = 0;
    std::string line;
    while (std::getline(f, line)) line_count++;
    CHECK(line_count == 3);  // header + 2 alpha rows
}

TEST_CASE("dQ_leiden >= dQ_fixed majority of karate trials", "[experiment]") {
    Graph g = make_karate();
    ExperimentConfig cfg;
    cfg.alphas   = {0.5, 0.6, 0.7, 0.8, 0.9};
    cfg.n_trials = 10;
    cfg.run_leiden_resample = true;

    auto results = run_experiment(g, "karate", cfg);
    int pass = 0;
    for (auto& r : results)
        if (r.modularity_leiden_change >= r.modularity_fixed_change) pass++;
    double rate = static_cast<double>(pass) / results.size();
    CHECK(rate >= 0.8);
}

TEST_CASE("preserved_intra + preserved_inter == m_sparse", "[experiment]") {
    Graph g = make_karate();
    ExperimentConfig cfg;
    cfg.alphas   = {0.6};
    cfg.n_trials = 5;
    cfg.run_leiden_resample = false;

    auto results = run_experiment(g, "karate", cfg);
    for (auto& r : results)
        CHECK(r.preserved_intra + r.preserved_inter == r.m_sparse);
}

TEST_CASE("parallel matches serial (OpenMP)", "[experiment][parallel]") {
    // Run with 1 thread then with multiple threads — results must be bit-identical
    // because the seed scheme is deterministic per (ai, t).
    Graph g = make_karate();
    ExperimentConfig cfg;
    cfg.alphas            = {0.5, 0.6, 0.7, 0.8, 0.9};
    cfg.n_trials          = 6;
    cfg.base_seed         = 123;
    cfg.run_leiden_resample = false;  // deterministic; Leiden is randomised

#ifdef _OPENMP
    omp_set_num_threads(1);
    auto serial = run_experiment(g, "karate", cfg);

    omp_set_num_threads(4);
    auto parallel = run_experiment(g, "karate", cfg);

    REQUIRE(serial.size() == parallel.size());
    for (int i = 0; i < (int)serial.size(); i++) {
        CHECK(serial[i].seed     == parallel[i].seed);
        CHECK(serial[i].m_sparse == parallel[i].m_sparse);
        CHECK_THAT(serial[i].F_observed,
                   WithinAbs(parallel[i].F_observed, 1e-15));
        CHECK_THAT(serial[i].modularity_fixed_change,
                   WithinAbs(parallel[i].modularity_fixed_change, 1e-15));
        CHECK_THAT(serial[i].dQ_reconstructed,
                   WithinAbs(parallel[i].dQ_reconstructed, 1e-15));
    }
    // restore default
    omp_set_num_threads(omp_get_max_threads());
#else
    // Without OpenMP the loop is single-threaded; test still validates seed scheme
    auto r1 = run_experiment(g, "karate", cfg);
    auto r2 = run_experiment(g, "karate", cfg);
    REQUIRE(r1.size() == r2.size());
    for (int i = 0; i < (int)r1.size(); i++)
        CHECK_THAT(r1[i].modularity_fixed_change,
                   WithinAbs(r2[i].modularity_fixed_change, 1e-15));
#endif
}
