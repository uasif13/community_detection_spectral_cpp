#pragma once
#include "graph.hpp"
#include <cstdint>
#include <string>
#include <vector>

struct ExperimentConfig {
    std::vector<double> alphas      = {0.5, 0.6, 0.7, 0.8, 0.9};
    int                 n_trials    = 10;
    uint64_t            base_seed   = 42;
    bool                run_leiden_resample = true;
};

struct TrialResult {
    // Config
    double   retention;
    int      trial;
    uint64_t seed;
    // Graph stats
    int m, n1, n2;
    // DSpar structural metrics (constant per dataset)
    double mu_intra, mu_inter, delta;
    // Original graph modularity
    double F_original, G_original, modularity_fixed_original;
    // Sparse graph (fixed membership)
    int    m_sparse, preserved_intra, preserved_inter;
    double F_observed, G_sparse_observed;
    double intra_rate, inter_rate, ratio_observed;
    double modularity_fixed_sparse;
    // Derived changes
    double F_improvement_observed;  // ΔF
    double dG_observed;             // ΔG
    double modularity_fixed_change; // ΔQ_fixed
    double dQ_reconstructed;        // ΔF − ΔG (must == modularity_fixed_change)
    // Leiden re-optimisation
    double modularity_leiden_sparse;
    double modularity_leiden_change;
    // Verification
    double identity_error;          // |dQ_fixed − dQ_reconstructed|
};

// Run the complete experiment for one dataset.
std::vector<TrialResult> run_experiment(
    const Graph&            g,
    const std::string&      dataset_name,
    const ExperimentConfig& cfg);

// Write raw CSV (one row per trial) matching Python column order exactly.
void write_raw_csv(
    const std::vector<TrialResult>& results,
    const std::string&              path);

// Write summary CSV (mean/std per retention level).
void write_summary_csv(
    const std::vector<TrialResult>& results,
    const std::string&              path);
