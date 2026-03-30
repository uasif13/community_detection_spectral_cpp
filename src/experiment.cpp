#include "experiment.hpp"
#include "dspar.hpp"
#include "sampler.hpp"
#include "modularity.hpp"
#include "leiden.hpp"
#include "csv_writer.hpp"
#include "rng.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

// ── internal helpers ──────────────────────────────────────────────────────────

static TrialResult build_trial_result(
    double              alpha,
    int                 trial,
    uint64_t            seed,
    const DsparScores&  scores,
    const Graph&        g,
    const ModResult&    orig,
    const ModResult&    sp,
    double              leiden_q,
    double              modularity_fixed_original)
{
    TrialResult r;
    r.retention  = alpha;
    r.trial      = trial;
    r.seed       = seed;

    r.m  = g.m;
    r.n1 = scores.n1;
    r.n2 = scores.n2;

    r.mu_intra = scores.mu_intra;
    r.mu_inter = scores.mu_inter;
    r.delta    = scores.delta;

    r.F_original               = orig.F;
    r.G_original               = orig.G;
    r.modularity_fixed_original = modularity_fixed_original;

    r.m_sparse        = sp.m_sparse;
    r.preserved_intra = sp.preserved_intra;
    r.preserved_inter = sp.preserved_inter;
    r.F_observed      = sp.F;
    r.G_sparse_observed = sp.G;
    r.intra_rate      = sp.intra_rate;
    r.inter_rate      = sp.inter_rate;
    r.ratio_observed  = sp.ratio_obs;
    r.modularity_fixed_sparse = sp.Q;

    r.F_improvement_observed  = sp.F - orig.F;          // ΔF
    r.dG_observed             = sp.G - orig.G;          // ΔG
    r.modularity_fixed_change = sp.Q - modularity_fixed_original;  // ΔQ_fixed
    r.dQ_reconstructed        = r.F_improvement_observed - r.dG_observed;
    r.identity_error          = std::abs(r.modularity_fixed_change - r.dQ_reconstructed);

    r.modularity_leiden_sparse = leiden_q;
    r.modularity_leiden_change = leiden_q - modularity_fixed_original;

    return r;
}

// ── public API ────────────────────────────────────────────────────────────────

std::vector<TrialResult> run_experiment(
    const Graph&            g,
    const std::string&      /*dataset_name*/,
    const ExperimentConfig& cfg)
{
    // ── 1. One-time setup ────────────────────────────────────────────────────
    DsparScores scores = compute_scores(g);

    LeidenResult leiden_orig = run_leiden(g);
    const std::vector<int>& membership_fixed = leiden_orig.membership;

    // Original-graph modularity with fixed membership
    ModResult orig = compute_modularity(
        g, nullptr, membership_fixed, scores.n1, scores.n2);
    // Note: scores.n1/n2 are 0 until compute_structural_metrics is called.
    // Call it now — it also fills mu_intra, mu_inter, delta.
    compute_structural_metrics(scores, g, membership_fixed);
    // Recompute orig with correct n1/n2
    orig = compute_modularity(g, nullptr, membership_fixed, scores.n1, scores.n2);

    const double Q_fixed_original = orig.Q;

    // ── 2. Per-alpha, per-trial loop ─────────────────────────────────────────
    const int n_alphas = static_cast<int>(cfg.alphas.size());
    std::vector<TrialResult> results;
    results.reserve(n_alphas * cfg.n_trials);

    std::vector<bool> keep_buf;
    keep_buf.reserve(g.m);

    for (int ai = 0; ai < n_alphas; ai++) {
        double alpha = cfg.alphas[ai];
        for (int t = 0; t < cfg.n_trials; t++) {
            uint64_t seed = cfg.base_seed
                + static_cast<uint64_t>(t)  * 1000003ULL
                + static_cast<uint64_t>(ai) * 7919ULL;

            Xoshiro256pp rng(seed);
            sample_edges_inplace(scores, alpha, rng, keep_buf);

            ModResult sp = compute_modularity(
                g, &keep_buf, membership_fixed, scores.n1, scores.n2);

            double leiden_q = 0.0;
            if (cfg.run_leiden_resample) {
                LeidenResult lr = run_leiden(g, &keep_buf);
                leiden_q = lr.modularity;
            }

            results.push_back(build_trial_result(
                alpha, t, seed, scores, g, orig, sp, leiden_q, Q_fixed_original));
        }
    }

    return results;
}

// ── CSV output ────────────────────────────────────────────────────────────────

void write_raw_csv(
    const std::vector<TrialResult>& results,
    const std::string&              path)
{
    CsvWriter w(path);
    w.write_header({
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
    });

    for (auto& r : results) {
        w.col(r.retention)
         .col(static_cast<long long>(r.seed))
         .col(r.m).col(r.n1).col(r.n2)
         .col(r.mu_intra).col(r.mu_inter).col(r.delta)
         .col(r.F_original).col(r.G_original).col(r.modularity_fixed_original)
         .col(r.m_sparse).col(r.F_observed).col(r.G_sparse_observed)
         .col(r.dG_observed).col(r.F_improvement_observed)
         .col(r.modularity_fixed_sparse).col(r.modularity_fixed_change)
         .col(r.modularity_leiden_sparse).col(r.modularity_leiden_change)
         .col(r.dQ_reconstructed)
         .col(r.preserved_intra).col(r.preserved_inter)
         .col(r.intra_rate).col(r.inter_rate).col(r.ratio_observed);
        w.end_row();
    }
}

void write_summary_csv(
    const std::vector<TrialResult>& results,
    const std::string&              path)
{
    // Collect unique retention values in order
    std::vector<double> alphas;
    for (auto& r : results) {
        if (std::find(alphas.begin(), alphas.end(), r.retention) == alphas.end())
            alphas.push_back(r.retention);
    }
    std::sort(alphas.begin(), alphas.end());

    CsvWriter w(path);
    w.write_header({
        "retention",
        "modularity_fixed_change_mean", "modularity_fixed_change_std",
        "modularity_leiden_change_mean", "modularity_leiden_change_std",
        "F_improvement_observed_mean", "F_improvement_observed_std",
        "dG_observed_mean", "dG_observed_std",
        "delta_mean"
    });

    for (double alpha : alphas) {
        // Gather rows for this alpha
        std::vector<double> dq, dl, df, dg, delta;
        for (auto& r : results) {
            if (r.retention != alpha) continue;
            dq.push_back(r.modularity_fixed_change);
            dl.push_back(r.modularity_leiden_change);
            df.push_back(r.F_improvement_observed);
            dg.push_back(r.dG_observed);
            delta.push_back(r.delta);
        }
        int n = static_cast<int>(dq.size());
        if (n == 0) continue;

        auto mean_std = [&](const std::vector<double>& v) -> std::pair<double,double> {
            double mu = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
            double var = 0.0;
            for (double x : v) var += (x - mu) * (x - mu);
            // ddof=1 to match Python default numpy.std with ddof=1 (sample std)
            double std = (v.size() > 1) ? std::sqrt(var / (v.size() - 1)) : 0.0;
            return {mu, std};
        };

        auto [dq_mu, dq_sd] = mean_std(dq);
        auto [dl_mu, dl_sd] = mean_std(dl);
        auto [df_mu, df_sd] = mean_std(df);
        auto [dg_mu, dg_sd] = mean_std(dg);
        double delta_mu = std::accumulate(delta.begin(), delta.end(), 0.0) / delta.size();

        w.col(alpha)
         .col(dq_mu).col(dq_sd)
         .col(dl_mu).col(dl_sd)
         .col(df_mu).col(df_sd)
         .col(dg_mu).col(dg_sd)
         .col(delta_mu);
        w.end_row();
    }
}
