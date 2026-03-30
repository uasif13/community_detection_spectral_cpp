#include "config.hpp"
#include "downloader.hpp"
#include "experiment.hpp"
#include "graph.hpp"

#include <chrono>
#include <filesystem>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

// ── CLI parsing ───────────────────────────────────────────────────────────────

struct CliArgs {
    std::string config_path  = "configs/exp1.toml";
    std::string output_dir   = "";          // empty = use TOML value
    std::string datasets_filter = "";       // empty = run all
    int64_t     seed_override  = -1;        // -1 = use TOML value
    int         trials_override = -1;       // -1 = use TOML value
    bool        resume   = false;
    bool        dry_run  = false;
    bool        no_leiden = false;
    bool        help     = false;
};

static void print_help(const char* argv0) {
    std::cout <<
        "Usage: " << argv0 << " [options]\n\n"
        "Options:\n"
        "  --config PATH         TOML config file [default: configs/exp1.toml]\n"
        "  --datasets D1,D2,...  Run only these datasets (comma-separated)\n"
        "  --seed N              Override base_seed\n"
        "  --output DIR          Override output_dir\n"
        "  --trials N            Override n_trials\n"
        "  --resume              Skip datasets that already have output files\n"
        "  --dry-run             Print config and exit\n"
        "  --no-leiden           Skip Leiden re-optimisation (faster)\n"
        "  --help                Show this help\n";
}

static CliArgs parse_args(int argc, char** argv) {
    CliArgs args;
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc)
                throw std::runtime_error("Expected value after " + a);
            return argv[++i];
        };
        if      (a == "--config")   args.config_path     = next();
        else if (a == "--datasets") args.datasets_filter = next();
        else if (a == "--seed")     args.seed_override   = std::stoll(next());
        else if (a == "--output")   args.output_dir      = next();
        else if (a == "--trials")   args.trials_override = std::stoi(next());
        else if (a == "--resume")   args.resume   = true;
        else if (a == "--dry-run")  args.dry_run  = true;
        else if (a == "--no-leiden")args.no_leiden = true;
        else if (a == "--help")     args.help     = true;
        else throw std::runtime_error("Unknown option: " + a);
    }
    return args;
}

// ── helpers ───────────────────────────────────────────────────────────────────

static std::set<std::string> parse_dataset_filter(const std::string& s) {
    std::set<std::string> result;
    if (s.empty()) return result;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ','))
        if (!tok.empty()) result.insert(tok);
    return result;
}

static void print_config(const AppConfig& cfg) {
    std::cout << "=== Configuration ===\n";
    std::cout << "  output_dir : " << cfg.output_dir << "\n";
    std::cout << "  n_trials   : " << cfg.exp.n_trials << "\n";
    std::cout << "  base_seed  : " << cfg.exp.base_seed << "\n";
    std::cout << "  run_leiden : " << (cfg.exp.run_leiden_resample ? "yes" : "no") << "\n";
    std::cout << "  alphas     : [";
    for (int i = 0; i < (int)cfg.exp.alphas.size(); i++) {
        if (i) std::cout << ", ";
        std::cout << cfg.exp.alphas[i];
    }
    std::cout << "]\n";
    std::cout << "  datasets   :\n";
    for (auto& d : cfg.datasets)
        std::cout << "    " << d.name << "  (" << d.url << ")\n";
}

// ── main ──────────────────────────────────────────────────────────────────────

int main(int argc, char** argv) {
    CliArgs cli;
    try {
        cli = parse_args(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Argument error: " << e.what() << "\n";
        print_help(argv[0]);
        return 1;
    }

    if (cli.help) { print_help(argv[0]); return 0; }

    // Load config
    AppConfig cfg;
    try {
        cfg = load_config(cli.config_path);
    } catch (const std::exception& e) {
        std::cerr << "Config error: " << e.what() << "\n";
        return 1;
    }

    // Apply CLI overrides
    if (!cli.output_dir.empty())   cfg.output_dir          = cli.output_dir;
    if (cli.seed_override >= 0)    cfg.exp.base_seed        = static_cast<uint64_t>(cli.seed_override);
    if (cli.trials_override > 0)   cfg.exp.n_trials         = cli.trials_override;
    if (cli.no_leiden)             cfg.exp.run_leiden_resample = false;

    if (cli.dry_run) { print_config(cfg); return 0; }

    // Parse dataset filter
    auto filter = parse_dataset_filter(cli.datasets_filter);

    // Create output directory
    try {
        fs::create_directories(cfg.output_dir);
        fs::create_directories("data");
    } catch (const std::exception& e) {
        std::cerr << "Cannot create directories: " << e.what() << "\n";
        return 1;
    }

    // Run each dataset
    int n_done = 0, n_skipped = 0, n_failed = 0;
    for (auto& ds : cfg.datasets) {
        if (!filter.empty() && !filter.count(ds.name)) continue;

        std::string raw_csv = cfg.output_dir + "/" + ds.name
                              + "_theoretical_validation_FIXED.csv";
        std::string sum_csv = cfg.output_dir + "/" + ds.name
                              + "_summary_FIXED.csv";

        if (cli.resume && fs::exists(raw_csv) && fs::exists(sum_csv)) {
            std::cout << "[" << ds.name << "] Skipping (--resume, output exists)\n";
            n_skipped++;
            continue;
        }

        auto wall_start = std::chrono::steady_clock::now();
        std::cout << "[" << ds.name << "] Downloading/checking dataset...\n";

        std::string data_path;
        try {
            data_path = ensure_dataset(ds.name, ds.url, "data");
        } catch (const std::exception& e) {
            std::cerr << "[" << ds.name << "] FAILED (download): " << e.what() << "\n";
            n_failed++;
            continue;
        }

        std::cout << "[" << ds.name << "] Loading graph...\n";
        Graph g;
        try {
            g = load_snap(data_path);
        } catch (const std::exception& e) {
            std::cerr << "[" << ds.name << "] FAILED (load): " << e.what() << "\n";
            n_failed++;
            continue;
        }
        std::cout << "  n=" << g.n << " m=" << g.m << "\n";

        std::cout << "[" << ds.name << "] Running experiment ("
                  << cfg.exp.alphas.size() << " alphas × "
                  << cfg.exp.n_trials << " trials)...\n";
        std::vector<TrialResult> results;
        try {
            results = run_experiment(g, ds.name, cfg.exp);
        } catch (const std::exception& e) {
            std::cerr << "[" << ds.name << "] FAILED (experiment): " << e.what() << "\n";
            n_failed++;
            continue;
        }

        // Verify identity
        double max_err = 0.0;
        for (auto& r : results)
            if (r.identity_error > max_err) max_err = r.identity_error;
        std::cout << "  max identity error: " << max_err << "\n";

        write_raw_csv(results, raw_csv);
        write_summary_csv(results, sum_csv);

        auto wall_end = std::chrono::steady_clock::now();
        double secs = std::chrono::duration<double>(wall_end - wall_start).count();
        std::cout << "[" << ds.name << "] Done in " << secs << "s. Wrote:\n"
                  << "  " << raw_csv << "\n"
                  << "  " << sum_csv << "\n";
        n_done++;
    }

    std::cout << "\nFinished: " << n_done << " done, "
              << n_skipped << " skipped, "
              << n_failed  << " failed.\n";
    return (n_failed > 0) ? 1 : 0;
}
