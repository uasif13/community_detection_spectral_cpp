#include "config.hpp"
#include <toml.hpp>
#include <stdexcept>

AppConfig load_config(const std::string& path) {
    auto data = toml::parse(path);

    AppConfig cfg;

    // [experiment]
    const auto& exp = toml::find(data, "experiment");
    cfg.exp.n_trials  = toml::find<int>   (exp, "n_trials");
    cfg.exp.base_seed = static_cast<uint64_t>(toml::find<int64_t>(exp, "base_seed"));
    cfg.exp.run_leiden_resample = toml::find<bool>(exp, "run_leiden_resample");
    cfg.output_dir = toml::find<std::string>(exp, "output_dir");

    auto alphas_raw = toml::find<std::vector<double>>(exp, "alphas");
    cfg.exp.alphas = alphas_raw;

    // [[datasets]]
    const auto& datasets = toml::find(data, "datasets");
    for (const auto& d : datasets.as_array()) {
        DatasetEntry entry;
        entry.name = toml::find<std::string>(d, "name");
        entry.url  = toml::find<std::string>(d, "url");
        cfg.datasets.push_back(entry);
    }

    return cfg;
}
