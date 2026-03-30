#pragma once
#include "experiment.hpp"
#include <string>
#include <vector>

struct DatasetEntry {
    std::string name;
    std::string url;
};

struct AppConfig {
    ExperimentConfig exp;
    std::string      output_dir  = "results/exp1";
    std::string      data_dir    = "data";
    std::vector<DatasetEntry> datasets;
};

// Load AppConfig from a TOML file.
AppConfig load_config(const std::string& path);
