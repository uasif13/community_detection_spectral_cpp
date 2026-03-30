#pragma once
#include <string>

// Download url and decompress .gz to local_path (plain .txt).
// Returns true on success. Prints progress to stderr.
bool download_and_decompress(
    const std::string& url,
    const std::string& local_path);

// Ensure dataset exists locally; download+decompress if not.
// Returns path to the decompressed .txt file.
std::string ensure_dataset(
    const std::string& name,
    const std::string& url,
    const std::string& data_dir = "data");
