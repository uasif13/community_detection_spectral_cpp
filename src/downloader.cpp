#include "downloader.hpp"
#include <curl/curl.h>
#include <zlib.h>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace fs = std::filesystem;

// ── libcurl write callback: appends to a byte vector ─────────────────────────

static size_t curl_write_cb(void* ptr, size_t size, size_t nmemb, void* userdata) {
    auto* buf = static_cast<std::vector<char>*>(userdata);
    size_t total = size * nmemb;
    buf->insert(buf->end(), static_cast<char*>(ptr),
                             static_cast<char*>(ptr) + total);
    return total;
}

// ── download URL into memory ──────────────────────────────────────────────────

static std::vector<char> curl_download(const std::string& url) {
    CURL* curl = curl_easy_init();
    if (!curl) throw std::runtime_error("curl_easy_init failed");

    std::vector<char> buf;
    curl_easy_setopt(curl, CURLOPT_URL,           url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, curl_write_cb);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA,     &buf);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_USERAGENT,     "dspar-exp1/1.0");
    curl_easy_setopt(curl, CURLOPT_FAILONERROR,   1L);

    CURLcode res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);

    if (res != CURLE_OK)
        throw std::runtime_error(
            std::string("curl download failed: ") + curl_easy_strerror(res));
    return buf;
}

// ── zlib decompress gz bytes to file ─────────────────────────────────────────

static void gz_decompress(
    const std::vector<char>& compressed,
    const std::string&       out_path)
{
    // zlib inflateInit2 with windowBits=31 handles gzip format
    z_stream zs{};
    if (inflateInit2(&zs, 31) != Z_OK)
        throw std::runtime_error("inflateInit2 failed");

    zs.next_in  = reinterpret_cast<Bytef*>(
                    const_cast<char*>(compressed.data()));
    zs.avail_in = static_cast<uInt>(compressed.size());

    std::ofstream out(out_path, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot write to: " + out_path);

    constexpr size_t CHUNK = 65536;
    std::vector<char> outbuf(CHUNK);
    int ret;
    do {
        zs.next_out  = reinterpret_cast<Bytef*>(outbuf.data());
        zs.avail_out = static_cast<uInt>(CHUNK);
        ret = inflate(&zs, Z_NO_FLUSH);
        if (ret == Z_STREAM_ERROR || ret == Z_DATA_ERROR || ret == Z_MEM_ERROR) {
            inflateEnd(&zs);
            throw std::runtime_error("zlib inflate error");
        }
        size_t have = CHUNK - zs.avail_out;
        out.write(outbuf.data(), static_cast<std::streamsize>(have));
    } while (ret != Z_STREAM_END);

    inflateEnd(&zs);
}

// ── public API ────────────────────────────────────────────────────────────────

bool download_and_decompress(
    const std::string& url,
    const std::string& local_path)
{
    try {
        std::cerr << "  Downloading " << url << " ...\n";
        auto compressed = curl_download(url);
        std::cerr << "  Decompressing (" << compressed.size()
                  << " bytes) -> " << local_path << " ...\n";
        gz_decompress(compressed, local_path);
        return true;
    } catch (const std::exception& e) {
        std::cerr << "  ERROR: " << e.what() << '\n';
        return false;
    }
}

std::string ensure_dataset(
    const std::string& name,
    const std::string& url,
    const std::string& data_dir)
{
    fs::create_directories(data_dir);
    std::string local_path = data_dir + "/" + name + ".txt";
    if (fs::exists(local_path)) {
        std::cerr << "[" << name << "] Using cached " << local_path << "\n";
        return local_path;
    }
    if (!download_and_decompress(url, local_path))
        throw std::runtime_error("Failed to download dataset: " + name);
    return local_path;
}
