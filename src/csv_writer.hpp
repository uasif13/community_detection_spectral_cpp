#pragma once
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

// Simple header-only CSV writer.
// Call write_header() once, then write_row() for each row.
class CsvWriter {
public:
    explicit CsvWriter(const std::string& path)
        : out_(path), header_written_(false)
    {
        if (!out_.is_open())
            throw std::runtime_error("CsvWriter: cannot open " + path);
        out_ << std::fixed << std::setprecision(15);
    }

    // Write comma-separated header names
    void write_header(std::initializer_list<std::string> cols) {
        bool first = true;
        for (auto& c : cols) {
            if (!first) out_ << ',';
            out_ << c;
            first = false;
        }
        out_ << '\n';
        header_written_ = true;
    }

    // Stream-style row builder — call col() for each value, then end_row()
    CsvWriter& col(int v)    { write_sep(); out_ << v;  return *this; }
    CsvWriter& col(double v) { write_sep(); out_ << v;  return *this; }
    CsvWriter& col(long long v) { write_sep(); out_ << v; return *this; }
    CsvWriter& col(uint64_t v)  { write_sep(); out_ << v; return *this; }

    void end_row() { out_ << '\n'; first_in_row_ = true; }

private:
    std::ofstream out_;
    bool header_written_;
    bool first_in_row_ = true;

    void write_sep() {
        if (!first_in_row_) out_ << ',';
        first_in_row_ = false;
    }
};
