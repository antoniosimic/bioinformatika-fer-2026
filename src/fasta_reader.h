// FASTA file reader — header-only utility.
// Reads a .fna / .fasta file and returns the concatenated nucleotide
// sequence as a single uppercase string. Header lines (starting with '>')
// and whitespace are stripped automatically.
// Author: Jakov

#pragma once
#include <cctype>
#include <fstream>
#include <stdexcept>
#include <string>

class FastaReader {
 public:
    // Reads the file at `path` and returns the full concatenated sequence.
    // Throws std::runtime_error if the file cannot be opened.
    static std::string Read(const std::string& path) {
        std::ifstream file(path);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open FASTA file: " + path);
        }

        std::string sequence;
        std::string line;

        while (std::getline(file, line)) {
            // Skip header lines and empty lines
            if (line.empty() || line[0] == '>') continue;

            // Append only valid nucleotide characters, uppercased
            for (char c : line) {
                c = std::toupper(c);
                if (c == 'A' || c == 'T' || c == 'G' || c == 'C') {
                    sequence += c;
                }
                // 'N' and other ambiguous bases are silently skipped
            }
        }

        if (sequence.empty()) {
            throw std::runtime_error("No sequence data found in: " + path);
        }

        return sequence;
    }

    // Returns basic stats about the sequence without storing it
    struct Stats {
        size_t length;
        size_t count_a;
        size_t count_t;
        size_t count_g;
        size_t count_c;
    };
};