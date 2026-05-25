// Author: Antonio Šimić
// Benchmark our Bamboo Filter on synthetic random DNA and (optionally) on
// the E. coli K-12 MG1655 reference genome. For each combination of
// (length or genome, k), extract every k-mer via sliding window, insert
// them into the filter, then measure:
//   - insert time
//   - lookup time on inserted k-mers (no false negatives expected)
//   - lookup time + false positive rate on random non-member k-mers
//   - memory usage and bits-per-item
// Output is one CSV row per measurement; progress goes to stderr.

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "bamboo_filter.h"
#include "fasta_reader.h"
#include "kmer.h"

class Timer {
 public:
    Timer() : start_(std::chrono::high_resolution_clock::now()) {}
    double ElapsedMs() const {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(end - start_).count();
    }
 private:
    std::chrono::high_resolution_clock::time_point start_;
};

static std::vector<size_t> ParseList(const std::string& s) {
    std::vector<size_t> out;
    size_t start = 0;
    while (start < s.size()) {
        size_t end = s.find(',', start);
        if (end == std::string::npos) end = s.size();
        out.push_back(std::stoull(s.substr(start, end - start)));
        start = end + 1;
    }
    return out;
}

static void PrintUsage(const char* prog) {
    std::cerr << "Usage: " << prog
              << " --output <csv_file>"
              << " [--fasta <genome.fna>]"
              << " [--synthetic-lengths L1,L2,...]"
              << " [--k K1,K2,...]\n";
}

// One (filter, source, length, k) run — used by both the synthetic loop and
// the E. coli loop so we don't duplicate the timing + CSV-writing code.
static void RunOne(const std::vector<std::string>& kmers,
                   size_t k, size_t seq_length,
                   const std::string& data_source,
                   std::ofstream& csv) {
    BambooFilter bf([&]{
    size_t n = kmers.size();
    if (n < 2000)   return size_t(16);        // mali n: 16 segmenata
    if (n < 20000)  return std::max(size_t(16), n / 512);  // srednji n: gustoća ~0.12
    return std::max(size_t(16), n / 1024);    // veliki n: gustoća ~0.25
    }());

    Timer t_ins;
    for (const auto& km : kmers) bf.Insert(km);
    double insert_ms = t_ins.ElapsedMs();

    Timer t_pos;
    size_t found = 0;
    for (const auto& km : kmers) if (bf.Lookup(km)) found++;
    double lookup_pos_ms = t_pos.ElapsedMs();

    size_t num_negs = std::min<size_t>(kmers.size(), 100000);
    auto negs = KmerUtils::GenerateTrueNegatives(k, num_negs, kmers);

    Timer t_neg;
    size_t fp = 0;
    for (const auto& km : negs) if (bf.Lookup(km)) fp++;
    double lookup_neg_ms = t_neg.ElapsedMs();
    double fpr = negs.empty()
        ? 0.0
        : static_cast<double>(fp) / negs.size();

    double bits = kmers.empty()
        ? 0.0
        : static_cast<double>(bf.MemoryUsage()) * 8.0 / kmers.size();

    csv << "Bamboo," << data_source << "," << k << "," << seq_length << ","
        << kmers.size() << "," << insert_ms << ","
        << lookup_pos_ms << "," << lookup_neg_ms << ","
        << fpr << "," << bf.MemoryUsage() << ","
        << bits << "\n";
    csv.flush();

    std::cerr << "  inserted=" << kmers.size()
              << " found=" << found
              << " fpr=" << fpr
              << " mem=" << bf.MemoryUsage()
              << " bits/item=" << bits << "\n";
}

int main(int argc, char* argv[]) {
    std::string output_file;
    std::string fasta_path;
    std::vector<size_t> lengths = {1000, 10000, 100000, 1000000};
    std::vector<size_t> ks = {10, 20, 50, 100, 200};

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--output" && i + 1 < argc) {
            output_file = argv[++i];
        } else if (arg == "--fasta" && i + 1 < argc) {
            fasta_path = argv[++i];
        } else if (arg == "--synthetic-lengths" && i + 1 < argc) {
            lengths = ParseList(argv[++i]);
        } else if (arg == "--k" && i + 1 < argc) {
            ks = ParseList(argv[++i]);
        } else if (arg == "--help") {
            PrintUsage(argv[0]);
            return 0;
        }
    }

    if (output_file.empty()) {
        PrintUsage(argv[0]);
        return 1;
    }

    std::ofstream csv(output_file);
    if (!csv.is_open()) {
        std::cerr << "Cannot open output file: " << output_file << "\n";
        return 1;
    }

    csv << "filter,data_source,k,seq_length,num_inserted,"
        << "insert_ms,lookup_pos_ms,lookup_neg_ms,"
        << "false_positive_rate,memory_bytes,bits_per_item\n";

    // 1. Synthetic random DNA.
    for (size_t length : lengths) {
        std::string sequence = KmerUtils::GenerateRandom(length, 1)[0];
        for (size_t k : ks) {
            if (k >= length) {
                std::cerr << "  Skip k=" << k << " for length=" << length
                          << " (k >= length)\n";
                continue;
            }
            std::cerr << "Synthetic length=" << length
                      << " k=" << k << "\n";
            auto kmers = KmerUtils::Extract(sequence, k);
            RunOne(kmers, k, length, "synthetic", csv);
        }
    }

    // 2. E. coli genome (if --fasta was supplied).
    if (!fasta_path.empty()) {
        std::cerr << "Reading FASTA: " << fasta_path << "\n";
        std::string genome;
        try {
            genome = FastaReader::Read(fasta_path);
        } catch (const std::exception& ex) {
            std::cerr << "FASTA read failed: " << ex.what() << "\n";
            return 1;
        }
        std::cerr << "Genome length: " << genome.size() << " bp\n";

        for (size_t k : ks) {
            if (k >= genome.size()) continue;
            std::cerr << "E. coli k=" << k << "\n";
            auto kmers = KmerUtils::Extract(genome, k);
            RunOne(kmers, k, genome.size(), "ecoli", csv);
        }
    }

    std::cerr << "Results written to: " << output_file << "\n";
    return 0;
}
