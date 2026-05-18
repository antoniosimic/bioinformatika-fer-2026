// K-mer extraction and random DNA generation utilities — header-only.
// Used by the benchmark and test programs to prepare query sets.
// Author: Jakov

#pragma once
#include <algorithm>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

class KmerUtils {
 public:
    // Extracts all k-mers from sequence via sliding window.
    // Returns a vector of strings, each of length k.
    static std::vector<std::string> Extract(const std::string& sequence,
                                            size_t k) {
        std::vector<std::string> kmers;
        if (sequence.size() < k) return kmers;

        kmers.reserve(sequence.size() - k + 1);
        for (size_t i = 0; i <= sequence.size() - k; i++) {
            kmers.push_back(sequence.substr(i, k));
        }
        return kmers;
    }

    // Generates n random DNA strings of length k.
    // Used as the negative query set — these are probably not in the genome
    // so any hit is a false positive (with very high probability).
    static std::vector<std::string> GenerateRandom(size_t k,
                                                   size_t n,
                                                   uint32_t seed = 123) {
        std::mt19937 rng(seed);
        std::uniform_int_distribution<> dist(0, 3);
        const char bases[] = "ACGT";

        std::vector<std::string> kmers;
        kmers.reserve(n);
        for (size_t i = 0; i < n; i++) {
            std::string km(k, ' ');
            for (size_t j = 0; j < k; j++) km[j] = bases[dist(rng)];
            kmers.push_back(km);
        }
        return kmers;
    }

    // Generates n random k-mers guaranteed NOT to be in the given set.
    // Slower than GenerateRandom but gives exact FPR measurements.
    static std::vector<std::string> GenerateTrueNegatives(
            size_t k,
            size_t n,
            const std::vector<std::string>& existing,
            uint32_t seed = 456) {

        // Build a lookup set for fast membership check
        std::unordered_set<std::string> present(existing.begin(),
                                                existing.end());
        std::mt19937 rng(seed);
        std::uniform_int_distribution<> dist(0, 3);
        const char bases[] = "ACGT";

        std::vector<std::string> negatives;
        negatives.reserve(n);
        while (negatives.size() < n) {
            std::string km(k, ' ');
            for (size_t j = 0; j < k; j++) km[j] = bases[dist(rng)];
            if (present.find(km) == present.end()) {
                negatives.push_back(km);
            }
        }
        return negatives;
    }
};