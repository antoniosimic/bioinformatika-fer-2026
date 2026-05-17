// Bloom Filter test — insert real k-mers, query with random ones, measure FPR
// Author: Jakov

#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <cassert>
#include "bloom_filter.h"

// Generates a random DNA string of length k
std::string RandomKmer(size_t k, std::mt19937& rng) {
    const std::string bases = "ACGT";
    std::uniform_int_distribution<> dist(0, 3);
    std::string result(k, ' ');
    for (size_t i = 0; i < k; i++) result[i] = bases[dist(rng)];
    return result;
}

// Extracts all k-mers from a sequence via sliding window
std::vector<std::string> ExtractKmers(const std::string& seq, size_t k) {
    std::vector<std::string> kmers;
    if (seq.size() < k) return kmers;
    for (size_t i = 0; i <= seq.size() - k; i++) {
        kmers.push_back(seq.substr(i, k));
    }
    return kmers;
}

void TestFPR(size_t k, size_t seq_len, double target_fpr) {
    std::mt19937 rng(42); // fixed seed — reproducible results

    // Build a fake "genome" of random DNA
    std::string genome = RandomKmer(seq_len, rng);

    // Extract all k-mers and insert into filter
    auto kmers = ExtractKmers(genome, k);
    BloomFilter bf(kmers.size(), target_fpr);
    for (const auto& km : kmers) bf.Insert(km);

    // Sanity check — every inserted k-mer must be found (no false negatives)
    for (const auto& km : kmers) {
        assert(bf.Lookup(km) && "False negative detected — this should never happen!");
    }

    // Measure false positive rate with random queries
    size_t queries = 10000;
    size_t false_positives = 0;
    for (size_t i = 0; i < queries; i++) {
        std::string q = RandomKmer(k, rng);
        // Only count as FP if the k-mer was not actually inserted
        if (bf.Lookup(q)) {
            bool actually_present = false;
            for (const auto& km : kmers) {
                if (km == q) { actually_present = true; break; }
            }
            if (!actually_present) false_positives++;
        }
    }

    double measured_fpr = (double)false_positives / (double)queries;

    std::cout << "k=" << k
              << "  seq_len=" << seq_len
              << "  inserted=" << kmers.size()
              << "  bits=" << bf.NumBits()
              << "  target_fpr=" << target_fpr
              << "  measured_fpr=" << measured_fpr
              << "\n";
}

int main() {
    std::cout << "=== Bloom Filter Tests ===\n\n";

    // 1. Basic correctness
    BloomFilter bf(100, 0.01);
    bf.Insert("ACGTACGT");
    assert(bf.Lookup("ACGTACGT") && "Inserted element not found!");
    assert(!bf.Lookup("ZZZZZZZZ") && "Impossible element returned true!");
    std::cout << "Basic correctness: PASSED\n\n";

    // 2. FPR measurement for different k values
    std::cout << "False positive rate measurements:\n";
    TestFPR(10,  1000,  0.01);
    TestFPR(20,  1000,  0.01);
    TestFPR(50,  5000,  0.01);
    TestFPR(100, 10000, 0.01);

    // 3. Show the no-delete limitation
    std::cout << "\nDelete limitation demo:\n";
    BloomFilter bf2(10, 0.01);
    bf2.Insert("AAAA");
    bf2.Insert("CCCC");
    // Can't delete — if we could zero out bits, we'd corrupt other entries
    std::cout << "Lookup AAAA after insert: "
              << (bf2.Lookup("AAAA") ? "FOUND" : "NOT FOUND") << "\n";
    std::cout << "No delete supported — this is why Cuckoo Filter exists.\n";

    return 0;
}


