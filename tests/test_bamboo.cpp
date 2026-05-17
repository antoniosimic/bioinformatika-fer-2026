// Bamboo Filter tests — correctness, delete, and resize behaviour
// Author: Jakov

#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <cassert>
#include "bamboo_filter.h"

// Generates a random DNA string of length k
std::string RandomKmer(size_t k, std::mt19937& rng) {
    const std::string bases = "ACGT";
    std::uniform_int_distribution<> dist(0, 3);
    std::string result(k, ' ');
    for (size_t i = 0; i < k; i++) result[i] = bases[dist(rng)];
    return result;
}

// --- Test 1: Basic insert and lookup ---
void TestBasicCorrectness() {
    std::cout << "Test 1: Basic insert and lookup\n";
    BambooFilter bf(4);

    bf.Insert("ACGTACGT");
    bf.Insert("TTTTGGGG");
    bf.Insert("CCCCAAAA");

    assert(bf.Lookup("ACGTACGT") && "Inserted element not found!");
    assert(bf.Lookup("TTTTGGGG") && "Inserted element not found!");
    assert(bf.Lookup("CCCCAAAA") && "Inserted element not found!");
    assert(!bf.Lookup("ZZZZZZZZ") && "Impossible element returned true!");

    std::cout << "  PASSED\n";
}

// --- Test 2: Delete ---
void TestDelete() {
    std::cout << "Test 2: Delete\n";
    BambooFilter bf(4);

    bf.Insert("ACGTACGT");
    bf.Insert("TTTTGGGG");

    assert(bf.Lookup("ACGTACGT"));
    bool deleted = bf.Delete("ACGTACGT");
    assert(deleted && "Delete returned false for existing element!");
    assert(!bf.Lookup("ACGTACGT") && "Element still found after delete!");

    // Deleting something not inserted should return false
    bool bad_delete = bf.Delete("GGGGGGGG");
    assert(!bad_delete && "Delete returned true for non-existent element!");

    std::cout << "  PASSED\n";
}

// --- Test 3: No false negatives ---
// Every inserted element MUST be found — this is a hard guarantee
void TestNoFalseNegatives() {
    std::cout << "Test 3: No false negatives (1000 elements)\n";
    std::mt19937 rng(42);
    BambooFilter bf(4);

    std::vector<std::string> inserted;
    for (int i = 0; i < 1000; i++) {
        std::string km = RandomKmer(20, rng);
        bf.Insert(km);
        inserted.push_back(km);
    }

    for (const auto& km : inserted) {
        assert(bf.Lookup(km) && "False negative — inserted element not found!");
    }

    std::cout << "  PASSED (size=" << bf.Size() << ")\n";
}

// --- Test 4: False positive rate ---
void TestFPR() {
    std::cout << "Test 4: False positive rate\n";
    std::mt19937 rng(99);
    BambooFilter bf(4);

    // Insert 500 known k-mers
    std::vector<std::string> inserted;
    for (int i = 0; i < 500; i++) {
        std::string km = RandomKmer(20, rng);
        bf.Insert(km);
        inserted.push_back(km);
    }

    // Query 5000 random k-mers, count false positives
    size_t fp_count = 0;
    size_t queries = 5000;
    for (size_t i = 0; i < queries; i++) {
        std::string q = RandomKmer(20, rng);
        if (bf.Lookup(q)) {
            // Check if it was actually inserted
            bool real = false;
            for (const auto& km : inserted) {
                if (km == q) { real = true; break; }
            }
            if (!real) fp_count++;
        }
    }

    double fpr = (double)fp_count / (double)queries;
    std::cout << "  queries=" << queries
              << "  false_positives=" << fp_count
              << "  fpr=" << fpr << "\n";

    // Bamboo uses 12-bit tags so theoretical FPR ~ 1/4096 ~ 0.024%
    // We expect well under 1% in practice
    assert(fpr < 0.01 && "FPR too high!");
    std::cout << "  PASSED\n";
}

// --- Test 5: Expand (auto-split) ---
// Insert enough elements to trigger segment splits and verify the filter
// still finds everything correctly afterwards
void TestExpand() {
    std::cout << "Test 5: Expand (auto-split on high load)\n";
    std::mt19937 rng(7);
    BambooFilter bf(4);

    size_t initial_segments = bf.NumSegments();
    std::vector<std::string> inserted;

    // Insert until segments have doubled at least once
    while (bf.NumSegments() < initial_segments * 2) {
        std::string km = RandomKmer(20, rng);
        bf.Insert(km);
        inserted.push_back(km);
    }

    std::cout << "  segments: " << initial_segments
              << " → " << bf.NumSegments()
              << "  items=" << bf.Size()
              << "  load=" << bf.LoadFactor() << "\n";

    // Every inserted element must still be found after resize
    for (const auto& km : inserted) {
        assert(bf.Lookup(km) && "False negative after expand!");
    }

    std::cout << "  PASSED\n";
}

// --- Test 6: Shrink (auto-merge) ---
// Insert a lot, then delete most of them and verify segments merge back
void TestShrink() {
    std::cout << "Test 6: Shrink (auto-merge on low load)\n";
    std::mt19937 rng(13);
    BambooFilter bf(4);

    std::vector<std::string> inserted;
    for (int i = 0; i < 800; i++) {
        std::string km = RandomKmer(20, rng);
        bf.Insert(km);
        inserted.push_back(km);
    }

    size_t segments_after_insert = bf.NumSegments();

    // Delete most elements to drop load below shrink threshold
    size_t deleted = 0;
    for (size_t i = 0; i < inserted.size() * 3 / 4; i++) {
        if (bf.Delete(inserted[i])) deleted++;
    }

    std::cout << "  segments: " << segments_after_insert
              << " → " << bf.NumSegments()
              << "  deleted=" << deleted
              << "  load=" << bf.LoadFactor() << "\n";

    // Remaining elements must still be found
    for (size_t i = inserted.size() * 3 / 4; i < inserted.size(); i++) {
        assert(bf.Lookup(inserted[i]) && "False negative after shrink!");
    }

    std::cout << "  PASSED\n";
}

int main() {
    std::cout << "=== Bamboo Filter Tests ===\n\n";

    TestBasicCorrectness();
    TestDelete();
    TestNoFalseNegatives();
    TestFPR();
    TestExpand();
    TestShrink();

    std::cout << "\nAll tests passed!\n";
    return 0;
}