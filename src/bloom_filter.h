#pragma once
#include <cstdint>
#include <string>
#include <vector>

// Bloom Filter implementation using Kirsch-Mitzenmacher trick:
// generates k hash functions from only 2 base hashes.
// Supports insert and lookup — NO delete (fundamental limitation).
// Author: Jakov

class BloomFilter {
public:
    // n = expected number of elements
    // fpr = desired false positive rate (e.g. 0.01 = 1%)
    BloomFilter(size_t n, double fpr);

    // Insert a k-mer string into the filter
    void Insert(const std::string& key);

    // Returns true if key was probably inserted, false if definitely not
    bool Lookup(const std::string& key) const;

    // Returns current number of bits set to 1
    size_t CountSetBits() const;

    // Returns total number of bits in the filter
    size_t NumBits() const { return m_; }

private:
    size_t m_;               // total number of bits
    size_t k_;               // number of hash functions
    std::vector<bool> bits_; // the bit array

    // Returns the i-th hash position for key
    // Uses Kirsch-Mitzenmacher: h_i(x) = h1(x) + i * h2(x)
    size_t GetBitIndex(const std::string& key, size_t i) const;
};