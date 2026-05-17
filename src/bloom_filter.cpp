#include "bloom_filter.h"
#include "hash.h"
#include <cmath>

// Author: Jakov

// Optimal m (bits) and k (hash functions) from n and fpr:
//   m = -n * ln(fpr) / (ln2)^2
//   k =  m/n * ln2
BloomFilter::BloomFilter(size_t n, double fpr) {
    double ln2 = std::log(2.0);
    m_ = (size_t)std::ceil(-1.0 * (double)n * std::log(fpr) / (ln2 * ln2));
    k_ = (size_t)std::ceil((double)m_ / (double)n * ln2);
    bits_.assign(m_, false);
}

size_t BloomFilter::GetBitIndex(const std::string& key, size_t i) const {
    HashTag h = Hash64(key);
    // Kirsch-Mitzenmacher: combine two hashes to simulate k independent ones
    return (h.h1 + i * h.h2) % m_;
}

void BloomFilter::Insert(const std::string& key) {
    for (size_t i = 0; i < k_; i++) {
        bits_[GetBitIndex(key, i)] = true;
    }
}

bool BloomFilter::Lookup(const std::string& key) const {
    for (size_t i = 0; i < k_; i++) {
        if (!bits_[GetBitIndex(key, i)]) return false;
    }
    return true;
}

size_t BloomFilter::CountSetBits() const {
    size_t count = 0;
    for (bool b : bits_) count += b;
    return count;
}