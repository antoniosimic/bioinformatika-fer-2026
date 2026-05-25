#pragma once
#include <cstdint>
#include <string>
#include "murmurhash3.h"

// Hashes a k-mer string into a pair of 64-bit values derived from a
// single MurmurHash3_x64_128 call (the two halves of the 128-bit output).
// Author: Jakov Malić

struct HashTag {
    uint64_t h1;  // used for segment + bucket index (and tag bits)
    uint64_t h2;  // spare (kept for symmetry / future use)
};

// Returns a HashTag for the given key string.
HashTag Hash64(const std::string& key);