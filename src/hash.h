#pragma once
#include <cstdint>
#include <string>
#include "murmurhash3.h"

// Hashes a k-mer string into two independent 64-bit values.
// Uses two different seeds so the outputs are independent.
// Author: Jakov
struct HashTag {
    uint64_t h1;  // used for bucket index
    uint64_t h2;  // used for fingerprint
};

// Returns a HashTag for the given key string
HashTag Hash64(const std::string& key);