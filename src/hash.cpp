#include "hash.h"

// Author: Jakov

HashTag Hash64(const std::string& key) {
    // Murmur3_x64_128 already produces a full 128-bit mix in a single call;
    // its two 64-bit halves are independent enough for the Kirsch-Mitzenmacher
    // trick (Bloom) and for h1/h2 reuse (Cuckoo, Bamboo). Halves the per-key
    // hash cost vs. the previous two-seed version — material for small k.
    uint64_t out[2];
    MurmurHash3_x64_128(key.data(), (int)key.size(), 0x9747b28c, out);
    return HashTag{ out[0], out[1] };
}