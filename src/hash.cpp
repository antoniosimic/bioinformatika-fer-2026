#include "hash.h"

// Author: Jakov Malić

// Hashes `key` with MurmurHash3_x64_128 and returns both 64-bit halves of
// the 128-bit output as a HashTag. One call produces two values that are
// independent enough for our needs (h1 drives segment/bucket/tag bits;
// h2 is reserved). MurmurHash3 is chosen for speed and good distribution
// on short DNA strings; the fixed seed makes runs reproducible.
HashTag Hash64(const std::string& key) {
    uint64_t out[2];
    MurmurHash3_x64_128(key.data(), (int)key.size(), 0x9747b28c, out);
    return HashTag{ out[0], out[1] };
}