#include "hash.h"

// Author: Jakov

HashTag Hash64(const std::string& key) {
    uint64_t out1[2];
    uint64_t out2[2];

    // Two calls with different seeds → two independent hash values
    MurmurHash3_x64_128(key.data(), (int)key.size(), 0x9747b28c, out1);
    MurmurHash3_x64_128(key.data(), (int)key.size(), 0x85ebca6b, out2);

    return HashTag{ out1[0], out2[0] };
}