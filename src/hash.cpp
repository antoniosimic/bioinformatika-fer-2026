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

uint8_t Fingerprint(const std::string& key) {
    HashTag h = Hash64(key);
    // Clamp to 8 bits, ensure never 0 (0 = empty slot in buckets)
    uint8_t fp = (uint8_t)(h.h2 & 0xFF);
    return fp == 0 ? 1 : fp;
}