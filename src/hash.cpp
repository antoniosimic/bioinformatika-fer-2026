#include "hash.h"

// Author: Jakov

HashTag Hash64(const std::string& key) {
    uint64_t out[2];
    MurmurHash3_x64_128(key.data(), (int)key.size(), 0x9747b28c, out);
    return HashTag{ out[0], out[1] };
}