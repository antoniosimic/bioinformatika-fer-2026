// MurmurHash3 by Austin Appleby — public domain
// https://github.com/aappleby/smhasher
#pragma once
#include <cstdint>

void MurmurHash3_x86_32(const void* key, int len, uint32_t seed, void* out);
void MurmurHash3_x64_128(const void* key, int len, uint32_t seed, void* out);