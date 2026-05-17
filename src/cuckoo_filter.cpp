// Author: Antonio Šimić

#include "cuckoo_filter.h"

#include <cstdlib>
#include <utility>

#include "hash.h"

static size_t UpperPowerOfTwo(size_t n) {
    size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

CuckooFilter::CuckooFilter(size_t expected_items) : num_items_(0) {
    // Aim for ~50% load at expected_items so cuckoo evictions stay rare.
    // 4 slots per bucket × num_buckets = total capacity.
    size_t target_buckets = (expected_items + 1) / 2;
    if (target_buckets < 4) target_buckets = 4;
    num_buckets_ = UpperPowerOfTwo(target_buckets);
    bucket_mask_ = num_buckets_ - 1;
    buckets_.assign(num_buckets_ * kTagsPerBucket, 0);
}

size_t CuckooFilter::TagHash(uint16_t fp) const {
    // Small integer mixer; the masked result must be nonzero so that the
    // alternate bucket is distinct from the primary.
    uint32_t h = static_cast<uint32_t>(fp) * 0x5bd1e995u;
    h ^= h >> 15;
    size_t r = h & bucket_mask_;
    return r == 0 ? 1 : r;
}

bool CuckooFilter::InsertIntoBucket(size_t b, uint16_t fp) {
    size_t base = b * kTagsPerBucket;
    for (int i = 0; i < kTagsPerBucket; i++) {
        if (buckets_[base + i] == 0) {
            buckets_[base + i] = fp;
            return true;
        }
    }
    return false;
}

bool CuckooFilter::BucketContains(size_t b, uint16_t fp) const {
    size_t base = b * kTagsPerBucket;
    for (int i = 0; i < kTagsPerBucket; i++) {
        if (buckets_[base + i] == fp) return true;
    }
    return false;
}

bool CuckooFilter::DeleteFromBucket(size_t b, uint16_t fp) {
    size_t base = b * kTagsPerBucket;
    for (int i = 0; i < kTagsPerBucket; i++) {
        if (buckets_[base + i] == fp) {
            buckets_[base + i] = 0;
            return true;
        }
    }
    return false;
}

bool CuckooFilter::Insert(const std::string& key) {
    HashTag h = Hash64(key);
    // 12-bit fingerprint from the lower bits of the second hash; bump 0 to 1
    // to preserve the "0 == empty" invariant in buckets_.
    uint16_t fp = static_cast<uint16_t>(h.h2 & 0xFFFu);
    if (fp == 0) fp = 1;
    size_t i1 = static_cast<size_t>(h.h1) & bucket_mask_;
    size_t i2 = i1 ^ TagHash(fp);

    if (InsertIntoBucket(i1, fp) || InsertIntoBucket(i2, fp)) {
        num_items_++;
        return true;
    }

    // Both candidate buckets are full — start cuckoo eviction at a random
    // one of them. On each kick, swap the carried fp with a random slot
    // and chase the displaced tag to its alternate bucket.
    size_t i = (std::rand() & 1) ? i1 : i2;
    for (int n = 0; n < kMaxKicks; n++) {
        size_t base = i * kTagsPerBucket;
        int slot = std::rand() % kTagsPerBucket;
        std::swap(fp, buckets_[base + slot]);
        i = i ^ TagHash(fp);
        if (InsertIntoBucket(i, fp)) {
            num_items_++;
            return true;
        }
    }
    return false;  // filter is effectively full
}

bool CuckooFilter::Lookup(const std::string& key) const {
    HashTag h = Hash64(key);
    uint16_t fp = static_cast<uint16_t>(h.h2 & 0xFFFu);
    if (fp == 0) fp = 1;
    size_t i1 = static_cast<size_t>(h.h1) & bucket_mask_;
    size_t i2 = i1 ^ TagHash(fp);
    return BucketContains(i1, fp) || BucketContains(i2, fp);
}

bool CuckooFilter::Delete(const std::string& key) {
    HashTag h = Hash64(key);
    uint16_t fp = static_cast<uint16_t>(h.h2 & 0xFFFu);
    if (fp == 0) fp = 1;
    size_t i1 = static_cast<size_t>(h.h1) & bucket_mask_;
    size_t i2 = i1 ^ TagHash(fp);
    if (DeleteFromBucket(i1, fp) || DeleteFromBucket(i2, fp)) {
        num_items_--;
        return true;
    }
    return false;
}
