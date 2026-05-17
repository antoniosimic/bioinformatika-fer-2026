// Author: Antonio Šimić

#include "bamboo_filter.h"

#include <cstdlib>
#include <utility>

#include "hash.h"

static size_t UpperPowerOfTwo(size_t n) {
    size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

static int Log2Power2(size_t n) {
    int r = 0;
    while ((static_cast<size_t>(1) << r) < n) r++;
    return r;
}

BambooFilter::BambooFilter(size_t initial_num_segments)
    : initial_num_segments_(
          UpperPowerOfTwo(initial_num_segments < 4 ? 4
                                                  : initial_num_segments)),
      initial_seg_bits_(Log2Power2(initial_num_segments_)),
      level_(0),
      split_pointer_(0),
      num_items_(0) {
    segments_.reserve(initial_num_segments_);
    for (size_t i = 0; i < initial_num_segments_; i++) {
        segments_.emplace_back(kBucketsPerSegment * kTagsPerBucket, 0);
    }
}

size_t BambooFilter::SegmentIndex(uint64_t hv) const {
    // At level L, the segment index uses log2(N0) + L bits of the hash,
    // starting just above the bucket bits.
    size_t base = BaseSegments();
    size_t s = static_cast<size_t>(hv >> kBucketBits) & (base - 1);
    // If this segment has already been split in the current round, the
    // item belongs to either s or its buddy (s + base). Rehash with the
    // expanded modulus to find out. (split_pointer_ is always 0 in this
    // commit; the branch is dead code until expansion is added.)
    if (s < split_pointer_) {
        s = static_cast<size_t>(hv >> kBucketBits) & ((base << 1) - 1);
    }
    return s;
}

uint16_t BambooFilter::MakeTag(uint64_t hv) const {
    // Tag occupies bits just above the segment bits at level 0. By placing
    // it there, the tag's low bits coincide with future segment-index high
    // bits — which is exactly what SplitSegment will need to decide how to
    // partition tags between source and buddy segment.
    int shift = kBucketBits + initial_seg_bits_;
    uint16_t tag = static_cast<uint16_t>(
        (hv >> shift) & ((1u << kTagBits) - 1));
    return tag == 0 ? 1 : tag;
}

size_t BambooFilter::AltBucket(size_t b, uint16_t tag) const {
    uint32_t h = static_cast<uint32_t>(tag) * 0x5bd1e995u;
    h ^= h >> 15;
    size_t r = h & (kBucketsPerSegment - 1);
    if (r == 0) r = 1;
    return b ^ r;
}

bool BambooFilter::InsertIntoBucket(size_t seg, size_t b, uint16_t tag) {
    auto& slots = segments_[seg];
    size_t base = b * kTagsPerBucket;
    for (int i = 0; i < kTagsPerBucket; i++) {
        if (slots[base + i] == 0) {
            slots[base + i] = tag;
            return true;
        }
    }
    return false;
}

bool BambooFilter::BucketContains(size_t seg, size_t b, uint16_t tag) const {
    const auto& slots = segments_[seg];
    size_t base = b * kTagsPerBucket;
    for (int i = 0; i < kTagsPerBucket; i++) {
        if (slots[base + i] == tag) return true;
    }
    return false;
}

bool BambooFilter::Insert(const std::string& key) {
    HashTag h = Hash64(key);
    uint16_t tag = MakeTag(h.h1);
    size_t s = SegmentIndex(h.h1);
    size_t b1 = BucketIndex(h.h1);
    size_t b2 = AltBucket(b1, tag);

    if (InsertIntoBucket(s, b1, tag) ||
        InsertIntoBucket(s, b2, tag)) {
        num_items_++;
        return true;
    }

    // Both candidate buckets are full — evict within this segment.
    size_t b = (std::rand() & 1) ? b1 : b2;
    for (int n = 0; n < kMaxKicks; n++) {
        size_t base = b * kTagsPerBucket;
        int slot = std::rand() % kTagsPerBucket;
        std::swap(tag, segments_[s][base + slot]);
        b = AltBucket(b, tag);
        if (InsertIntoBucket(s, b, tag)) {
            num_items_++;
            return true;
        }
    }
    return false;  // segment full and we have no expansion yet
}

bool BambooFilter::Lookup(const std::string& key) const {
    HashTag h = Hash64(key);
    uint16_t tag = MakeTag(h.h1);
    size_t s = SegmentIndex(h.h1);
    size_t b1 = BucketIndex(h.h1);
    size_t b2 = AltBucket(b1, tag);
    return BucketContains(s, b1, tag) || BucketContains(s, b2, tag);
}
