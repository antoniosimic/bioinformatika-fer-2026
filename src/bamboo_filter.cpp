// Author: Antonio Šimić

#include "bamboo_filter.h"

#include <array>
#include <cstdint>
#include <cstdlib>
#include <utility>

#include "hash.h"

namespace {
// Fast per-thread pseudo-random generator used only in the cuckoo eviction
// path. xorshift64 is a few XOR/shift instructions vs. std::rand()'s LCG +
// internal state, and being thread_local it avoids any library-level
// synchronisation. Quality is more than enough for picking a random bucket
// / slot during eviction — we are not using it for anything cryptographic.
inline uint64_t FastRand() {
    thread_local uint64_t state = 0x9E3779B97F4A7C15ULL;
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return state;
}
}  // namespace

constexpr uint16_t BambooFilter::kEmpty;

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
        segments_.emplace_back(kBucketsPerSegment * kTagsPerBucket, kEmpty);
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
    // partition tags between source and buddy segment. We MUST NOT remap
    // any value here (e.g. bump 0→1) because that would break the bit-level
    // correspondence with the hash; the empty-slot collision is handled by
    // using kEmpty (0xFFFF) as the empty sentinel instead of 0.
    int shift = kBucketBits + initial_seg_bits_;
    return static_cast<uint16_t>(
        (hv >> shift) & ((1u << kTagBits) - 1));
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
        if (slots[base + i] == kEmpty) {
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

bool BambooFilter::DeleteFromBucket(size_t seg, size_t b, uint16_t tag) {
    auto& slots = segments_[seg];
    size_t base = b * kTagsPerBucket;
    for (int i = 0; i < kTagsPerBucket; i++) {
        if (slots[base + i] == tag) {
            slots[base + i] = kEmpty;
            return true;
        }
    }
    return false;
}

bool BambooFilter::TryInsertOnce(uint16_t tag, size_t s, size_t b1) {
    size_t b2 = AltBucket(b1, tag);
    if (InsertIntoBucket(s, b1, tag) ||
        InsertIntoBucket(s, b2, tag)) {
        return true;
    }
    // Cuckoo eviction within segment s. We snapshot the segment first so
    // that if the chain fails after kMaxKicks we can revert — otherwise
    // the last carried tag would be silently dropped, producing a real
    // false negative for a previously-inserted item.
    //
    // Snapshot lives on the stack (std::array, fixed at the compile-time
    // segment capacity) to avoid a per-eviction heap alloc/free, which the
    // old std::vector path paid every time we kicked.
    constexpr size_t kSegmentSlots = kBucketsPerSegment * kTagsPerBucket;
    std::array<uint16_t, kSegmentSlots> backup;
    const auto& seg_data = segments_[s];
    for (size_t i = 0; i < kSegmentSlots; i++) backup[i] = seg_data[i];

    size_t b = (FastRand() & 1) ? b1 : b2;
    for (int n = 0; n < kMaxKicks; n++) {
        size_t base = b * kTagsPerBucket;
        // kTagsPerBucket is a power of two so we mask instead of dividing —
        // a single AND vs. std::rand's modulo on a global state.
        int slot = static_cast<int>(FastRand() & (kTagsPerBucket - 1));
        std::swap(tag, segments_[s][base + slot]);
        b = AltBucket(b, tag);
        if (InsertIntoBucket(s, b, tag)) {
            return true;
        }
    }
    auto& seg_restore = segments_[s];
    for (size_t i = 0; i < kSegmentSlots; i++) seg_restore[i] = backup[i];
    return false;
}

bool BambooFilter::Insert(const std::string& key) {
    HashTag h = Hash64(key);
    uint16_t tag = MakeTag(h.h1);

    // Dedupe: as a set filter (per the Bamboo / Cuckoo Filter contract),
    // inserting the same key twice is a no-op. We have to check before
    // splitting because duplicate keys share the same hash → same tag →
    // same tag bit at every level, so SplitSegment can never distribute
    // them across the source/buddy pair. Without this guard, a stream of
    // duplicate inserts (e.g. dup k-mers from a small alphabet) packs the
    // same bucket until it overflows, then loops splitting forever until
    // the level cap is hit, blowing memory to the maximum.
    {
        size_t s = SegmentIndex(h.h1);
        size_t b1 = BucketIndex(h.h1);
        size_t b2 = AltBucket(b1, tag);
        if (BucketContains(s, b1, tag) || BucketContains(s, b2, tag)) {
            return true;
        }
    }

    // Try; if the target segment is full, split and retry. SplitSegment
    // only splits segments_[split_pointer_], so if our target s is well
    // ahead of split_pointer_ we may need many splits before s itself
    // gets touched. Loop until either we succeed or SplitSegment is a
    // no-op (level exhausted — tag bits all used up as routing bits).
    while (true) {
        size_t s = SegmentIndex(h.h1);
        size_t b1 = BucketIndex(h.h1);
        if (TryInsertOnce(tag, s, b1)) {
            num_items_++;
            // Pre-emptively split when the overall load gets uncomfortable
            // so we don't have to chase failures on the next insertion.
            if (LoadFactor() > kExpandThreshold) {
                SplitSegment();
            }
            return true;
        }
        size_t prev_sp = split_pointer_;
        size_t prev_level = level_;
        SplitSegment();
        if (split_pointer_ == prev_sp && level_ == prev_level) {
            return false;  // can't split further; filter is genuinely full
        }
    }
}

void BambooFilter::SplitSegment() {
    // Safety: we can't split past the bits available in the tag.
    if (level_ >= kTagBits) return;

    size_t p = split_pointer_;
    size_t new_seg = segments_.size();
    segments_.emplace_back(kBucketsPerSegment * kTagsPerBucket, kEmpty);

    // Move tags whose `level_`-th bit is 1 from source segment p into the
    // new buddy segment, keeping the bucket-within-segment position. The
    // bit choice works because MakeTag was set up so tag bit L equals the
    // hash bit that becomes the new high bit of the segment index after
    // the L-th doubling round.
    auto& src = segments_[p];
    auto& dst = segments_[new_seg];
    for (size_t b = 0; b < kBucketsPerSegment; b++) {
        size_t base = b * kTagsPerBucket;
        int dst_slot = 0;
        for (int slot = 0; slot < kTagsPerBucket; slot++) {
            uint16_t tag = src[base + slot];
            if (tag != kEmpty && ((tag >> level_) & 1)) {
                dst[base + dst_slot++] = tag;
                src[base + slot] = kEmpty;
            }
        }
    }

    split_pointer_++;
    if (split_pointer_ >= BaseSegments()) {
        // Completed a full doubling round — promote to the next level.
        level_++;
        split_pointer_ = 0;
    }
}

bool BambooFilter::Delete(const std::string& key) {
    HashTag h = Hash64(key);
    uint16_t tag = MakeTag(h.h1);
    size_t s = SegmentIndex(h.h1);
    size_t b1 = BucketIndex(h.h1);
    size_t b2 = AltBucket(b1, tag);
    if (DeleteFromBucket(s, b1, tag) ||
        DeleteFromBucket(s, b2, tag)) {
        num_items_--;
        if (LoadFactor() < kShrinkThreshold) {
            MergeSegment();
        }
        return true;
    }
    return false;
}

void BambooFilter::MergeSegment() {
    // Never shrink below the initial allocation.
    if (segments_.size() <= initial_num_segments_) return;

    // The last segment is always the most recently split-off buddy. Its
    // parent (target of the merge) is the segment we split FROM in that
    // round. With BaseSegments = N0 << level_ and a fresh split_pointer_
    // value, the parent index is (split_pointer_ - 1) at the current level,
    // wrapping back into the previous level when split_pointer_ == 0.
    bool will_decrement_level = (split_pointer_ == 0);
    if (will_decrement_level && level_ == 0) return;  // already at base

    size_t source = segments_.size() - 1;
    size_t target;
    if (will_decrement_level) {
        // The buddy we're absorbing is from the previous round.
        target = (BaseSegments() / 2) - 1;
    } else {
        target = split_pointer_ - 1;
    }

    // Refuse to merge if the combined load is uncomfortably close to a
    // single segment's capacity — we'd just expand again immediately.
    size_t src_count = 0, dst_count = 0;
    for (auto t : segments_[source]) if (t != kEmpty) src_count++;
    for (auto t : segments_[target]) if (t != kEmpty) dst_count++;
    size_t cap = kBucketsPerSegment * kTagsPerBucket;
    if (src_count + dst_count > static_cast<size_t>(cap * 0.95)) return;

    // Move tags back into the parent, preserving bucket-within-segment.
    for (size_t b = 0; b < kBucketsPerSegment; b++) {
        size_t base = b * kTagsPerBucket;
        for (int slot = 0; slot < kTagsPerBucket; slot++) {
            uint16_t tag = segments_[source][base + slot];
            if (tag == kEmpty) continue;
            TryInsertOnce(tag, target, b);
            segments_[source][base + slot] = kEmpty;
        }
    }

    segments_.pop_back();
    if (will_decrement_level) {
        level_--;
        split_pointer_ = BaseSegments() - 1;
    } else {
        split_pointer_--;
    }
}

bool BambooFilter::Lookup(const std::string& key) const {
    HashTag h = Hash64(key);
    uint16_t tag = MakeTag(h.h1);
    size_t s = SegmentIndex(h.h1);
    size_t b1 = BucketIndex(h.h1);
    size_t b2 = AltBucket(b1, tag);
    return BucketContains(s, b1, tag) || BucketContains(s, b2, tag);
}
