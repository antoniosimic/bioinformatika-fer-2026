// Author: Antonio Šimić
// Bamboo Filter (Wang et al. 2022) — a Cuckoo-style AMQ filter that supports
// smooth, incremental resizing. Insert auto-splits segments as the filter
// fills; Delete auto-merges them when it empties. Both operations touch a
// single segment at a time, so resizing is amortised constant per insert
// (no stop-the-world rebuild).
//
// Storage is organised as a vector of fixed-size segments. Each segment
// holds `kBucketsPerSegment` buckets, each bucket holds `kTagsPerBucket`
// 12-bit tags stored in uint16_t slots (slot value 0 == empty).
//
// Hash decomposition (using the first 64-bit hash from hash.h):
//   bits [0, 6)                 → bucket index within segment
//   bits [6, 6 + log2(N0_seg))  → segment index at level 0
//   bits [6 + log2(N0_seg), …)  → 12-bit tag (overlaps with segment-index
//                                 high bits — see SplitSegment in the next
//                                 commit for why)

#ifndef BAMBOO_FILTER_H_
#define BAMBOO_FILTER_H_

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

class BambooFilter {
 public:
    static constexpr int kBucketsPerSegment = 64;   // 2^6
    static constexpr int kBucketBits = 6;
    static constexpr int kTagsPerBucket = 4;
    static constexpr int kTagBits = 12;

    // initial_num_segments is rounded up to a power of two. Each segment
    // is preallocated to hold kBucketsPerSegment * kTagsPerBucket tags.
    explicit BambooFilter(size_t initial_num_segments = 4);

    bool Insert(const std::string& key);
    bool Lookup(const std::string& key) const;
    bool Delete(const std::string& key);

    // Introspection (used by the resize-demo plot in a later commit).
    size_t Size() const { return num_items_; }
    size_t NumSegments() const { return segments_.size(); }
    size_t NumBuckets() const {
        return segments_.size() * kBucketsPerSegment;
    }
    size_t Level() const { return level_; }
    size_t SplitPointer() const { return split_pointer_; }
    double LoadFactor() const {
        return num_items_ /
            static_cast<double>(NumBuckets() * kTagsPerBucket);
    }
    size_t MemoryUsage() const {
        return segments_.size() *
            kBucketsPerSegment * kTagsPerBucket * sizeof(uint16_t);
    }

 private:
    static constexpr int kMaxKicks = 500;
    // Auto-split when the filter exceeds this overall load factor,
    // auto-merge when it drops below this one.
    static constexpr double kExpandThreshold = 0.9;
    static constexpr double kShrinkThreshold = 0.4;

    // Number of segments before any splits in the current round.
    size_t BaseSegments() const {
        return initial_num_segments_ << level_;
    }

    // Hash decomposition.
    size_t SegmentIndex(uint64_t hv) const;
    size_t BucketIndex(uint64_t hv) const {
        return static_cast<size_t>(hv) & (kBucketsPerSegment - 1);
    }
    uint16_t MakeTag(uint64_t hv) const;

    // Alternate bucket within the same segment. The XOR is masked to
    // kBucketsPerSegment so it stays valid, and forced nonzero so the
    // alt is distinct from the primary.
    size_t AltBucket(size_t b, uint16_t tag) const;

    bool InsertIntoBucket(size_t seg, size_t b, uint16_t tag);
    bool BucketContains(size_t seg, size_t b, uint16_t tag) const;
    bool DeleteFromBucket(size_t seg, size_t b, uint16_t tag);

    // One attempt to place `tag` into segment `s`, starting at bucket b1
    // and falling back to the alt bucket / cuckoo eviction. Returns false
    // if the segment is full (caller can then split and retry).
    bool TryInsertOnce(uint16_t tag, size_t s, size_t b1);

    // Append a buddy segment for segments_[split_pointer_] and move tags
    // whose `level_`-th bit is 1 into it. Advances split_pointer_; if it
    // wraps past BaseSegments(), promotes the filter to the next level.
    void SplitSegment();

    // Reverse of SplitSegment: pour the tags from the last segment back
    // into its buddy and remove the last segment. Aborts safely if the
    // combined load would no longer fit, or if we'd shrink below the
    // initial segment count.
    void MergeSegment();

    // segments_[s] is a flat vector of kBucketsPerSegment * kTagsPerBucket
    // tags. Bucket b within that segment occupies indices [b*4, b*4 + 4).
    std::vector<std::vector<uint16_t>> segments_;
    size_t initial_num_segments_;   // N0_seg, power of two
    int initial_seg_bits_;          // log2(initial_num_segments_)
    size_t level_;                  // L: number of completed doublings
    size_t split_pointer_;          // p: next segment to split (unused
                                    //    until the next commit)
    size_t num_items_;
};

#endif  // BAMBOO_FILTER_H_
