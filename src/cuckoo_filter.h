// Author: Antonio Šimić
// Cuckoo Filter (Fan et al. 2014) — fixed-size AMQ with 12-bit fingerprints
// stored in 4-slot buckets. Supports insert / lookup / delete via partial-key
// cuckoo hashing. Slot value 0 is reserved as the "empty" marker; the
// fingerprint computation bumps a zero result to 1 so it can never collide
// with the empty sentinel.

#ifndef CUCKOO_FILTER_H_
#define CUCKOO_FILTER_H_

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

class CuckooFilter {
 public:
    // Sizes the table for ~half-full at expected_items inserts; bucket
    // count is rounded up to a power of two so we can use a bit mask.
    explicit CuckooFilter(size_t expected_items);

    // Returns false if the filter is too full (after kMaxKicks evictions).
    bool Insert(const std::string& key);
    bool Lookup(const std::string& key) const;
    bool Delete(const std::string& key);

    size_t Size() const { return num_items_; }
    size_t NumBuckets() const { return num_buckets_; }
    size_t MemoryUsage() const { return buckets_.size() * sizeof(uint16_t); }

 private:
    static constexpr int kTagsPerBucket = 4;
    static constexpr int kMaxKicks = 500;

    // Mixing function used to compute the alternate bucket from a tag.
    // Must be nonzero so that AltIndex(i, fp) != i.
    size_t TagHash(uint16_t fp) const;

    bool InsertIntoBucket(size_t b, uint16_t fp);
    bool BucketContains(size_t b, uint16_t fp) const;
    bool DeleteFromBucket(size_t b, uint16_t fp);

    // Flat storage: bucket b lives at indices [b*kTagsPerBucket, ...).
    // We store each 12-bit tag in a uint16_t for clarity (4 unused bits per
    // slot) rather than bit-packing — much easier to read than the bit
    // arithmetic in the reference impl.
    std::vector<uint16_t> buckets_;
    size_t num_buckets_;   // power of two
    size_t bucket_mask_;   // num_buckets_ - 1
    size_t num_items_;
};

#endif  // CUCKOO_FILTER_H_
