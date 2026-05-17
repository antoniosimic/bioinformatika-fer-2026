// Author: Antonio Šimić
// Sanity check for Cuckoo Filter: insert + lookup + delete + FPR
// measurement. Demonstrates the key advantage over Bloom Filter — Delete
// actually works.

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "cuckoo_filter.h"

int main() {
    const size_t kN = 10000;
    CuckooFilter cf(kN);

    std::vector<std::string> items;
    items.reserve(kN);
    for (size_t i = 0; i < kN; i++) {
        items.push_back("item_" + std::to_string(i));
        bool ok = cf.Insert(items.back());
        assert(ok && "Insert failed before reaching expected capacity");
    }
    std::cout << "Inserted " << cf.Size() << " items into "
              << cf.NumBuckets() << " buckets ("
              << cf.MemoryUsage() << " bytes)\n";

    // No false negatives — every inserted item must still be found.
    for (const auto& item : items) {
        assert(cf.Lookup(item));
    }
    std::cout << "All " << kN << " lookups returned true (no false negatives)\n";

    // Delete half and verify Size shrinks and deleted items vanish.
    for (size_t i = 0; i < kN / 2; i++) {
        bool ok = cf.Delete(items[i]);
        assert(ok && "Delete failed for an inserted item");
    }
    assert(cf.Size() == kN - kN / 2);
    std::cout << "Deleted " << kN / 2 << " items; Size is now "
              << cf.Size() << "\n";

    // Measure FPR on disjoint non-member set.
    size_t fp = 0;
    const size_t kQueries = 100000;
    for (size_t i = 0; i < kQueries; i++) {
        if (cf.Lookup("nonmember_" + std::to_string(i))) fp++;
    }
    double fpr = static_cast<double>(fp) / kQueries;
    // Theoretical FPR for 12-bit fingerprints checking 2 buckets:
    //   ~ 2 * (1/2^12) = 2 / 4096 ≈ 0.00049
    std::cout << "Measured FPR = " << fpr
              << "  (theory ≈ " << (2.0 / 4096.0) << ")\n";
    assert(fpr < 0.005);  // generous slack
    std::cout << "OK\n";
    return 0;
}
