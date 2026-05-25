# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

C++ implementation of the **Bamboo Filter** (Wang et al. 2022) — a resizable approximate membership query (AMQ) filter — applied to k-mer indexing over the *E. coli* K-12 MG1655 genome. Built for the FER *Bioinformatika 1* course (2025/2026), task 2.

Authors: Antonio Šimić (filter core, benchmark harness), Jakov Malić (hash, kmer utils, fasta reader, tests).

## Build

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

Release build uses `-O3 -march=native -flto` (non-portable; fine for local benchmarking).

To build the optional Wang et al. reference comparison (Linux + AVX2 only):
```bash
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_REFERENCE=ON
make -j$(nproc)
```

## Common commands

All binaries are built into `build/`. Run from `build/`:

**Tests:**
```bash
./test_bamboo
```

**Benchmark — synthetic data only:**
```bash
./benchmark_runner --output ../results/results.csv \
  --synthetic-lengths 1000,10000,100000,1000000 \
  --k 10,20,50,100,200
```

**Benchmark — synthetic + E. coli genome:**
```bash
./benchmark_runner --output ../results/results.csv \
  --fasta ../data/GCA_000005845_2_ASM584v2_genomic.fna \
  --synthetic-lengths 1000,10000,100000,1000000 \
  --k 10,20,50,100,200
```

**Download genome (idempotent):**
```bash
./scripts/fetch_data.sh
```

**Plot results (from repo root):**
```bash
python3 scripts/plot.py --mine results/results.csv --output results/
# With reference data:
python3 scripts/plot.py --mine results/results_mine.csv --ref results/results_ref.csv --output results/
```

## Architecture

### Library (`src/`, compiled as `filter_lib`)

- **`bamboo_filter.h/.cpp`** — core AMQ filter. Storage is `std::vector<std::vector<uint16_t>>` (vector of segments, each a flat array of `kBucketsPerSegment × kTagsPerBucket` 16-bit slots). Tag value `0xFFFF` is the empty-slot sentinel (not 0, to avoid corrupting bit 0 which drives level-0 split partitioning). Key operations:
  - `Insert` — cuckoo-style placement; triggers `SplitSegment()` when load > 90%.
  - `Delete` — removes one matching tag; triggers `MergeSegment()` when load < 40%.
  - `SplitSegment()` — appends a buddy segment, redistributes tags whose level-bit is 1, advances `split_pointer_`; promotes `level_` when pointer wraps.
  - `MergeSegment()` — reverse: pours last segment back into its buddy, removes it. Aborts if combined load wouldn't fit.
  - Hash decomposition: bits `[0,6)` → bucket index; bits `[6, 6+log2(N0_seg))` → segment index; remaining bits → 12-bit tag.

- **`hash.h/.cpp`** — `Hash64(key)` returns `HashTag{h1, h2}` via MurmurHash3 with two seeds. `h1` drives bucket placement; `h2` is available for fingerprint construction.

- **`murmurhash3.h/.cpp`** — verbatim MurmurHash3 (Appleby, public domain).

- **`kmer.h`** — header-only `KmerUtils`: `Extract(seq, k)` sliding window, `GenerateRandom(k, n)`, `GenerateTrueNegatives(k, n, existing)`.

- **`fasta_reader.h`** — header-only `FastaReader::Read(path)`: strips headers and ambiguous bases (`N`), returns uppercase ACGT-only sequence. Throws `std::runtime_error` on failure.

### Benchmark (`benchmark/benchmark_main.cpp`)

CLI tool. Initial segment count is sized from the k-mer count: `n < 2000 → 16 segs`, `n < 20000 → max(16, n/512)`, `else max(16, n/1024)`. Measures insert time, positive lookup time, negative lookup time + FPR (using `GenerateTrueNegatives`), memory, and bits-per-item. Outputs one CSV row per (data source, k, length) tuple.

### Tests (`tests/test_bamboo.cpp`)

Six standalone tests: basic correctness, delete, no false negatives (1000 elements), FPR < 1%, expand (auto-split), shrink (auto-merge).

### Plotting (`scripts/plot.py`)

Takes `--mine` CSV (and optionally `--ref` for Wang et al. reference CSV) and produces up to 5 PNGs: insert time, FPR, memory (bits/item), E. coli lookup time, and a ratio plot (ours vs reference). Uses pandas + matplotlib.

## Key invariants

- The empty-slot sentinel is `0xFFFF`, not `0`. Do not change this without updating `SplitSegment` and `MergeSegment`.
- `initial_num_segments_` is always a power of two (enforced in the constructor).
- Tags are 12-bit values stored in `uint16_t`; bit 0 of the tag is used by `SplitSegment` to decide which buddy segment a tag belongs to after a split.
- `kMaxKicks = 500` caps cuckoo eviction cycles; `Insert` returns `false` if this is exceeded (filter is full/pathological).
