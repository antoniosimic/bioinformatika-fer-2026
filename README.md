# Bamboo Filter

C++ implementation of the Bamboo Filter (Wang et al. 2022) — a resizable
approximate membership query (AMQ) filter that extends Cuckoo Filters with
smooth, incremental expansion and shrinkage.

Built for the FER course *Bioinformatika 1* (2025./2026.), task (2).

Subject page link: https://www.fer.unizg.hr/predmet/bio1

**Goal.** Search for random k-mers (k ∈ {10, 20, 50, 100, 200}) in the
*E. coli* K-12 MG1655 genome and in synthetic DNA using our own Bamboo
Filter implementation.

**Authors:** Antonio Šimić, Jakov Malić

---

## Repository structure

bioinformatika-fer-2026/
├── src/                  C++ library (Bamboo Filter + utilities)
├── tests/                Unit tests (Bamboo)
├── benchmark/            Benchmark harness + CSV output
├── reference/            Wang et al. reference impl wrapper (Linux only)
├── scripts/              Python (Matplotlib) plotting
└── data/                 E. coli K-12 MG1655 genome (.fna)

---

## Requirements

| Tool | Version |
|------|---------|
| g++ | ≥ 9 |
| CMake | ≥ 3.14 |
| Python 3 | ≥ 3.8 (plotting only) |

Install on Ubuntu / WSL:

```bash
sudo apt update && sudo apt install -y g++ cmake build-essential
pip3 install -r requirements.txt
```

---

## Build

```bash
git clone https://github.com/AntoniSimic/bioinformatika-fer-2026.git
cd bioinformatika-fer-2026
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

---

## Download genome data

The *E. coli* genome is not stored in the repository. Fetch it from NCBI:

```
./scripts/fetch_data.sh
```

This downloads ~1.4 MB (gzipped, ~4.5 MB unpacked) into `data/`.
Source: NCBI Assembly [GCA_000005845.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000005845.2/).

The script is idempotent — running it again when the file already exists is a no-op.

---

## Run tests

```bash
./test_bamboo
```

Should print `PASSED` for every test case.

## Run benchmark

Synthetic data only:

```bash
./benchmark_runner \
  --output ../results/results.csv \
  --synthetic-lengths 1000,10000,100000,1000000 \
  --k 10,20,50,100,200
```

Synthetic + E. coli genome:

```bash
./benchmark_runner \
  --output ../results/results.csv \
  --fasta ../data/GCA_000005845_2_ASM584v2_genomic.fna \
  --synthetic-lengths 1000,10000,100000,1000000 \
  --k 10,20,50,100,200
```

Output is one CSV row per (data source, k, sequence length) combination.
CSV columns: `filter, data_source, k, seq_length, num_inserted,
insert_ms, lookup_pos_ms, lookup_neg_ms, false_positive_rate,
memory_bytes, bits_per_item`

---

## Plot results

```bash
python3 scripts/plot.py \
  --mine results/results_mine.csv \
  --ref  results/results_ref.csv \
  --output results/
```

The `--ref` argument is optional — if omitted, only our results are plotted.

Generates four PNG files in `results/`:

| File | What it shows |
|------|--------------|
| `plot_time.png`         | Insert and lookup time vs k-mer size |
| `plot_fpr.png`          | False positive rate vs k-mer size |
| `plot_memory.png`       | Memory efficiency — bits per item |
| `plot_memory_total.png` | Total memory used by the filter [KB], vs k |
| `plot_ratio.png`        | Ours / reference ratio summary (when `--ref` is given) |

---

## Algorithm overview

A **Bamboo Filter** stores a compact fingerprint (12 bits) for each
inserted element rather than the element itself. Lookups check two
candidate buckets (4 slots each) using partial-key cuckoo hashing —
identical to a Cuckoo Filter — but additionally support **smooth
incremental resizing**:

- **Expand** — when load exceeds 90%, a single segment is split and its
  fingerprints are redistributed. No full rebuild required.
- **Shrink** — when load drops below 40%, the last segment is merged back
  into its parent.

Unlike a standard Cuckoo Filter which requires a full rebuild to resize,
Bamboo touches only a single segment per split or merge operation. This
makes it suitable for workloads where the number of elements is not known
in advance — such as streaming k-mer indexing over genomes of unknown size.

---

## Genome data

The benchmark uses the *E. coli* K-12 MG1655 complete genome, downloaded
on demand by `scripts/fetch_data.sh` (see [Download genome data](#download-genome-data)).

- **File:** `data/GCA_000005845_2_ASM584v2_genomic.fna` (created by fetch script)
- **Source:** NCBI Assembly — accession [GCA_000005845.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000005845.2/)
- **Length:** 4,641,652 bp

---

## References

- H. Wang et al. *Bamboo Filters: Make Resizing Smooth.* ICDE 2022.
  doi:[10.1109/ICDE53745.2022.00078](https://doi.org/10.1109/ICDE53745.2022.00078)
- H. Wang et al. *Bamboo Filters: Make Resizing Smooth and Adaptive.*
  IEEE/ACM Trans. Networking, 2024.
  doi:[10.1109/TNET.2024.3403997](https://doi.org/10.1109/TNET.2024.3403997)
- B. Fan et al. *Cuckoo Filter: Practically Better Than Bloom.*
  CoNEXT 2014.
  doi:[10.1145/2674005.2674994](https://doi.org/10.1145/2674005.2674994)
- A. Appleby. *MurmurHash3.* Public domain.
  [github.com/aappleby/smhasher](https://github.com/aappleby/smhasher)
- NCBI. *Escherichia coli str. K-12 substr. MG1655 genome assembly ASM584v2.*
  [https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000005845.2/](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000005845.2/)
- Reference implementation: [github.com/wanghanchengchn/bamboofilters](https://github.com/wanghanchengchn/bamboofilters)