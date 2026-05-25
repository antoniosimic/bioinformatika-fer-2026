#!/bin/bash
# Downloads E. coli K-12 MG1655 reference genome from NCBI.
# Output: data/GCA_000005845_2_ASM584v2_genomic.fna
# Author: Jakov Malić
set -e

DATA_DIR="$(dirname "$0")/../data"
mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.fna.gz"
OUT="GCA_000005845_2_ASM584v2_genomic.fna"

if [ -f "$OUT" ]; then
    echo "Already downloaded: $OUT"
    exit 0
fi

echo "Downloading E. coli K-12 MG1655 genome from NCBI..."
wget -q --show-progress "$URL" -O "${OUT}.gz"
gunzip "${OUT}.gz"
echo "Downloaded: $DATA_DIR/$OUT ($(wc -c < "$OUT") bytes)"