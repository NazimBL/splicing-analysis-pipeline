#!/bin/bash
# Build a STAR genome index from a FASTA + GTF, timed with /usr/bin/time -v.
# Usage: build_star_index.sh <STAR binary> <fasta> <gtf> <out_dir> [threads] [sjdbOverhang]
# The same script builds the 2.7.11b production index and the 2.7.2a benchmark
# index (Parabricks only loads the 2.7.2a format); pass the matching STAR binary.
set -euo pipefail
STAR=$1 FASTA=$2 GTF=$3 OUT=$4 THREADS=${5:-32} OVERHANG=${6:-100}
mkdir -p "$OUT"
"$STAR" --version > "$OUT/star_version.txt"
md5sum "$FASTA" "$GTF" > "$OUT/inputs.md5"
/usr/bin/time -v -o "$OUT/build_time.txt" \
  "$STAR" --runMode genomeGenerate --runThreadN "$THREADS" \
    --genomeDir "$OUT" --genomeFastaFiles "$FASTA" --sjdbGTFfile "$GTF" \
    --sjdbOverhang "$OVERHANG" --limitGenomeGenerateRAM 64000000000 \
    --outFileNamePrefix "$OUT/"
grep -E 'Elapsed|Maximum resident|Percent of CPU' "$OUT/build_time.txt"
