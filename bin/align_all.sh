#!/bin/bash
# Interim driver until the Nextflow pipeline exists: FastQC(raw) -> fastp -> STAR
# for every sample in a samplesheet, N samples at a time. Same per-sample scripts
# and settings the pipeline will use (docs/DECISIONS.md, docs/walkthrough_s001.md).
#
# Resumable: a sample whose BAM index and Log.final.out exist is skipped.
# Trimmed FASTQs are deleted once the sample's BAM is indexed.
#
# Usage: align_all.sh <samplesheet.csv> <raw_dir> <star_index> <outdir> [n_parallel] [star_threads]
set -euo pipefail
SHEET=$1 RAW=$2 IDX=$3 OUT=$4 NPAR=${5:-3} STHR=${6:-14}
BIN=$(cd "$(dirname "$0")" && pwd)
ADAPTER=CTGTCTCTTATACACATCT   # Nextera transposase (FastQC on s001, 2026-10-15)
mkdir -p "$OUT"/{fastqc_raw,fastp,star,trim,logs}

one () {
  local s=$1 r1=$RAW/$2 r2=$RAW/$3 st=$OUT/star/$1
  if [ -s "$st/$s.Aligned.sortedByCoord.out.bam.bai" ] && [ -s "$st/$s.Log.final.out" ]; then
    echo "[$(date +%T)] skip $s (done)"; return 0
  fi
  echo "[$(date +%T)] start $s"
  fastqc -q -t 2 -o "$OUT/fastqc_raw" "$r1" "$r2" > "$OUT/logs/$s.fastqc.log" 2>&1 &
  local fq=$!
  /usr/bin/time -v -o "$OUT/fastp/$s.time.txt" fastp -i "$r1" -I "$r2" \
    -o "$OUT/trim/${s}_R1.fq.gz" -O "$OUT/trim/${s}_R2.fq.gz" \
    --adapter_sequence $ADAPTER --adapter_sequence_r2 $ADAPTER \
    --disable_quality_filtering --length_required 20 --trim_poly_g -w 8 \
    -j "$OUT/fastp/$s.fastp.json" -h "$OUT/fastp/$s.fastp.html" 2> "$OUT/fastp/$s.fastp.log" \
    || { echo "FAIL $s fastp"; return 1; }
  "$BIN/star_align.sh" "$IDX" "$s" "$OUT/trim/${s}_R1.fq.gz" "$OUT/trim/${s}_R2.fq.gz" "$st" "$STHR" \
    > "$OUT/logs/$s.star.log" 2>&1 || { echo "FAIL $s star (see logs/$s.star.log)"; return 1; }
  wait $fq || echo "warn $s fastqc"
  rm -f "$OUT/trim/${s}_R1.fq.gz" "$OUT/trim/${s}_R2.fq.gz"
  echo "[$(date +%T)] done $s  $(grep 'Uniquely mapped reads %' "$st/$s.Log.final.out" | cut -f2)  $(grep Elapsed "$st/$s.time.txt" | awk '{print $NF}')"
}
export -f one
export RAW IDX OUT STHR BIN ADAPTER

echo "[$(date)] align_all: $(($(wc -l < "$SHEET") - 1)) samples, $NPAR parallel x $STHR STAR threads -> $OUT"
SECONDS=0
python3 -I -c 'import csv,sys
for r in csv.DictReader(open(sys.argv[1])): print(r["sample"], r["fastq_1"], r["fastq_2"])' "$SHEET" \
  | xargs -P "$NPAR" -L 1 bash -c 'one "$@"' _ || true
n=$(ls "$OUT"/star/*/*.bam.bai 2>/dev/null | wc -l)
echo "[$(date)] align_all finished: $n BAMs, wall ${SECONDS}s"
