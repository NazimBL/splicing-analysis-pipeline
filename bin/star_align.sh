#!/bin/bash
# Align one paired-end sample with STAR, timed with /usr/bin/time -v.
#
# One parameter set for every sample, layout and aligner:
#  - 2-pass (--twopassMode Basic): pass 1 finds this sample's novel junctions,
#    pass 2 realigns with them added, so reads across new junctions are placed.
#  - ENCODE long-RNA filters (STAR manual 3.3.2), minus --outFilterType BySJout,
#    which Parabricks rna_fq2bam does not offer; the GPU run must match.
#  - XS strand tag (--outSAMstrandField intronMotif) for unstranded libraries.
#  - MD tag, read group = sample, and GeneCounts (used to infer strandedness).
#
# Usage: star_align.sh <index_dir> <sample> <R1.fq.gz> <R2.fq.gz> <out_dir> [threads]
set -euo pipefail
IDX=$1 SAMPLE=$2 R1=$3 R2=$4 OUT=$5 THREADS=${6:-16}
STAR=${STAR:-STAR}
mkdir -p "$OUT"
"$STAR" --version > "$OUT/star_version.txt"

/usr/bin/time -v -o "$OUT/${SAMPLE}.time.txt" \
"$STAR" --runThreadN "$THREADS" \
  --genomeDir "$IDX" \
  --readFilesIn "$R1" "$R2" --readFilesCommand zcat \
  --outFileNamePrefix "$OUT/${SAMPLE}." \
  --twopassMode Basic \
  --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 30000000000 \
  --outSAMattributes NH HI AS nM MD \
  --outSAMattrRGline ID:"$SAMPLE" SM:"$SAMPLE" LB:"$SAMPLE" PL:ILLUMINA \
  --outSAMstrandField intronMotif \
  --quantMode GeneCounts \
  --outFilterMultimapNmax 20 \
  --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
  --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
  --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000

samtools index -@ 4 "$OUT/${SAMPLE}.Aligned.sortedByCoord.out.bam"
rm -rf "$OUT/${SAMPLE}._STARgenome" "$OUT/${SAMPLE}._STARpass1" "$OUT/${SAMPLE}._STARtmp"
grep -E 'Elapsed|Maximum resident' "$OUT/${SAMPLE}.time.txt"
