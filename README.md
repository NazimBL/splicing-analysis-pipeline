# Splicing analysis pipeline

An alternative-splicing pipeline for bulk RNA-seq: STAR for alignment, rMATS for splicing
events. It is set up for four macaque (*Macaca mulatta*) viral-infection datasets (SIV,
SARS-CoV-2, H5N1, mpox; H5N1 includes cynomolgus). It works on any dataset through config
files, with no script edits per dataset.

> **Status: under construction.** The alignment step and the environment are in place. The
> rMATS steps and the Nextflow workflow are not written yet. See [Status](#status).

## Why this exists

An earlier pipeline for these datasets was hand-edited for each dataset and silently lost
reads. Rebuilding it, we keep to these rules. Some are already in the scripts; the checks for the
rest are listed in [Status](#status).

- **Read length** is measured on the raw reads and passed to rMATS. A wrong `--readLength`
  makes rMATS silently drop reads.
- **Clipped reads** are kept with `--allow-clipping`. Multi-mapped reads (NH != 1) are always
  dropped by rMATS.
- **Strandedness** is inferred from STAR gene counts, not hardcoded.
- **Layouts are not mixed.** Single-end and paired-end samples are never compared in one
  contrast. Both layouts use the same STAR settings.
- **Paired animals** are analysed as paired (`--paired-stats`), not as independent groups.
- **Event IDs stay comparable.** rMATS runs once per virus x tissue on a shared event set,
  not once per comparison.
- **Output paths come from the contrast table**, never typed by hand, so one comparison
  cannot overwrite another.

## Pipeline

```mermaid
flowchart LR
    FASTQ[FASTQ] --> TRIM["fastp (adapter + poly-G trimming)"]
    TRIM --> STAR["STAR 2-pass alignment (sorted BAM)"]
    STAR --> PREP["rMATS --task prep (per BAM)"]
    PREP --> POST["rMATS --task post (per virus x tissue)"]
    POST --> STAT["rMATS --task stat (per contrast)"]
```

FastQC runs on the raw and trimmed reads for QC.

## Repository layout

| Path | Contents |
|------|----------|
| `envs/` | `splice.yaml` (the environment, with pinned key tools) and `splice.lock.yaml` (exact versions) |
| `bin/` | Command-line scripts: STAR index build, per-sample alignment, and the sample driver |
| `datasets/<virus>/` | One folder per dataset: `samplesheet.csv`, `contrasts.csv`, `params.yaml` |

## Requirements

- Linux, conda or mamba.
- Java 17 or newer, and Nextflow 26.04.6 once the workflow exists.

Create the environment from the pinned lock file:

```bash
conda env create -f envs/splice.lock.yaml   # exact versions
conda activate splice
```

Key tools and versions: STAR 2.7.11b, rMATS-turbo 4.4.0 with PAIRADISE 1.26.0 (needed by
`--paired-stats`), fastp 1.3.7, FastQC, samtools, MultiQC, seqkit, pigz.

## Reference

- Ensembl 110, Mmul_10 toplevel FASTA, and the **full** Ensembl GTF. Do not use a filtered
  GTF (for example the Cell Ranger one): it drops transcripts that splicing analysis needs.
- STAR index built with STAR 2.7.11b and `--sjdbOverhang 100`. An index from a STAR version
  older than 2.7.4a will not load.

Build the index (paths are arguments; nothing is hardcoded):

```bash
bin/build_star_index.sh <STAR> <genome.fa> <annotation.gtf> <index_dir> [threads] [sjdbOverhang]
```

## Running a dataset

Each dataset is described by three files in `datasets/<virus>/`:

- `samplesheet.csv`: sample, animal, tissue, timepoint, and the FASTQ file names.
- `contrasts.csv`: the two groups to compare, and whether the comparison is paired by animal.
- `params.yaml`: input and output locations.

Align all samples of a dataset (resumable; a finished sample is skipped):

```bash
bin/align_all.sh <samplesheet.csv> <raw_fastq_dir> <star_index> <outdir> [n_parallel] [star_threads]
```

For each sample this runs FastQC on the raw reads, fastp, then `bin/star_align.sh` (STAR
2-pass, with the same settings for every sample). It writes BAMs and STAR logs under
`<outdir>/star/`.

## Status

| Step | State |
|------|-------|
| Environment (`envs/`) | done |
| STAR index build (`bin/build_star_index.sh`) | done |
| Trimming, QC, alignment (`bin/align_all.sh`) | done, run by hand per dataset |
| rMATS prep / post / stat | not yet written |
| Nextflow workflow | not yet written |
| Read-length, strandedness and layout checks | not yet written (`bin/star_align.sh` writes the GeneCounts they need) |
| Single-end alignment | not yet written (`bin/star_align.sh` is paired-end only) |
| GPU alignment (NVIDIA Parabricks) | planned, to compare speed and concordance with the CPU run |
| Variant calling, then AlphaGenome | later |

## Contributing

Open an issue or a pull request. Please do not commit data or results: FASTQ, BAM, rMATS
output tables and figures are excluded by `.gitignore`.
