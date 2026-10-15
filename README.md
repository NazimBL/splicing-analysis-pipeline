#  Splicing-Analysis  
A Modular RNA-seq Alternative Splicing Pipeline 

This repository provides a **modular**, and **reproducible** pipeline for performing alternative splicing analysis on RNA-seq data using **rMATS-turbo**, followed by **post-processing**, **gene annotation**, and **transcript-level interpretation**.

This workflow was originally developed for the splicing analysis of multiple viral infection models in *Macaca mulatta* (rhesus macaque), but the scripts can adapted to *any species or dataset.

## What we are doing
The pipeline runs STAR alignment and rMATS splicing analysis on bulk RNA-seq for four macaque
viral-infection datasets (SIV, SARS-CoV-2, H5N1, mpox; H5N1 includes cynomolgus). It works on
any dataset through config files. It replaces an earlier PhD pipeline (git history at commit
d7a1d7b), which was hand-edited for each dataset and silently lost reads.

### Order of work
1. **CPU pipeline (current):** STAR + rMATS, correct and documented, with runtime measured per step.
2. **GPU:** the same pipeline with NVIDIA Parabricks for alignment; compare speed and concordance.
3. **Later:** variant calling, then AlphaGenome.

### Goals
- One command per dataset. Nobody edits scripts to run a new dataset.
- Anyone can rerun it: pinned environment, configs in the repo, no hardcoded paths.
- Data and results are never committed: FASTQ, BAM, rMATS tables, figures.

### Environment
- conda env `splice`, built from `envs/splice.yaml`; exact versions in `envs/splice.lock.yaml`:
  STAR 2.7.11b, rMATS 4.4.0 + PAIRADISE 1.26.0 (needed by `--paired-stats`),
  Trimmomatic 0.39, samtools, FastQC, MultiQC, seqkit, pigz.
- Nextflow 26.04.6 (pinned with NXF_VER), Java >= 17.
- Reference: Ensembl 110 Mmul_10, toplevel FASTA plus the full GTF. Do not use a Cell Ranger
  filtered GTF, which drops transcripts that splicing analysis needs.
- STAR index built with STAR 2.7.11b and `--sjdbOverhang 100`. An index from STAR older than
  2.7.4a will not load.

### Design
One config folder per dataset, `datasets/<virus>/`:
- `samplesheet.csv`: sample, animal, tissue, timepoint, fastq paths.
- `contrasts.csv`: which two groups to compare, and whether the comparison is paired by animal.
- `params.yaml`.

The pipeline works out the rest itself and refuses to run on inconsistent input:
- Single- vs paired-end: from whether a second FASTQ exists. A contrast must not mix layouts.
- Read length: measured on the raw reads and passed to rMATS together with
  `--variable-read-length`. Warn if the groups in a contrast differ.
- Strandedness: inferred from STAR `--quantMode GeneCounts` (ReadsPerGene column 4 vs
  column 3). Stop if the samples disagree.

Always:
- rMATS runs with `--allow-clipping`.
- Use `--paired-stats` wherever the same animals appear in both groups, with b1 and b2
  ordered by animal.
- One STAR index and one parameter set for every layout: `--twopassMode Basic`, read
  groups, MD tags, GeneCounts.
- rMATS runs `--task prep` once per BAM, then one `--task post --statoff` per virus x tissue
  (a shared event set), then `--task stat` for each contrast. This keeps event IDs comparable
  across timepoints.
- Output paths are generated from `contrasts.csv`, never typed by hand.

### Bugs from the earlier pipeline that must not return
- Copy-pasted output folder names, which let one comparison overwrite another.
- `--readLength` not matching the real reads. rMATS silently drops every read whose length
  differs, and also drops soft-clipped reads unless `--allow-clipping` is set (it always drops
  multi-mapped reads, NH != 1).
- Strandedness hardcoded to `fr-unstranded`.
- Animals sampled at several timepoints analysed as unpaired groups.
- Single-end and paired-end data aligned with different STAR settings.
- A separate rMATS run per comparison, so event IDs didn't match across comparisons.
- Three filter scripts with three different thresholds, and no read-coverage floor.

---
#  Features

### End-to-end RNA-seq splicing pipeline
1. Download raw FASTQ data  
2. Perform trimming + QC  
3. Build STAR index  
4. Align reads (Pair-Ended (PE) or Single-Ended (SE))  
5. Run rMATS for differential splicing  
6. Post-process events (filter, annotate, rank)  
7. Map events to transcript-level consequences  

### HPC-friendly (SLURM)
- All steps come with ready-to-run **SLURM scripts**  
- Modular design → you can run steps independently  
- Works with any genome FASTA + GTF  , just need to plug the reference genome file and the appropriate path

### downstream analysis tools
- Filter significant SE/RI/A5SS/A3SS/MXE events  
- Map Ensembl IDs → gene symbols  
- Rank genes by splicing burden  
- Summarize transcript-level impact  
 
---
# Pipeline Overview
```mermaid
flowchart LR
    FASTQ[FASTQ] --> QC["Trimming + QC (Trimmomatic, FastQC, MultiQC)"]
    QC --> STAR["STAR Alignment (produces sorted BAM files)"]
    STAR --> RMATS["rMATS Turbo (SE, RI, A3SS, A5SS, MXE)"]
    RMATS --> POST["Post-Splicing Analysis (filtering, gene mapping, event ranking)"]
    POST --> TX["Transcript-Level Mapping (isoform consequences)"]
```


---
# Quick start guide

### 1️⃣ Download FASTQ Files
Edit your SRA accession list inside:

```
scripts/download_sra.sh
```

Run:
```bash
sbatch scripts/download_sra.sh
```
**Outputs →** `data/raw/`

### 2️⃣ Trim Reads + Run QC
Run:
```bash
sbatch scripts/trim_qc.sh
```

This script performs:
- Trimmomatic trimming  
- FastQC on raw + trimmed reads  
- (Optional) MultiQC summary  

**Outputs →**
- `data/trimmed/`
- `data/qc/`

### 3️⃣ Determine Read Length
Accurate read length is required for STAR and rMATS.

Run:
```bash
sbatch scripts/readLength.sh
sbatch scripts/readLengthMode.sh
```
Use the **mode read length** reported.

### 4️⃣ Build STAR Genome Index
Place reference files in `reference/`:

```
reference/genome.fasta
reference/annotation.gtf
```

Then run:
```bash
sbatch scripts/build_star_index.sh
```

**Outputs →** `reference/star_index/`

### 5️⃣ Align Reads (PE or SE)

**Paired-end alignment:**
```bash
sbatch scripts/run_star_PE.sh
```

**Single-end alignment:**
```bash
sbatch scripts/align_star_SE.sh
```

**Outputs →** sorted BAM files in `data/bam/`

### 6️⃣ Define Experimental Groups for rMATS
Create two text files listing BAMs:

```
baseline.txt
peak.txt
```
(One BAM per line or comma-separated depending on script, I got stuck on checking for bugs for an hour+ because of this smh.)

### 7️⃣ Run rMATS
Run:
```bash
sbatch scripts/rMATS_run.sh
```

**Outputs →**
```
rmats_out/
    ├── SE.MATS.JCEC.txt
    ├── RI.MATS.JCEC.txt
    ├── A5SS.MATS.JCEC.txt
    ├── A3SS.MATS.JCEC.txt
    └── MXE.MATS.JCEC.txt
```

### 8️⃣ Post-Splicing Analysis (Python Tools)

#### 8.1 Filter Events
```bash
python post-splicing-analysis/filter_rmats_events.py \
  --input rmats_out/SE.MATS.JCEC.txt \
  --fdr 0.05 \
  --delta_psi 0.1 \
  --output SE.filtered.tsv
```
Applies:
- FDR threshold  
- |ΔPSI| threshold  
- Optional read filters  

#### 8.2 Map Gene Names from GTF
```bash
python post-splicing-analysis/gene_name_map.py \
  --input SE.filtered.tsv \
  --gtf reference/annotation.gtf \
  --output SE.filtered.named.tsv
```

Adds:
- `gene_symbol`
- `gene_description` (if available)

#### 8.3 Rank Genes by Number of Splicing Events
```bash
python post-splicing-analysis/rank_genes_by_event_count.py \
  --input SE.filtered.named.tsv \
  --output gene_splicing_rank.tsv
```

Outputs:
- Genes sorted by splicing burden  
- Useful for identifying key splicing-regulated genes  

### 9️⃣ Transcript-Level Impact Analysis
Scripts in `transcripts-mapping/` allow you to:

- Map exon-level events to transcript annotations  
- Identify isoform gain/loss  
- Detect events affecting:
  - CDS
  - 5′ UTR
  - 3′ UTR
  - Internal exons  
- Understand functional impact:
  - Isoform switching  
  - Coding sequence disruptions  
  - Alternative protein generation  

---
# Requirments

### Core Tools
- STAR ≥ 2.7  
- rMATS-turbo ≥ 4.1.2  
- Trimmomatic ≥ 0.39  
- FastQC ≥ 0.11.9  
- MultiQC  
- samtools  
- Python ≥ 3.8  
for Python dependencies nothing crazy just pandas, numpy etc..

---
# Contributions
To contribute, open an issue, submit a pull request, request a new feature ! 
Looking forward to you bioinformaticians feedback !

Cheers.



