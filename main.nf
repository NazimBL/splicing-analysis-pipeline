#!/usr/bin/env nextflow
// Splicing analysis pipeline. Currently: FastQC (raw) -> fastp + STAR ->
// strandedness check -> MultiQC. rMATS steps follow.
//
//   nextflow run main.nf -params-file datasets/<virus>/params.yaml \
//       --raw_dir <fastq_dir> --star_index <index_dir> --outdir <outdir>

include { FASTQC       } from './modules/fastqc'
include { TRIM_ALIGN   } from './modules/trim_align'
include { STRANDEDNESS } from './modules/strandedness'
include { MULTIQC      } from './modules/multiqc'

workflow {
    ['dataset', 'samplesheet', 'adapter', 'raw_dir', 'star_index', 'outdir'].each { p ->
        if (!params[p]) error "Missing --${p} (see README, 'Running a dataset')"
    }

    // One item per sample: [meta, R1, R2]. meta is the samplesheet row without
    // the file names (sample, animal, species, timepoint, ...), used later to
    // group and pair samples for rMATS.
    samples = channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = row.findAll { k, _v -> !(k in ['fastq_1', 'fastq_2']) }
            tuple(meta,
                  file("${params.raw_dir}/${row.fastq_1}", checkIfExists: true),
                  file("${params.raw_dir}/${row.fastq_2}", checkIfExists: true))
        }

    FASTQC(samples)
    TRIM_ALIGN(samples, params.adapter, file(params.star_index, type: 'dir', checkIfExists: true))

    STRANDEDNESS(TRIM_ALIGN.out.gene_counts.collect())

    MULTIQC(
        FASTQC.out.reports
            .mix(TRIM_ALIGN.out.fastp_json, TRIM_ALIGN.out.star_log, TRIM_ALIGN.out.gene_counts)
            .collect()
    )
}
