// FastQC on the raw reads (QC only; nothing downstream depends on it).
process FASTQC {
    tag "${meta.sample}"
    publishDir "${params.outdir}/align/fastqc_raw", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    path "*_fastqc.{zip,html}", emit: reports

    script:
    """
    fastqc -q -t ${task.cpus} ${r1} ${r2}
    """
}
