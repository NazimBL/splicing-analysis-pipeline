// One MultiQC report over FastQC (raw), fastp and STAR.
process MULTIQC {
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path reports

    output:
    path "${params.dataset}_multiqc.html"
    path "${params.dataset}_multiqc_data"

    script:
    """
    multiqc -f -n ${params.dataset}_multiqc .
    """
}
