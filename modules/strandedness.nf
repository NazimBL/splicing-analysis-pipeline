// One rMATS --libType for the whole dataset, inferred from STAR GeneCounts.
// Fails the run if samples disagree (bin/infer_strandedness.py).
process STRANDEDNESS {
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path gene_counts

    output:
    path 'strandedness.tsv'
    path 'lib_type.txt', emit: lib_type

    script:
    """
    infer_strandedness.py strandedness.tsv ${gene_counts} > lib_type.txt
    """
}
