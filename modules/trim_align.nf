// fastp adapter trimming, then STAR 2-pass (bin/star_align.sh), in one task.
// The trimmed reads are as large as the raw ones (~235 GB for H5N1), so they
// stay inside the task directory and are deleted as soon as STAR is done,
// instead of piling up in work/ as a separate process's output.
process TRIM_ALIGN {
    tag "${meta.sample}"
    // fastp files -> align/fastp/, the sample's STAR folder -> align/star/<sample>/.
    // Files inside the STAR folder are also declared on their own (for later
    // steps); null skips publishing them a second time.
    publishDir "${params.outdir}/align", mode: 'copy', saveAs: { f ->
        f.contains('/')            ? null
        : f.endsWith('.fastp_time.txt') ? "fastp/${f.replace('.fastp_time.txt', '.time.txt')}"
        : f.contains('.fastp.')    ? "fastp/${f}"
        : "star/${f}"
    }

    input:
    tuple val(meta), path(r1), path(r2)
    val adapter
    path star_index

    output:
    tuple val(meta), path("${meta.sample}/${meta.sample}.Aligned.sortedByCoord.out.bam"),
          path("${meta.sample}/${meta.sample}.Aligned.sortedByCoord.out.bam.bai"), emit: bam
    path "${meta.sample}/${meta.sample}.ReadsPerGene.out.tab", emit: gene_counts
    path "${meta.sample}/${meta.sample}.Log.final.out",        emit: star_log
    path "${meta.sample}.fastp.json",                          emit: fastp_json
    path "${meta.sample}.fastp.{html,log}"
    path "${meta.sample}.fastp_time.txt"
    path "${meta.sample}"

    script:
    def s = meta.sample
    // fastp -w 8 as in the 2026-10-15 run; its output does not depend on threads.
    """
    /usr/bin/time -v -o ${s}.fastp_time.txt fastp -i ${r1} -I ${r2} \\
        -o ${s}_R1.trim.fq.gz -O ${s}_R2.trim.fq.gz \\
        --adapter_sequence ${adapter} --adapter_sequence_r2 ${adapter} \\
        --disable_quality_filtering --length_required 20 --trim_poly_g -w 8 \\
        -j ${s}.fastp.json -h ${s}.fastp.html 2> ${s}.fastp.log

    star_align.sh ${star_index} ${s} ${s}_R1.trim.fq.gz ${s}_R2.trim.fq.gz ${s} ${task.cpus}

    rm -f ${s}_R1.trim.fq.gz ${s}_R2.trim.fq.gz
    """
}
