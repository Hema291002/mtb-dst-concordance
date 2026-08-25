#!/usr/bin/env nextflow

/*
 * M. tuberculosis DST concordance pipeline
 * https://github.com/Hema291002/mtb-dst-concordance
 */

include { FASTP          } from './modules/local/fastp'
include { BWAMEM2_ALIGN  } from './modules/local/bwamem2_align'
include { SAMTOOLS_MARKDUP } from './modules/local/samtools_markdup'
include { BCFTOOLS_CALL     } from './modules/local/bcftools_call'

workflow {

    // Fail early and clearly rather than partway through a run.
    if (!params.input) { error "No samplesheet given. Use --input or -profile test" }
    if (!params.fasta) { error "No reference given. Use --fasta" }

    // Read the samplesheet into a channel of [meta, [read1, read2]].
    // checkIfExists makes a missing FASTQ an error at channel construction,
    // not a confusing failure inside a task twenty minutes later.
    ch_reads = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [ id: row.sample ]
            def r1 = file(row.fastq_1, checkIfExists: true)
            def r2 = file(row.fastq_2, checkIfExists: true)
            tuple(meta, [r1, r2])
        }

    FASTP(ch_reads)

    // The reference and its bwa-mem2 index are staged into every alignment
    // task. They are single values, not per-sample, so they go into value
    // channels: a value channel can be consumed repeatedly, whereas a queue
    // channel is drained after one use and the second sample would hang.
    ch_fasta = Channel.value(file(params.fasta, checkIfExists: true))
    ch_index = Channel.value(
        file("${params.fasta}.{0123,amb,ann,bwt.2bit.64,pac}", checkIfExists: true)
    )

    BWAMEM2_ALIGN(FASTP.out.reads, ch_fasta, ch_index)

    SAMTOOLS_MARKDUP(BWAMEM2_ALIGN.out.bam)

    ch_fai = Channel.value(file("${params.fasta}.fai", checkIfExists: true))

    BCFTOOLS_CALL(SAMTOOLS_MARKDUP.out.bam, ch_fasta, ch_fai)

    BCFTOOLS_CALL.out.vcf.view { meta, vcf, csi -> "called: ${meta.id}" }
}
