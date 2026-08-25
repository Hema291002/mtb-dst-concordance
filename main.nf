#!/usr/bin/env nextflow

/*
 * M. tuberculosis DST concordance pipeline
 * https://github.com/Hema291002/mtb-dst-concordance
 */

include { FASTP } from './modules/local/fastp'

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

    FASTP.out.reads.view { meta, reads -> "trimmed: ${meta.id} -> ${reads*.name}" }
}

