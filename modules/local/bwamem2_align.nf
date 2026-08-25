process BWAMEM2_ALIGN {
    tag "$meta.id"

    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.21"
    container "biocontainers/bwa-mem2:2.2.1--he513fc3_0"

    input:
    tuple val(meta), path(reads)
    path  fasta
    path  index          // the bwa-mem2 index files

    output:
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam
    path "versions.yml"                                            , emit: versions

    script:
    def rg = "@RG\\\\tID:${meta.id}\\\\tSM:${meta.id}\\\\tLB:${meta.id}\\\\tPL:ILLUMINA"
    """
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R "${rg}" \\
        ${fasta} \\
        ${reads[0]} ${reads[1]} \\
      | samtools sort -@ 2 -m 512M -o ${meta.id}.sorted.bam -

    samtools index ${meta.id}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | tail -1)
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
