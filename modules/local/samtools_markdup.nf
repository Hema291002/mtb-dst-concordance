process SAMTOOLS_MARKDUP {
    tag "$meta.id"

    conda "bioconda::samtools=1.21"
    container "biocontainers/samtools:1.21--h50ea8bc_0"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.markdup.bam"), path("*.markdup.bam.bai"), emit: bam
    tuple val(meta), path("*.markdup.stats")                         , emit: stats
    path "versions.yml"                                              , emit: versions

    script:
    """
    # Four stages, and each is required:
    #   collate  regroups mates adjacently; the BAM is coordinate-sorted,
    #            which separates mates that mapped far apart
    #   fixmate  -m adds mate score tags. WITHOUT -m, markdup cannot choose a
    #            representative and duplicate marking silently does nothing.
    #            This is the commonest mistake in this workflow.
    #   sort     restores coordinate order, which markdup requires
    #   markdup  flags duplicates. Marked, never removed: flagging is
    #            reversible and lets downstream tools decide for themselves.
    #
    # -u between stages keeps the intermediate uncompressed, avoiding a
    # pointless compress/decompress cycle through each pipe.

    samtools collate -@ ${task.cpus} -O -u ${bam} tmp_collate_${meta.id} \\
      | samtools fixmate -@ ${task.cpus} -m -u - - \\
      | samtools sort -@ ${task.cpus} -m 512M -u - \\
      | samtools markdup -@ ${task.cpus} \\
            -f ${meta.id}.markdup.stats \\
            - ${meta.id}.markdup.bam

    samtools index ${meta.id}.markdup.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
