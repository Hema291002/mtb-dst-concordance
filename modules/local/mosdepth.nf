process MOSDEPTH {
    tag "$meta.id"

    conda "bioconda::mosdepth=0.3.10"
    container "biocontainers/mosdepth:0.3.10--h4e814b3_0"

    input:
    tuple val(meta), path(bam), path(bai)
    path  targets

    output:
    tuple val(meta), path("*.regions.bed.gz")    , emit: regions
    tuple val(meta), path("*.thresholds.bed.gz") , emit: thresholds
    tuple val(meta), path("*.mosdepth.summary.txt"), emit: summary
    path "versions.yml"                          , emit: versions

    script:
    """
    # Coverage over the resistance genes is what makes a wild-type call
    # evidence rather than absence of data. If katG is not covered, "no
    # mutation in katG" means nothing.
    #
    # --thresholds reports, per region, how many bases were covered at each
    # depth. That catches a zero-coverage gap inside an otherwise
    # well-covered gene, which a mean would hide. This is how the 3,275 bp
    # ethA deletion in ERR8975061 was found: no variant caller emits a record
    # for absent sequence.
    #
    # -n skips the per-base output, which is large and not needed here.

    mosdepth \\
        -t ${task.cpus} \\
        -n \\
        --by ${targets} \\
        --thresholds 1,10,30,50 \\
        ${meta.id} \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version | sed 's/mosdepth //')
    END_VERSIONS
    """
}
