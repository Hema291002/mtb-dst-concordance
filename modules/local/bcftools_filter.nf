process BCFTOOLS_FILTER {
    tag "$meta.id"

    conda "bioconda::bcftools=1.21"
    container "biocontainers/bcftools:1.21--h8b25389_0"

    input:
    tuple val(meta), path(vcf), path(csi)
    path  mask

    output:
    tuple val(meta), path("*.filt.vcf.gz"), path("*.filt.vcf.gz.csi"), emit: vcf
    path "versions.yml"                                              , emit: versions

    script:
    """
    # Thresholds are parameters, not literals, so their effect on the result
    # can be tested without editing code. Justification for each is in
    # docs/filter_rationale.md; they were chosen by inspecting the observed
    # distributions rather than imported from another study.
    #
    # The allele-fraction filter is the haploid-specific one. In a pure culture
    # of a one-chromosome organism, a real variant sits near 1.0. Intermediate
    # values mean mixed infection or mapping error, not heterozygosity.
    #
    # -T ^mask excludes repetitive regions where reads map ambiguously and
    # generate variants that reflect placement failure rather than sequence
    # difference.

    bcftools view \\
        -T ^${mask} \\
        -i 'QUAL>=${params.min_qual} && FORMAT/DP>=${params.min_dp} && MQ>=${params.min_mq} && (FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1])>=${params.min_af}' \\
        -Oz -o ${meta.id}.filt.vcf.gz \\
        ${vcf}

    bcftools index ${meta.id}.filt.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
