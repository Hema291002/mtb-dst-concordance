process BCFTOOLS_CALL {
    tag "$meta.id"

    conda "bioconda::bcftools=1.21"
    container "biocontainers/bcftools:1.21--h8b25389_0"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.raw.vcf.gz"), path("*.raw.vcf.gz.csi"), emit: vcf
    path "versions.yml"                                            , emit: versions

    script:
    """
    # --ploidy 1 is not cosmetic. M. tuberculosis is haploid: one chromosome,
    # one true base per position. A diploid model reserves prior probability
    # for a heterozygous state that cannot exist, distorting the likelihood.
    # Enforcing haploidy also makes disagreement between reads meaningful:
    # sequencing error, mapping error, or genuinely mixed infection.
    #
    # No BQSR. It requires a trusted known-sites resource, and none exists for
    # M. tuberculosis. Applied without one it would treat every true variant
    # as error, depressing quality exactly where variation is.
    #
    # -q filters reads by mapping quality, the read-level counterpart to the
    # region-level repeat mask applied downstream.

    bcftools mpileup \\
        -f ${fasta} \\
        -a AD,DP,SP \\
        -q ${params.min_mapq} \\
        -Q 20 \\
        -d ${params.max_depth} \\
        -Ou ${bam} \\
      | bcftools call -mv --ploidy 1 -Oz -o ${meta.id}.raw.vcf.gz

    bcftools index ${meta.id}.raw.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
