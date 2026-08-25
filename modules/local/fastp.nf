process FASTP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::fastp=1.3.6"
    container "biocontainers/fastp:1.3.6--heae3180_0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trim.fastq.gz"), emit: reads
    tuple val(meta), path("*.fastp.json")   , emit: json
    tuple val(meta), path("*.fastp.html")   , emit: html
    path "versions.yml"                     , emit: versions

    script:
    """
    fastp \\
        -i ${reads[0]} \\
        -I ${reads[1]} \\
        -o ${meta.id}_1.trim.fastq.gz \\
        -O ${meta.id}_2.trim.fastq.gz \\
        --detect_adapter_for_pe \\
        --cut_tail \\
        --cut_tail_window_size 4 \\
        --cut_tail_mean_quality ${params.cut_mean_q} \\
        --length_required ${params.min_length} \\
        --n_base_limit ${params.n_base_limit} \\
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e 's/fastp //g')
    END_VERSIONS
    """
}
