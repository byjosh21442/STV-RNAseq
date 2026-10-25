process SALMON_QUANT {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${meta.id}"), emit: results
    tuple val(meta), path("${meta.id}/lib_format_counts.json"), emit: lib_format
    tuple val("${task.process}"), val('salmon'), eval('salmon --version | sed "s/salmon //"'), topic: versions

    script:
    if (!meta.single_end) {
        error("SALMON_QUANT: sample ${meta.id} is not single-end; this workflow only supports single-end reads")
    }
    // -l must precede -r (Salmon docs); libtype comes from the samplesheet strandedness or --salmon_libtype.
    """
    salmon quant \\
        --threads ${task.cpus} \\
        --index ${index} \\
        --libType ${meta.libtype} \\
        -r ${reads} \\
        ${task.ext.args ?: ''} \\
        -o ${meta.id}
    """

    stub:
    """
    mkdir ${meta.id}
    printf 'Name\\tLength\\tEffectiveLength\\tTPM\\tNumReads\\n' > ${meta.id}/quant.sf
    printf '{"expected_format": "${meta.libtype}"}\\n' > ${meta.id}/lib_format_counts.json
    """
}
