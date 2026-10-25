process FASTQC {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(reads, stageAs: 'input/*')

    output:
    tuple val(meta), path("${meta.id}_fastqc.html"), emit: html
    tuple val(meta), path("${meta.id}_fastqc.zip"), emit: zip
    tuple val("${task.process}"), val('fastqc'), eval('fastqc --version | sed "s/FastQC v//"'), topic: versions

    script:
    // FastQC heap = threads x --memory; keep ~20% of the task memory for JVM overhead.
    def mem_mb = Math.max(100, Math.min(10000, (task.memory.toMega() * 0.8 / task.cpus).intValue()))
    """
    ln -s ${reads} ${meta.id}.fastq.gz
    fastqc --threads ${task.cpus} --memory ${mem_mb} ${task.ext.args ?: ''} ${meta.id}.fastq.gz
    """

    stub:
    """
    touch ${meta.id}_fastqc.html ${meta.id}_fastqc.zip
    """
}
