process SALMON_INDEX {
    label 'process_high'

    input:
    path transcript_fasta
    path genome_fasta   // optional decoys: pass [] for a transcript-only index

    output:
    path 'salmon_index', emit: index
    path 'decoys.txt', emit: decoys, optional: true
    path 'salmon_index.log', emit: log
    tuple val("${task.process}"), val('salmon'), eval('salmon --version | sed "s/salmon //"'), topic: versions

    script:
    // Full-genome decoys (Salmon docs: decoy-aware transcriptome; decoy records must come last in the FASTA).
    def decoy_cmds = genome_fasta
        ? """
          gunzip -cf ${genome_fasta} | grep '^>' | cut -d ' ' -f 1 | sed 's/^>//' > decoys.txt
          cat <(gunzip -cf ${transcript_fasta}) <(gunzip -cf ${genome_fasta}) > gentrome.fa
          """
        : """
          gunzip -cf ${transcript_fasta} > gentrome.fa
          """
    def decoy_arg = genome_fasta ? '--decoys decoys.txt' : ''
    """
    set -o pipefail
    ${decoy_cmds}
    salmon index --threads ${task.cpus} -t gentrome.fa ${decoy_arg} -i salmon_index ${task.ext.args ?: ''} 2>&1 | tee salmon_index.log
    rm -f gentrome.fa
    """

    stub:
    """
    mkdir salmon_index
    touch salmon_index/versionInfo.json salmon_index.log
    ${genome_fasta ? 'touch decoys.txt' : ''}
    """
}
