process MULTIQC {
    label 'process_single'

    input:
    path qc_files, stageAs: '?/*'
    val title

    output:
    path '*multiqc_report.html', emit: report
    path '*_data', emit: data
    tuple val("${task.process}"), val('multiqc'), eval('multiqc --version | sed "s/.*version //"'), topic: versions

    script:
    """
    multiqc --force --title "${title}" --filename multiqc_report.html ${task.ext.args ?: ''} .
    """

    stub:
    """
    mkdir multiqc_report_data
    touch multiqc_report.html multiqc_report_data/multiqc.log
    """
}
