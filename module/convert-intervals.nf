/*
    Nextflow module for converting IntervalList format to BED

    input:
        intervals: path to IntervalList format intervals
        interval_id: interval ID

    params:
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_gatk: string
*/
process convert_IntervalListToBed_GATK {
    container params.docker_image_gatk

    publishDir path: "${META.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "*-contig.bed",
               enabled: params.save_intermediate_files
    publishDir "${params.log_output_dir}/process-log",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${interval_id}/log${file(it).getName()}" }

    input:
    val(META)
    tuple val(interval_id), path(intervals)

    output:
    tuple val(interval_id), path(output_filename), emit: interval_bed
    path ".command.*"

    script:
    output_filename = "${file(intervals).baseName}.bed"
    """
    set -euo pipefail
    gatk IntervalListToBed \
        --INPUT ${intervals} \
        --OUTPUT ${output_filename}
    """
}
