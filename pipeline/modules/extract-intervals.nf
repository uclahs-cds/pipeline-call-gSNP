/*
    Nextflow module for extracting genome intervals from reference dictionary

    input:
        reference_dict: path to .dict associated with reference genome

    params:
        params.docker_image_validate: string
        params.output_dir: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
*/
process extract_GenomeIntervals {
    container params.docker_image_validate

    publishDir path: "${params.output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "genomic_intervals.list",
               enabled: params.save_intermediate_files
    publishDir path: "${params.log_output_dir}/process-log",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
    path reference_dict

    output:
    path "genomic_intervals.list", emit: genomic_intervals
    path ".command.*"

    script:
    """
    set -euo pipefail
    grep -e "^@SQ" ${reference_dict} | \
    cut -f 2 | \
    sed -e 's/SN://g' \
    > genomic_intervals.list
    """
}