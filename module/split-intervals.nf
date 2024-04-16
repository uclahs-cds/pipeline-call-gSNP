/*
    Nextflow module for splitting input intervals

    input:
        intervals: path to set of target intervals to split
        reference: path to reference genome fasta file
        reference_index: path to index for reference fasta
        reference_dict: path to dictionary for reference fasta

    params:
        params.is_targeted: bool.
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_gatk: string
*/
process run_SplitIntervals_GATK {
    container params.docker_image_gatk

    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "*-contig.interval_list",
               enabled: params.save_intermediate_files
    publishDir "${params.log_output_dir}/process-log",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
    path intervals
    path reference
    path reference_index
    path reference_dict

    output:
    path "*-contig.interval_list", emit: interval_list
    path ".command.*"

    script:
    subdivision_mode = params.is_targeted ? "--subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION" : ""
    """
    set -euo pipefail
    gatk SplitIntervals \
        -R ${reference} \
        -L "\$(realpath ${intervals})" \
        --scatter-count ${params.scatter_count} \
        ${subdivision_mode} \
        ${params.split_intervals_extra_args} \
        -O ./

    for interval in `ls *-scattered.interval_list`
    do
        mv \$interval `echo \$interval | cut -d '-' -f 1`-contig.interval_list
    done
    """
}
