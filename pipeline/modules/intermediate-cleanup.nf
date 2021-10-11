/*
    Nextflow module for removing intermediate files. Follows symlinks to remove original files.

    input:
        file_to_remove: path to file to be removed
        merge_sams_completion_signal: val to indicate that merge SAMs process has completed

    params:
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_validate: string
*/
process remove_intermediate_files {
    container params.docker_image_validate
    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process}/${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path(file_to_remove)
    val(merge_sams_completion_signal)

    output:
    path(".command.*")

    when:
    !params.save_intermediate_files

    script:
    """
    set -euo pipefail

    if [[ -L ${file_to_remove} ]]
    then
        rm `readlink -f ${file_to_remove}`
    fi
    
    rm ${file_to_remove}
    """
}
