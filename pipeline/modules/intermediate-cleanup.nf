/*
    Nextflow module for removing intermediate files. Follows symlinks to remove original files.

    input:
        file_to_remove: path to file to be removed
        ready_for_deletion_signal: val to indicate that the file can be removed. For example, 
            if multiple processes need to use the file to be deleted, the signal allows for a
            way to wait until all of those processes have completed before deleting the file.

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
      saveAs: { (task.process.split(':').size() > 1)  ?
          "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" :
          "${task.process}/${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}"
          }

    input:
    path(file_to_remove)
    val(ready_for_deletion_signal)

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
