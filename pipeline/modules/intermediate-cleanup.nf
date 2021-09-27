/*
    Removes intermediate files as they're no longer needed.
    Follows symlinks to remove the original file and also the symlink.
*/
process remove_intermediate_files {
    container params.docker_image_validate
    publishDir path: "${params.log_output_dir}/process-log"
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path(file_to_remove)

    output:
    path(".command.*")

    when:
    !params.save_intermediate_files

    script:
    """
    if [[ -L "${file_to_remove}" ]]
    then
        rm "$(readlink -f ${file_to_remove})"
    fi

    rm ${file_to_remove}
    """
}