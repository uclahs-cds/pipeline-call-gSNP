process run_validate_PipeVal {
    container params.docker_image_validate

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process}/${task.process}-${task.index}/log${file(it).getName()}" }

    input:
    path file_to_validate

    output:
    path(".command.*")
    path("input_validation.txt"), emit: val_file

    """
    set -euo pipefail
    python -m validate -t file-input ${file_to_validate} > 'input_validation.txt'
    """
}

process calculate_sha512 {
    container params.docker_image_validate
    publishDir path: "${params.output_dir}/output",
      mode: "copy",
      pattern: "*.sha512"

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process}/${task.process}-${task.index}/log${file(it).getName()}" }

    input:
    path(file_for_calc)

    output:
    path(".command.*")
    path("*.sha512"), emit: sha512_sum

    script:
    """
    set -euo pipefail
    sha512sum ${file_for_calc} > ${file_for_calc}.sha512
    """
}
