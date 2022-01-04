/*
    Nextflow module for validating files

    input:
        file_to_validate: path to file to validate
        
    params:
        params.log_output_dir: string(path)
        params.docker_image_validate: string
*/
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

/*
  Nextflow module for calculating SHA512 checksum

  input:
      file_for_calc: path to file for whichc to calculate checksum

  params:
      params.output_dir: string(path)
      params.log_output_dir: string(path)
      params.docker_image_validate: string
*/
process calculate_sha512 {
    container params.docker_image_validate
    publishDir path: "${params.output_dir}/output",
      mode: "copy",
      pattern: "*.sha512",
      saveAs: { filename -> (filename.endsWith(".bai.sha512") && !filename.endsWith(".bam.bai.sha512")) ? "${file(file(filename).baseName).baseName}.bam.bai.sha512" : "${filename}"}

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

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
