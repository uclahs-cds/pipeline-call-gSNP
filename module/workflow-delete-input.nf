nextflow.enable.dsl=2

include { remove_intermediate_files } from '../external/nextflow-module/modules/common/intermediate_file_removal/main.nf' addParams(
    options: [
        save_intermediate_files: !params.metapipeline_delete_input_bams || params.save_intermediate_files,
        output_dir: params.output_dir,
        log_output_dir: params.log_input_deletion
        ]
    )

/*
    Nextflow process to check if given input can be deleted.

    input:
        file_to_check: path to file to be checked for deletion
*/
process check_deletion_status {
    container params.docker_image_validate
    containerOptions "--volume ${params.final_metapipeline_output_dir.split('/')[0..1].join('/')}:${params.final_metapipeline_output_dir.split('/')[0..1].join('/')}"

    publishDir path: "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}-${id}/log${file(it).getName()}" }

    input:
    path(file_to_check)

    output:
    path(file_to_check), emit: file_to_delete
    path(".command.*")

    when:
    params.metapipeline_delete_input_bams

    debug true

    script:
    """
    set -euo pipefail

    FILE_NAME=`basename ${file_to_check}`
    until ls ${params.final_metapipeline_output_dir}/\$FILE_NAME &> /dev/null
    do
        sleep 30
    done

    FINAL_PATH_TO_CHECK=`ls ${params.final_metapipeline_output_dir}/\$FILE_NAME`

    EXPECTED_SIZE=`stat --printf="%s" \$(readlink ${file_to_check})`
    COPIED_SIZE=`stat --printf="%s" \$FINAL_PATH_TO_CHECK`

    while [ \$EXPECTED_SIZE != \$COPIED_SIZE ]
    do
        sleep 30
        COPIED_SIZE=`stat --printf="%s" \$FINAL_PATH_TO_CHECK`
    done
    """
}

workflow delete_input {
    take:
    files_to_delete

    main:
    files_to_delete.view()
    check_deletion_status(files_to_delete)

    remove_intermediate_files(
        check_deletion_status.out.file_to_delete,
        "ready_to_delete"
        )
}
