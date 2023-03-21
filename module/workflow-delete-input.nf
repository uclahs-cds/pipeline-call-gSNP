nextflow.enable.dsl=2

include { remove_intermediate_files } from '../external/nextflow-module/modules/common/intermediate_file_removal/main.nf' addParams(
    options: [
        save_intermediate_files: !params.metapipeline_delete_input_bams || params.save_intermediate_files,
        output_dir: params.output_dir,
        log_output_dir: "${params.log_output_dir}/process-log/multi_sample_targeted"
        ]
    )

/*
    Nextflow process to check if given input can be deleted.

    input:
        file_to_check: path to file to be checked for deletion
*/
process check_deletion_status {
    container params.docker_image_validate

    input:
    path(file_to_check)

    output:
    path(file_to_check), emit: file_to_delete

    when:
    params.metapipeline_delete_input_bams

    script:
    """
    FILE_NAME=`basename ${file_to_check}`
    FINAL_PATH_TO_CHECK=`ls ${params.final_metapipeline_output_dir}/output/align-DNA-*/*/BWA-MEM2-*/output/\$FILE_NAME`

    EXPECTED_SIZE=`stat --printf="%s" ${file_to_check}`
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
    check_deletion_status(files_to_delete)

    remove_intermediate_files(
        check_deletion_status.out.file_to_delete,
        "ready_to_delete"
        )
}
