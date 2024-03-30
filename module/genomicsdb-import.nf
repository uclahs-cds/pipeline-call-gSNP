include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

/*
    Nextflow module for importing GVCFs into GenomicsDB for joint genotyping with GATK
*/
process run_GenomicsDBImport_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
        mode: "copy",
        enabled: params.save_intermediate_files,
        pattern: '*genomicsdb'

    publishDir path: "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/${task.process.split(':')[-1]}-${interval_id}/log${file(it).getName()}" }

    input:
    tuple path(gvcfs), path(gvcf_indices), path(interval_path), val(interval_id)

    output:
    path(".command.*")
    tuple path(output_filename), path(interval_path), val(interval_id), emit: genomicsdb

    script:
    output_filename = generate_standard_filename(
        "GATK-${params.gatk_version}",
        params.dataset_id,
        params.patient_id,
        [
            'additional_information': "${interval_id}.genomicsdb"
        ]
    )
    gvcf_input_str = gvcfs.collect{ "--variant '${it}'" }.join(' ')
    interval_str = "--intervals ${interval_path}"
    interval_padding = params.is_targeted ? "--interval-padding 100" : ""
    """
    set -euo pipefail

    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m" \
        GenomicsDBImport \
        ${gvcf_input_str} \
        --genomicsdb-workspace-path ${output_filename} \
        --verbosity INFO \
        ${interval_str} \
        ${interval_padding}
    """
}
