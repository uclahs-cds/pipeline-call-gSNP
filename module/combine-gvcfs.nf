include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

/*
    Nextflow module for merging GVCFs for joint genotyping with GATK
*/
process run_CombineGVCFs_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
        mode: "copy",
        enabled: params.save_intermediate_files,
        pattern: '*g.vcf.gz*'
    publishDir path: "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/${task.process.split(':')[-1]}-${interval_id}/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    tuple path(gvcfs), path(gvcf_indices), path(interval_path), val(interval_id)

    output:
    path(".command.*")
    tuple path(output_filename), path("${output_filename}.tbi"), path(interval_path), val(interval_id), emit: combined_gvcf

    script:
    output_filename = generate_standard_filename(
        "GATK-${params.gatk_version}",
        params.dataset_id,
        params.patient_id,
        [
            'additional_information': "${interval_id}.g.vcf.gz"
        ]
    )
    gvcf_input_str = gvcfs.collect{ "--variant '${it}'" }.join(' ')
    interval_str = "--intervals ${interval_path}"
    interval_padding = params.is_targeted ? "--interval-padding 100" : ""
    """
    set -euo pipefail

    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m" \
        CombineGVCFs \
        --reference ${reference_fasta} \
        ${gvcf_input_str} \
        --output ${output_filename} \
        --create-output-variant-index true \
        --verbosity INFO \
        ${interval_str} \
        ${interval_padding}
    """
}
