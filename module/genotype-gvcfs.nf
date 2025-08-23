include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

/*
    Nextflow module for joint genotyping merged GVCFs with GATK
*/
process run_GenotypeGVCFs_GATK {
    container params.docker_image_gatk
    publishDir path: "${META.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
        mode: "copy",
        enabled: params.save_intermediate_files,
        pattern: '*.vcf*'

    publishDir path: "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/${task.process.split(':')[-1]}-${interval_id}/log${file(it).getName()}" }

    input:
    val(META)
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(dbsnp_bundle)
    path(dbsnp_bundle_index)
    tuple path(combined_gvcf), path(combined_gvcf_index), path(interval_path), val(interval_id)

    output:
    path(".command.*")
    tuple path(output_filename), path("${output_filename}.tbi"), emit: vcfs

    script:
    output_filename = generate_standard_filename(
        "GATK-${params.gatk_version}",
        params.dataset_id,
        params.patient_id,
        [
            'additional_information': "${interval_id}.vcf.gz"
        ]
    )
    interval_str = "--intervals ${interval_path}"
    interval_padding = params.is_targeted ? "--interval-padding 100" : ""
    """
    set -euo pipefail

    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m" \
        GenotypeGVCFs \
        --variant ${combined_gvcf} \
        --reference ${reference_fasta} \
        --verbosity INFO \
        --output ${output_filename} \
        --dbsnp ${dbsnp_bundle} \
        --standard-min-confidence-threshold-for-calling 50 \
        ${interval_str} \
        ${interval_padding}
    """
}
