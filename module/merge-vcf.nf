include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

/*
    Nextflow module for merging input VCFs

    input:
        vcfs: VCFs to merge
        vcf_indices: VCF indices
        vcf_type: string to indicate whether calling in germline mode
        id: Identifier for merged VCF
        
    params:
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.docker_image_picard: string
        params.gatk_command_mem_diff: float(memory)
*/
process run_MergeVcfs_Picard {
    container params.docker_image_picard

    publishDir path: "${params.output_dir_base}/output",
      mode: "copy",
      pattern: "*.vcf*"

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}-${id}/log${file(it).getName()}" }

    input:
    tuple path(vcfs), path(vcf_indices), val(vcf_type), val(id)

    output:
    path(".command.*")
    tuple val(id), path("*.vcf{,.gz}"), path("*.vcf.{idx,gz.tbi}"), emit: merged_vcf

    script:
    all_vcfs = vcfs.collect{ "-INPUT '${it}'" }.join(' ')
    output_filename_base = generate_standard_filename("GATK-${params.gatk_version}", params.dataset_id, id, [:])
    output_filename = (vcf_type == "GVCF") ?
        "${output_filename_base}.g.vcf.gz" :
        "${output_filename_base}.vcf.gz"
    """
    set -euo pipefail

    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -Djava.io.tmpdir=${workDir} \
      -jar /usr/local/share/picard-slim-2.26.10-0/picard.jar MergeVcfs \
      ${all_vcfs} \
      -OUTPUT ${output_filename} \
      -VALIDATION_STRINGENCY LENIENT
    """
}
