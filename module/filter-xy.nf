include { generate_standard_filename; sanitize_string } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

/*
    Nextflow module for filtering chrX and chrY variant calls based on sample sex

    input:
        sample_id: identifier for sample
        sample_vcf: path to VCF to filter
        sample_vcf_tbi: path to index of VCF to filter

    params:
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.docker_image_hail: string
        params.sample_sex: string
        params.par_bed: string(path)
*/

process filter_XY_Hail {
    container params.docker_image_hail

    publishDir path: "${params.output_dir_base}/output",
      mode: "copy",
      pattern: '*.vcf.bgz*'

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: {
        "${task.process.replace(':', '/')}-${sample_id}/log${file(it).getName()}"
        }

    input:
    tuple val(sample_id), path(recalibrated_vcf), path(recalibrated_vcf_tbi)
    path(reference_fasta)
    path(reference_index)
    path(par_bed)
    path(script_dir)

    output:
    path(".command.*")
    tuple path("${output_filename}.vcf.bgz"), path("${output_filename}.vcf.bgz.tbi"), emit: xy_filtered_vqsr

    script:
    output_filename = generate_standard_filename(
        "Hail-${params.hail_version}",
        params.dataset_id,
        sample_id,
        [additional_tools:["GATK-${params.gatk_version}"]]
        )
    """
    set -euo pipefail

    zgrep "##source=" ${recalibrated_vcf} > ./vcf_source.txt

    python ${script_dir}/filter_xy_call.py \
        --input_vcf ${recalibrated_vcf} \
        --vcf_source_file ./vcf_source.txt \
        --sample_sex ${params.sample_sex} \
        --par_bed ${par_bed} \
        --ref-genome ${reference_fasta} \
        --ref-index ${reference_index} \
        --output_name ${output_filename} \
        --output_dir .
    """
}
