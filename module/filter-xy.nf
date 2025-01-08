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

process filter_XY {
    container params.docker_image_hail

    publishDir path: "${params.output_dir_base}/output",
      mode: "copy",
      pattern: '*.vcf.bgz*',
      saveAs: {
        "${output_filename}_${sanitize_string(file(it).getName().replace("${sample_id}_", ""))}"
        }

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: {
        "${task.process.replace(':', '/')}/${task.process.split(':')[-1]}-${sample_id}-${interval_id}/log${file(it).getName()}"
        }

    input:
    tuple val(sample_id), path(recalibrated_vcf), path(recalibrated_vcf_tbi)
    path(par_bed)
    path(script_dir)

    output:
    path(".command.*")
    path("${output_filename}_XY_filtered.vcf.bgz")
    path("${output_filename}_XY_filtered.vcf.bgz.tbi")

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
        --sample_name ${output_filename} \
        --input_vcf ${recalibrated_vcf} \
        --vcf_source_file ./vcf_source.txt \
        --sample_sex ${params.sample_sex} \
        --par_bed ${par_bed} \
        --genome_build ${params.genome_build} \
        --output_dir .
    """
}
