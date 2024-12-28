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
    publishDir path:
    publishDir path:

    input:
    val sample_id
    val sample_sex
    tuple path(recalibrated_vcf), path(recalibrated_vcf_tbi)

    output:

    script:
    """
    set -euo pipefail

    zgrep "##source=" ${recalibrated_vcf} > ./vcf_source.txt

    python ${script_dir}/filter_xy_call.py \
        --sample_name ${sample_id} \
        --input_vcf ${recalibrated_vcf} \
        --vcf_source ./vcf_source.txt \
        --sample_sex ${params.sample_sex} \
        --par_bed ${params.par_bed} \
        --output_dir .
    """
}