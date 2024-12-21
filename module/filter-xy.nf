include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

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

    output:

    script:
    """
    set -euo pipefail

    """
}