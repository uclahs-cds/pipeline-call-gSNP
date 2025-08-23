include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'
/*
    Nextflow module for filtering GATK variant calls

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        sample_id: identifier for sample
        sample_vcf: path to VCF to filter
        sample_vcf_tbi: path to index of VCF to filter

    params:
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.docker_image_gatkfilter: string
        params.is_NT_paired: bool. Indicator of whether input has normal and tumour samples
*/
process filter_gSNP_GATK {
    container params.docker_image_gatkfilter
    publishDir path: "${META.output_dir_base}/output",
      mode: "copy",
      pattern: "${output_filename}*"

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
    val(META)
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    tuple val(sample_id), path(sample_vcf), path(sample_vcf_tbi)

    output:
    path(".command.*")
    tuple path("${output_filename}_snv.vcf.gz"),
          path("${output_filename}_snv.vcf.gz.tbi"),
          path("${output_filename}_indel.vcf.gz"),
          path("${output_filename}_indel.vcf.gz.tbi"), emit: germline_filtered
    tuple path("${output_filename}_variant-class-count.tsv"),
          path("${output_filename}_genotype-count.tsv"), emit: germline_filtered_tsv

    script:
    identifier_opts = (params.is_NT_paired)
        ? params.samples_to_process.collect{ "--${it.sample_type.replace('tumor', 'tumour')} ${it.id}" }.join(' ')
        : params.samples_to_process.collect{ "--tumour ${it.id} --normal ${it.id}" }.join(' ')

    output_filename = generate_standard_filename(
        "GATK-${params.gatk_version}",
        params.dataset_id,
        sample_id,
        [:]
    )
    """
    set -euo pipefail
    /src/NGS-Tools-GATK/bin/filter_GATK_SNV_calls.pl \
        --input ${sample_vcf} \
        --sample ${sample_id} \
        --ref ${reference_fasta} \
        ${identifier_opts} \
        --filter_somatic Y \
        --filter_ambiguous Y \
        --split_calls Y \
        --output_dir `pwd`

    mv filtered_germline_snv_${sample_id}_nosomatic.vcf.gz ${output_filename}_snv.vcf.gz
    mv filtered_germline_snv_${sample_id}_nosomatic.vcf.gz.tbi ${output_filename}_snv.vcf.gz.tbi
    mv filtered_germline_indel_${sample_id}_nosomatic.vcf.gz ${output_filename}_indel.vcf.gz
    mv filtered_germline_indel_${sample_id}_nosomatic.vcf.gz.tbi ${output_filename}_indel.vcf.gz.tbi
    mv filtered_germline_variant_class_count_${sample_id}.tsv ${output_filename}_variant-class-count.tsv
    mv filtered_germline_genotype_count_${sample_id}.tsv ${output_filename}_genotype-count.tsv
    """
}
