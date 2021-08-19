process run_VariantRecalibratorINDEL_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir}/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: "*_output_indel.*"

    publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "run_VariantRecalibratorINDEL_GATK/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(bundle_mills_and_1000g_gold_standards_vcf_gz)
    path(bundle_mills_and_1000g_gold_standards_vcf_gz_tbi)
    tuple val(sample_id), val(normal_id), val(tumour_id)
    path(sample_vcf)
    path(sample_vcf_tbi)


    output:
    path(".command.*")
    path("${sample_id}_output_indel.plots.R")
    path("${sample_id}_output_indel.plots.R.pdf")
    tuple path(sample_vcf),
          path(sample_vcf_tbi),
          path("${sample_id}_output_indel.recal"),
          path("${sample_id}_output_indel.recal.idx"),
          path("${sample_id}_output_indel.tranches"),
          val(sample_id), emit: indel_recal

    script:
    variable_mode_options = params.is_NT_paired ? "--use-annotation MQRankSum --use-annotation ReadPosRankSum" : ""
    """
    set -euo pipefail

    gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
            VariantRecalibrator \
            --variant ${sample_vcf} \
            --reference ${reference_fasta} \
            --resource:mills,known=true,training=true,truth=true,prior=12.0 ${bundle_mills_and_1000g_gold_standards_vcf_gz} \
            --use-annotation DP \
            --use-annotation FS \
            ${variable_mode_options} \
            --mode INDEL \
            --tranches-file ${sample_id}_output_indel.tranches \
            --output ${sample_id}_output_indel.recal \
            --truth-sensitivity-tranche 100.0 \
            --truth-sensitivity-tranche 99.9 \
            --truth-sensitivity-tranche 99.0 \
            --truth-sensitivity-tranche 90.0 \
            --max-gaussians 2 \
            --max-attempts 5 \
            --rscript-file ${sample_id}_output_indel.plots.R
    """
}

process run_VariantRecalibratorSNP_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir}/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: "*_output_snp.*"

    publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "run_VariantRecalibratorSNP_GATK/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(bundle_v0_dbsnp138_vcf_gz)
    path(bundle_v0_dbsnp138_vcf_gz_tbi)
    path(bundle_hapmap_3p3_vcf_gz)
    path(bundle_hapmap_3p3_vcf_gz_tbi)
    path(bundle_omni_1000g_2p5_vcf_gz)
    path(bundle_omni_1000g_2p5_vcf_gz_tbi)
    path(bundle_phase1_1000g_snps_high_conf_vcf_gz)
    path(bundle_phase1_1000g_snps_high_conf_vcf_gz_tbi)
    tuple val(sample_id), val(normal_id), val(tumour_id)
    path(sample_vcf)
    path(sample_vcf_tbi)


    output:
    path(".command.*")
    path("${sample_id}_output_snp.plots.R")
    path("${sample_id}_output_snp.plots.R.pdf")
    tuple path(sample_vcf),
          path(sample_vcf_tbi),
          path("${sample_id}_output_snp.recal"),
          path("${sample_id}_output_snp.recal.idx"),
          path("${sample_id}_output_snp.tranches"),
          val(sample_id), emit: snp_recal

    script:
    variable_mode_options = params.is_NT_paired ? "--use-annotation MQRankSum --use-annotation ReadPosRankSum" : ""
    """
    set -euo pipefail

    gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
        VariantRecalibrator \
        --variant ${sample_vcf} \
        --reference ${reference_fasta} \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${bundle_hapmap_3p3_vcf_gz} \
        --resource:omni,known=false,training=true,truth=true,prior=12.0 ${bundle_omni_1000g_2p5_vcf_gz} \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${bundle_phase1_1000g_snps_high_conf_vcf_gz} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${bundle_v0_dbsnp138_vcf_gz} \
        --use-annotation DP \
        --use-annotation QD \
        --use-annotation FS \
        --use-annotation MQ \
        ${variable_mode_options} \
        --mode SNP \
        --tranches-file ${sample_id}_output_snp.tranches \
        --output ${sample_id}_output_snp.recal \
        --truth-sensitivity-tranche 100.0 \
        --truth-sensitivity-tranche 99.9 \
        --truth-sensitivity-tranche 99.0 \
        --truth-sensitivity-tranche 90.0 \
        --max-gaussians 4 \
        --max-attempts 5 \
        --rscript-file ${sample_id}_output_snp.plots.R
    """
}

process run_ApplyVQSR_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir}/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: "*_merged_recalibrated_${suffix}.vcf.gz{,.tbi}"

    publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "run_ApplyVQSR_GATK/log${file(it).getName()}" }

    input:
    val(mode)
    val(suffix)
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    tuple(path(sample_vcf), path(sample_vcf_tbi), path(recal_file), path(recal_index_file), path(tranches_file), val(sample_id))


    output:
    path(".command.*")
    path("${sample_id}_merged_recalibrated_${suffix}.vcf.gz"), emit: vcf
    path("${sample_id}_merged_recalibrated_${suffix}.vcf.gz.tbi"), emit: vcf_index

    script:
    """
    set -euo pipefail
    gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
           ApplyVQSR \
           --variant ${sample_vcf} \
           --reference ${reference_fasta} \
           --mode ${mode} \
           --ts-filter-level 99 \
           --tranches-file ${tranches_file} \
           --recal-file ${recal_file} \
           --output ${sample_id}_merged_recalibrated_${suffix}.vcf.gz
    """
}

process filter_gSNP_GATK {
    container params.docker_image_gatkfilter
    publishDir path: "${params.output_dir}/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: "filtered_germline_*"

    publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "filter_gSNP_GATK/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(sample_vcf)
    path(sample_vcf_tbi)
    tuple val(sample_id), val(normal_id), val(tumour_id)

    output:
    path(".command.*")
    tuple path("filtered_germline_snv_${sample_id}_nosomatic.vcf.gz"),
          path("filtered_germline_snv_${sample_id}_nosomatic.vcf.gz.tbi"),
          path("filtered_germline_indel_${sample_id}_nosomatic.vcf.gz"),
          path("filtered_germline_indel_${sample_id}_nosomatic.vcf.gz.tbi"), emit: germline_filtered
    tuple path("filtered_germline_variant_class_count_${sample_id}.tsv"),
          path("filtered_germline_genotype_count_${sample_id}.tsv"), emit: germline_filtered_tsv

    script:
    tumour_option = params.is_NT_paired ? "--tumour ${tumour_id}" : "--tumour ${normal_id}"
    """
    set -euo pipefail
    /src/NGS-Tools-GATK/bin/filter_GATK_SNV_calls.pl \
        --input ${sample_vcf} \
        --sample ${sample_id} \
        --ref ${reference_fasta} \
        --normal ${normal_id} \
        ${tumour_option} \
        --filter_somatic Y \
        --filter_ambiguous Y \
        --split_calls Y \
        --output_dir `pwd`
    """
}

workflow recalibrate_snps {
  take:
  merge_identifiers
  sample_vcf
  sample_vcf_tbi

  main:
  run_VariantRecalibratorSNP_GATK(
      params.reference_fasta,
      "${params.reference_fasta}.fai",
      params.reference_dict,
      params.bundle_v0_dbsnp138_vcf_gz,
      "${params.bundle_v0_dbsnp138_vcf_gz}.tbi",
      params.bundle_hapmap_3p3_vcf_gz,
      "${params.bundle_hapmap_3p3_vcf_gz}.tbi",
      params.bundle_omni_1000g_2p5_vcf_gz,
      "${params.bundle_omni_1000g_2p5_vcf_gz}.tbi",
      params.bundle_phase1_1000g_snps_high_conf_vcf_gz,
      "${params.bundle_phase1_1000g_snps_high_conf_vcf_gz}.tbi",
      merge_identifiers,
      sample_vcf,
      sample_vcf_tbi
  )

  run_ApplyVQSR_GATK(
      'SNP',
      'SNP',
      params.reference_fasta,
      "${params.reference_fasta}.fai",
      params.reference_dict,
      run_VariantRecalibratorSNP_GATK.out.snp_recal
  )

  emit:
  vcf = run_ApplyVQSR_GATK.out.vcf
  vcf_index = run_ApplyVQSR_GATK.out.vcf_index
}

workflow recalibrate_indels {
  take:
  merge_identifiers
  sample_vcf
  sample_vcf_tbi

  main:
  run_VariantRecalibratorINDEL_GATK(
      params.reference_fasta,
      "${params.reference_fasta}.fai",
      params.reference_dict,
      params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz,
      "${params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz}.tbi",
      merge_identifiers,
      sample_vcf,
      sample_vcf_tbi
  )

  run_ApplyVQSR_GATK(
      'INDEL',
      'SNP_AND_INDEL',
      params.reference_fasta,
      "${params.reference_fasta}.fai",
      params.reference_dict,
      run_VariantRecalibratorINDEL_GATK.out.indel_recal
  )

  emit:
  vcf = run_ApplyVQSR_GATK.out.vcf
  vcf_index = run_ApplyVQSR_GATK.out.vcf_index
}
