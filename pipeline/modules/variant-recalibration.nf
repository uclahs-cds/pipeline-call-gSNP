process run_VariantRecalibratorINDEL_GATK {
    container params.docker_image_gatk
    publishDir path: params.output_dir,
      mode: "copy",
      pattern: "cohort_indel*"

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
    tuple(path(cohort_vcf), path(cohort_vcf_tbi))


    output:
    path(".command.*")
    path("cohort_indel.plots.R")
    path("cohort_indel.plots.R.pdf")
    tuple path(cohort_vcf),
          path(cohort_vcf_tbi),
          path("cohort_indel.recal"),
          path("cohort_indel.recal.idx"),
          path("cohort_indel.tranches"), emit: indel_recal

    script:
    """
    set -euo pipefail
    gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/tmp" \
           VariantRecalibrator \
           --variant ${cohort_vcf} \
           --reference ${reference_fasta} \
           --resource:mills,known=true,training=true,truth=true,prior=12.0 ${bundle_mills_and_1000g_gold_standards_vcf_gz} \
           --use-annotation DP \
           --use-annotation FS \
           --use-annotation QD \
           --use-annotation MQ \
           --use-annotation MQRankSum \
           --use-annotation ReadPosRankSum \
           --mode INDEL \
           --tranches-file cohort_indel.tranches \
           --output cohort_indel.recal \
           --truth-sensitivity-tranche 100.0 \
           --truth-sensitivity-tranche 99.9 \
           --truth-sensitivity-tranche 99.0 \
           --truth-sensitivity-tranche 90.0 \
           --max-gaussians 4 \
           --max-attempts 5 \
           --rscript-file cohort_indel.plots.R
    """
}

process run_VariantRecalibratorSNP_GATK {
    container params.docker_image_gatk
    publishDir path: params.output_dir,
      mode: "copy",
      pattern: "cohort_snp*"

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
    tuple(path(cohort_vcf), path(cohort_vcf_tbi))


    output:
    path(".command.*")
    path("cohort_snp.plots.R")
    path("cohort_snp.plots.R.pdf")
    path("cohort_snp.tranches.pdf")
    tuple path(cohort_vcf),
          path(cohort_vcf_tbi),
          path("cohort_snp.recal"),
          path("cohort_snp.recal.idx"),
          path("cohort_snp.tranches"), emit: snp_recal


    script:
    """
    set -euo pipefail
    gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/tmp" \
           VariantRecalibrator \
           --variant ${cohort_vcf} \
           --reference ${reference_fasta} \
           --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${bundle_hapmap_3p3_vcf_gz} \
           --resource:omni,known=false,training=true,truth=true,prior=12.0 ${bundle_omni_1000g_2p5_vcf_gz} \
           --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${bundle_phase1_1000g_snps_high_conf_vcf_gz} \
           --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${bundle_v0_dbsnp138_vcf_gz} \
           --use-annotation DP \
           --use-annotation FS \
           --use-annotation QD \
           --use-annotation MQ \
           --use-annotation MQRankSum \
           --use-annotation ReadPosRankSum \
           --mode SNP \
           --tranches-file cohort_snp.tranches \
           --output cohort_snp.recal \
           --truth-sensitivity-tranche 100.0 \
           --truth-sensitivity-tranche 99.9 \
           --truth-sensitivity-tranche 99.0 \
           --truth-sensitivity-tranche 90.0 \
           --max-gaussians 4 \
           --max-attempts 5 \
           --rscript-file cohort_snp.plots.R
    """
}

process run_ApplyVQSR_GATK {
    container params.docker_image_gatk
    publishDir path: params.output_dir,
      mode: "copy",
      pattern: "regenotyped_merged_recalibrated_${suffix}.vcf.gz{,.tbi}"

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
    tuple(path(cohort_vcf), path(cohort_vcf_tbi), path(recal_file), path(recal_index_file), path(tranches_file))


    output:
    path(".command.*")
    path("regenotyped_merged_recalibrated_${suffix}.vcf.gz{,.tbi}"), emit: vcf

    script:
    """
    set -euo pipefail
    gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/tmp" \
           ApplyVQSR \
           --variant ${cohort_vcf} \
           --reference ${reference_fasta} \
           --mode ${mode} \
           --truth-sensitivity-filter-level 99 \
           --tranches-file ${tranches_file} \
           --recal-file ${recal_file} \
           --output regenotyped_merged_recalibrated_${suffix}.vcf.gz
    """
}

process run_vcfstats_RTG {
    container params.docker_image_rtg
    publishDir path: params.output_dir,
      mode: "copy",
      pattern: "*.txt"

    publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "run_vcfstats_RTG/log${file(it).getName()}" }

    input:
    tuple(path(cohort_vcf), path(cohort_vcf_tbi))

    output:
    path(".command.*")
    path("*.txt")

    script:
    """
    set -euo pipefail
    rtg vcfstats \
      ${cohort_vcf} \
      > ${cohort_vcf.baseName.substring(0, cohort_vcf.baseName.indexOf('.'))}_vcfstats.txt
    """
}

workflow recalibrate_snps {
  take: cohort_vcf_channel
  main:
  run_VariantRecalibratorSNP_GATK(
      params.reference_fasta,
      "${params.reference_fasta}.fai",
      params.reference_dict, params.bundle_v0_dbsnp138_vcf_gz,
      "${params.bundle_v0_dbsnp138_vcf_gz}.tbi",
      params.bundle_hapmap_3p3_vcf_gz,
      "${params.bundle_hapmap_3p3_vcf_gz}.tbi",
      params.bundle_omni_1000g_2p5_vcf_gz,
      "${params.bundle_omni_1000g_2p5_vcf_gz}.tbi",
      params.bundle_phase1_1000g_snps_high_conf_vcf_gz,
      "${params.bundle_phase1_1000g_snps_high_conf_vcf_gz}.tbi",
      cohort_vcf_channel
  )

  run_ApplyVQSR_GATK(
      'SNP',
      'SNP',
      params.reference_fasta,
      "${params.reference_fasta}.fai",
      params.reference_dict,
      run_VariantRecalibratorSNP_GATK.out.snp_recal
  )

  run_vcfstats_RTG(
      run_ApplyVQSR_GATK.out.vcf
  )

  emit:
  vcf = run_ApplyVQSR_GATK.out.vcf
}

workflow recalibrate_indels {
  take: cohort_vcf_channel
  main:
  run_VariantRecalibratorINDEL_GATK(
      params.reference_fasta,
      "${params.reference_fasta}.fai",
      params.reference_dict,
      params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz,
      "${params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz}.tbi",
      cohort_vcf_channel
  )

  run_ApplyVQSR_GATK(
      'INDEL',
      'SNP_AND_INDEL',
      params.reference_fasta,
      "${params.reference_fasta}.fai",
      params.reference_dict,
      run_VariantRecalibratorINDEL_GATK.out.indel_recal
  )

  run_vcfstats_RTG(
      run_ApplyVQSR_GATK.out.vcf
  )

  emit:
  vcf = run_ApplyVQSR_GATK.out.vcf
}
