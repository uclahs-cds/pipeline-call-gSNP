include { generate_standard_filename } from '../external/nextflow-module/modules/common/generate_standardized_filename/main.nf'
/*
    Nextflow module for generating INDEL variant recalibration

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        bundle_mills_and_1000g_gold_standards_vcf_gz: path to standard Mills and 1000 genomes variants
        bundle_mills_and_1000g_gold_standards_vcf_gz_tbi: path to index file for Mills and 1000g variants
        (sample_id, normal_id, tumour_id): tuples of string identifiers for the samples
        sample_vcf: path to VCF to recalibrate
        sample_vcf_tbi: path to index of VCF to recalibrate
        
    params:
        params.output_dir: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_gatk: string
        params.is_NT_paired: bool. Indicator of whether input has normal and tumour samples
*/
process run_VariantRecalibratorINDEL_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir}/intermediate/${task.process.replace(':', '/')}",
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: "*_output-indel.{recal*,tranches}"

    publishDir path: "${params.output_dir}/QC/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: "*_output-indel.{plots*}"

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(bundle_mills_and_1000g_gold_standards_vcf_gz)
    path(bundle_mills_and_1000g_gold_standards_vcf_gz_tbi)
    val(sample_id)
    path(sample_vcf)
    path(sample_vcf_tbi)


    output:
    path(".command.*")
    path("${output_filename}_output-indel.plots.R")
    path("${output_filename}_output-indel.plots.R.pdf")
    tuple path(sample_vcf),
          path(sample_vcf_tbi),
          path("${output_filename}_output-indel.recal"),
          path("${output_filename}_output-indel.recal.idx"),
          path("${output_filename}_output-indel.tranches"),
          val(sample_id), emit: indel_recal

    script:
    variable_mode_options = params.is_NT_paired ? "--use-annotation MQRankSum --use-annotation ReadPosRankSum" : ""
    output_filename = generate_standard_filename(
        "GATK-${params.gatk_version}",
        params.dataset_id,
        sample_id,
        [:]
    )
    """
    set -euo pipefail

    gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir}" \
            VariantRecalibrator \
            --variant ${sample_vcf} \
            --reference ${reference_fasta} \
            --resource:mills,known=true,training=true,truth=true,prior=12.0 ${bundle_mills_and_1000g_gold_standards_vcf_gz} \
            --use-annotation DP \
            --use-annotation FS \
            ${variable_mode_options} \
            --mode INDEL \
            --tranches-file ${output_filename}_output-indel.tranches \
            --output ${output_filename}_output-indel.recal \
            --truth-sensitivity-tranche 100.0 \
            --truth-sensitivity-tranche 99.9 \
            --truth-sensitivity-tranche 99.0 \
            --truth-sensitivity-tranche 90.0 \
            --max-gaussians 2 \
            --max-attempts 5 \
            --rscript-file ${output_filename}_output-indel.plots.R
    """
}

/*
    Nextflow module for generating SNP variant recalibration

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        bundle_hapmap_3p3_vcf_gz: path to hapmap variants
        bundle_hapmap_3p3_vcf_gz_tbi: path to index of hapmap variants
        bundle_omni_1000g_2p5_vcf_gz: path to omni variants
        bundle_omni_1000g_2p5_vcf_gz_tbi: path to index of omni variants
        bundle_phase1_1000g_snps_high_conf_vcf_gz: path to high confidence SNPs
        bundle_phase1_1000g_snps_high_conf_vcf_gz_tbi: path to index of high confidence SNPs
        (sample_id, normal_id, tumour_id): tuples of string identifiers for the samples
        sample_vcf: path to VCF to recalibrate
        sample_vcf_tbi: path to index of VCF to recalibrate
        
    params:
        params.output_dir: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_gatk: string
        params.is_NT_paired: bool. Indicator of whether input has normal and tumour samples
*/
process run_VariantRecalibratorSNP_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir}/intermediate/${task.process.replace(':', '/')}",
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: "*_output-snp.{recal*,tranches}"

    publishDir path: "${params.output_dir}/QC/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: "*_output-snp.{plots*,tranches.pdf}"

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

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
    val(sample_id)
    path(sample_vcf)
    path(sample_vcf_tbi)


    output:
    path(".command.*")
    path("${output_filename}_output-snp.plots.R")
    path("${output_filename}_output-snp.plots.R.pdf")
    path("${output_filename}_output-snp.tranches.pdf")
    tuple path(sample_vcf),
          path(sample_vcf_tbi),
          path("${output_filename}_output-snp.recal"),
          path("${output_filename}_output-snp.recal.idx"),
          path("${output_filename}_output-snp.tranches"),
          val(sample_id), emit: snp_recal

    script:
    variable_mode_options = params.is_NT_paired ? "--use-annotation MQRankSum --use-annotation ReadPosRankSum" : ""
    output_filename = generate_standard_filename(
        "GATK-${params.gatk_version}",
        params.dataset_id,
        sample_id,
        [:]
    )
    """
    set -euo pipefail

    gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir}" \
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
        --tranches-file ${output_filename}_output-snp.tranches \
        --output ${output_filename}_output-snp.recal \
        --truth-sensitivity-tranche 100.0 \
        --truth-sensitivity-tranche 99.9 \
        --truth-sensitivity-tranche 99.0 \
        --truth-sensitivity-tranche 90.0 \
        --max-gaussians 4 \
        --max-attempts 5 \
        --rscript-file ${output_filename}_output-snp.plots.R
    """
}

/*
    Nextflow module for applying variant recalibration

    input:
        mode: string to indicate SNP or INDEL mode
        suffix: string to indicate which of SNP and INDELs are recalibrated
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        (sample_vcf, sample_vcf_tbi, recal_file, recal_index_file, tranches_file, sample_id):
          tuple of input VCF and index, recalibration files, and sample ID
        
    params:
        params.output_dir: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_gatk: string
*/
process run_ApplyVQSR_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir}/intermediate/${task.process.replace(':', '/')}/",
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: "*_merged_recalibrated_SNP.vcf.gz{,.tbi}"

    publishDir path: "${params.output_dir}/output/",
      mode: "copy",
      pattern: "*_merged_recalibrated_SNP_AND_INDEL.vcf.gz{,.tbi}"

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
    val(mode)
    val(suffix)
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    tuple(path(sample_vcf), path(sample_vcf_tbi), path(recal_file), path(recal_index_file), path(tranches_file), val(sample_id))


    output:
    path(".command.*")
    path(output_filename), emit: vcf
    path("${output_filename}.tbi"), emit: vcf_index
    val(sample_id), emit: associated_id

    script:
    output_filename = generate_standard_filename(
        "GATK-${params.gatk_version}",
        params.dataset_id,
        sample_id,
        [
            'additional_information': "merged_recalibrated_${suffix}.vcf.gz"
        ]
    )
    """
    set -euo pipefail
    gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir}" \
           ApplyVQSR \
           --variant ${sample_vcf} \
           --reference ${reference_fasta} \
           --mode ${mode} \
           --ts-filter-level 99 \
           --tranches-file ${tranches_file} \
           --recal-file ${recal_file} \
           --output ${output_filename}
    """
}

/*
    Nextflow module for filtering GATK variant calls

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        sample_vcf: path to VCF to filter
        sample_vcf_tbi: path to index of VCF to filter
        (sample_id, normal_id, tumour_id): tuples of string identifiers for the samples
        
    params:
        params.output_dir: string(path)
        params.log_output_dir: string(path)
        params.docker_image_gatkfilter: string
        params.is_NT_paired: bool. Indicator of whether input has normal and tumour samples
*/
process filter_gSNP_GATK {
    container params.docker_image_gatkfilter
    publishDir path: "${params.output_dir}/output",
      mode: "copy",
      pattern: "filtered_germline_*"

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(sample_vcf)
    path(sample_vcf_tbi)
    val(sample_id)
    val(identifiers_options)

    output:
    path(".command.*")
    tuple path("filtered_germline_snv_${sample_id}_nosomatic.vcf.gz"),
          path("filtered_germline_snv_${sample_id}_nosomatic.vcf.gz.tbi"),
          path("filtered_germline_indel_${sample_id}_nosomatic.vcf.gz"),
          path("filtered_germline_indel_${sample_id}_nosomatic.vcf.gz.tbi"), emit: germline_filtered
    tuple path("filtered_germline_variant_class_count_${sample_id}.tsv"),
          path("filtered_germline_genotype_count_${sample_id}.tsv"), emit: germline_filtered_tsv

    script:
    identifier_opts = identifiers_options.collect{ "$it" }.join(' ')
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
    """
}

/*
    Nextflow module for generating VCF stats

    input:
        (cohort_vcf, cohort_vcf_tbi): tuple of paths to VCF and index

    params:
        params.output_dir: string(path)
        params.log_output_dir: string(path)
        params.docker_image_rtg: string
*/
process run_vcfstats_RTG {
    container params.docker_image_rtg
    publishDir path: "${params.output_dir}/output",
      mode: "copy",
      pattern: "*.txt"

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

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
  take:
  merge_identifiers
  sample_vcf
  sample_vcf_tbi

  main:
  run_VariantRecalibratorSNP_GATK(
      params.reference_fasta,
      "${params.reference_fasta}.fai",
      "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
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
      "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
      run_VariantRecalibratorSNP_GATK.out.snp_recal
  )

  emit:
  vcf = run_ApplyVQSR_GATK.out.vcf
  vcf_index = run_ApplyVQSR_GATK.out.vcf_index
  associated_id = run_ApplyVQSR_GATK.out.associated_id
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
      "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
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
      "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
      run_VariantRecalibratorINDEL_GATK.out.indel_recal
  )

  emit:
  vcf = run_ApplyVQSR_GATK.out.vcf
  vcf_index = run_ApplyVQSR_GATK.out.vcf_index
  associated_id = run_ApplyVQSR_GATK.out.associated_id
}
