nextflow.enable.dsl=2

include {
    run_BaseRecalibrator_GATK
    run_ApplyBQSR_GATK
    } from './base-recalibration.nf'

workflow recalibrate_base {
    take:
    realigned_bam
    realigned_bam_index
    associated_interval
    includes_unmapped
    bqsr_generator_identifiers
    intervals

    main:
    run_BaseRecalibrator_GATK(
      params.reference_fasta,
      "${params.reference_fasta}.fai",
      "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
      params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz,
      "${params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz}.tbi",
      params.bundle_known_indels_vcf_gz,
      "${params.bundle_known_indels_vcf_gz}.tbi",
      params.bundle_v0_dbsnp138_vcf_gz,
      "${params.bundle_v0_dbsnp138_vcf_gz}.tbi",
      intervals,
      realigned_bam.collect(),
      realigned_bam_index.collect(),
      bqsr_generator_identifiers
      )

    // Extract the normal and tumour IDs
    bqsr_ids = bqsr_generator_identifiers
      .map{ it ->
        it[1]
        }
      .flatten()
      .unique()
      .map{ it ->
        [it, 'normal']
        }
      .mix(
        bqsr_generator_identifiers
          .map{ it ->
            it[2..-1]
            }
          .flatten()
          .unique()
          .filter{ it != 'NA' }
          .map{ it ->
            [it, 'tumour']
            }
        )

    // Generate the input channels
    counter_bam = 0
    bams_with_idx = realigned_bam
      .map{
        [counter_bam = counter_bam + 1, it]
        }

    counter_index = 0
    bam_index_with_idx = realigned_bam_index
      .map{
        [counter_index = counter_index + 1, it]
        }

    counter_interval = 0
    intervals_with_idx = associated_interval
      .map{
        [counter_interval = counter_interval + 1, it]
        }

    counter_unmapped = 0
    unmapped_with_idx = includes_unmapped
      .map{
        [counter_unmapped = counter_unmapped + 1, it]
        }

    apply_bqsr_ich = bams_with_idx
      .join(bam_index_with_idx, by: 0)
      .join(intervals_with_idx, by: 0)
      .join(unmapped_with_idx, by: 0)
      .map{ it ->
        it[1..-1]
        }
      .combine(bqsr_ids)

    run_ApplyBQSR_GATK(
      params.reference_fasta,
      "${params.reference_fasta}.fai",
      "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
      run_BaseRecalibrator_GATK.out.recalibration_table,
      apply_bqsr_ich
      )

    // Split the normal and tumour channels
    run_ApplyBQSR_GATK.out.apply_bqsr_och
      .branch{
        normal: it[1] == 'normal'
        tumour: it[1] == 'tumour'
        }
      .set{ apply_bqsr_split_och }

    apply_bqsr_split_och.normal
      .map{it ->
        [it[2].getFileName(), it]
        }.set{ normal_with_key }

    apply_bqsr_split_och.tumour
      .map{it ->
        [it[2].getFileName(), it]
        }.set{ tumour_with_key }
      
    // Split up the output channels
    if (params.is_NT_paired) {
      normal_with_key
        .join(tumour_with_key, by: 0) // Join on interval
        .view()
        .multiMap{it ->
          assoc_intervals_och: it[1][2]
          normal_bam_och: it[1][3]
          normal_bam_index_och: it[1][4]
          tumour_bam_och: it[2][3]
          tumour_bam_index_och: it[2][4]
          }.set{ bqsr_och }
    } else {
      normal_with_key
        .join(tumour_with_key, by: 0, remainder: true) // Join on interval
        .multiMap{it ->
          assoc_intervals_och: it[1][2]
          normal_bam_och: it[1][3]
          normal_bam_index_och: it[1][4]
          tumour_bam_och: "placeholder_bam"
          tumour_bam_index_och: "placeholder_index"
          }.set{ bqsr_och }
    }

    // Filter the deletion channel
    run_ApplyBQSR_GATK.out.deletion_och
      .filter{ it[0] == 'normal' }
      .multiMap{it ->
        bam_deletion_och: it[1]
        bam_index_deletion_och: it[2]
        }
      .set{ filtered_deletion_och }

    emit:
    recalibrated_normal_bam = bqsr_och.normal_bam_och
    recalibrated_normal_bam_index = bqsr_och.normal_bam_index_och
    recalibrated_tumour_bam = bqsr_och.tumour_bam_och
    recalibrated_tumour_bam_index = bqsr_och.tumour_bam_index_och
    associated_interval = bqsr_och.assoc_intervals_och
    bam_for_deletion = filtered_deletion_och.bam_deletion_och
    bam_index_for_deletion = filtered_deletion_och.bam_index_deletion_och
}
