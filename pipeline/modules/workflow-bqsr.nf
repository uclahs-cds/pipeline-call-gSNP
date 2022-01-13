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
    // run_BaseRecalibrator_GATK(
    //   params.reference_fasta,
    //   "${params.reference_fasta}.fai",
    //   "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
    //   params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz,
    //   "${params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz}.tbi",
    //   params.bundle_known_indels_vcf_gz,
    //   "${params.bundle_known_indels_vcf_gz}.tbi",
    //   params.bundle_v0_dbsnp138_vcf_gz,
    //   "${params.bundle_v0_dbsnp138_vcf_gz}.tbi",
    //   intervals,
    //   realigned_bam.collect(),
    //   realigned_bam_index.collect(),
    //   bqsr_generator_identifiers
    //   )

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
      "/hot/users/yashpatel/pipeline-call-gSNP/work/recal_table_multi.grp",
      apply_bqsr_ich
      )

    // run_ApplyBQSR_GATK(
    //   params.reference_fasta,
    //   "${params.reference_fasta}.fai",
    //   "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
    //   run_BaseRecalibrator_GATK.out.recalibration_table,
    //   apply_bqsr_ich
    //   )

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
        }
        .groupTuple()
        .set{ tumour_with_key }
      
    // Split up the output channels
    if (params.is_NT_paired) {
      normal_with_key
        .join(tumour_with_key, by: 0) // Join on interval
        .multiMap{it ->
          assoc_intervals_och: it[1][2]
          normal_bam_och: it[1][0,3,2]
          normal_bam_index_och: it[1][0,4]
          tumour_bam_raw: it[2]
          tumour_bam_index_raw: it[2]
          }.set{ bqsr_och }

          bqsr_och.tumour_bam_raw.map{ it ->
            mapped_it_bam = []
            s_bam = it.size
            while (!(s_bam instanceof Integer)) {
              s_bam = s_bam.size
              }
            for(i_bam = 0; i_bam < s_bam; i_bam = i_bam + 1) {
              mapped_it_bam = mapped_it_bam + [it[i_bam][0,3,2]]
              }
            mapped_it_bam
            }.set{ tumour_bam_och }

          bqsr_och.tumour_bam_index_raw.map{ it ->
            mapped_it_bai = []
            s_bai = it.size
            while (!(s_bai instanceof Integer)) {
              s_bai = s_bai.size
              }
            // print(s)
            // print(it)
            for(i_bai = 0; i_bai < s_bai; i_bai = i_bai + 1) {
              // print('after')
              // print(s)
              // print(it[i])
              // print('end')
              mapped_it_bai = mapped_it_bai + [it[i_bai][0,4]]
              }
            mapped_it_bai
            }.set{ tumour_bam_index_och }
    } else {
      normal_with_key
        .join(tumour_with_key, by: 0, remainder: true) // Join on interval
        .multiMap{it ->
          assoc_intervals_och: it[1][2]
          normal_bam_och: it[1][0,3,2]
          normal_bam_index_och: it[1][0,4]
          tumour_bam_raw: "placeholder_bam"
          tumour_bam_index_raw: "placeholder_index"
          }.set{ bqsr_och }

          tumour_bam_och = bqsr_och.tumour_bam_raw
          tumour_bam_index_och = bqsr_och.tumour_bam_index_raw
    }

    bqsr_och.normal_bam_och.view{"nchanbam: $it"}
    bqsr_och.normal_bam_index_och.view{"nchanbai: $it"}
    tumour_bam_och.view{"tchanbam: $it"}
    tumour_bam_index_och.view{"tchanbai: $it"}

    // // Filter the deletion channel
    // run_ApplyBQSR_GATK.out.deletion_och
    //   .filter{ it[0] == 'normal' }
    //   .multiMap{it ->
    //     bam_deletion_och: it[1]
    //     bam_index_deletion_och: it[2]
    //     }
    //   .set{ filtered_deletion_och }

    // emit:
    // recalibrated_normal_bam = bqsr_och.normal_bam_och
    // recalibrated_normal_bam_index = bqsr_och.normal_bam_index_och
    // recalibrated_tumour_bam = tumour_bam_och
    // recalibrated_tumour_bam_index = tumour_bam_index_och
    // associated_interval = bqsr_och.assoc_intervals_och
    // bam_for_deletion = filtered_deletion_och.bam_deletion_och
    // bam_index_for_deletion = filtered_deletion_och.bam_index_deletion_och
}
