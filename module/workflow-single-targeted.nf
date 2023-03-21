nextflow.enable.dsl=2

include { calculate_sha512 } from './validation.nf'
include {
    run_SplitIntervals_GATK as run_SplitIntervals_GATK_targeted
    run_HaplotypeCallerVCF_GATK
    run_HaplotypeCallerGVCF_GATK
    run_MergeVcfs_Picard as run_MergeVcfs_Picard_VCF
    run_MergeVcfs_Picard as run_MergeVcfs_Picard_GVCF
    } from './genotype-processes.nf'
include {
    recalibrate_snps
    recalibrate_indels
    filter_gSNP_GATK
    } from './variant-recalibration.nf'
include { realign_indels } from './indel-realignment.nf'
include { recalibrate_base } from './workflow-bqsr.nf'
include { merge_and_deduplicate as run_MergeSamFiles_Picard } from './workflow-merge-dedup.nf' addParams(
    log_output_dir_deletion: "${params.log_output_dir}/process-log/single_sample_targeted"
    )
include {
    calculate_contamination_normal
    run_DepthOfCoverage_GATK
    } from './summary-processes.nf'
include {
    remove_intermediate_files as remove_realigned_bams
    remove_intermediate_files as remove_recalibrated_bams
    remove_intermediate_files as remove_reheadered_bams
    } from '../external/nextflow-module/modules/common/intermediate_file_removal/main.nf' addParams(
        options: [
            save_intermediate_files: params.save_intermediate_files,
            output_dir: params.output_dir,
            log_output_dir: "${params.log_output_dir}/process-log/single_sample_targeted"
            ]
        )
include { delete_input } from './workflow-delete-input.nf'

workflow single_sample_targeted {
    take:
    intervals
    split_intervals
    ir_input
    ir_input_no_interval
    identifiers
    identifier_sample

    main:
    identifiers
        .map{ it ->
            it[1]
            }
        .flatten()
        .unique()
        .set{ normal_identifier }

    run_SplitIntervals_GATK_targeted(
        params.intervals,
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        "targeted-intervals",
        true
        )

    split_targeted_intervals = run_SplitIntervals_GATK_targeted.out.interval_list.flatten()

    realign_indels(
        ir_input,
        ir_input_no_interval
        )

    ir_input
        .map{ it -> [it[3], it[5]]}
        .flatten()
        .filter{ !it.endsWith('NO_FILE.bam') }
        .unique()
        .set{ input_files }

    realign_indels.out.includes_unmapped
        .collect()
        .set{ ir_complete_signal }

    input_files
        .combine(ir_complete_signal)
        .map{ it[0] }
        .set{ input_files_to_delete }

    delete_input(input_files_to_delete)

    base_recal_intervals = params.intervals

    recalibrate_base(
        realign_indels.out.realigned_bam,
        realign_indels.out.realigned_bam_index,
        realign_indels.out.associated_interval,
        realign_indels.out.includes_unmapped,
        identifiers,
        identifier_sample,
        base_recal_intervals
        )

    remove_realigned_bams(
        recalibrate_base.out.bam_for_deletion.mix(recalibrate_base.out.bam_index_for_deletion),
        "decoy signal"
        )

    // Extract the BAMs and indices separately
    normal_bam_ch = recalibrate_base.out.recalibrated_normal_bam.map{ it -> it[1] }
    normal_bam_index_ch = recalibrate_base.out.recalibrated_normal_bam_index.map{ it -> it[1] }
    hc_interval = recalibrate_base.out.associated_interval

    run_MergeSamFiles_Picard(
        normal_bam_ch.collect(),
        normal_identifier
        )

    hc_interval = split_targeted_intervals
    summary_intervals = split_targeted_intervals

    calculate_contamination_normal(
        run_MergeSamFiles_Picard.out.merged_bam,
        run_MergeSamFiles_Picard.out.merged_bam_index,
        summary_intervals,
        run_MergeSamFiles_Picard.out.associated_id
        )

    run_DepthOfCoverage_GATK(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        summary_intervals.collect(),
        run_MergeSamFiles_Picard.out.merged_bam,
        run_MergeSamFiles_Picard.out.merged_bam_index,
        "normal",
        run_MergeSamFiles_Picard.out.associated_id
        )

    // Prepare input for VCF calling
    // Replicate the merged normal BAM and BAI for every split interval
    run_MergeSamFiles_Picard.out.merged_bam
        .collect()
        .combine(hc_interval)
        .map{ it ->
            it[0]
            }
        .set{ hc_vcf_bams_ich }

    run_MergeSamFiles_Picard.out.merged_bam_index
        .collect()
        .combine(hc_interval)
        .map{ it ->
            it[0]
            }
        .set{ hc_vcf_bais_ich }

    // Replicate the sample identifier for every interval
    identifier_sample.combine(hc_vcf_bams_ich)
        .map{ it ->
            it[0]
            }
        .combine(hc_interval)
        .map{ it ->
            it[0]
            }
        .set{ hc_vcf_ids_ich }

    run_HaplotypeCallerVCF_GATK(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        params.bundle_v0_dbsnp138_vcf_gz,
        "${params.bundle_v0_dbsnp138_vcf_gz}.tbi",
        hc_vcf_ids_ich,
        hc_vcf_bams_ich,
        hc_vcf_bais_ich,
        hc_interval
        )

    // Prepare input for GVCF calling
    // Add indices for joining inputs for GVCF calling
    hc_bam_counter = 0
    run_MergeSamFiles_Picard.out.merged_bam
        .map{ it ->
            [hc_bam_counter = hc_bam_counter + 1, it]
            }
        .set{ hc_gvcf_bams }

    hc_bai_counter = 0
    run_MergeSamFiles_Picard.out.merged_bam_index
        .map{ it ->
            [hc_bai_counter = hc_bai_counter + 1, it]
            }
        .set{ hc_gvcf_bais }
    
    hc_interval_counter = 0
    run_MergeSamFiles_Picard.out.associated_id
        .map{ it ->
            [hc_interval_counter = hc_interval_counter + 1, it]
            }
        .set{ hc_gvcf_ids }

    // Join and remove join index
    // Replicate for each split interval
    hc_gvcf_bams
        .join(hc_gvcf_bais, by: 0)
        .join(hc_gvcf_ids, by: 0)
        .map{ it ->
            it[1..-1]
            }
        .combine(hc_interval)
        .map{ it ->
            [it[0], it[1], it[3], it[2]] // Re-order to match the input order
            }
        .set{ gvcf_caller_ich }

    run_HaplotypeCallerGVCF_GATK(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        params.bundle_v0_dbsnp138_vcf_gz,
        "${params.bundle_v0_dbsnp138_vcf_gz}.tbi",
        gvcf_caller_ich
        )

    reheadered_bams_to_delete = recalibrate_base.out.recalibrated_normal_bam.map{ it -> it[1] }.mix(
        recalibrate_base.out.recalibrated_normal_bam_index.map{ it -> it[1] }
        )

    reheadered_deletion_signal = run_MergeSamFiles_Picard.out.merging_complete_signal.last()

    remove_reheadered_bams(
        reheadered_bams_to_delete,
        reheadered_deletion_signal
        )

    // Prep input for merging VCFs
    // Keep only the id and VCF and group by id
    run_HaplotypeCallerVCF_GATK.out.vcfs
        .map{ it ->
            it[0,1]
            }
        .groupTuple(by: 0)
        .multiMap{ it ->
            vcfs: it[1].flatten()
            id: it[0]
            }
        .set{ merge_vcf_ich }
    
    run_MergeVcfs_Picard_VCF(
        merge_vcf_ich.vcfs,
        "VCF",
        merge_vcf_ich.id
        )

    // Prep input for merging GVCFs
    // Keep only the id and GVCF and group by id
    run_HaplotypeCallerGVCF_GATK.out.gvcfs
        .map{ it ->
            it[0,1]
            }
        .groupTuple(by: 0)
        .multiMap{ it ->
            gvcfs: it[1].flatten()
            id: it[0]
            }
        .set{ merge_gvcf_ich }

    run_MergeVcfs_Picard_GVCF(
        merge_gvcf_ich.gvcfs,
        "GVCF",
        merge_gvcf_ich.id
        )

    recalibrate_snps(
        run_MergeVcfs_Picard_VCF.out.associated_id,
        run_MergeVcfs_Picard_VCF.out.vcf,
        run_MergeVcfs_Picard_VCF.out.vcf_index
        )

    recalibrate_indels(
        recalibrate_snps.out.associated_id,
        recalibrate_snps.out.vcf,
        recalibrate_snps.out.vcf_index
        )

    // Prep identifiers for filtering script
    // For filtering, the normal and tumour identifiers are required
    // Generate the options to pass as inputs
    identifiers
        .map{ it ->
            it[1]
            }
        .flatten()
        .unique()
        .map{ it ->
            "--normal $it --tumour $it"
            }
        .set{ filter_identifiers_ich }

    filter_gSNP_GATK(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        recalibrate_indels.out.vcf,
        recalibrate_indels.out.vcf_index,
        recalibrate_indels.out.associated_id,
        filter_identifiers_ich.collect()
        )

    files_for_sha512 = run_MergeVcfs_Picard_GVCF.out.vcf.flatten().mix(
        run_MergeVcfs_Picard_GVCF.out.vcf_index.flatten(),
        recalibrate_indels.out.vcf.flatten(),
        recalibrate_indels.out.vcf_index.flatten(),
        filter_gSNP_GATK.out.germline_filtered.flatten(),
        run_MergeSamFiles_Picard.out.merged_bam.flatten(),
        run_MergeSamFiles_Picard.out.merged_bam_index.flatten()
        )

    calculate_sha512(files_for_sha512)
}