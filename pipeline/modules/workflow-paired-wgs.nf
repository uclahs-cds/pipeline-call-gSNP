nextflow.enable.dsl=2

include { calculate_sha512 } from './validation.nf'
include {
    run_HaplotypeCallerVCF_GATK
    run_HaplotypeCallerGVCF_GATK as run_HaplotypeCallerGVCF_GATK_normal
    run_HaplotypeCallerGVCF_GATK as run_HaplotypeCallerGVCF_GATK_tumour
    run_MergeVcfs_Picard as run_MergeVcfs_Picard_VCF
    run_MergeVcfs_Picard as run_MergeVcfs_Picard_normal_GVCF
    run_MergeVcfs_Picard as run_MergeVcfs_Picard_tumour_GVCF
    } from './genotype-processes.nf'
include {
    recalibrate_snps
    recalibrate_indels
    filter_gSNP_GATK
    } from './variant-recalibration.nf'
include { realign_indels } from './indel-realignment.nf'
include { recalibrate_base } from './workflow-bqsr.nf'
include { reheader_interval_bams } from './workflow-reheader.nf'
include {
    run_MergeSamFiles_Picard as run_MergeSamFiles_Picard_normal
    run_MergeSamFiles_Picard as run_MergeSamFiles_Picard_tumour
    } from './bam-processing.nf'
include {
    calculate_contamination_normal
    calculate_contamination_tumour
    run_DepthOfCoverage_GATK as run_DepthOfCoverage_GATK_normal
    run_DepthOfCoverage_GATK as run_DepthOfCoverage_GATK_tumour
    } from './summary-processes.nf'
include {
    remove_intermediate_files as remove_realigned_bams
    remove_intermediate_files as remove_recalibrated_bams
    remove_intermediate_files as remove_reheadered_bams
    } from './intermediate-cleanup.nf'

workflow paired_sample_wgs {
    take:
    intervals
    split_intervals
    ir_input
    ir_input_no_interval
    identifiers

    main:
    realign_indels(
        ir_input,
        ir_input_no_interval
        )

    base_recal_intervals = intervals

    recalibrate_base(
        realign_indels.out.realigned_bam,
        realign_indels.out.realigned_bam_index,
        realign_indels.out.associated_interval,
        realign_indels.out.includes_unmapped,
        identifiers,
        base_recal_intervals
        )

    remove_realigned_bams(
        recalibrate_base.out.bam_for_deletion.mix(recalibrate_base.out.bam_index_for_deletion),
        recalibrate_base.out.recalibrated_normal_bam.collect().mix(recalibrate_base.out.recalibrated_tumour_bam.collect()) // Let BQSR finish before deletion
        )

    reheader_interval_bams(
        identifiers,
        recalibrate_base.out.recalibrated_normal_bam,
        recalibrate_base.out.recalibrated_normal_bam_index,
        recalibrate_base.out.recalibrated_tumour_bam,
        recalibrate_base.out.recalibrated_tumour_bam_index,
        recalibrate_base.out.associated_interval
        )

    normal_bam_ch = reheader_interval_bams.out.reheadered_normal_bam
    normal_bam_index_ch = reheader_interval_bams.out.reheadered_normal_bam_index
    tumour_bam_ch = reheader_interval_bams.out.reheadered_tumour_bam
    tumour_bam_index_ch = reheader_interval_bams.out.reheadered_tumour_bam_index
    hc_interval = reheader_interval_bams.out.associated_interval

    recalibrated_bams_to_delete = reheader_interval_bams.out.normal_bam_for_deletion.mix(
        reheader_interval_bams.out.normal_bam_index_for_deletion,
        reheader_interval_bams.out.tumour_bam_for_deletion,
        reheader_interval_bams.out.tumour_bam_index_for_deletion
        )

    remove_recalibrated_bams(
        recalibrated_bams_to_delete,
        "mergesams_complete" // Decoy signal to let these files be deleted
        )

    run_MergeSamFiles_Picard_normal(
        normal_bam_ch.collect(),
        "normal",
        identifiers
        )

    run_MergeSamFiles_Picard_tumour(
        tumour_bam_ch.collect(),
        "tumour",
        identifiers
        )

    merged_tumour_bam = run_MergeSamFiles_Picard_tumour.out.merged_bam
    merged_tumour_bam_index = run_MergeSamFiles_Picard_tumour.out.merged_bam_index

    summary_intervals = split_intervals

    calculate_contamination_normal(
        run_MergeSamFiles_Picard_normal.out.merged_bam,
        run_MergeSamFiles_Picard_normal.out.merged_bam_index,
        summary_intervals,
        identifiers
        )

    calculate_contamination_tumour(
        merged_tumour_bam,
        merged_tumour_bam_index,
        summary_intervals,
        identifiers,
        calculate_contamination_normal.out.pileupsummaries
        )

    run_DepthOfCoverage_GATK_normal(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        summary_intervals.collect(),
        run_MergeSamFiles_Picard_normal.out.merged_bam,
        run_MergeSamFiles_Picard_normal.out.merged_bam_index,
        "normal",
        identifiers
        )

    run_DepthOfCoverage_GATK_tumour(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        summary_intervals.collect(),
        merged_tumour_bam,
        merged_tumour_bam_index,
        "tumour",
        identifiers
        )

    run_HaplotypeCallerVCF_GATK(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        params.bundle_v0_dbsnp138_vcf_gz,
        "${params.bundle_v0_dbsnp138_vcf_gz}.tbi",
        identifiers,
        normal_bam_ch,
        normal_bam_index_ch,
        tumour_bam_ch,
        tumour_bam_index_ch,
        hc_interval
        )

    run_HaplotypeCallerGVCF_GATK_normal(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        params.bundle_v0_dbsnp138_vcf_gz,
        "${params.bundle_v0_dbsnp138_vcf_gz}.tbi",
        identifiers,
        normal_bam_ch,
        normal_bam_index_ch,
        hc_interval,
        "normal"
        )

    run_HaplotypeCallerGVCF_GATK_tumour(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        params.bundle_v0_dbsnp138_vcf_gz,
        "${params.bundle_v0_dbsnp138_vcf_gz}.tbi",
        identifiers,
        tumour_bam_ch,
        tumour_bam_index_ch,
        hc_interval,
        "tumour"
        )

    hc_completion_signal = run_HaplotypeCallerVCF_GATK.out.vcf.collect().mix(
        run_HaplotypeCallerGVCF_GATK_normal.out.gvcf.collect(),
        run_HaplotypeCallerGVCF_GATK_tumour.out.gvcf.collect()
        )
        .collect()

    reheadered_bams_to_delete = reheader_interval_bams.out.reheadered_normal_bam.mix(
        reheader_interval_bams.out.reheadered_normal_bam_index,
        reheader_interval_bams.out.reheadered_tumour_bam,
        reheader_interval_bams.out.reheadered_tumour_bam_index
        )

    reheadered_deletion_signal = run_MergeSamFiles_Picard_normal.out.merged_bam.mix(
        merged_tumour_bam,
        hc_completion_signal
        )
        .collect()

    remove_reheadered_bams(
        reheadered_bams_to_delete,
        reheadered_deletion_signal
        )
    
    run_MergeVcfs_Picard_VCF(
        run_HaplotypeCallerVCF_GATK.out.vcf.collect(),
        "VCF",
        "-",
        identifiers
        )

    run_MergeVcfs_Picard_normal_GVCF(
        run_HaplotypeCallerGVCF_GATK_normal.out.gvcf.collect(),
        "GVCF",
        "normal",
        identifiers
        )

    run_MergeVcfs_Picard_tumour_GVCF(
        run_HaplotypeCallerGVCF_GATK_tumour.out.gvcf.collect(),
        "GVCF",
        "tumour",
        identifiers
        )

    recalibrate_snps(
        identifiers,
        run_MergeVcfs_Picard_VCF.out.vcf,
        run_MergeVcfs_Picard_VCF.out.vcf_index
        )

    recalibrate_indels(
        identifiers,
        recalibrate_snps.out.vcf,
        recalibrate_snps.out.vcf_index
        )

    filter_gSNP_GATK(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        recalibrate_indels.out.vcf,
        recalibrate_indels.out.vcf_index,
        identifiers
        )

    files_for_sha512 = run_MergeVcfs_Picard_normal_GVCF.out.vcf.flatten().mix(
        run_MergeVcfs_Picard_normal_GVCF.out.vcf_index.flatten(),
        filter_gSNP_GATK.out.germline_filtered.flatten(),
        run_MergeSamFiles_Picard_normal.out.merged_bam.flatten(),
        run_MergeSamFiles_Picard_normal.out.merged_bam_index.flatten(),
        run_MergeVcfs_Picard_tumour_GVCF.out.vcf.flatten(),
        run_MergeVcfs_Picard_tumour_GVCF.out.vcf_index.flatten(),
        run_MergeSamFiles_Picard_tumour.out.merged_bam.flatten(),
        run_MergeSamFiles_Picard_tumour.out.merged_bam_index.flatten()
        )

    calculate_sha512(files_for_sha512)
}