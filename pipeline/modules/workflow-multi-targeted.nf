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
include { flatten_samples } from './functions.nf'

workflow multi_sample_targeted {
    take:
    intervals
    split_intervals
    ir_input
    ir_input_no_interval
    identifiers
    identifier_sample

    main:
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
        recalibrate_base.out.recalibrated_normal_bam.collect().mix(recalibrate_base.out.recalibrated_tumour_bam.collect()).collect() // Let BQSR finish before deletion
        )

    reheader_interval_bams(
        recalibrate_base.out.recalibrated_normal_bam,
        recalibrate_base.out.recalibrated_normal_bam_index,
        recalibrate_base.out.recalibrated_tumour_bam,
        recalibrate_base.out.recalibrated_tumour_bam_index
        )

    recalibrated_bams_to_delete = reheader_interval_bams.out.normal_bam_for_deletion.mix(
        reheader_interval_bams.out.normal_bam_index_for_deletion,
        reheader_interval_bams.out.tumour_bam_for_deletion,
        reheader_interval_bams.out.tumour_bam_index_for_deletion
        )

    remove_recalibrated_bams(
        recalibrated_bams_to_delete,
        "mergesams_complete" // Decoy signal to let these files be deleted
        )

    // Prep the input for merging normal BAMs
    // Get only the BAMs
    reheader_interval_bams.out.reheadered_normal_bam
        .map{ it -> it[0]}
        .set{ normal_merge_bams_ich }
    
    // Get only the ids
    reheader_interval_bams.out.reheadered_normal_bam
        .map{it -> it[-1]}
        .unique()
        .set{ normal_merge_id_ich }

    run_MergeSamFiles_Picard_normal(
        normal_merge_bams_ich.collect(),
        normal_merge_id_ich
        )

    // Prep the input for merging tumour BAMs
    // Flatten the input by 1 level, keeping BAM-id associations
    // Extract the BAM and id from each
    flatten_samples(reheader_interval_bams.out.reheadered_tumour_bam)
    
    flatten_samples.out.och
        .map{ it ->
            it[-1,0]
            }
        .groupTuple()
        .multiMap{ it ->
            bams_ich: it[1]
            id_ich: it[0]
            }
        .set{ tumour_merge_ich }

    run_MergeSamFiles_Picard_tumour(
        tumour_merge_ich.bams_ich,
        tumour_merge_ich.id_ich
        )

    hc_interval = split_targeted_intervals
    summary_intervals = split_targeted_intervals

    calculate_contamination_normal(
        run_MergeSamFiles_Picard_normal.out.merged_bam,
        run_MergeSamFiles_Picard_normal.out.merged_bam_index,
        summary_intervals,
        run_MergeSamFiles_Picard_normal.out.associated_id
        )

    calculate_contamination_tumour(
        run_MergeSamFiles_Picard_tumour.out.merged_bam,
        run_MergeSamFiles_Picard_tumour.out.merged_bam_index,
        summary_intervals,
        run_MergeSamFiles_Picard_tumour.out.associated_id,
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
        run_MergeSamFiles_Picard_normal.out.associated_id
        )

    run_DepthOfCoverage_GATK_tumour(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        summary_intervals.collect(),
        run_MergeSamFiles_Picard_tumour.out.merged_bam,
        run_MergeSamFiles_Picard_tumour.out.merged_bam_index,
        "tumour",
        run_MergeSamFiles_Picard_tumour.out.associated_id
        )

    // Prep input for VCF calling
    // Replicate the merged normal + tumour BAMs and BAIs for every split interval
    run_MergeSamFiles_Picard_normal.out.merged_bam
        .concat(run_MergeSamFiles_Picard_tumour.out.merged_bam)
        .collect()
        .combine(hc_interval)
        .map{ it ->
            it[0..-2]
            }
        .set{ hc_vcf_bams_ich }

    run_MergeSamFiles_Picard_normal.out.merged_bam_index
        .concat(run_MergeSamFiles_Picard_tumour.out.merged_bam_index)
        .collect()
        .combine(hc_interval)
        .map{ it ->
            it[0..-2]
            }
        .set{ hc_vcf_bais_ich }

    // Replicate the sample identifier for every interval
    identifier_sample.combine(hc_vcf_bams_ich)
        .map{it ->
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

    // Prep input for GVCF calling
    // Add indices for joining inputs for GVCF calling
    hc_bam_counter = 0
    run_MergeSamFiles_Picard_normal.out.merged_bam
        .concat(run_MergeSamFiles_Picard_tumour.out.merged_bam)
        .map{ it ->
            [hc_bam_counter = hc_bam_counter + 1, it]
            }
        .set{ hc_gvcf_bams }

    hc_bai_counter = 0
    run_MergeSamFiles_Picard_normal.out.merged_bam_index
        .concat(run_MergeSamFiles_Picard_tumour.out.merged_bam_index)
        .map{ it ->
            [hc_bai_counter = hc_bai_counter + 1, it]
            }
        .set{ hc_gvcf_bais }
    
    hc_ids_counter = 0
    run_MergeSamFiles_Picard_normal.out.associated_id
        .concat(run_MergeSamFiles_Picard_tumour.out.associated_id)
        .map{ it ->
            [hc_ids_counter = hc_ids_counter + 1, it]
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

    hc_completion_signal = run_HaplotypeCallerVCF_GATK.out.vcfs.collect().mix(
        run_HaplotypeCallerGVCF_GATK.out.gvcfs.collect()
        )
        .collect()

    reheadered_bams_to_delete = reheader_interval_bams.out.reheadered_normal_bam.flatten()
        .mix(
            reheader_interval_bams.out.reheadered_tumour_bam.flatten()
            )
        .filter { it.toString().endsWith('.bam') || it.toString().endsWith('.bai') }
        .unique()

    reheadered_deletion_signal = run_MergeSamFiles_Picard_normal.out.merged_bam.mix(
        run_MergeSamFiles_Picard_tumour.out.merged_bam,
        hc_completion_signal
        )
        .collect()

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
            "--normal $it"
            }
        .mix(
            identifiers
                .map{ it ->
                    it[2..-1]
                    }
                .flatten()
                .unique()
                .filter{ it != 'NA' }
                .map{ it ->
                    "--tumour $it"
                    }
            )
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
        filter_gSNP_GATK.out.germline_filtered.flatten(),
        run_MergeSamFiles_Picard_normal.out.merged_bam.flatten(),
        run_MergeSamFiles_Picard_normal.out.merged_bam_index.flatten(),
        run_MergeSamFiles_Picard_tumour.out.merged_bam.flatten(),
        run_MergeSamFiles_Picard_tumour.out.merged_bam_index.flatten()
        )

    calculate_sha512(files_for_sha512)
}
