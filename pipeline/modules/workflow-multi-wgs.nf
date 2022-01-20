nextflow.enable.dsl=2

include { calculate_sha512 } from './validation.nf'
include {
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

workflow multi_sample_wgs {
    take:
    intervals
    split_intervals
    ir_input
    ir_input_no_interval
    identifiers
    identifier_sample

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
        identifier_sample,
        base_recal_intervals
        )

    remove_realigned_bams(
        recalibrate_base.out.bam_for_deletion.mix(recalibrate_base.out.bam_index_for_deletion),
        recalibrate_base.out.recalibrated_normal_bam.collect().mix(recalibrate_base.out.recalibrated_tumour_bam.collect()) // Let BQSR finish before deletion
        )

    reheader_interval_bams(
        recalibrate_base.out.recalibrated_normal_bam,
        recalibrate_base.out.recalibrated_normal_bam_index,
        recalibrate_base.out.recalibrated_tumour_bam,
        recalibrate_base.out.recalibrated_tumour_bam_index
        )

    // normal_bam_ch = reheader_interval_bams.out.reheadered_normal_bam
    // normal_bam_index_ch = reheader_interval_bams.out.reheadered_normal_bam_index
    // tumour_bam_ch = reheader_interval_bams.out.reheadered_tumour_bam
    // tumour_bam_index_ch = reheader_interval_bams.out.reheadered_tumour_bam_index
    // hc_interval = reheader_interval_bams.out.associated_interval

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
    reheader_interval_bams.out.reheadered_normal_bam
        .map{ it -> it[0]}
        .set{ normal_merge_bams_ich }
    
    reheader_interval_bams.out.reheadered_normal_bam
        .map{it -> it[-1]}
        .unique()
        .set{ normal_merge_id_ich }

    run_MergeSamFiles_Picard_normal(
        normal_merge_bams_ich.collect(),
        normal_merge_id_ich
        )

    // Prep the input for merging tumour BAMs
    reheader_interval_bams.out.reheadered_tumour_bam
        .map{ it ->
            tumour_it = []
            s_tum = it.size
            while (!(s_tum instanceof Integer)) {
                s_tum = s_tum.size
                }
            for(i_tum = 0; i_tum < s_tum; i_tum = i_tum + 1) {
                tumour_it = tumour_it + [it[i_tum][-1] + 'my_merging_separator' + it[i_tum][0]]
                }
            tumour_it
            }
        .flatten()
        .map{ it ->
            it.split('my_merging_separator')[0..-1]
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

    // merged_tumour_bam = run_MergeSamFiles_Picard_tumour.out.merged_bam
    // merged_tumour_bam_index = run_MergeSamFiles_Picard_tumour.out.merged_bam_index

    summary_intervals = split_intervals

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
    reheader_interval_bams.out.reheadered_normal_bam
        .map{ it ->
            [it[2].getFileName().toString(), (it[0] + 'my_vcfcalling_separator' + it[1] + 'my_vcfcalling_separator' + it[2]).toString()]
            }
        .set{ normal_for_join }

    reheader_interval_bams.out.reheadered_tumour_bam
        .map{ it ->
            tumour_vcf = []
            s_vcf = it.size
            while (!(s_vcf instanceof Integer)) {
                s_vcf = s_vcf.size
                }
            for(i_vcf = 0; i_vcf < s_vcf; i_vcf = i_vcf + 1) {
                tumour_vcf = tumour_vcf + [it[i_vcf][2].getFileName().toString() + 'my_vcfcalling_separator' + it[i_vcf][0] + 'my_vcfcalling_separator' + it[i_vcf][1] + 'my_vcfcalling_separator' + it[i_vcf][2]]
                }
            tumour_vcf
            }
        .flatten()
        .map{ it ->
            split_it = it.split('my_vcfcalling_separator')[0..-1]
            [split_it[0], split_it[1] + 'my_vcfcalling_separator' + split_it[2] + 'my_vcfcalling_separator' + split_it[3]]
            }
        .set{ tumour_for_join }

    normal_for_join
        .mix(tumour_for_join)
        .groupTuple(by: 0)
        .map{ it -> it[1].flatten()}
        .set{ grouped_bams }

    grouped_bams
        .map{ it ->
            mapped_bams = []
            s_map = it.size
            while (!(s_map instanceof Integer)) {
                s_map = s_map.size
                }
            for(i_map = 0; i_map < s_map; i_map = i_map + 1) {
                mapped_bams = mapped_bams + it[i_map].split('my_vcfcalling_separator')[0]
                }
            mapped_bams
            }
        .set{ vcf_caller_bams_ich }

    grouped_bams
        .map{ it ->
            mapped_bais = []
            s_map_bai = it.size
            while (!(s_map_bai instanceof Integer)) {
                s_map_bai = s_map_bai.size
                }
            for(i_map_bai = 0; i_map_bai < s_map_bai; i_map_bai = i_map_bai + 1) {
                mapped_bais = mapped_bais + it[i_map_bai].split('my_vcfcalling_separator')[1]
                }
            mapped_bais
            }
        .set{ vcf_caller_bais_ich }
    
    grouped_bams
        .map{ it ->
            it[0].split('my_vcfcalling_separator')[2]
            }
        .set{ vcf_caller_intervals_ich }
    
    identifier_sample.combine(vcf_caller_bams_ich)
        .map{it ->
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
        vcf_caller_bams_ich,
        vcf_caller_bais_ich,
        vcf_caller_intervals_ich
        )

    // Prep input for GVCF calling
    reheader_interval_bams.out.reheadered_normal_bam
        .map{ it ->
            it.join('my_gvcfcalling_separator')
            }
        .set{ gvcf_normal_intermediate }

    reheader_interval_bams.out.reheadered_tumour_bam
        .map{ it ->
            tumour_gvcf = []
            s_gvcf = it.size
            while (!(s_gvcf instanceof Integer)) {
                s_gvcf = s_gvcf.size
                }
            for(i_gvcf = 0; i_gvcf < s_gvcf; i_gvcf = i_gvcf + 1) {
                tumour_gvcf = tumour_gvcf + [it[i_gvcf].join('my_gvcfcalling_separator')]
                }
            tumour_gvcf
            }
        .flatten()
        .set{ gvcf_tumour_intermediate }

    gvcf_normal_intermediate
        .mix(gvcf_tumour_intermediate)
        .map{ it ->
            it.split('my_gvcfcalling_separator')[0..-1]
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