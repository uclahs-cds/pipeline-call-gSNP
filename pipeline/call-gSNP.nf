#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """\
======================================
C A L L - G S N P  P I P E L I N E
======================================
Boutros Lab

Current Configuration:

      - pipeline:
          name: ${workflow.manifest.name}
          version: ${workflow.manifest.version}

      - input:
          input_csv: ${params.input_csv}
          intervals: ${params.intervals}
          bundle_v0_dbsnp138_vcf_gz: ${params.bundle_v0_dbsnp138_vcf_gz}
          bundle_mills_and_1000g_gold_standard_indels_vcf_gz: ${params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz}
          bundle_known_indels_vcf_gz: ${params.bundle_known_indels_vcf_gz}
          bundle_v0_dbsnp138_vcf_gz: ${params.bundle_v0_dbsnp138_vcf_gz}
          bundle_hapmap_3p3_vcf_gz: ${params.bundle_hapmap_3p3_vcf_gz}
          bundle_omni_1000g_2p5_vcf_gz: ${params.bundle_omni_1000g_2p5_vcf_gz}
          bundle_phase1_1000g_snps_high_conf_vcf_gz: ${params.bundle_phase1_1000g_snps_high_conf_vcf_gz}
          bundle_contest_hapmap_3p3_vcf_gz: ${params.bundle_contest_hapmap_3p3_vcf_gz}

       - output: 
          output: ${params.output_dir}
          temp: ${params.temp_dir}

      Tools Used:
          tool GATK: ${params.docker_image_gatk}
          tool validate: ${params.docker_image_validate}
          tool Picard: ${params.docker_image_picard}
          tool SAMtools: ${params.docker_image_samtools}
          tool GATK3: ${params.docker_image_gatk3}
          tool GATK Filter: ${params.docker_image_gatkfilter}

      Extra parameters:
          ${params}

------------------------------------
Starting workflow...
------------------------------------
        """

include { run_validate_PipeVal; calculate_sha512 } from './modules/validation.nf'
include { run_SplitIntervals_GATK; run_HaplotypeCallerVCF_GATK; run_HaplotypeCallerGVCF_GATK as run_HaplotypeCallerGVCF_GATK_normal; run_HaplotypeCallerGVCF_GATK as run_HaplotypeCallerGVCF_GATK_tumour; run_MergeVcfs_Picard as run_MergeVcfs_Picard_VCF; run_MergeVcfs_Picard as run_MergeVcfs_Picard_normal_GVCF; run_MergeVcfs_Picard as run_MergeVcfs_Picard_tumour_GVCF } from './modules/genotype-processes.nf'
include { recalibrate_snps; recalibrate_indels; filter_gSNP_GATK } from './modules/variant-recalibration.nf'
include { realign_indels } from './modules/indel-realignment.nf'
include { recalibrate_base } from './modules/base-recalibration.nf'
include { reheader_interval_bams; run_MergeSamFiles_Picard as run_MergeSamFiles_Picard_normal; run_MergeSamFiles_Picard as run_MergeSamFiles_Picard_tumour } from './modules/bam-processing.nf'
include { calculate_contamination_normal; calculate_contamination_tumour; run_DepthOfCoverage_GATK as run_DepthOfCoverage_GATK_normal; run_DepthOfCoverage_GATK as run_DepthOfCoverage_GATK_tumour } from './modules/summary-processes.nf'
include { remove_intermediate_files as remove_realigned_bams; remove_intermediate_files as remove_recalibrated_bams; remove_intermediate_files as remove_reheadered_bams } from './modules/intermediate-cleanup.nf'

// Returns the index file for the given bam or vcf
def indexFile(bam_or_vcf) {
  if(bam_or_vcf.endsWith('.bam')) {
    return "${bam_or_vcf}.bai"
  }
  else if(bam_or_vcf.endsWith('vcf.gz')) {
    return "${bam_or_vcf}.tbi"
  }
  else {
    throw new Exception("Index file for ${bam_or_vcf} file type not supported. Use .bam or .vcf.gz files.")
  }
}

if (params.is_NT_paired) {
    Channel
        .fromPath(params.input_csv, checkIfExists: true)
        .splitCsv(header:true)
        .map{
            [normal_index: indexFile(it.normal_BAM)] + [tumour_index: indexFile(it.tumour_BAM)] + it
            }
        .set { input_ch_input_csv }

    // Validation channel    
    input_ch_input_csv.flatMap{it -> [it.normal_BAM, it.normal_index, it.tumour_BAM, it.tumour_index]}.set{input_validation}

} else {
    Channel
        .fromPath(params.input_csv, checkIfExists: true)
        .splitCsv(header:true)
        .map{
            // Add filler values for tumour sample if in single sample mode
            [normal_index: indexFile(it.normal_BAM)] + [tumour_id: 'NA'] + [tumour_BAM: '/scratch/NA.bam'] + [tumour_index: '/scratch/NA.bam.bai'] + it
            }
        .set { input_ch_input_csv }

    // Validation channel    
    input_ch_input_csv.flatMap{it -> [it.normal_BAM, it.normal_index]}.set{input_validation}
}

identifiers = input_ch_input_csv.map{it -> [it.sample_id, it.normal_id, it.tumour_id]}.collect()
identifiers.set{ merge_identifiers }
identifiers.set{ recal_snp_identifiers }
identifiers.set{ recal_indels_identifiers }
identifiers.set{ filter_gSNP_identifiers }
identifiers.set{ bqsr_generator_identifiers }
identifiers.set{ hc_identifiers }
identifiers.set{ bam_reheadering_identifiers }
identifiers.set{ merge_bams_identifiers }
identifiers.set{ contamination_identifiers }
identifiers.set{ doc_identifiers }

workflow {
    run_validate_PipeVal(input_validation)
    // Collect and store input validation output
    run_validate_PipeVal.out.val_file.collectFile(
      name: 'input_validation.txt',
      storeDir: "${params.output_dir}/validation"
      )

    if (params.intervals) {
      intervals = params.intervals
    } else {
      intervals = "${projectDir}/config/hg38_decoy_chromosomes_canonical.list"
    }

    run_SplitIntervals_GATK(
      intervals,
      params.reference_fasta,
      "${params.reference_fasta}.fai",
      params.reference_dict
    )

    split_intervals = run_SplitIntervals_GATK.out.interval_list.flatten()


    if (params.is_targeted) {
      // Cross the input files with all the exome targets
      ir_input = input_ch_input_csv.combine(Channel.of(params.intervals))
          .map{ input_csv,interval -> [input_csv.sample_id, input_csv.normal_id, input_csv.tumour_id, input_csv.normal_BAM, input_csv.normal_index, input_csv.tumour_BAM, input_csv.tumour_index, interval] }
      ir_input_no_interval = input_ch_input_csv.combine(Channel.of(params.intervals))
          .map{ input_csv,interval -> [input_csv.sample_id, input_csv.normal_id, input_csv.tumour_id, input_csv.normal_BAM, input_csv.normal_index, input_csv.tumour_BAM, input_csv.tumour_index] }
    } else {
      // Cross the input files with all the chr list
      ir_input = input_ch_input_csv.combine(split_intervals)
          .map{ input_csv,interval -> [input_csv.sample_id, input_csv.normal_id, input_csv.tumour_id, input_csv.normal_BAM, input_csv.normal_index, input_csv.tumour_BAM, input_csv.tumour_index, interval] }
      ir_input_no_interval = input_ch_input_csv.combine(split_intervals)
          .map{ input_csv,interval -> [input_csv.sample_id, input_csv.normal_id, input_csv.tumour_id, input_csv.normal_BAM, input_csv.normal_index, input_csv.tumour_BAM, input_csv.tumour_index] }
    }

    realign_indels(
        ir_input,
        ir_input_no_interval
        )

    recalibrate_base(
      realign_indels.out.realigned_bam,
      realign_indels.out.realigned_bam_index,
      realign_indels.out.associated_interval,
      realign_indels.out.includes_unmapped,
      bqsr_generator_identifiers
      )

    remove_realigned_bams(
      recalibrate_base.out.bam_for_deletion.mix(recalibrate_base.out.bam_index_for_deletion),
      "mergesams_complete" // Decoy signal to let these files be deleted
      )

    // Reheader interval-level bams in NT paired mode
    if (params.is_NT_paired) {
      reheader_interval_bams(
        bam_reheadering_identifiers,
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
    } else {
      // Generate decoy tumour bam and index channels for single sample mode
      normal_bam_ch = recalibrate_base.out.recalibrated_normal_bam
      normal_bam_index_ch = recalibrate_base.out.recalibrated_normal_bam_index
      tumour_bam_ch = Channel.of(1..params.scatter_count).map{"/scratch/placeholder_${it}.txt"}
      tumour_bam_index_ch = Channel.of(1..params.scatter_count).map{"/scratch/placeholder_${it}_index.txt"}
      hc_interval = recalibrate_base.out.associated_interval
    }

    run_MergeSamFiles_Picard_normal(
      normal_bam_ch.collect(),
      "normal",
      merge_bams_identifiers
      )
    
    if (params.is_NT_paired) {
      run_MergeSamFiles_Picard_tumour(
        tumour_bam_ch.collect(),
        "tumour",
        merge_bams_identifiers
        )

      merged_tumour_bam = run_MergeSamFiles_Picard_tumour.out.merged_bam
      merged_tumour_bam_index = run_MergeSamFiles_Picard_tumour.out.merged_bam_index
    } else {
      merged_tumour_bam = Channel.of("/scratch/placeholder.txt")
      merged_tumour_bam_index = Channel.of("/scratch/placeholder_index.txt")
    }

    calculate_contamination_normal(
      run_MergeSamFiles_Picard_normal.out.merged_bam,
      run_MergeSamFiles_Picard_normal.out.merged_bam_index,
      split_intervals,
      contamination_identifiers
      )

    if (params.is_NT_paired) {
      calculate_contamination_tumour(
        merged_tumour_bam,
        merged_tumour_bam_index,
        split_intervals,
        contamination_identifiers,
        calculate_contamination_normal.out.pileupsummaries
        )
    }

    run_DepthOfCoverage_GATK_normal(
      params.reference_fasta,
      "${params.reference_fasta}.fai",
      params.reference_dict,
      split_intervals.collect(),
      run_MergeSamFiles_Picard_normal.out.merged_bam,
      run_MergeSamFiles_Picard_normal.out.merged_bam_index,
      "normal",
      doc_identifiers
      )
    
    if (params.is_NT_paired) {
      run_DepthOfCoverage_GATK_tumour(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        params.reference_dict,
        split_intervals.collect(),
        merged_tumour_bam,
        merged_tumour_bam_index,
        "tumour",
        doc_identifiers
        )
    }

    if (params.is_targeted) {
      normal_bam_ch = run_MergeSamFiles_Picard_normal.out.merged_bam
      normal_bam_index_ch = run_MergeSamFiles_Picard_normal.out.merged_bam_index
      tumour_bam_ch = merged_tumour_bam
      tumour_bam_index_ch = merged_tumour_bam_index
      hc_interval = split_intervals
    }

    run_HaplotypeCallerVCF_GATK(
      params.reference_fasta,
      "${params.reference_fasta}.fai",
      params.reference_dict,
      params.bundle_v0_dbsnp138_vcf_gz,
      "${params.bundle_v0_dbsnp138_vcf_gz}.tbi",
      hc_identifiers,
      normal_bam_ch,
      normal_bam_index_ch,
      tumour_bam_ch,
      tumour_bam_index_ch,
      hc_interval
      )

    run_HaplotypeCallerGVCF_GATK_normal(
      params.reference_fasta,
      "${params.reference_fasta}.fai",
      params.reference_dict,
      params.bundle_v0_dbsnp138_vcf_gz,
      "${params.bundle_v0_dbsnp138_vcf_gz}.tbi",
      hc_identifiers,
      normal_bam_ch,
      normal_bam_index_ch,
      hc_interval,
      "normal"
      )

    if (params.is_NT_paired) {
      run_HaplotypeCallerGVCF_GATK_tumour(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        params.reference_dict,
        params.bundle_v0_dbsnp138_vcf_gz,
        "${params.bundle_v0_dbsnp138_vcf_gz}.tbi",
        hc_identifiers,
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
    } else {
      hc_completion_signal = run_HaplotypeCallerVCF_GATK.out.vcf.collect().mix(
        run_HaplotypeCallerGVCF_GATK_normal.out.gvcf.collect()
        )
        .collect()
    }

    reheadered_bams_to_delete = run_HaplotypeCallerVCF_GATK.out.normal_bam_for_deletion.mix(
      run_HaplotypeCallerVCF_GATK.out.normal_bam_index_for_deletion,
      run_HaplotypeCallerVCF_GATK.out.tumour_bam_for_deletion,
      run_HaplotypeCallerVCF_GATK.out.tumour_bam_index_for_deletion
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
      merge_identifiers
      )

    run_MergeVcfs_Picard_normal_GVCF(
      run_HaplotypeCallerGVCF_GATK_normal.out.gvcf.collect(),
      "GVCF",
      "normal",
      merge_identifiers
      )

    if (params.is_NT_paired) {
      run_MergeVcfs_Picard_tumour_GVCF(
        run_HaplotypeCallerGVCF_GATK_tumour.out.gvcf.collect(),
        "GVCF",
        "tumour",
        merge_identifiers
        )
    }

    recalibrate_snps(
      recal_snp_identifiers,
      run_MergeVcfs_Picard_VCF.out.vcf,
      run_MergeVcfs_Picard_VCF.out.vcf_index
      )

    recalibrate_indels(
      recal_indels_identifiers,
      recalibrate_snps.out.vcf,
      recalibrate_snps.out.vcf_index
      )

    filter_gSNP_GATK(
      params.reference_fasta,
      "${params.reference_fasta}.fai",
      params.reference_dict,
      recalibrate_indels.out.vcf,
      recalibrate_indels.out.vcf_index,
      filter_gSNP_identifiers
      )

    files_for_sha512 = run_MergeVcfs_Picard_normal_GVCF.out.vcf.flatten().mix(
      run_MergeVcfs_Picard_normal_GVCF.out.vcf_index.flatten(),
      filter_gSNP_GATK.out.germline_filtered.flatten(),
      run_MergeSamFiles_Picard_normal.out.merged_bam.flatten(),
      run_MergeSamFiles_Picard_normal.out.merged_bam_index.flatten()
      )

    if (params.is_NT_paired) {
      files_for_sha512 = files_for_sha512.mix(
        run_MergeVcfs_Picard_tumour_GVCF.out.vcf.flatten(),
        run_MergeVcfs_Picard_tumour_GVCF.out.vcf_index.flatten(),
        run_MergeSamFiles_Picard_tumour.out.merged_bam.flatten(),
        run_MergeSamFiles_Picard_tumour.out.merged_bam_index.flatten()
        )
    }

    calculate_sha512(files_for_sha512)
}
