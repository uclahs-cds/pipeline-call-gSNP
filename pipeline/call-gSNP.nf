#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """\
======================================
C A L L - G S N P  P I P E L I N E
======================================
Boutros Lab

Current Configuration:

      - pipeline:
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
          tool RTG: ${params.docker_image_rtg}
          tool Picard: ${params.docker_image_picard}

      Extra parameters:
          ${params}

------------------------------------
Starting workflow...
------------------------------------
        """

include { run_validate } from './modules/validation'
include { run_GenomicsDBImport_GATK; run_SplitIntervals_GATK; run_HaplotypeCaller_GATK; run_GenotypeGVCFs_GATK; run_SortVcf_GATK; run_MergeVcfs_Picard } from './modules/joint-genotype-processes'
include { recalibrate_snps; recalibrate_indels} from './modules/variant-recalibration'

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

merge_identifiers = input_ch_input_csv.map{it -> [it.sample_id, it.normal_id, it.tumour_id]}.collect()

workflow {
    run_validate(input_validation)

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

    hc_input = input_ch_input_csv.combine(split_intervals) // Cross the input files with all the chr list
        .map{ input_csv,interval -> [input_csv.sample_id, input_csv.normal_id, input_csv.tumour_id, input_csv.normal_BAM, input_csv.normal_index, input_csv.tumour_BAM, input_csv.tumour_index, interval] }

    run_HaplotypeCaller_GATK(
      params.reference_fasta,
      "${params.reference_fasta}.fai",
      params.reference_dict,
      params.bundle_v0_dbsnp138_vcf_gz,
      "${params.bundle_v0_dbsnp138_vcf_gz}.tbi",
      hc_input
      )

    run_MergeVcfs_Picard(
      run_HaplotypeCaller_GATK.out.vcf.collect(),
      run_HaplotypeCaller_GATK.out.gvcf_normal.collect(),
      run_HaplotypeCaller_GATK.out.gvcf_tumour.collect().ifEmpty("/scratch/placeholder.txt"),
      merge_identifiers
      )

}
