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

include { run_validate_PipeVal } from './modules/validation.nf'
include { run_SplitIntervals_GATK } from './modules/genotype-processes.nf'
include { extract_GenomeIntervals } from './modules/extract-intervals.nf'
include { single_sample_wgs } from './modules/workflow-single-wgs.nf'
include { single_sample_targeted } from './modules/workflow-single-targeted.nf'
include { paired_sample_wgs } from './modules/workflow-paired-wgs.nf'
include { paired_sample_targeted } from './modules/workflow-paired-targeted.nf'

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

workflow {
    run_validate_PipeVal(input_validation)
    // Collect and store input validation output
    run_validate_PipeVal.out.val_file.collectFile(
      name: 'input_validation.txt',
      storeDir: "${params.output_dir}/validation"
      )

    extract_GenomeIntervals(
      "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict"
      )
    
    intervals = extract_GenomeIntervals.out.genomic_intervals

    run_SplitIntervals_GATK(
      intervals,
      params.reference_fasta,
      "${params.reference_fasta}.fai",
      "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
      "genome-intervals",
      false
    )

    split_intervals = run_SplitIntervals_GATK.out.interval_list.flatten()

    ir_input = input_ch_input_csv.combine(split_intervals)
        .map{ input_csv,interval -> [input_csv.sample_id, input_csv.normal_id, input_csv.tumour_id, input_csv.normal_BAM, input_csv.normal_index, input_csv.tumour_BAM, input_csv.tumour_index, interval] }
    ir_input_no_interval = input_ch_input_csv.combine(split_intervals)
        .map{ input_csv,interval -> [input_csv.sample_id, input_csv.normal_id, input_csv.tumour_id, input_csv.normal_BAM, input_csv.normal_index, input_csv.tumour_BAM, input_csv.tumour_index] }

    if (params.is_NT_paired) {
        if (params.is_targeted) {
            paired_sample_targeted(
                intervals,
                split_intervals,
                ir_input,
                ir_input_no_interval,
                identifiers
                )
        } else {
            paired_sample_wgs(
                intervals,
                split_intervals,
                ir_input,
                ir_input_no_interval,
                identifiers
                )
        }
    } else {
        if (params.is_targeted) {
            single_sample_targeted(
                intervals,
                split_intervals,
                ir_input,
                ir_input_no_interval,
                identifiers
                )
        } else {
            single_sample_wgs(
                intervals,
                split_intervals,
                ir_input,
                ir_input_no_interval,
                identifiers
                )
        }
    }
}
