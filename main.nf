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
          input_csv: ${(params.containsKey("input")) ? "YAML input used" : params.input_csv}
          aligner: ${params.aligner}
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

include { run_validate_PipeVal } from './external/nextflow-module/modules/PipeVal/validate/main.nf' addParams(
    options: [
        docker_image_version: params.pipeval_version,
        main_process: "./" //Save logs in <log_dir>/process-log/run_validate_PipeVal
        ]
    )
include { run_SplitIntervals_GATK } from './module/genotype-processes.nf'
include { extract_GenomeIntervals } from './external/nextflow-module/modules/common/extract_genome_intervals/main.nf' addParams(
    options: [
        save_intermediate_files: params.save_intermediate_files
        ]
    )
include { single_sample_wgs } from './module/workflow-single-wgs.nf'
include { single_sample_targeted } from './module/workflow-single-targeted.nf'
include { multi_sample_wgs } from './module/workflow-multi-wgs.nf'
include { multi_sample_targeted } from './module/workflow-multi-targeted.nf'

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

// Inputs to channel
Channel
    .from( params.input.BAM.tumour )
    .map{
        [normal_index: indexFile(it.normal_BAM)] + [tumour_index: indexFile(it.tumour_BAM)] + it
    }
    .set { input_ch_input_csv }

// Set validation channel
if (params.is_NT_paired) {
    input_ch_input_csv.flatMap{it -> [it.normal_BAM, it.normal_index, it.tumour_BAM, it.tumour_index]}.unique().map{it -> ['file-input', it]}.set{input_validation}
} else {
    input_ch_input_csv.flatMap{it -> [it.normal_BAM, it.normal_index]}.unique().map{it -> ['file-input', it]}.set{input_validation}
}

// Gather the inputs into a single emission
input_ch_input_csv.multiMap{ it ->
  sample_id: it.sample_id
  normal_id: it.normal_id
  normal_BAM: it.normal_BAM
  normal_index: it.normal_index
  tumour_id: it.tumour_id
  tumour_BAM: it.tumour_BAM
  tumour_index: it.tumour_index
  }
  .set{ branched_input }

// For the sample ID and the normal input fields, keep only the first since they should be identical for the sample
branched_input.sample_id.first().collect().map{ it -> [sample_id: it] }.set{formatted_sample_id}
branched_input.normal_id.first().collect().map{ it -> [normal_id: it] }.set{formatted_normal_id}
branched_input.normal_BAM.first().collect().map{ it -> [normal_BAM: it] }.set{formatted_normal_BAM}
branched_input.normal_index.first().collect().map{ it -> [normal_index: it] }.set{formatted_normal_index}

// For the tumour input fields, keep all lines from csv
branched_input.tumour_id.collect().map{ it -> [tumour_id: it] }.set{formatted_tumour_id}
branched_input.tumour_BAM.collect().map{ it -> [tumour_BAM: it] }.set{formatted_tumour_BAM}
branched_input.tumour_index.collect().map{ it -> [tumour_index: it] }.set{formatted_tumour_index}

// Mix the formatted channels and gather into single emission
formatted_sample_id.mix(
  formatted_normal_id,
  formatted_normal_BAM,
  formatted_normal_index,
  formatted_tumour_id,
  formatted_tumour_BAM,
  formatted_tumour_index
  )
  .reduce{ a, b -> a + b }
  .set{ input_csv_formatted_ich }

// Create the identifier channels
identifiers = input_csv_formatted_ich.map{it -> it.sample_id + it.normal_id + it.tumour_id}.collect()
identifier_sample = input_csv_formatted_ich.map{it -> it.sample_id}.flatten()

workflow {
    run_validate_PipeVal(input_validation)
    // Collect and store input validation output
    run_validate_PipeVal.out.validation_result.collectFile(
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

    ir_input = input_csv_formatted_ich.combine(split_intervals)
        .map{ input_csv,interval -> [input_csv.sample_id, input_csv.normal_id, input_csv.tumour_id, input_csv.normal_BAM, input_csv.normal_index, input_csv.tumour_BAM, input_csv.tumour_index, interval] }
    ir_input_no_interval = input_csv_formatted_ich.combine(split_intervals)
        .map{ input_csv,interval -> [input_csv.sample_id, input_csv.normal_id, input_csv.tumour_id, input_csv.normal_BAM, input_csv.normal_index, input_csv.tumour_BAM, input_csv.tumour_index] }

    if (params.is_NT_paired) {
        if (params.is_targeted) {
            multi_sample_targeted(
                intervals,
                split_intervals,
                ir_input,
                ir_input_no_interval,
                identifiers,
                identifier_sample
                )
        } else {
            multi_sample_wgs(
                intervals,
                split_intervals,
                ir_input,
                ir_input_no_interval,
                identifiers,
                identifier_sample
                )
        }
    } else {
        if (params.is_targeted) {
            single_sample_targeted(
                intervals,
                split_intervals,
                ir_input,
                ir_input_no_interval,
                identifiers,
                identifier_sample
                )
        } else {
            single_sample_wgs(
                intervals,
                split_intervals,
                ir_input,
                ir_input_no_interval,
                identifiers,
                identifier_sample
                )
        }
    }
}
