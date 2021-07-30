#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """\
======================================
R E G E N O T Y P E - G S N P  P I P E L I N E
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

       - options:
          option scatter_count: ${params.scatter_count}
          option split_intervals_extra_args = ${params.split_intervals_extra_args}

      Tools Used:
          tool GATK: ${params.docker_image_gatk}
          tool validate: ${params.docker_image_validate}
          tool RTG: ${params.docker_image_rtg}

      Extra parameters:
          ${params}

------------------------------------
Starting workflow...
------------------------------------
        """

include { run_validate } from './modules/validation'
include { run_GenomicsDBImport_GATK; run_SplitIntervals_GATK; run_HaplotypeCaller_GATK; run_GenotypeGVCFs_GATK; run_SortVcf_GATK; run_MergeVcfs_GATK } from './modules/joint-genotype-processes'
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

Channel
 .fromPath(params.input_csv, checkIfExists: true)
 .splitCsv(header:true)
 .map{
    // Add the index file as a new key in the input_csv map
    // If the index file was explicitly added as a column then it will not be overwritten.
    [index_file: indexFile(it.input_file)] + it
  }
 .tap{ input_csv }
 .branch {
   bams: it.input_file.endsWith('.bam')
   gvcfs: it.input_file.endsWith('vcf.gz')
   other: true
 }
 .set { input_ch_input_csv }

// Create channel for validation
input_csv
  .flatMap{[it.input_file, it.index_file]}
  .set{ input_validation }

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

    gvcf_input = input_ch_input_csv.bams.combine(split_intervals) // Cross the input files with all the chr list
        .map{ input_csv,interval -> [input_csv.sample, input_csv.input_file, input_csv.index_file, interval] }

    run_HaplotypeCaller_GATK(
      params.reference_fasta,
      "${params.reference_fasta}.fai",
      params.reference_dict,
      params.bundle_v0_dbsnp138_vcf_gz,
      "${params.bundle_v0_dbsnp138_vcf_gz}.tbi",
      gvcf_input
      )

    chr_gvcfs_input = run_HaplotypeCaller_GATK.out.gvcf
        .mix(
          input_ch_input_csv.gvcfs
          // Cross the g.vcfs with all the chromosomes
          .combine(split_intervals)
          // Use all of the g.vcfs that are not split by chromosome (e.g. the field missing in csv)
          // If the g.vcf is already split by chromosome then
          //   only keep if the chromosomes match
          .filter{ input_csv,chr -> ! input_csv.chr || input_csv.chr == chr }
          .map{ input_csv,chr -> [chr, input_csv.input_file, input_csv.index_file] }
        )
        // Aggregate by chr
        .groupTuple()

    run_GenomicsDBImport_GATK(
      params.reference_fasta,
      "${params.reference_fasta}.fai",
      params.reference_dict,
      chr_gvcfs_input
    )


    run_GenotypeGVCFs_GATK(
      params.reference_fasta,
      "${params.reference_fasta}.fai",
      params.reference_dict,
      params.bundle_v0_dbsnp138_vcf_gz,
      "${params.bundle_v0_dbsnp138_vcf_gz}.tbi",
      run_GenomicsDBImport_GATK.out.gdb
    )

    // Sort the vcfs.
    run_SortVcf_GATK(run_GenotypeGVCFs_GATK.out.vcf)

    // Merge the sorted vcfs
    run_MergeVcfs_GATK(run_SortVcf_GATK.out.vcf.collect())

    // Build model for variant recalibration for SNPs
    // Separated SNPs and INDELs from GATK best practices
    recalibrate_snps(run_MergeVcfs_GATK.out.vcf)

    // Build model for variant recalibration for INDELs from the SNP recalibrated vcf
    recalibrate_indels(recalibrate_snps.out.vcf)
}
