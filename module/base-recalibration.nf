include { generate_standard_filename } from '../external/nextflow-module/modules/common/generate_standardized_filename/main.nf'
/*
    Nextflow module for generating base recalibration table

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        bundle_mills_and_1000g_gold_standards_vcf_gz: path to standard Mills and 1000 genomes variants
        bundle_mills_and_1000g_gold_standards_vcf_gz_tbi: path to index file for Mills and 1000g variants
        bundle_known_indels_vcf_gz: path to set of known indels
        bundle_known_indels_vcf_gz_tbi: path to index of known indels VCF
        bundle_v0_dbsnp138_vcf_gz: path to dbSNP variants
        bundle_v0_dbsnp138_vcf_gz_tbi: path to index of dbSNP variants
        all_intervals: list or tuple of paths to all split intervals
        indelrealigned_bams: list or tuple of paths to indel realigned BAMs
        indelrealigned_bams_bai: list or tuple of paths to indel realigned BAM indices
        (sample_id, normal_id, tumour_id):  tuples of string identifiers for the samples
        
    params:
        params.output_dir: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_gatk: string
        params.is_targeted: bool. Indicator of whether in targeted exome mode or in WGS mode
        params.gatk_command_mem_diff: float(memory)
*/
process run_BaseRecalibrator_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir}/intermediate/${task.process.replace(':', '/')}",
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: "*.grp"

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}-${id}/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(bundle_mills_and_1000g_gold_standards_vcf_gz)
    path(bundle_mills_and_1000g_gold_standards_vcf_gz_tbi)
    path(bundle_known_indels_vcf_gz)
    path(bundle_known_indels_vcf_gz_tbi)
    path(bundle_v0_dbsnp138_vcf_gz)
    path(bundle_v0_dbsnp138_vcf_gz_tbi)
    path(all_target_intervals)
    path(indelrealigned_bams)
    path(indelrealigned_bams_bai)
    val(id)

    output:
    path(".command.*")
    path("${id}_recalibration_table.grp"), emit: recalibration_table

    script:
    all_ir_bams = indelrealigned_bams.collect{ "--input '$it'" }.join(' ')
    targeted_options = params.is_targeted ? "--intervals ${all_target_intervals} --interval-padding 100" : ""
    """
    set -euo pipefail
    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir}" \
        BaseRecalibrator \
        ${all_ir_bams} \
        --reference ${reference_fasta} \
        --verbosity INFO \
        --known-sites ${bundle_mills_and_1000g_gold_standards_vcf_gz} \
        --known-sites ${bundle_known_indels_vcf_gz} \
        --known-sites ${bundle_v0_dbsnp138_vcf_gz} \
        --output ${id}_recalibration_table.grp \
        ${targeted_options} \
        --read-filter SampleReadFilter \
        --sample ${id}
    """
}

/*
    Nextflow module for recalibrating base quality scores in BAM file

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        recalibration_table: path to base recalibration table
        indelrealigned_bam: list or tuple of paths to indel realigned BAMs
        indelrealigned_bam_index: list or tuple of paths to indel realigned BAM indices
        interval: path to specific intervals file associated with input BAM
        includes_unmapped: boolean to indicate if unmapped reads are included in input BAM
        (sample_id, normal_id, tumour_id):  tuples of string identifiers for the samples

    params:
        params.output_dir: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_gatk: string
        params.gatk_command_mem_diff: float(memory)
        params.is_targeted: bool. Indicator of whether in targeted exome mode or in WGS mode
        params.is_emit_original_quals: bool. Indicator of whether to keep original base quality scores
        params.is_NT_paired: bool. Indicator of whether the input samples include both normal and tumour
*/
process run_ApplyBQSR_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir}/intermediate/${task.process.replace(':', '/')}",
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: "*_recalibrated-*",
      saveAs: { filename -> (file(filename).getExtension() == "bai") ? "${file(filename).baseName}.bam.bai" : "${filename}" }

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}-${interval_id}/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(recalibration_table)
    tuple path(indelrealigned_bam),
          path(indelrealigned_bam_index),
          path(interval),
          val(includes_unmapped),
          val(id)

    output:
    path(".command.*")
    tuple path("*_recalibrated-${interval_id}.bam"),
          path("*_recalibrated-${interval_id}.bai"), emit: apply_bqsr_och
    tuple path(indelrealigned_bam),
          path(indelrealigned_bam_index), emit: deletion_och

    script:
    // Get split interval number to serve as task ID
    interval_id = interval.baseName.split('-')[0]
    unmapped_interval_option = (includes_unmapped) ? "--intervals unmapped" : ""
    combined_interval_options = "--intervals ${interval} ${unmapped_interval_option}"
    all_commands = id.collect{
      """
      gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \\
        ApplyBQSR \\
        --input ${indelrealigned_bam} \\
        --bqsr-recal-file ${it.split("-")[0..-2].join("-")}_recalibration_table.grp \\
        --reference ${reference_fasta} \\
        --read-filter SampleReadFilter \\
        ${combined_interval_options} \\
        --emit-original-quals ${params.is_emit_original_quals} \\
        --output ${generate_standard_filename(params.aligner, params.dataset_id, it.split("-")[0..-2].join("-"), ['additional_tools': "GATK-${params.gatk_version}"])}_recalibrated-${interval_id}.bam \\
        --sample ${it.split("-")[0..-2].join("-")}
      """
      }
      .join("\n")
    """
    set -euo pipefail
    ${all_commands}
    """
}
