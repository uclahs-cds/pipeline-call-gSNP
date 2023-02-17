include { generate_standard_filename } from '../external/nextflow-module/modules/common/generate_standardized_filename/main.nf'
/*
    Nextflow module for generating realignment targets

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        bundle_mills_and_1000g_gold_standards_vcf_gz: path to standard Mills and 1000 genomes variants
        bundle_mills_and_1000g_gold_standards_vcf_gz_tbi: path to index file for Mills and 1000g variants
        bundle_known_indels_vcf_gz: path to set of known indels
        bundle_known_indels_vcf_gz_tbi: path to index of known indels VCF
        (sample_id, normal_id, tumour_id, bam, bam_index, bam_tumour, bam_index_tumour, interval):  
          tuples of string identifiers for the samples, input BAM and index files, and target interval
        
    params:
        params.output_dir: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_gatk3: string
        params.is_targeted: bool. Indicator of whether in targeted exome mode or in WGS mode
        params.is_NT_paired: bool. Indicator of whether input has normal and tumour samples
        params.gatk_command_mem_diff: float(memory)
*/
process run_RealignerTargetCreator_GATK {
    container params.docker_image_gatk3
    publishDir path: "${params.output_dir}/intermediate/${task.process.replace(':', '/')}",
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: "*.intervals"

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}-${interval_id}/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(bundle_mills_and_1000g_gold_standards_vcf_gz)
    path(bundle_mills_and_1000g_gold_standards_vcf_gz_tbi)
    path(bundle_known_indels_vcf_gz)
    path(bundle_known_indels_vcf_gz_tbi)
    tuple val(sample_id), val(normal_id), val(tumour_id), path(bam), path(bam_index), path(bam_tumour), path(bam_index_tumour), path(interval)

    output:
    path(".command.*")
    path(interval), emit: scatter_intervals
    path("${sample_id[0]}_RTC_${interval_id}.intervals"), emit: intervals_RTC

    script:
    // Get split interval number to serve as task ID
    interval_id = interval.baseName.split('-')[0]
    tumour_bams = bam_tumour.collect{ "--input_file '$it'" }.join(' ')
    bam_input_str = params.is_NT_paired ? "--input_file ${bam} ${tumour_bams}" : "--input_file ${bam}"
    interval_padding = params.is_targeted ? "--interval_padding 100" : ""
    targeted_interval_params = params.is_targeted ? "--intervals ${params.intervals} --interval_set_rule INTERSECTION" : ""
    """
    set -euo pipefail
    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir} \
        -jar /GenomeAnalysisTK.jar \
        --analysis_type RealignerTargetCreator \
        ${bam_input_str} \
        --reference_sequence ${reference_fasta} \
        --known ${bundle_mills_and_1000g_gold_standards_vcf_gz} \
        --known ${bundle_known_indels_vcf_gz} \
        --intervals ${interval} \
        --out ${sample_id[0]}_RTC_${interval_id}.intervals \
        --allow_potentially_misencoded_quality_scores \
        --num_threads 2 \
        ${interval_padding} \
        ${targeted_interval_params} || touch ${sample_id[0]}_RTC_${interval_id}.intervals
    """
}

/*
    Nextflow module for realigning indels

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        bundle_mills_and_1000g_gold_standards_vcf_gz: path to standard Mills and 1000 genomes variants
        bundle_mills_and_1000g_gold_standards_vcf_gz_tbi: path to index file for Mills and 1000g variants
        bundle_known_indels_vcf_gz: path to set of known indels
        bundle_known_indels_vcf_gz_tbi: path to index of known indels VCF
        (sample_id, normal_id, tumour_id, bam, bam_index, bam_tumour, bam_index_tumour, interval): 
          tuples of string identifiers for the samples, input BAM and index files, and target interval
        target_intervals_RTC: path to realignment target intervals
        scatter_intervals: path to intervals being operated on
        
    params:
        params.output_dir: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_gatk3: string
        params.is_targeted: bool. Indicator of whether in targeted exome mode or in WGS mode
        params.is_NT_paired: bool. Indicator of whether input has normal and tumour samples
        params.gatk_command_mem_diff: float(memory)
*/
process run_IndelRealigner_GATK {
    container params.docker_image_gatk3
    publishDir path: "${params.output_dir}/intermediate/${task.process.replace(':', '/')}",
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: "*indelrealigned*",
      saveAs: { filename -> (file(filename).getExtension() == "bai") ? "${file(filename).baseName}.bam.bai" : "${filename}" }

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}-${interval_id}/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(bundle_mills_and_1000g_gold_standards_vcf_gz)
    path(bundle_mills_and_1000g_gold_standards_vcf_gz_tbi)
    path(bundle_known_indels_vcf_gz)
    path(bundle_known_indels_vcf_gz_tbi)
    tuple val(sample_id), val(normal_id), val(tumour_id), path(bam), path(bam_index), path(bam_tumour), path(bam_index_tumour)
    path(target_intervals_RTC)
    path(scatter_intervals)

    output:
    path(".command.*")
    path(scatter_intervals), emit: associated_interval
    val(has_unmapped), emit: includes_unmapped
    path("${output_filename}.bam"), emit: realigned_indels_bam
    path("${output_filename}.bai"), emit: realigned_indels_bam_index

    script:
    // Get split interval number to serve as task ID
    interval_id = scatter_intervals.baseName.split('-')[0]
    tumour_bams = bam_tumour.collect{ "--input_file '$it'" }.join(' ')
    bam_input_str = params.is_NT_paired ? "--input_file ${bam} ${tumour_bams}" : "--input_file ${bam}"
    unmapped_interval_option = (interval_id == '0000') ? "--intervals unmapped" : ""
    has_unmapped = (interval_id == '0000') ? true : false
    combined_interval_options = "--intervals ${scatter_intervals} ${unmapped_interval_option}"
    output_filename = generate_standard_filename(
        "GATK-${params.gatk_version}",
        params.dataset_id,
        sample_id[0],
        [
            'additional_information': "indelrealigned_${interval_id}"
        ]
    )
    """
    set -euo pipefail
    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir} \
        -jar /GenomeAnalysisTK.jar \
        --analysis_type IndelRealigner \
        ${bam_input_str} \
        --reference_sequence ${reference_fasta} \
        --bam_compression ${params.gatk_ir_compression} \
        --knownAlleles ${bundle_mills_and_1000g_gold_standards_vcf_gz} \
        --knownAlleles ${bundle_known_indels_vcf_gz} \
        --allow_potentially_misencoded_quality_scores \
        --targetIntervals ${target_intervals_RTC} \
        --out ${output_filename}.bam \
        ${combined_interval_options}
    """
}

workflow realign_indels {
    take:
    ir_input
    ir_input_no_interval

    main:
    run_RealignerTargetCreator_GATK(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz,
        "${params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz}.tbi",
        params.bundle_known_indels_vcf_gz,
        "${params.bundle_known_indels_vcf_gz}.tbi",
        ir_input
        )

    run_IndelRealigner_GATK(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz,
        "${params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz}.tbi",
        params.bundle_known_indels_vcf_gz,
        "${params.bundle_known_indels_vcf_gz}.tbi",
        ir_input_no_interval,
        run_RealignerTargetCreator_GATK.out.intervals_RTC,
        run_RealignerTargetCreator_GATK.out.scatter_intervals
        )

    emit:
    associated_interval = run_IndelRealigner_GATK.out.associated_interval
    includes_unmapped = run_IndelRealigner_GATK.out.includes_unmapped
    realigned_bam = run_IndelRealigner_GATK.out.realigned_indels_bam
    realigned_bam_index = run_IndelRealigner_GATK.out.realigned_indels_bam_index
}
