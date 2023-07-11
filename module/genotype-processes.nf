include { generate_standard_filename } from '../external/nextflow-module/modules/common/generate_standardized_filename/main.nf'
/*
    Nextflow module for splitting input intervals into multiple intervals for parallelization

    input:
        intervals: path to set of target intervals to split
        reference: path to reference genome fasta file
        reference_index: path to index for reference fasta
        reference_dict: path to dictionary for reference fasta
        out_folder_name: name of output folder for storing scattered intervals
        targeted_intervals: bool. to indicate whether the intervals are specifically targeted

    params:
        params.output_dir: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_gatk: string
        params.scatter_count: integer. Number of intervals to split into
        params.split_intervals_extra_args: string. Additional arguments for splitting intervals
*/
process run_SplitIntervals_GATK {
    container params.docker_image_gatk

    publishDir path: "${params.output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "${out_folder_name}/*-scattered.interval_list",
               enabled: params.save_intermediate_files
    publishDir "${params.log_output_dir}/process-log",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
    path intervals
    path reference
    path reference_index
    path reference_dict
    val out_folder_name
    val targeted_intervals

    output:
    path "${out_folder_name}/*-scattered.interval_list", emit: interval_list
    path ".command.*"

    script:
    subdivision_mode = targeted_intervals ? "--subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION" : ""
    """
    set -euo pipefail
    gatk SplitIntervals \
        -R $reference \
        -L $intervals \
        --scatter-count ${params.scatter_count} \
        ${subdivision_mode} \
        ${params.split_intervals_extra_args} \
        -O ${out_folder_name}
    """
}

/*
    Nextflow module for calling haplotypes in VCF mode

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        dbsnp_bundle: path to dbSNP variants
        dbsnp_bundle_index: path to index of dbSNP variants
        (sample_id, normal_id, tumour_id):  tuples of string identifiers for the samples
        bam: path to normal BAM for calling
        bam_index: path to index of normal BAM
        bam_tumour: path to tumour BAM for calling
        bam_tumour_index: path to index of tumour BAM
        interval: path to specific intervals for calling
        
    params:
        params.output_dir: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_gatk: string
        params.is_targeted: bool. Indicator of whether in targeted exome mode or in WGS mode
        params.is_NT_paired: bool. Indicator of whether input has normal and tumour samples
        params.gatk_command_mem_diff: float(memory)
*/
process run_HaplotypeCallerVCF_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir}/intermediate/${task.process.replace(':', '/')}",
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: '*.vcf*'

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}/${task.process.split(':')[-1]}-${interval_id}/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(dbsnp_bundle)
    path(dbsnp_bundle_index)
    val(sample_id)
    path(bams)
    path(bams_index)
    path(interval)


    output:
    path(".command.*")
    tuple val(sample_id), path(output_filename), path("${output_filename}.idx"), emit: vcfs
    path(bams), emit: bams_for_deletion
    path(bams_index), emit: bams_index_for_deletion

    script:
    // Get split interval number to serve as task ID
    interval_id = interval.baseName.split('-')[0]
    output_filename = generate_standard_filename(
        "GATK-${params.gatk_version}",
        params.dataset_id,
        sample_id,
        [
            'additional_information': "${interval_id}.vcf"
        ]
    )
    interval_str = "--intervals ${interval}"
    bam_input_str = bams.collect{ "--input '$it'" }.join(' ')
    interval_padding = params.is_targeted ? "--interval-padding 100" : ""

    """
    set -euo pipefail

    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir}" \
        HaplotypeCaller \
        ${bam_input_str} \
        --output ${output_filename} \
        --reference ${reference_fasta} \
        --verbosity INFO \
        --output-mode EMIT_VARIANTS_ONLY \
        --dbsnp ${dbsnp_bundle} \
        --sample-ploidy 2 \
        --standard-min-confidence-threshold-for-calling 50 \
        ${interval_str} \
        ${interval_padding}
    """
}


/*
    Nextflow module for calling haplotypes in GVCF mode

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        dbsnp_bundle: path to dbSNP variants
        dbsnp_bundle_index: path to index of dbSNP variants
        (sample_id, normal_id, tumour_id):  tuples of string identifiers for the samples
        bam: path to BAM for calling
        bam_index: path to index of BAM
        interval: path to specific intervals for calling
        sample_type: val to indicate whether processing normal or tumour sample
        
    params:
        params.output_dir: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_gatk: string
        params.is_targeted: bool. Indicator of whether in targeted exome mode or in WGS mode
        params.is_NT_paired: bool. Indicator of whether input has normal and tumour samples
        params.gatk_command_mem_diff: float(memory)
*/
process run_HaplotypeCallerGVCF_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir}/intermediate/${task.process.replace(':', '/')}",
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: '*.vcf*'

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}/${task.process.split(':')[-1]}-${id}-${interval_id}/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(dbsnp_bundle)
    path(dbsnp_bundle_index)
    tuple path(bam), path(bam_index), path(interval), val(id)


    output:
    path(".command.*")
    tuple val(id), path(output_filename), path("${output_filename}.tbi"), emit: gvcfs
    path(bam), emit: bam_for_deletion
    path(bam_index), emit: bam_index_for_deletion

    script:
    // Get split interval number to serve as task ID
    interval_id = interval.baseName.split('-')[0]
    output_filename = generate_standard_filename(
        "GATK-${params.gatk_version}",
        params.dataset_id,
        id,
        [
            'additional_information': "${interval_id}_raw_variants.g.vcf.gz"
        ]
    )
    interval_str = "--intervals ${interval}"
    interval_padding = params.is_targeted ? "--interval-padding 100" : ""
    output_mode = params.emit_all_confident_sites ? "EMIT_ALL_CONFIDENT_SITES" : "EMIT_VARIANTS_ONLY"
    """
    set -euo pipefail

    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir}" \
        HaplotypeCaller \
        --input ${bam} \
        --output ${output_filename} \
        --reference ${reference_fasta} \
        --verbosity INFO \
        --output-mode ${output_mode} \
        --emit-ref-confidence GVCF \
        --dbsnp ${dbsnp_bundle} \
        --sample-ploidy 2 \
        ${interval_str} \
        ${interval_padding}
    """
}

/*
    Nextflow module for merging input VCFs

    input:
        vcfs: list or tuples of inputs VCFs to merge
        vcf_type: string to indicate whether calling in germline mode
        sample_type: string to indicate whether merging normal or tumour files
        (sample_id, normal_id, tumour_id):  tuples of string identifiers for the samples
        
    params:
        params.output_dir: string(path)
        params.log_output_dir: string(path)
        params.docker_image_picard: string
        params.gatk_command_mem_diff: float(memory)
*/
process run_MergeVcfs_Picard {
    container params.docker_image_picard

    publishDir path: "${params.output_dir}/output",
      mode: "copy",
      pattern: "*.vcf*"

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}-${id}/log${file(it).getName()}" }

    input:
    path(vcfs)
    val(vcf_type)
    val(id)

    output:
    path(".command.*")
    val(id), emit: associated_id
    path("*.vcf{,.gz}"), emit: vcf
    path("*.vcf.{idx,gz.tbi}"), emit: vcf_index

    script:
    all_vcfs = vcfs.collect{ "-INPUT '$it'" }.join(' ')
    output_filename_base = generate_standard_filename("GATK-${params.gatk_version}", params.dataset_id, id, [:])
    output_filename = (vcf_type == "GVCF") ?
        "${output_filename_base}.g.vcf.gz" :
        "${output_filename_base}.vcf"
    """
    set -euo pipefail

    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -Djava.io.tmpdir=${workDir} \
      -jar /usr/local/share/picard-slim-2.26.10-0/picard.jar MergeVcfs \
      ${all_vcfs} \
      -OUTPUT ${output_filename} \
      -VALIDATION_STRINGENCY LENIENT
    """
}
