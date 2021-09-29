process run_SplitIntervals_GATK {
    container params.docker_image_gatk

    publishDir path: "${params.output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "interval-files/*-scattered.interval_list",
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

    output:
    path 'interval-files/*-scattered.interval_list', emit: interval_list
    path ".command.*"

    """
    set -euo pipefail
    gatk SplitIntervals \
        -R $reference \
        -L $intervals \
        --scatter-count ${params.scatter_count} \
        ${params.split_intervals_extra_args} \
        -O interval-files
    """
}

process run_HaplotypeCaller_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir}/intermediate/${task.process.replace(':', '/')}",
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: '*.vcf*'

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}/${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(dbsnp_bundle)
    path(dbsnp_bundle_index)
    tuple val(sample_id), val(normal_id), val(tumour_id)
    path(bam)
    path(bam_index)
    path(bam_tumour)
    path(bam_index_tumour)
    path(interval)


    output:
      path(".command.*")
      path("${sample_id}_${task.index}.vcf"), emit: vcf
      path("${sample_id}_${task.index}.vcf.idx"), emit: vcf_index
      path("${normal_id}_${task.index}_raw_variants.g.vcf.gz"), emit: gvcf_normal
      path("${normal_id}_${task.index}_raw_variants.g.vcf.gz.tbi"), emit: gvcf_normal_index
      path("${tumour_id}_${task.index}_raw_variants.g.vcf.gz"), emit: gvcf_tumour optional true
      path("${tumour_id}_${task.index}_raw_variants.g.vcf.gz.tbi"), emit: gvcf_tumour_index optional true
      path(bam), emit: normal_bam_for_deletion
      path(bam_index), emit: normal_bam_index_for_deletion
      path(bam_tumour), emit: tumour_bam_for_deletion optional true
      path(bam_index_tumour), emit: tumour_bam_index_for_deletion optional true

    script:
        out_filename_normal = "${normal_id}_${task.index}_raw_variants.g.vcf.gz"
        out_filename_tumour = "${tumour_id}_${task.index}_raw_variants.g.vcf.gz"
        out_filename_vcf = "${sample_id}_${task.index}.vcf"
        interval_str = "--intervals ${interval}"
        bam_input_str = params.is_NT_paired ? "--input ${bam} --input ${bam_tumour}" : "--input ${bam}"
        interval_padding = params.is_targeted ? "--interval-padding 100" : ""

    """
    set -euo pipefail

    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
        HaplotypeCaller \
        ${bam_input_str} \
        --output ${out_filename_vcf} \
        --reference ${reference_fasta} \
        --verbosity INFO \
        --output-mode EMIT_VARIANTS_ONLY \
        --dbsnp ${dbsnp_bundle} \
        --sample-ploidy 2 \
        --standard-min-confidence-threshold-for-calling 50 \
        ${interval_str} \
        ${interval_padding}

    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
        HaplotypeCaller \
        --input ${bam} \
        --output ${out_filename_normal} \
        --reference ${reference_fasta} \
        --verbosity INFO \
        --output-mode EMIT_VARIANTS_ONLY \
        --emit-ref-confidence GVCF \
        --dbsnp ${dbsnp_bundle} \
        --sample-ploidy 2 \
        ${interval_str} \
        ${interval_padding}

    if ${params.is_NT_paired}
    then
      gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
        HaplotypeCaller \
        --input ${bam_tumour} \
        --output ${out_filename_tumour} \
        --reference ${reference_fasta} \
        --verbosity INFO \
        --output-mode EMIT_VARIANTS_ONLY \
        --emit-ref-confidence GVCF \
        --dbsnp ${dbsnp_bundle} \
        --sample-ploidy 2 \
        ${interval_str} \
        ${interval_padding}
    fi
    """
}

process run_MergeVcfs_Picard {
    container params.docker_image_picard

    publishDir path: "${params.output_dir}/output",
      mode: "copy",
      pattern: "*.vcf*"

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
    path(vcfs)
    val(vcf_type)
    val(sample_type)
    tuple val(sample_id), val(normal_id), val(tumour_id)

    output:
    path(".command.*")
    path("*.vcf{,.gz}"), emit: vcf
    path("*.vcf.{idx,gz.tbi}"), emit: vcf_index

    script:
    all_vcfs = vcfs.collect{ "-INPUT '$it'" }.join(' ')
    output_gvcf_id = (sample_type == "normal") ? "${normal_id}" : "${tumour_id}"
    output_filename = (vcf_type == "GVCF") ? "${output_gvcf_id}_merged_raw_variants.g.vcf.gz" : "${sample_id}_merged_raw.vcf"

    """
    set -euo pipefail

    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -Djava.io.tmpdir=/scratch \
      -jar /picard-tools/picard.jar MergeVcfs \
      ${all_vcfs} \
      -OUTPUT ${output_filename} \
      -VALIDATION_STRINGENCY LENIENT
    """
}
