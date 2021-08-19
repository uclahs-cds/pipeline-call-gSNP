process run_SplitIntervals_GATK {
    container params.docker_image_gatk

    publishDir params.output_dir,
               mode: "copy",
               pattern: "interval-files/*-scattered.interval_list",
               enabled: params.save_intermediate_files
    publishDir params.log_output_dir,
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

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
    publishDir path: "${params.output_dir}/${task.process.replace(':', '/')}",
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: '*.vcf*'

    publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

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

    script:
        out_filename_normal = "${normal_id}_${task.index}_raw_variants.g.vcf.gz"
        out_filename_tumour = "${tumour_id}_${task.index}_raw_variants.g.vcf.gz"
        out_filename_vcf = "${sample_id}_${task.index}.vcf"
        interval_str = "--intervals ${interval}"
        bam_input_str = params.is_NT_paired ? "--input ${bam} --input ${bam_tumour}" : "--input ${bam}"

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
        ${interval_str}

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
        --standard-min-confidence-threshold-for-calling 30 \
        ${interval_str}

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
        --standard-min-confidence-threshold-for-calling 30 \
        ${interval_str}
    fi
    """
}

process run_MergeVcfs_Picard {
    container params.docker_image_picard

    publishDir path: "${params.output_dir}/${task.process.replace(':', '/')}",
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: "*.vcf*"

    publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

    input:
    path(vcfs)
    path(gvcfs_normal)
    path(gvcfs_tumour)
    tuple val(sample_id), val(normal_id), val(tumour_id)

    output:
    path(".command.*")
    path("${sample_id}_merged_raw.vcf"), emit: vcf
    path("${sample_id}_merged_raw.vcf.idx"), emit: vcf_index
    path("${normal_id}_merged_raw_variants.g.vcf.gz"), emit: gvcf_normal
    path("${normal_id}_merged_raw_variants.g.vcf.gz.tbi"), emit: gvcf_normal_index
    path("${tumour_id}_merged_raw_variants.g.vcf.gz"), emit: gvcf_tumour optional true
    path("${tumour_id}_merged_raw_variants.g.vcf.gz.tbi"), emit: gvcf_tumour_index optional true

    script:
    vcf_args = vcfs.collect{ "-INPUT '$it'" }.join(' ')
    gvcf_args_normal = gvcfs_normal.collect{ "-INPUT '$it'" }.join(' ')
    gvcf_args_tumour = gvcfs_tumour.collect{ "-INPUT '$it'" }.join(' ')

    """
    set -euo pipefail

    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -Djava.io.tmpdir=/scratch \
      -jar /picard-tools/picard.jar MergeVcfs \
      ${vcf_args} \
      -OUTPUT ${sample_id}_merged_raw.vcf \
      -VALIDATION_STRINGENCY LENIENT

    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -Djava.io.tmpdir=/scratch \
      -jar /picard-tools/picard.jar MergeVcfs \
      ${gvcf_args_normal} \
      -OUTPUT ${normal_id}_merged_raw_variants.g.vcf.gz \
      -VALIDATION_STRINGENCY LENIENT
    
    if ${params.is_NT_paired}
    then
      java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -Djava.io.tmpdir=/scratch \
        -jar /picard-tools/picard.jar MergeVcfs \
        ${gvcf_args_tumour} \
        -OUTPUT ${tumour_id}_merged_raw_variants.g.vcf.gz \
        -VALIDATION_STRINGENCY LENIENT
    fi
    """
}
