process run_BaseRecalibrator_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir}/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: "*.grp"

    publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "run_BaseRecalibrator_GATK/log${file(it).getName()}" }

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
    path(indelrealigned_bams)
    path(indelrealigned_bams_bai)
    tuple val(sample_id), val(normal_id), val(tumour_id)

    output:
    path(".command.*")
    path("${sample_id}_recalibration_table.grp"), emit: recalibration_table

    script:
    all_ir_bams = indelrealigned_bams.collect{ "--input '$it'" }.join(' ')
    """
    set -euo pipefail
    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
        BaseRecalibrator \
        ${all_ir_bams} \
        --reference ${reference_fasta} \
        --verbosity INFO \
        --known-sites ${bundle_mills_and_1000g_gold_standards_vcf_gz} \
        --known-sites ${bundle_known_indels_vcf_gz} \
        --known-sites ${bundle_v0_dbsnp138_vcf_gz} \
        --output ${sample_id}_recalibration_table.grp
    """
}

process run_ApplyBQSR_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir}/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: "*_recalibrated_*"

    publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(recalibration_table)
    path(indelrealigned_bam)
    path(indelrealigned_bam_index)
    tuple val(sample_id), val(normal_id), val(tumour_id), path(bam), path(bam_index), path(bam_tumour), path(bam_index_tumour), path(interval)

    output:
    path(".command.*")
    path(interval), emit: assoc_interval
    path("${normal_id}_recalibrated_${task.index}.bam"), emit: recalibrated_normal_bam
    path("${normal_id}_recalibrated_${task.index}.bai"), emit: recalibrated_normal_bam_index
    path("${tumour_id}_recalibrated_${task.index}.bam"), emit: recalibrated_tumour_bam optional true
    path("${tumour_id}_recalibrated_${task.index}.bam"), emit: recalibrated_tumour_bam_index optional true

    script:
    """
    set -euo pipefail
    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
        ApplyBQSR \
        --input ${indelrealigned_bam} \
        --bqsr-recal-file ${recalibration_table} \
        --reference ${reference_fasta} \
        --output ${normal_id}_recalibrated_${task.index}.bam \
        --read-filter SampleReadFilter \
        --sample ${normal_id} \
        --intervals ${interval}

    if ${params.is_NT_paired}
    then
        gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
            ApplyBQSR \
            --input ${indelrealigned_bam} \
            --bqsr-recal-file ${recalibration_table} \
            --reference ${reference_fasta} \
            --output ${tumour_id}_recalibrated_${task.index}.bam \
            --read-filter SampleReadFilter \
            --sample ${tumour_id} \
            --intervals ${interval}
    fi
    """
}

