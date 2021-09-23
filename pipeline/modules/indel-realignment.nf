process run_RealignerTargetCreator_GATK {
    container params.docker_image_gatk3
    publishDir path: "${params.output_dir}/${params.docker_image_gatk3.split("/")[1].replace(':', '-')}/intermediate/${task.process.replace(':', '/')}",
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: "*.intervals"

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

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
    path("${sample_id}_RTC_${task.index}.intervals"), emit: intervals_RTC

    script:
    bam_input_str = params.is_NT_paired ? "--input_file ${bam} --input_file ${bam_tumour}" : "--input_file ${bam}"
    interval_padding = params.is_targeted ? "--interval_padding 100" : ""
    """
    set -euo pipefail
    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch \
        -jar /GenomeAnalysisTK.jar \
        --analysis_type RealignerTargetCreator \
        ${bam_input_str} \
        --reference_sequence ${reference_fasta} \
        --known ${bundle_mills_and_1000g_gold_standards_vcf_gz} \
        --known ${bundle_known_indels_vcf_gz} \
        --intervals ${interval} \
        --out ${sample_id}_RTC_${task.index}.intervals \
        --allow_potentially_misencoded_quality_scores \
        --num_threads 2 \
        ${interval_padding}
    """
}

process run_IndelRealigner_GATK {
    container params.docker_image_gatk3
    publishDir path: "${params.output_dir}/${params.docker_image_gatk3.split("/")[1].replace(':', '-')}/intermediate/${task.process.replace(':', '/')}",
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: "*_indelrealigned_*"

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

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
    path("${sample_id}_indelrealigned_${task.index}.bam"), emit: realigned_indels_bam
    path("${sample_id}_indelrealigned_${task.index}.bai"), emit: realigned_indels_bam_index

    script:
    bam_input_str = params.is_NT_paired ? "--input_file ${bam} --input_file ${bam_tumour}" : "--input_file ${bam}"
    unmapped_interval_option = (task.index == 1) ? "--intervals unmapped" : ""
    has_unmapped = (task.index == 1) ? true : false
    combined_interval_options = (params.is_targeted) ? "" : "--intervals ${scatter_intervals} ${unmapped_interval_option}"
    """
    set -euo pipefail
    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch \
        -jar /GenomeAnalysisTK.jar \
        --analysis_type IndelRealigner \
        ${bam_input_str} \
        --reference_sequence ${reference_fasta} \
        --bam_compression 0 \
        --knownAlleles ${bundle_mills_and_1000g_gold_standards_vcf_gz} \
        --knownAlleles ${bundle_known_indels_vcf_gz} \
        --allow_potentially_misencoded_quality_scores \
        --targetIntervals ${target_intervals_RTC} \
        --out ${sample_id}_indelrealigned_${task.index}.bam \
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
        params.reference_dict,
        params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz,
        "${params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz}.tbi",
        params.bundle_known_indels_vcf_gz,
        "${params.bundle_known_indels_vcf_gz}.tbi",
        ir_input
        )

    run_IndelRealigner_GATK(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        params.reference_dict,
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
