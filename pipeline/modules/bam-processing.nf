process run_reheader_SAMtools {
    container params.docker_image_samtools
    publishDir path: "${params.output_dir}/${task.process.replace(':', '/')}",
        mode: "copy",
        enabled: params.save_intermediate_files,
        pattern: "*_reheadered_*"

    publishDir path: params.log_output_dir,
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    tuple val(sample_id), val(normal_id), val(tumour_id)
    path(normal_bam)
    path(normal_bam_index)
    path(tumour_bam)
    path(tumour_bam_index)
    path(interval)

    output:
    path(".command.*")
    path("${normal_id}_recalibrated_reheadered_${task.index}.bam"), emit: normal_bam_reheadered
    path("${tumour_id}_recalibrated_reheadered_${task.index}.bam"), emit: tumour_bam_reheadered
    path(interval), emit: associated_interval

    script:
    """
    set -euo pipefail

    samtools view -H ${normal_bam} | \
    grep -v -P "^@RG.*SM:${tumour_id}" | \
    samtools reheader - ${normal_bam} \
        > ${normal_id}_recalibrated_reheadered_${task.index}.bam

    samtools view -H ${tumour_bam} | \
    grep -v -P "^@RG.*SM:${normal_id}" | \
    samtools reheader - ${tumour_bam} \
        > ${tumour_id}_recalibrated_reheadered_${task.index}.bam
    """
}

process run_BuildBamIndex_Picard {
    container params.docker_image_picard
    publishDir path: "${params.output_dir}/${task.process.replace(':', '/')}",
        mode: "copy",
        enabled: params.save_intermediate_files,
        pattern: "*_reheadered_*"

    publishDir path: params.log_output_dir,
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path(normal_bam)
    path(tumour_bam)
    path(interval)

    output:
    path(".command.*")
    path(normal_bam), emit: normal_bam_reheadered
    path("${normal_bam}.bai"), emit: normal_bam_reheadered_index
    path(tumour_bam), emit: tumour_bam_reheadered
    path("${tumour_bam}.bai"), emit: tumour_bam_reheadered_index
    path(interval), emit: associated_interval

    script:
    """
    set -euo pipefail
    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -Djava.io.tmpdir=/scratch \
        -jar /picard-tools/picard.jar BuildBamIndex \
        -VALIDATION_STRINGENCY LENIENT \
        -INPUT ${normal_bam} \
        -OUTPUT ${normal_bam}.bai \

    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -Djava.io.tmpdir=/scratch \
        -jar /picard-tools/picard.jar BuildBamIndex \
        -VALIDATION_STRINGENCY LENIENT \
        -INPUT ${tumour_bam} \
        -OUTPUT ${tumour_bam}.bai \
    """
}

process run_MergeSamFiles_Picard {
    container params.docker_image_picard
    publishDir path: "${params.output_dir}/${task.process.replace(':', '/')}",
        mode: "copy",
        pattern: "*_merged*"

    publishDir path: params.log_output_dir,
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path(bams)
    val(sample_type)
    tuple val(sample_id), val(normal_id), val(tumour_id)

    output:
    path(".command.*")
    path("${output_id}_realigned_recalibrated_merged.bam"), emit: merged_bam
    path("${output_id}_realigned_recalibrated_merged.bai"), emit: merged_bam_index

    script:
    all_bams = bams.collect{ "-INPUT '$it'" }.join(' ')
    output_id = (sample_type == "normal") ? "${normal_id}" : "${tumour_id}"
    """
    set -euo pipefail
    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -Djava.io.tmpdir=/scratch \
        -jar /picard-tools/picard.jar MergeSamFiles \
        ${all_bams} \
        -OUTPUT ${output_id}_realigned_recalibrated_merged.bam \
        -CREATE_INDEX true \
        -SORT_ORDER coordinate \
        -ASSUME_SORTED false \
        -USE_THREADING false \
        -VALIDATION_STRINGENCY LENIENT
    """
}

workflow reheader_interval_bams {
    take:
    identifiers
    normal_bams
    normal_bams_index
    tumour_bams
    tumour_bams_index
    intervals

    main:
    run_reheader_SAMtools(
        identifiers,
        normal_bams,
        normal_bams_index,
        tumour_bams,
        tumour_bams_index,
        intervals
        )

    run_BuildBamIndex_Picard(
        run_reheader_SAMtools.out.normal_bam_reheadered,
        run_reheader_SAMtools.out.tumour_bam_reheadered,
        run_reheader_SAMtools.out.associated_interval
        )

    emit:
    reheadered_normal_bam = run_BuildBamIndex_Picard.out.normal_bam_reheadered
    reheadered_normal_bam_index = run_BuildBamIndex_Picard.out.normal_bam_reheadered_index
    reheadered_tumour_bam = run_BuildBamIndex_Picard.out.tumour_bam_reheadered
    reheadered_tumour_bam_index = run_BuildBamIndex_Picard.out.tumour_bam_reheadered_index
    associated_interval = run_BuildBamIndex_Picard.out.associated_interval
}
