/*
    Nextflow module for reheadering BAM files

    input:
        (sample_id, normal_id, tumour_id):  tuples of string identifiers for the samples
        normal_bam: path to input normal BAM
        normal_bam_index: path to input normal BAM index
        tumour_bam: path to input tumour BAM
        tumour_bam_index: path to input tumour BAM index
        interval: path to associated intervals file to propagate

    params:
        params.output_dir: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_samtools: string
*/
process run_reheader_SAMtools {
    container params.docker_image_samtools
    publishDir path: "${params.output_dir}/intermediate/${task.process.replace(':', '/')}",
        mode: "copy",
        enabled: params.save_intermediate_files,
        pattern: "*_reheadered_*"

    publishDir path: "${params.log_output_dir}/process-log",
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
    path(normal_bam), emit: normal_bam_for_deletion
    path(tumour_bam), emit: tumour_bam_for_deletion
    path(normal_bam_index), emit: normal_bam_index_for_deletion
    path(tumour_bam_index), emit: tumour_bam_index_for_deletion

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

/*
    Nextflow module for indexing BAM files

    input:
        normal_bam: path to input normal BAM
        tumour_bam_index: path to input tumour BAM index
        interval: path to associated intervals file to propagate

    params:
        params.output_dir: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_picard: string
        params.gatk_command_mem_diff: float(memory)
*/
process run_BuildBamIndex_Picard {
    container params.docker_image_picard
    publishDir path: "${params.output_dir}/intermediate/${task.process.replace(':', '/')}",
        mode: "copy",
        enabled: params.save_intermediate_files,
        pattern: "*_reheadered_*"

    publishDir path: "${params.log_output_dir}/process-log",
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

/*
    Nextflow module for merging BAM files

    input:
        bams: list or tuple of paths to BAMs to be merged
        sample_type: string. Indicator of normal or tumour samples
        (sample_id, normal_id, tumour_id):  tuples of string identifiers for the samples

    params:
        params.output_dir: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_picard: string
        params.gatk_command_mem_diff: float(memory)
*/
process run_MergeSamFiles_Picard {
    container params.docker_image_picard
    publishDir path: "${params.output_dir}/output",
        mode: "copy",
        pattern: "*_merged*",
        saveAs: { filename -> (file(filename).getExtension() == "bai") ? "${file(filename).baseName}.bam.bai" : "${filename}" }

    publishDir path: "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

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
    normal_bam_for_deletion = run_reheader_SAMtools.out.normal_bam_for_deletion
    tumour_bam_for_deletion = run_reheader_SAMtools.out.tumour_bam_for_deletion
    normal_bam_index_for_deletion = run_reheader_SAMtools.out.normal_bam_index_for_deletion
    tumour_bam_index_for_deletion = run_reheader_SAMtools.out.tumour_bam_index_for_deletion
}
