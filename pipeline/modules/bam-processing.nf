/*
    Nextflow module for reheadering BAM files

    input:
        (sample_id, normal_id, tumour_id):  tuples of string identifiers for the samples
        bam: path to input BAM
        bam_index: path to input BAM index
        interval: path to associated intervals file to propagate
        bam_type: indicator of whether processing normal or tumour BAM

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
    path(bam)
    path(bam_index)
    path(interval)
    val(bam_type)

    output:
    path(".command.*")
    path("${keep_id}_recalibrated_reheadered_${task.index}.bam"), emit: bam_reheadered
    path(interval), emit: associated_interval
    path(bam), emit: bam_for_deletion
    path(bam_index), emit: bam_index_for_deletion

    script:
    keep_id = (bam_type == 'normal') ? normal_id : tumour_id
    remove_id = (bam_type == 'normal') ? tumour_id : normal_id
    """
    set -euo pipefail

    samtools reheader \
        -c 'grep -v -P "^@RG.*SM:${remove_id}"' \
        ${bam} \
        > ${keep_id}_recalibrated_reheadered_${task.index}.bam
    """
}

/*
    Nextflow module for indexing BAM files

    input:
        bam: path to input normal BAM
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
    path(bam)
    path(interval)

    output:
    path(".command.*")
    tuple path(bam), path("${bam}.bai"), path(interval), emit: indexed_out

    script:
    """
    set -euo pipefail
    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -Djava.io.tmpdir=/scratch \
        -jar /usr/local/share/picard-slim-2.26.8-0/picard.jar BuildBamIndex \
        -VALIDATION_STRINGENCY LENIENT \
        -INPUT ${bam} \
        -OUTPUT ${bam}.bai \
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
        -jar /usr/local/share/picard-slim-2.26.8-0/picard.jar MergeSamFiles \
        ${all_bams} \
        -OUTPUT ${output_id}_realigned_recalibrated_merged.bam \
        -CREATE_INDEX true \
        -SORT_ORDER coordinate \
        -ASSUME_SORTED false \
        -USE_THREADING true \
        -VALIDATION_STRINGENCY LENIENT
    """
}