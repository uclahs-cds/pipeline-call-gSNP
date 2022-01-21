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
    tuple val(id), path(bam), path(interval)
    tuple val(index_id), path(bam_index)

    output:
    path(".command.*")
    tuple val(id), path("${id}_recalibrated_reheadered_${task.index}.bam"), emit: bam_reheadered
    path(interval), emit: associated_interval
    path(bam), emit: bam_for_deletion
    path(bam_index), emit: bam_index_for_deletion

    script:
    """
    set -euo pipefail

    samtools reheader \
        -c 'sed "/^@RG/! s/.*/keep&/" | sed "/^@RG/ s/.*SM:${id}\$/keep&/" | sed "/^@RG/ s/.*SM:${id}\t/keep&/" | grep "^keep" | sed "s/^keep//g"' \
        ${bam} \
        > ${id}_recalibrated_reheadered_${task.index}.bam
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
        params.docker_image_samtools: string
*/
process run_index_SAMtools {
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
    tuple val(id), path(bam)
    path(interval)

    output:
    path(".command.*")
    tuple path(bam), path("${bam}.bai"), path(interval), val(id), emit: indexed_out

    script:
    """
    set -euo pipefail
    
    samtools index ${bam}
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
        saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path(bams)
    val(id)

    output:
    path(".command.*")
    path("${id}_realigned_recalibrated_merged.bam"), emit: merged_bam
    path("${id}_realigned_recalibrated_merged.bai"), emit: merged_bam_index
    val(id), emit: associated_id

    script:
    all_bams = bams.collect{ "-INPUT '$it'" }.join(' ')
    """
    set -euo pipefail
    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -Djava.io.tmpdir=/scratch \
        -jar /usr/local/share/picard-slim-2.26.8-0/picard.jar MergeSamFiles \
        ${all_bams} \
        -OUTPUT ${id}_realigned_recalibrated_merged.bam \
        -CREATE_INDEX true \
        -SORT_ORDER coordinate \
        -ASSUME_SORTED false \
        -USE_THREADING true \
        -VALIDATION_STRINGENCY LENIENT
    """
}