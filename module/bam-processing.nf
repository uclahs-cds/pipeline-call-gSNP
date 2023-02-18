include { generate_standard_filename } from '../external/nextflow-module/modules/common/generate_standardized_filename/main.nf'
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
        pattern: "*reheadered*"

    publishDir path: "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}-${id}-${interval_id}/log${file(it).getName()}" }

    input:
    tuple val(id), path(bam), path(interval)
    tuple val(index_id), path(bam_index)

    output:
    path(".command.*")
    tuple val(id), path(output_file_name), emit: bam_reheadered
    path(interval), emit: associated_interval
    path(bam), emit: bam_for_deletion
    path(bam_index), emit: bam_index_for_deletion

    script:
    // Get split interval number to serve as task ID
    interval_id = interval.baseName.split('-')[0]
    output_file_name = generate_standard_filename(
        params.aligner,
        params.dataset_id,
        id,
        [
            'additional_tools': ["GATK-${params.gatk_version}"],
            'additional_information': "recalibrated_reheadered_${interval_id}.bam"
        ]
    )
    id_generalized = id.replace('-', '[-,_]')
    """
    set -euo pipefail

    samtools reheader \
        -c 'sed "/^@RG/! s/.*/keep&/" | sed "/^@RG/ s/.*SM:${id_generalized}\$/keep&/" | sed "/^@RG/ s/.*SM:${id_generalized}\t/keep&/" | grep "^keep" | sed "s/^keep//g"' \
        ${bam} \
        > ${output_file_name}
    """
}

/*
    Nextflow module for removing duplicated identical records from interval processing.
    Due to parallelization via interval splitting, reads that overlap two intervals end up in
    the BAMs for both intervals. When merged, these records get duplicated, causing potential
    issues in downstream pipelines and analysis. Since records are sorted, uniq is able to de-deuplicate
    these records.

    See Issue #79 (https://github.com/uclahs-cds/pipeline-call-gSNP/issues/79)

    input:
        bam: path to input BAM
        id: ID of sample in BAM

    params:
        params.output_dir: string(path)
        params.log_output_dir: string(path)
        params.docker_image_samtools: string
*/
process deduplicate_records_SAMtools {
    container params.docker_image_samtools
    publishDir path: "${params.output_dir}/output",
        mode: "copy",
        pattern: "*merged-dedup*"

    publishDir path: "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}-${id}/log${file(it).getName()}" }

    input:
    path(bam)
    val(id)

    output:
    path(".command.*")
    tuple val(id), path(output_file_name), emit: merged_bam
    path(bam), emit: bam_for_deletion

    script:
    output_file_name = generate_standard_filename(
        params.aligner,
        params.dataset_id,
        id,
        [
            'additional_tools': ["GATK-${params.gatk_version}"],
            'additional_information': "realigned_recalibrated_merged_dedup.bam"
        ]
    )
    """
    samtools view \
        -h \
        ${bam} | \
        awk '(\$1 \$2 \$3 \$4 \$5 \$6 \$7 \$8 \$9 \$10 \$11)!=f_p && NR>1 {print f} {f=\$0} {f_p=(\$1 \$2 \$3 \$4 \$5 \$6 \$7 \$8 \$9 \$10 \$11)} END {print f}' | \
        samtools view \
        -b \
        -o ${output_file_name}
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
        enabled: params.save_intermediate_files && !params.is_dedup_bam,
        pattern: "*reheadered*"

    publishDir path: "${params.output_dir}/output",
        mode: "copy",
        enabled: params.is_dedup_bam,
        pattern: "*.bai"

    publishDir path: "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { params.is_dedup_bam ?
            "${task.process.replace(':', '/')}-${id}-dedup/log${file(it).getName()}" :
            "${task.process.replace(':', '/')}-${id}-${interval_id}/log${file(it).getName()}" }

    input:
    tuple val(id), path(bam)
    path(interval)

    output:
    path(".command.*")
    tuple path(bam), path("${bam}.bai"), path(interval), val(id), emit: indexed_out

    script:
    // Get split interval number to serve as task ID
    interval_id = interval.baseName.split('-')[0]
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
    publishDir path: "${params.output_dir}/intermediate/${task.process.replace(':', '/')}",
        mode: "copy",
        enabled: params.save_intermediate_files,
        pattern: "*merged*"

    publishDir path: "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}-${id}/log${file(it).getName()}" }

    input:
    path(bams)
    val(id)

    output:
    path(".command.*")
    path(output_file_name), emit: merged_bam
    val(id), emit: associated_id

    script:
    all_bams = bams.collect{ "-INPUT '$it'" }.join(' ')
    output_file_name = generate_standard_filename(
        params.aligner,
        params.dataset_id,
        id,
        [
            'additional_tools': ["GATK-${params.gatk_version}"],
            'additional_information': "realigned_recalibrated_merged.bam"
        ]
    )
    """
    set -euo pipefail
    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -Djava.io.tmpdir=${workDir} \
        -jar /usr/local/share/picard-slim-2.26.10-0/picard.jar MergeSamFiles \
        ${all_bams} \
        -OUTPUT ${output_file_name} \
        -SORT_ORDER coordinate \
        -ASSUME_SORTED false \
        -USE_THREADING true \
        -VALIDATION_STRINGENCY LENIENT
    """
}
