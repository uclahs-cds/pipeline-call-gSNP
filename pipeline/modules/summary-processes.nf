process run_GetPileupSummaries_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir}/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: '*.table'

    publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(bundle_contest_hapmap_3p3_vcf_gz)
    path(bundle_contest_hapmap_3p3_vcf_gz_tbi)
    path(all_intervals)
    path(normal_bam)
    path(normal_bam_index)
    path(tumour_bam)
    path(tumour_bam_index)
    tuple val(sample_id), val(normal_id), val(tumour_id)

    output:
    path(".command.*")
    path("${normal_id}_getpileupsummaries.table"), emit: normal_pileupsummaries
    path("${tumour_id}_getpileupsummaries.table"), emit: tumour_pileupsummaries optional true

    script:
    interval_options = all_intervals.collect{ "--intervals '$it'" }.join(' ')
    """
    set -euo pipefail
    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
        GetPileupSummaries \
        --input ${normal_bam} \
        --reference ${reference_fasta} \
        --variant ${bundle_contest_hapmap_3p3_vcf_gz} \
        ${interval_options} \
        --output ${normal_id}_getpileupsummaries.table


    if ${params.is_NT_paired}
    then
        gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
            GetPileupSummaries \
            --input ${tumour_bam} \
            --reference ${reference_fasta} \
            --variant ${bundle_contest_hapmap_3p3_vcf_gz} \
            ${interval_options} \
            --output ${tumour_id}_getpileupsummaries.table
    fi
    """
}

process run_CalculateContamination_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir}/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: '*.table'

    publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

    input:
    path(normal_pileupsummaries)
    path(tumour_pileupsummaries)
    tuple val(sample_id), val(normal_id), val(tumour_id)

    output:
    path(".command.*")
    path("${normal_id}_calculatecontamination_normal_alone.table"), emit: normal_contamination
    path("${tumour_id}_calculatecontamination_tumour_alone.table"), emit: tumour_contamination optional true
    path("${tumour_id}_calculatecontamination_tumour_with_matched_normal.table"), emit: tumour_normal_matched_contamination optional true

    script:
    """
    set -euo pipefail
    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
        CalculateContamination \
        --input ${normal_pileupsummaries} \
        --output ${normal_id}_calculatecontamination_normal_alone.table


    if ${params.is_NT_paired}
    then
        gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
            CalculateContamination \
            --input ${tumour_pileupsummaries} \
            --output ${tumour_id}_calculatecontamination_tumour_alone.table

        gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
            CalculateContamination \
            --input ${tumour_pileupsummaries} \
            --matched-normal ${normal_pileupsummaries} \
            --output ${tumour_id}_calculatecontamination_tumour_with_matched_normal.table
    fi
    """
}

process run_DepthOfCoverage_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir}/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: '*_DOC*'

    publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(all_intervals)
    path(normal_bam)
    path(normal_bam_index)
    path(tumour_bam)
    path(tumour_bam_index)
    tuple val(sample_id), val(normal_id), val(tumour_id)

    output:
    path(".command.*")
    path("*_DOC*")

    when:
    params.is_DOC_run

    script:
    interval_options = all_intervals.collect{ "--intervals '$it'" }.join(' ')
    """
    set -euo pipefail
    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
        DepthOfCoverage \
        --input ${normal_bam} \
        --output ${normal_id}_DOC \
        --output-format TABLE \
        --reference ${reference_fasta} \
        --omit-depth-output-at-each-base \
        --omit-interval-statistics \
        --omit-locus-table \
        --partition-type sample \
        --partition-type readgroup \
        --partition-type library \
        ${interval_options}


    if ${params.is_NT_paired}
    then
        gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
            DepthOfCoverage \
            --input ${tumour_bam} \
            --output ${tumour_id}_DOC \
            --output-format TABLE \
            --reference ${reference_fasta} \
            --omit-depth-output-at-each-base \
            --omit-interval-statistics \
            --omit-locus-table \
            --partition-type sample \
            --partition-type readgroup \
            --partition-type library \
            ${interval_options}
    fi
    """
}

workflow calculate_contamination {
    take:
    normal_bam
    normal_bam_index
    tumour_bam
    tumour_bam_index
    all_intervals
    identifiers

    main:
    run_GetPileupSummaries_GATK(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        params.reference_dict,
        params.bundle_contest_hapmap_3p3_vcf_gz,
        "${params.bundle_contest_hapmap_3p3_vcf_gz}.tbi",
        all_intervals.collect(),
        normal_bam,
        normal_bam_index,
        tumour_bam,
        tumour_bam_index,
        identifiers
        )

    run_CalculateContamination_GATK(
        run_GetPileupSummaries_GATK.out.normal_pileupsummaries,
        run_GetPileupSummaries_GATK.out.tumour_pileupsummaries.ifEmpty("/scratch/placeholder.txt"),
        identifiers
        )
}