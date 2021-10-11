/*
    Nextflow module for getting pileup summaries of BAMs

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        bundle_contest_hapmap_3p3_vcf_gz: path to contamination estimate variants
        bundle_contest_hapmap_3p3_vcf_gz_tbi: path to index of contamination estimate VCFs
        all_intervals: path to set of full target intervals
        normal_bam: path to normal BAM
        normal_bam_index: path to normal BAM index
        tumour_bam: path to tumour BAM
        tumour_bam_index: path to tumour BAM index
        (sample_id, normal_id, tumour_id): tuples of string identifiers for the samples
        
    params:
        params.output_dir: string(path)
        params.log_output_dir: string(path)
        params.docker_image_gatk: string
        params.is_NT_paired: bool. Indicator of whether input has normal and tumour samples
        params.gatk_command_mem_diff: float(memory)
*/
process run_GetPileupSummaries_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir}/QC/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: '*.table'

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

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

/*
    Nextflow module for calculating contamination

    input:
        normal_pileupsummaries: path to pileup summary for normal BAM
        tumour_pileupsummaries: path to pileup summary for tumour BAM
        (sample_id, normal_id, tumour_id): tuples of string identifiers for the samples
        
    params:
        params.output_dir: string(path)
        params.log_output_dir: string(path)
        params.docker_image_gatk: string
        params.is_NT_paired: bool. Indicator of whether input has normal and tumour samples
        params.gatk_command_mem_diff: float(memory)
*/
process run_CalculateContamination_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir}/QC/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: '*.table'

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

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

/*
    Nextflow module for calculating depth of coverage

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        all_intervals: path to set of full target intervals
        normal_bam: path to normal BAM
        normal_bam_index: path to normal BAM index
        tumour_bam: path to tumour BAM
        tumour_bam_index: path to tumour BAM index
        (sample_id, normal_id, tumour_id): tuples of string identifiers for the samples
        
    params:
        params.output_dir: string(path)
        params.log_output_dir: string(path)
        params.docker_image_gatk: string
        params.is_NT_paired: bool. Indicator of whether input has normal and tumour samples
        params.gatk_command_mem_diff: float(memory)
        params.is_DOC_run: bool. Indicator of whether to run depth of coverage process
*/
process run_DepthOfCoverage_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir}/QC/${task.process.replace(':', '/')}",
      mode: "copy",
      pattern: '*_DOC*'

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

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
