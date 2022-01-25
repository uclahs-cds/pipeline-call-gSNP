/*
    Nextflow module for getting pileup summaries of BAMs

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        bundle_contest_hapmap_3p3_vcf_gz: path to contamination estimate variants
        bundle_contest_hapmap_3p3_vcf_gz_tbi: path to index of contamination estimate VCFs
        all_intervals: path to set of full target intervals
        bam: path to normal BAM
        bam_index: path to normal BAM index
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
    val(sample_type)
    path(bam)
    path(bam_index)
    val(id)

    output:
    path(".command.*")
    path("*_getpileupsummaries.table"), emit: pileupsummaries

    script:
    interval_options = all_intervals.collect{ "--intervals '$it'" }.join(' ')
    output_filename = "${id}_getpileupsummaries.table"
    """
    set -euo pipefail
    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
        GetPileupSummaries \
        --input ${bam} \
        --reference ${reference_fasta} \
        --variant ${bundle_contest_hapmap_3p3_vcf_gz} \
        ${interval_options} \
        --output ${output_filename}
    """
}

/*
    Nextflow module for calculating contamination

    input:
        sample_type: string to indicate whether processing normal or tumour sample
        matched_normal_pileupsummaries: path to matched pileup summary for normal BAM
        pileupsummaries: path to pileup summary for BAM
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
    val(sample_type)
    path(matched_normal_pileupsummaries)
    path(pileupsummaries)
    val(id)

    output:
    path(".command.*")
    path("*_alone.table"), emit: contamination
    path("*_with_matched_normal.table"), emit: tumour_normal_matched_contamination optional true

    script:
    single_output_filename = "${id}_calculatecontamination_${sample_type}"
    calc_matched = (sample_type == "normal") ? false : true
    """
    set -euo pipefail
    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
        CalculateContamination \
        --input ${pileupsummaries} \
        --output ${single_output_filename}_alone.table

    if ${calc_matched}
    then
      gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
          CalculateContamination \
          --input ${pileupsummaries} \
          --matched-normal ${matched_normal_pileupsummaries} \
          --output ${single_output_filename}_with_matched_normal.table
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
        bam: path to normal BAM
        bam_index: path to normal BAM index
        sample_type: string to indicate whether processing normal or tumour sample
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
    path(bam)
    path(bam_index)
    val(sample_type)
    val(id)

    output:
    path(".command.*")
    path("*_DOC*")

    when:
    params.is_DOC_run

    script:
    interval_options = params.is_targeted ? "--intervals ${params.intervals}" : all_intervals.collect{ "--intervals '$it'" }.join(' ')
    """
    set -euo pipefail
    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
        DepthOfCoverage \
        --input ${bam} \
        --output ${id}_DOC \
        --output-format TABLE \
        --reference ${reference_fasta} \
        --omit-depth-output-at-each-base \
        --omit-interval-statistics \
        --omit-locus-table \
        --partition-type sample \
        --partition-type readgroup \
        --partition-type library \
        ${interval_options}
    """
}

workflow calculate_contamination_normal {
    take:
    bam
    bam_index
    all_intervals
    id

    main:
    run_GetPileupSummaries_GATK(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        params.bundle_contest_hapmap_3p3_vcf_gz,
        "${params.bundle_contest_hapmap_3p3_vcf_gz}.tbi",
        all_intervals.collect(),
        "normal",
        bam,
        bam_index,
        id
        )

    run_CalculateContamination_GATK(
        "normal",
        "/scratch/summaryplaceholder.txt", // Decoy since processing normal sample
        run_GetPileupSummaries_GATK.out.pileupsummaries,
        id
        )
    
    emit:
    pileupsummaries = run_GetPileupSummaries_GATK.out.pileupsummaries
}

workflow calculate_contamination_tumour {
    take:
    bam
    bam_index
    all_intervals
    id
    normal_pileupsummaries

    main:
    bam.map{ it ->
        "tumour"
        }
        .set{ type_ich }

    normal_pileupsummaries.combine(bam)
        .map{it ->
            it[0]
            }
        .set{ normal_summaries_ich }

    run_GetPileupSummaries_GATK(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        params.bundle_contest_hapmap_3p3_vcf_gz,
        "${params.bundle_contest_hapmap_3p3_vcf_gz}.tbi",
        all_intervals.collect(),
        type_ich,
        bam,
        bam_index,
        id
        )

    run_CalculateContamination_GATK(
        type_ich,
        normal_summaries_ich,
        run_GetPileupSummaries_GATK.out.pileupsummaries,
        id
        )
}
