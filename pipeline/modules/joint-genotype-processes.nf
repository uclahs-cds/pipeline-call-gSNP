process run_GenomicsDBImport_GATK {
    container params.docker_image_gatk

    publishDir path: params.output_dir,
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: "*_db"

    publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    tuple(path(interval), path(gvcfs), path(gvcf_indicies))

    output:
    path(".command.*")
    tuple path(interval), path("gvcfs_${task.index}_db"), emit: gdb

    script:
    variants = gvcfs.join(' --variant ')

    // --tmp-dir was added to help resolve potential memory issues
    // https://gatk.broadinstitute.org/hc/en-us/community/posts/360074224671-GenomicsDBImport-running-out-of-memory-

    """
    set -euo pipefail
    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
           GenomicsDBImport \
           --genomicsdb-workspace-path gvcfs_${task.index}_db \
           --intervals ${interval} \
           --reference ${reference_fasta} \
           --variant ${variants} \
           --max-num-intervals-to-import-in-parallel ${params.max_num_intervals_to_import_in_parallel} \
           --reader-threads ${params.reader_threads} \
           --batch-size ${params.batch_size} \
           ${params.genomics_db_import_extra_args} \
           --tmp-dir \$PWD
    """
}

process run_SplitIntervals_GATK {
    container params.docker_image_gatk

    publishDir params.output_dir,
               mode: "copy",
               pattern: "interval-files/*-scattered.interval_list",
               enabled: params.save_intermediate_files
    publishDir params.log_output_dir,
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

    input:
    path intervals
    path reference
    path reference_index
    path reference_dict

    output:
    path 'interval-files/*-scattered.interval_list', emit: interval_list
    path ".command.*"

    """
    set -euo pipefail
    gatk SplitIntervals \
        -R $reference \
        -L $intervals \
        --scatter-count ${params.scatter_count} \
        ${params.split_intervals_extra_args} \
        -O interval-files
    """
}

process run_HaplotypeCaller_GATK {
    container params.docker_image_gatk
    publishDir path:params.output_dir,
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: '*.vcf*'

    publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(dbsnp_bundle)
    path(dbsnp_bundle_index)
    tuple val(sample_id), val(normal_id), val(tumour_id), path(bam), path(bam_index), path(bam_tumour), path(bam_index_tumour), path(interval)


    output:
      path(".command.*")
      tuple val(interval), path("${normal_id}_{${task.index}_,}raw_variants.g.vcf.gz"), path("${normal_id}_{${task.index}_,}raw_variants.g.vcf.gz.tbi"), emit: gvcf_normal
      tuple val(interval), path("${tumour_id}_{${task.index}_,}raw_variants.g.vcf.gz"), path("${tumour_id}_{${task.index}_,}raw_variants.g.vcf.gz.tbi"), emit: gvcf_tumour optional true
      tuple val(interval), path("${sample_id}_{${task.index}_,}.vcf"), path("${sample_id}_{${task.index}_,}.vcf.idx"), emit: vcf

    script:
        out_filename_normal = "${normal_id}_${task.index}_raw_variants.g.vcf.gz"
        out_filename_tumour = "${tumour_id}_${task.index}_raw_variants.g.vcf.gz"
        out_filename_vcf = "${sample_id}_${task.index}.vcf"
        interval_str = "--intervals ${interval}"
        vcf_input_str = params.is_NT_paired ? "--input ${bam} --input ${bam_tumour}" : "--input ${bam}"

    """
    set -euo pipefail

    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
        HaplotypeCaller \
        ${vcf_input_str} \
        --output ${out_filename_vcf} \
        --reference ${reference_fasta} \
        --verbosity INFO \
        --output-mode EMIT_VARIANTS_ONLY \
        --dbsnp ${dbsnp_bundle} \
        --sample-ploidy 2 \
        --standard-min-confidence-threshold-for-calling 50 \
        ${interval_str}

    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
        HaplotypeCaller \
        --input ${bam} \
        --output ${out_filename_normal} \
        --reference ${reference_fasta} \
        --verbosity INFO \
        --output-mode EMIT_VARIANTS_ONLY \
        --emit-ref-confidence GVCF \
        --dbsnp ${dbsnp_bundle} \
        --sample-ploidy 2 \
        --standard-min-confidence-threshold-for-calling 30 \
        ${interval_str}

    if ${params.is_NT_paired}
    then
      gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
        HaplotypeCaller \
        --input ${bam_tumour} \
        --output ${out_filename_tumour} \
        --reference ${reference_fasta} \
        --verbosity INFO \
        --output-mode EMIT_VARIANTS_ONLY \
        --emit-ref-confidence GVCF \
        --dbsnp ${dbsnp_bundle} \
        --sample-ploidy 2 \
        --standard-min-confidence-threshold-for-calling 30 \
        ${interval_str}
    fi
    """
}

process run_GenotypeGVCFs_GATK {
    container params.docker_image_gatk

    publishDir path:params.output_dir,
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: '*.vcf*'

    publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(dbsnp_bundle)
    path(dbsnp_bundle_index)
    tuple(path(interval), path(gendb))

    output:
    path(".command.*")
    path("regenotype_cohort_joint_genotyping_${task.index}.vcf.gz"), emit: vcf
    path("regenotype_cohort_joint_genotyping_${task.index}.vcf.gz.tbi"), emit: vcf_index

    script:
    """
    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
        GenotypeGVCFs\
        --reference ${reference_fasta} \
        --variant gendb://${gendb} \
        --dbsnp ${dbsnp_bundle} \
        --output regenotype_cohort_joint_genotyping_${task.index}.vcf.gz
    """
}

process run_SortVcf_GATK {
    container params.docker_image_gatk

    publishDir path: params.output_dir,
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: "sorted_regenotype_cohort_joint_genotyping_*.vcf.gz{,.tbi}"

    publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

    input:
    path(vcf)

    output:
    path(".command.*")
    path("sorted_${vcf}"), emit: vcf
    path("sorted_${vcf}.tbi"), emit: vcf_index

    script:
    """
    set -euo pipefail
    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/scratch" \
	    SortVcf \
      --INPUT ${vcf} \
      --OUTPUT sorted_${vcf}
    """
}

process run_MergeVcfs_Picard {
    container params.docker_image_picard

    publishDir path: params.output_dir,
      mode: "copy",
      pattern: "*.vcf*"

    publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

    input:
    path(vcfs)
    path(gvcfs_normal)
    path(gvcfs_tumour)
    tuple val(sample_id), val(normal_id), val(tumour_id)

    output:
    path(".command.*")
    tuple path("${sample_id}_merged_raw.vcf"), path("${sample_id}_merged_raw.vcf.idx"), emit: vcf
    tuple path("${normal_id}_merged_raw_variants.g.vcf.gz"), path("${normal_id}_merged_raw_variants.g.vcf.gz.tbi"), emit: gvcf_normal
    tuple path("${tumour_id}_merged_raw_variants.g.vcf.gz"), path("${tumour_id}_merged_raw_variants.g.vcf.gz.tbi"), emit: gvcf_tumour optional true

    script:
    vcf_args = vcfs.collect{ "-INPUT '$it'" }.join(' ')
    gvcf_args_normal = gvcfs_normal.collect{ "-INPUT '$it'" }.join(' ')
    gvcf_args_tumour = gvcfs_tumour.collect{ "-INPUT '$it'" }.join(' ')

    """
    set -euo pipefail

    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -Djava.io.tmpdir=/scratch \
      -jar /picard-tools/picard.jar MergeVcfs \
      ${vcf_args} \
      -OUTPUT ${sample_id}_merged_raw.vcf \
      -SEQUENCE_DICTIONARY ${params.reference_dict} \
      -VALIDATION_STRINGENCY LENIENT

    java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -Djava.io.tmpdir=/scratch \
      -jar /picard-tools/picard.jar MergeVcfs \
      ${gvcfs_args_normal} \
      -OUTPUT ${normal_id}_merged_raw_variants.g.vcf.gz \
      -SEQUENCE_DICTIONARY ${params.reference_dict} \
      -VALIDATION_STRINGENCY LENIENT
    
    if ${params.is_NT_paired}
    then
      java -Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -Djava.io.tmpdir=/scratch \
        -jar /picard-tools/picard.jar MergeVcfs \
        ${gvcfs_args_tumour} \
        -OUTPUT ${tumour_id}_merged_raw_variants.g.vcf.gz \
        -SEQUENCE_DICTIONARY ${params.reference_dict} \
        -VALIDATION_STRINGENCY LENIENT
    fi
    """
}
