include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

/*
    Nextflow module for calling haplotypes in GVCF mode

    input:
        reference_fasta: path to reference genome fasta file
        reference_fasta_fai: path to index for reference fasta
        reference_fasta_dict: path to dictionary for reference fasta
        dbsnp_bundle: path to dbSNP variants
        dbsnp_bundle_index: path to index of dbSNP variants
        sample_id:  sample ID
        bam: path to BAM for calling
        bam_index: path to index of BAM
        interval: path to specific intervals for calling
        interval_id: interval ID
        
    params:
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_gatk: string
        params.is_targeted: bool. Indicator of whether in targeted exome mode or in WGS mode
        params.gatk_command_mem_diff: float(memory)
*/
process run_HaplotypeCallerGVCF_GATK {
    container params.docker_image_gatk
    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
      mode: "copy",
      enabled: params.save_intermediate_files,
      pattern: '*.vcf*'

    publishDir path: "${params.log_output_dir}/process-log",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process.replace(':', '/')}/${task.process.split(':')[-1]}-${sample_id}-${interval_id}/log${file(it).getName()}" }

    input:
    path(reference_fasta)
    path(reference_fasta_fai)
    path(reference_fasta_dict)
    path(dbsnp_bundle)
    path(dbsnp_bundle_index)
    tuple val(sample_id), path(bam), path(bam_index), path(interval_path), val(interval_id)


    output:
    path(".command.*")
    tuple val(sample_id), path(output_filename), path("${output_filename}.tbi"), path(interval_path), val(interval_id), emit: gvcfs

    script:
    output_filename = generate_standard_filename(
        "GATK-${params.gatk_version}",
        params.dataset_id,
        sample_id,
        [
            'additional_information': "${interval_id}_raw_variants.g.vcf.gz"
        ]
    )
    interval_str = "--intervals ${interval_path}"
    interval_padding = params.is_targeted ? "--interval-padding 100" : ""
    output_mode = params.emit_all_confident_sites ? "EMIT_ALL_CONFIDENT_SITES" : "EMIT_VARIANTS_ONLY"
    """
    set -euo pipefail

    gatk --java-options "-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=${workDir}" \
        HaplotypeCaller \
        --input ${bam} \
        --output ${output_filename} \
        --reference ${reference_fasta} \
        --verbosity INFO \
        --output-mode ${output_mode} \
        --emit-ref-confidence GVCF \
        --dbsnp ${dbsnp_bundle} \
        --sample-ploidy 2 \
        ${interval_str} \
        ${interval_padding}
    """
}
