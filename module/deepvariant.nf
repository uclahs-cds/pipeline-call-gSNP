include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

/*
    Nextflow module for running DeepVariant

    input:
        sample_id: sample ID
        bam: path to sample BAM
        bam_index: path to BAM index
        intervals: path to BED format intervals
        interval_id: interval ID

    params:
        params.output_dir_base: string(path)
        params.log_output_dir: string(path)
        params.save_intermediate_files: bool.
        params.docker_image_deepvariant: string
*/
process call_gSNP_DeepVariant {
    container params.docker_image_deepvariant

    tag "${sample_id}-${interval_id}"

    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "*.vcf.gz*",
               enabled: params.save_intermediate_files
    publishDir "${params.log_output_dir}/process-log",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
    tuple val(sample_id), path(bam), path(bam_index), path(intervals), val(interval_id)
    path(reference_fasta)
    path(reference_index)
    path(reference_dict)
    path(par_regions)

    output:
    tuple val(sample_id), path(vcf_filename), path("${vcf_filename}.tbi"), emit: vcf
    tuple val(sample_id), path(gvcf_filename), path("${gvcf_filename}.tbi"), emit: gvcf
    path ".command.*"

    script:
    output_filename_base = generate_standard_filename(
        "DeepVariant-${params.deepvariant_version}",
        params.dataset_id,
        sample_id,
        [
            'additional_information': interval_id
        ]
    )
    vcf_filename = "${output_filename_base}.vcf.gz"
    gvcf_filename = "${output_filename_base}.g.vcf.gz"
    model_type = (params.is_targeted) ? "WES" : "WGS"
    if (params.genetic_sex == 'XY') {
        haploid_args_raw = "--haploid_contigs=PREFIXX,PREFIXY --par_regions_bed=${par_regions}"
    } else {
        haploid_args_raw = ""
    }
    """
    set -euo pipefail

    mkdir log
    mkdir work

    if grep -q \$'\\tSN:chr1\\t' ${reference_dict}
    then
        haploid_args=`echo ${haploid_args_raw} | sed 's:PREFIX:chr:g'`
    else
        haploid_args=`echo ${haploid_args_raw} | sed 's:PREFIX::g'`
    fi

    /opt/deepvariant/bin/run_deepvariant \
        --model_type=${model_type} \
        --ref=${reference_faste} \
        --reads=${bam} \
        --output_vcf=${vcf_filename} \
        --output_gvcf=${gvcf_filename} \
        --num_shards=${task.cpus} \
        --logging_dir=log \
        --intermediate_results_dir=work \
        --regions=${intervals} \
        \$haploid_args
    """
}
