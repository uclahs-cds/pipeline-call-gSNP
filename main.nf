#!/usr/bin/env nextflow

nextflow.enable.dsl=2


log.info """\
======================================
C A L L - G S N P  P I P E L I N E
======================================
Boutros Lab

Current Configuration:

    - pipeline:
        name: ${workflow.manifest.name}
        version: ${workflow.manifest.version}

    - input:
        samples: ${params.samples_to_process}
        intervals: ${params.intervals}
        bundle_v0_dbsnp138_vcf_gz: ${params.bundle_v0_dbsnp138_vcf_gz}
        bundle_mills_and_1000g_gold_standard_indels_vcf_gz: ${params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz}
        bundle_v0_dbsnp138_vcf_gz: ${params.bundle_v0_dbsnp138_vcf_gz}
        bundle_hapmap_3p3_vcf_gz: ${params.bundle_hapmap_3p3_vcf_gz}
        bundle_omni_1000g_2p5_vcf_gz: ${params.bundle_omni_1000g_2p5_vcf_gz}
        bundle_phase1_1000g_snps_high_conf_vcf_gz: ${params.bundle_phase1_1000g_snps_high_conf_vcf_gz}
        algorithms: ${params.algorithm}

    - output:
        output: ${params.output_dir}
        output_dir_base: ${params.output_dir_base}
        log_output_dir: ${params.log_output_dir}

    Tools Used:
        tool GATK: ${params.docker_image_gatk}
        tool PipeVal: ${params.docker_image_pipeval}
        tool Picard: ${params.docker_image_picard}
        tool GATK_filter: ${params.docker_image_gatkfilter}
        tool DeepVariant: ${params.docker_image_deepvariant}

    Extra parameters:
        ${params}

------------------------------------
Starting workflow...
------------------------------------
"""

include { run_validate_PipeVal } from './external/pipeline-Nextflow-module/modules/PipeVal/validate/main.nf'
include { run_SplitIntervals_GATK } from './module/split-intervals.nf'
include { extract_GenomeIntervals } from './external/pipeline-Nextflow-module/modules/common/extract_genome_intervals/main.nf'
include { deepvariant } from './module/workflow-deepvariant.nf' addParams(
    output_dir_base: "${params.output_dir_root}/DeepVariant-${params.deepvariant_version}",
    current_caller: "DeepVariant-${params.deepvariant_version}"
)
include { haplotypecaller } from './module/workflow-haplotypecaller.nf' addParams(
    current_caller: "GATK-${params.gatk_version}"
)

// Returns the index file for the given bam or vcf
def indexFile(bam_or_vcf) {
    if (bam_or_vcf.endsWith('.bam')) {
        return "${bam_or_vcf}.bai"
    } else if (bam_or_vcf.endsWith('vcf.gz')) {
        return "${bam_or_vcf}.tbi"
    } else {
        throw new Exception("Index file for ${bam_or_vcf} file type not supported. Use .bam or .vcf.gz files.")
    }
}

workflow {
    /**
    *   Input channel processing
    */
    Channel.from(params.samples_to_process)
        .map{ sample -> ['index': indexFile(sample.path)] + sample }
        .set{ input_ch_samples_with_index }

    input_ch_samples_with_index
        .map{ sample -> [sample.path, sample.index] }
        .flatten()
        .set{ input_ch_validate }

    input_ch_samples_with_index
        .reduce( ['bams': [], 'indices': []] ){ a, b ->
            a.bams.add(b.path);
            a.indices.add(b.index);
            return a
        }
        .set{ input_ch_collected_files }

    script_dir_ch = Channel.fromPath(
        "$projectDir/script",
        checkIfExists: true
        )
        .collect()

    base_meta = Channel.value([
        'log_output_dir': params.log_output_dir,
        'output_dir': params.output_dir_base
    ])

    /**
    *   Input validation
    */
    run_validate_PipeVal(
        base_meta.combine(input_ch_validate)
    )

    run_validate_PipeVal.out.validation_result
        .collectFile(
            name: 'input_validation.txt',
            storeDir: "${params.output_dir_base}/validation"
        )

    /**
    *   Handle interval splitting based on targeted or WGS mode
    */
    intervals_to_split = Channel.empty()

    if (params.is_targeted) {
        intervals_to_split = Channel.from(params.intervals)
    } else {
        extract_GenomeIntervals(
            base_meta.map{ metadata -> [metadata, "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict"]
        )

        intervals_to_split = extract_GenomeIntervals.out.genomic_intervals.map{ genome_intervals -> genome_intervals[1] }
    }

    run_SplitIntervals_GATK(
        intervals_to_split,
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict"
    )

    run_SplitIntervals_GATK.out.interval_list
        .flatten()
        .map{ interval_path ->
            [
                'interval_id': file(interval_path).getName().replace('-contig.interval_list', ''),
                'interval_path': interval_path
            ]
        }
        .set{ input_ch_intervals }

    /**
    *   DeepVariant
    */
    if (params.algorithm.contains('DeepVariant')) {
        deepvariant(
            input_ch_samples_with_index,
            input_ch_intervals
        )
    }

    /**
    *   HaplotypeCaller
    */
    if (params.algorithm.contains('HaplotypeCaller')) {
        haplotypecaller(
            input_ch_samples_with_index,
            input_ch_intervals
        )
    }
}
