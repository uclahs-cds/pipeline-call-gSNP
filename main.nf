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

    - output:
        output: ${params.output_dir}
        output_dir_base: ${params.output_dir_base}
        log_output_dir: ${params.log_output_dir}

    Tools Used:
        tool GATK: ${params.docker_image_gatk}
        tool PipeVal: ${params.docker_image_pipeval}
        tool Picard: ${params.docker_image_picard}
        tool GATK_filter: ${params.docker_image_gatkfilter}

    Extra parameters:
        ${params}

------------------------------------
Starting workflow...
------------------------------------
"""

include { run_validate_PipeVal } from './external/pipeline-Nextflow-module/modules/PipeVal/validate/main.nf' addParams(
    options: [
        docker_image_version: params.pipeval_version,
        main_process: "./" //Save logs in <log_dir>/process-log/run_validate_PipeVal
        ]
    )
include { run_SplitIntervals_GATK } from './module/split-intervals.nf'
include { extract_GenomeIntervals } from './external/pipeline-Nextflow-module/modules/common/extract_genome_intervals/main.nf' addParams(
    options: [
        save_intermediate_files: params.save_intermediate_files,
        output_dir: params.output_dir_base
        ]
    )
include {
    run_HaplotypeCallerGVCF_GATK
    } from './module/haplotypecaller.nf'
include { run_CombineGVCFs_GATK } from './module/combine-gvcfs.nf'
include { run_GenotypeGVCFs_GATK } from './module/genotype-gvcfs.nf'
include {
    run_MergeVcfs_Picard as run_MergeVcfs_Picard_VCF
    run_MergeVcfs_Picard as run_MergeVcfs_Picard_GVCF
    } from './module/merge-vcf.nf'
include { recalibrate_variants } from './module/workflow-recalibrate-variants.nf'
include { filter_gSNP_GATK } from './module/filter-gsnp.nf'
include { filter_XY } from './module/filter-xy.nf'
include { calculate_sha512 } from './module/checksum.nf'

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
    // TO-DO: Add validation for input BAMs

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

    /**
    *   Input validation
    */
    run_validate_PipeVal(input_ch_validate)

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
        extract_GenomeIntervals("${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict")

        intervals_to_split = extract_GenomeIntervals.out.genomic_intervals
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
    *   Haplotype calling
    */

    input_ch_samples_with_index.combine(input_ch_intervals)
        .map{ it ->
            [
                it[0].id,
                it[0].path,
                it[0].index,
                it[1].interval_path,
                it[1].interval_id
            ]
        }
        .set{ input_ch_haplotypecallergvcf }

    run_HaplotypeCallerGVCF_GATK(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        params.bundle_v0_dbsnp138_vcf_gz,
        "${params.bundle_v0_dbsnp138_vcf_gz}.tbi",
        input_ch_haplotypecallergvcf
    )

    run_HaplotypeCallerGVCF_GATK.out.gvcfs
        .groupTuple(by: 4) // Group by interval ID
        .map{ it ->
            [
                it[1].flatten(), // GVCFs
                it[2].flatten(), // Indices
                it[3][0], // Interval path
                it[4] // Interval ID
            ]
        }
    .set { input_ch_combine_gvcfs }

    run_CombineGVCFs_GATK(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        input_ch_combine_gvcfs
        )

    run_GenotypeGVCFs_GATK(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        params.bundle_v0_dbsnp138_vcf_gz,
        "${params.bundle_v0_dbsnp138_vcf_gz}.tbi",
        run_CombineGVCFs_GATK.out.combined_gvcf
    )

    /**
    *   Merge VCFs
    */
    run_GenotypeGVCFs_GATK.out.vcfs
        .reduce( ['vcfs': [], 'indices': []] ){ a, b ->
            a.vcfs.add(b[0]);
            a.indices.add(b[1]);
            return a
        }
        .map{ it ->
            [
                it.vcfs,
                it.indices,
                'VCF',
                params.patient_id
            ]
        }
        .set{ input_ch_merge_vcfs }

    run_MergeVcfs_Picard_VCF(input_ch_merge_vcfs)

    run_HaplotypeCallerGVCF_GATK.out.gvcfs
        .groupTuple(by: 0) // Group by sample
        .map{ it ->
            [
                it[1].flatten(), // GVCFs
                it[2].flatten(), // Indices
                'GVCF',
                it[0] // Sample ID
            ]
        }
        .set{ input_ch_merge_gvcfs }

    run_MergeVcfs_Picard_GVCF(input_ch_merge_gvcfs)

    /**
    *   Recalibrate variants
    */
    recalibrate_variants(run_MergeVcfs_Picard_VCF.out.merged_vcf)

    /**
    *   Filter variants with Perl script
    */
    filter_gSNP_GATK(
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        recalibrate_variants.out.output_ch_recalibrated_variants
    )

    filter_xy_ch = recalibrate_variants.out.output_ch_recalibrated_variants
        .map { it -> [it[0], it[1], it[2]] }

    script_dir_ch = Channel.fromPath(
        "$projectDir/script",
        checkIfExists: true
        )
        .collect()

    filter_XY(
        filter_xy_ch,
        script_dir_ch
        )
    /**
    *   Calculate checksums for output files
    */
    run_MergeVcfs_Picard_VCF.out.merged_vcf
        .mix(run_MergeVcfs_Picard_GVCF.out.merged_vcf)
        .mix(recalibrate_variants.out.output_ch_recalibrated_variants)
        .map{ [it[1], it[2]] }
        .mix(filter_gSNP_GATK.out.germline_filtered)
        .flatten()
        .set{ input_ch_calculate_checksum }

    calculate_sha512(input_ch_calculate_checksum)
}
