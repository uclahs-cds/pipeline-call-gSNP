include {
    run_HaplotypeCallerGVCF_GATK
    } from './haplotypecaller.nf'
include { run_CombineGVCFs_GATK } from './combine-gvcfs.nf'
include { run_GenotypeGVCFs_GATK } from './genotype-gvcfs.nf'
include {
    run_MergeVcfs_Picard as run_MergeVcfs_Picard_VCF
    run_MergeVcfs_Picard as run_MergeVcfs_Picard_GVCF
    } from './merge-vcf.nf'
include { recalibrate_variants } from './workflow-recalibrate-variants.nf'
include { filter_gSNP_GATK } from './filter-gsnp.nf'
include { filter_XY_Hail } from './filter-xy.nf'
include { calculate_sha512 } from './checksum.nf'

workflow haplotypecaller {
    take:
    samples
    intervals

    main:
    /**
    *   Haplotype calling
    */
    samples.combine(intervals)
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

    filter_xy_output_ch = Channel.empty()
    if (params.genetic_sex != 'unknown') {
        filter_xy_ch = recalibrate_variants.out.output_ch_recalibrated_variants
            .map { it -> [it[0], it[1], it[2]] }

        script_dir_ch = Channel.fromPath(
            "$projectDir/script",
            checkIfExists: true
            )
            .collect()

        filter_XY_Hail(
            filter_xy_ch,
            params.reference_fasta,
            "${params.reference_fasta}.fai",
            params.par_bed,
            script_dir_ch
            )
        filter_xy_output_ch = filter_xy_output_ch.mix(filter_XY_Hail.out.xy_filtered_vqsr)
        }
    /**
    *   Calculate checksums for output files
    */
    run_MergeVcfs_Picard_VCF.out.merged_vcf
        .mix(run_MergeVcfs_Picard_GVCF.out.merged_vcf)
        .mix(recalibrate_variants.out.output_ch_recalibrated_variants)
        .map{ [it[1], it[2]] }
        .mix(filter_xy_output_ch)
        .mix(filter_gSNP_GATK.out.germline_filtered)
        .flatten()
        .set{ input_ch_calculate_checksum }

    calculate_sha512(input_ch_calculate_checksum)
}
