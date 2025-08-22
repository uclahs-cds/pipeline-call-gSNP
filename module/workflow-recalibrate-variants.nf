include {
    run_VariantRecalibratorSNP_GATK
    run_VariantRecalibratorINDEL_GATK
    run_ApplyVQSR_GATK as run_ApplyVQSR_GATK_SNP
    run_ApplyVQSR_GATK as run_ApplyVQSR_GATK_INDEL
} from './variant-recalibrator.nf'

workflow recalibrate_variants {
    take:
    workflow_meta
    input_ch_samples

    main:
    run_VariantRecalibratorSNP_GATK(
        workflow_meta,
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        params.bundle_v0_dbsnp138_vcf_gz,
        "${params.bundle_v0_dbsnp138_vcf_gz}.tbi",
        params.bundle_hapmap_3p3_vcf_gz,
        "${params.bundle_hapmap_3p3_vcf_gz}.tbi",
        params.bundle_omni_1000g_2p5_vcf_gz,
        "${params.bundle_omni_1000g_2p5_vcf_gz}.tbi",
        params.bundle_phase1_1000g_snps_high_conf_vcf_gz,
        "${params.bundle_phase1_1000g_snps_high_conf_vcf_gz}.tbi",
        input_ch_samples
    )

    run_ApplyVQSR_GATK_SNP(
        workflow_meta,
        'SNP',
        'SNP',
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        run_VariantRecalibratorSNP_GATK.out.snp_recal
    )

    run_VariantRecalibratorINDEL_GATK(
        workflow_meta,
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz,
        "${params.bundle_mills_and_1000g_gold_standard_indels_vcf_gz}.tbi",
        run_ApplyVQSR_GATK_SNP.out.output_ch_vqsr
    )

    run_ApplyVQSR_GATK_INDEL(
        workflow_meta,
        'INDEL',
        'SNP_AND_INDEL',
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        run_VariantRecalibratorINDEL_GATK.out.indel_recal
    )

    emit:
    output_ch_recalibrated_variants = run_ApplyVQSR_GATK_INDEL.out.output_ch_vqsr
}
