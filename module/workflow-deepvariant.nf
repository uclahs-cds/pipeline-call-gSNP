include { convert_IntervalListToBed_GATK } from './convert-intervals.nf'
include { call_gSNP_DeepVariant } from './deepvariant.nf'
include { run_MergeVcfs_Picard } from './merge-vcf.nf'
include { calculate_sha512 } from './checksum.nf'

workflow deepvariant {
    take:
    samples
    intervals

    main:
    /**
    *   Convert intervals to BED format
    */
    intervals.map{ interval_list ->
        [interval_list.interval_id, interval_list.interval_path]
    }
    .set{ input_ch_intervallisttobed }

    convert_IntervalListToBed_GATK(
        input_ch_intervallisttobed
    )

    samples.combine(
        convert_IntervalListToBed_GATK.out.interval_bed.map{ bed ->
            [
                'interval_id': bed[0],
                'interval_path': bed[1]
            ]
        }
    )
    .map{ combined ->
        [
            combined[0].id,
            combined[0].path,
            combined[0].index,
            combined[1].interval_path,
            combined[1].interval_id
        ]
    }
    .set{ input_ch_deepvariant }

    /**
    *   Run DeepVariant
    */
    call_gSNP_DeepVariant(
        input_ch_deepvariant,
        params.reference_fasta,
        "${params.reference_fasta}.fai",
        "${file(params.reference_fasta).parent}/${file(params.reference_fasta).baseName}.dict",
        params.par_bed
    )

    /**
    *   Merge VCFs and GVCFs
    */
    call_gSNP_DeepVariant.out.vcf
        .filter{ raw_vcf_out -> (raw_vcf_out[3] != "0") } // Filter out empty VCFs
        .groupTuple(by: 0) // Group by sample
        .map{ vcf_group ->
            [
                vcf_group[1].flatten(), // VCFs
                vcf_group[2].flatten(), // Indices
                'VCF',
                vcf_group[0] // Sample ID
            ]
        }
        .set{ input_ch_merge_vcfs }

    call_gSNP_DeepVariant.out.gvcf
        .filter{ raw_gvcf_out -> (raw_gvcf_out[3] != "0") } // Filter out empty GVCFs
        .groupTuple(by: 0) // Group by sample
        .map{ gvcf_group ->
            [
                gvcf_group[1].flatten(), // VCFs
                gvcf_group[2].flatten(), // Indices
                'GVCF',
                gvcf_group[0] // Sample ID
            ]
        }
        .set{ input_ch_merge_gvcfs }

    run_MergeVcfs_Picard(
        input_ch_merge_vcfs.mix(input_ch_merge_gvcfs)
    )

    /**
    *   Compute checksums
    */
    run_MergeVcfs_Picard.out.merged_vcf
        .map{ merged -> [merged[1], merged[2]] }
        .flatten()
        .set{ input_ch_calculate_checksum }

    calculate_sha512(input_ch_calculate_checksum)
}
