nextflow.enable.dsl=2

include {
    deduplicate_records_SAMtools
    run_index_SAMtools
    run_MergeSamFiles_Picard
    } from './bam-processing.nf' addParams(is_dedup_bam: true)

workflow merge_and_deduplicate {
    take:
    bams
    id

    main:
    run_MergeSamFiles_Picard(
        bams,
        id
        )

    deduplicate_records_SAMtools(
        run_MergeSamFiles_Picard.out.merged_bam,
        run_MergeSamFiles_Picard.out.associated_id
        )

    run_index_SAMtools(
        deduplicate_records_SAMtools.out.merged_bam,
        params.bundle_omni_1000g_2p5_vcf_gz // Decoy for interval, needs to be an actual file since it's also an output of the process
        )

    run_index_SAMtools.out.indexed_out
        .multiMap{ it ->
            bam: it[0]
            bai: it[1]
            id: it[3]
            }
        .set{ merged_dedup }

    emit:
    merged_bam = merged_dedup.bam
    merged_bam_index = merged_dedup.bai
    associated_id = merged_dedup.id
}