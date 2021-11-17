nextflow.enable.dsl=2

include { run_reheader_SAMtools as run_reheader_SAMtools_normal; run_reheader_SAMtools as run_reheader_SAMtools_tumour } from './bam-processing.nf'
include { run_BuildBamIndex_Picard as run_BuildBamIndex_Picard_normal; run_BuildBamIndex_Picard as run_BuildBamIndex_Picard_tumour } from './bam-processing.nf' addParams(
    is_output_bam: false
    )

workflow reheader_interval_bams {
    take:
    identifiers
    normal_bams
    normal_bams_index
    tumour_bams
    tumour_bams_index
    intervals

    main:
    run_reheader_SAMtools_normal(
        identifiers,
        normal_bams,
        normal_bams_index,
        intervals,
        'normal'
        )

    run_reheader_SAMtools_tumour(
        identifiers,
        tumour_bams,
        tumour_bams_index,
        intervals,
        'tumour'
    )

    run_BuildBamIndex_Picard_normal(
        run_reheader_SAMtools_normal.out.bam_reheadered,
        run_reheader_SAMtools_normal.out.associated_interval
        )
    
    run_BuildBamIndex_Picard_tumour(
        run_reheader_SAMtools_tumour.out.bam_reheadered,
        run_reheader_SAMtools_tumour.out.associated_interval
        )

    // Merge the reheadered, indexed output based on the associated interval
    normal_for_match = run_BuildBamIndex_Picard_normal.out.indexed_out
        .map{ it ->
            [it[0], it[1], it[2], it[2].getFileName()]
            }

    tumour_for_match = run_BuildBamIndex_Picard_tumour.out.indexed_out
        .map{ it ->
            [it[0], it[1], it[2], it[2].getFileName()]
            }

    normal_for_match
        .join(tumour_for_match, by: 3)
        .multiMap { it ->
            associated_interval: it[3]
            normal_bam_reheadered: it[1]
            normal_bam_reheadered_index: it[2]
            tumour_bam_reheadered: it[4]
            tumour_bam_reheadered_index: it[5]
            }
        .set { matched_output_channel }

    emit:
    reheadered_normal_bam = matched_output_channel.normal_bam_reheadered
    reheadered_normal_bam_index = matched_output_channel.normal_bam_reheadered_index
    reheadered_tumour_bam = matched_output_channel.tumour_bam_reheadered
    reheadered_tumour_bam_index = matched_output_channel.tumour_bam_reheadered_index
    associated_interval = matched_output_channel.associated_interval
    normal_bam_for_deletion = run_reheader_SAMtools_normal.out.bam_for_deletion
    tumour_bam_for_deletion = run_reheader_SAMtools_tumour.out.bam_for_deletion
    normal_bam_index_for_deletion = run_reheader_SAMtools_normal.out.bam_index_for_deletion
    tumour_bam_index_for_deletion = run_reheader_SAMtools_tumour.out.bam_index_for_deletion
}
