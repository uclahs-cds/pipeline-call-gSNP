nextflow.enable.dsl=2

include { run_reheader_SAMtools as run_reheader_SAMtools_normal; run_reheader_SAMtools as run_reheader_SAMtools_tumour } from './bam-processing.nf'
include { run_index_SAMtools as run_index_SAMtools_normal; run_index_SAMtools as run_index_SAMtools_tumour } from './bam-processing.nf'
include {
    flatten_samples as flatten_samples_bam_index
    flatten_samples as flatten_samples_bam
    } from './functions.nf'

workflow reheader_interval_bams {
    take:
    normal_bams
    normal_bams_index
    tumour_bams
    tumour_bams_index

    main:
    run_reheader_SAMtools_normal(
        normal_bams,
        normal_bams_index,
        )

    // Partially flatten the tumour inputs
    // Flatten by one level, keeping the associations between id, BAM, and interval intact
    flatten_samples_bam(tumour_bams)
    
    flatten_samples_bam.out.och
        .map{ it ->
            it[0,1,2]
            }
        .set{ tumour_bams_flattened_ich }

    // Repeat for BAIs
    flatten_samples_bam_index(tumour_bams_index)
    
    flatten_samples_bam_index.out.och
        .map{ it ->
            it[0,1]
            }
        .set{ tumour_bams_index_flattened_ich }

    run_reheader_SAMtools_tumour(
        tumour_bams_flattened_ich,
        tumour_bams_index_flattened_ich
    )

    run_index_SAMtools_normal(
        run_reheader_SAMtools_normal.out.bam_reheadered,
        run_reheader_SAMtools_normal.out.associated_interval
        )
    
    run_index_SAMtools_tumour(
        run_reheader_SAMtools_tumour.out.bam_reheadered,
        run_reheader_SAMtools_tumour.out.associated_interval
        )

    // Merge the reheadered, indexed output based on the associated interval
    // Each element has the format [BAM, index, interval, id]
    normal_for_match = run_index_SAMtools_normal.out.indexed_out
        .map{ it ->
            [it[2].getFileName(), it[0], it[1], it[2], it[3]]
            }

    tumour_for_match = run_index_SAMtools_tumour.out.indexed_out
        .map{ it ->
            [it[2].getFileName(), [it[0], it[1], it[2], it[3]]]
            }
        .groupTuple()

    // Join and remove the join index
    normal_for_match
        .join(tumour_for_match, by: 0)
        .multiMap { it ->
            normal_bam_reheadered: it[1,2,3,4]
            tumour_bam_reheadered: it[5]
            }
        .set { matched_output_channel }

    emit:
    reheadered_normal_bam = matched_output_channel.normal_bam_reheadered
    reheadered_tumour_bam = matched_output_channel.tumour_bam_reheadered
    normal_bam_for_deletion = run_reheader_SAMtools_normal.out.bam_for_deletion
    tumour_bam_for_deletion = run_reheader_SAMtools_tumour.out.bam_for_deletion
    normal_bam_index_for_deletion = run_reheader_SAMtools_normal.out.bam_index_for_deletion
    tumour_bam_index_for_deletion = run_reheader_SAMtools_tumour.out.bam_index_for_deletion
}
