nextflow.enable.dsl=2

include { run_reheader_SAMtools as run_reheader_SAMtools_normal; run_reheader_SAMtools as run_reheader_SAMtools_tumour } from './bam-processing.nf'
include { run_BuildBamIndex_Picard as run_BuildBamIndex_Picard_normal; run_BuildBamIndex_Picard as run_BuildBamIndex_Picard_tumour } from './bam-processing.nf'

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
    tumour_bams.map{ it ->
        combined_it = []
        s = it.size
        while (!(s instanceof Integer)) {
            s = s.size
            }
        for(i = 0; i < s; i = i + 1) {
            combined_it = combined_it + [it[i][0] + 'my_reheader_separator' + it[i][1] + 'my_reheader_separator' + it[i][2]]
            }
        combined_it
        }
        .flatten()
        .map{ it ->
            it.split('my_reheader_separator')[0,1,2]
            }
        .set{ tumour_bams_flattened_ich }

    // tumour_bams_flattened_all
    //     .multiMap{it ->
    //         tumour_bams_flattened: it[0,1]
    //         tumour_bams_flattened_intervals: it[3]
    //         }
    //     .set{ tumour_bam_ich }

    tumour_bams_index.map{ it ->
        combined_it = []
        s = it.size
        while (!(s instanceof Integer)) {
            s = s.size
            }
        for(i = 0; i < s; i = i + 1) {
            combined_it = combined_it + [it[i][0] + 'my_reheader_separator' + it[i][1]]
            }
        combined_it
        }
        .flatten()
        .map{ it ->
            it.split('my_reheader_separator')[0,1]
            }
        .set{ tumour_bams_index_flattened_ich }

    tumour_bams_flattened_ich.map{it -> ['bam', it]}.view()
    tumour_bams_index_flattened_ich.map{it -> ['index', it]}.view()

    run_reheader_SAMtools_tumour(
        tumour_bams_flattened_ich,
        tumour_bams_index_flattened_ich
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
            [it[2].getFileName(), it[0], it[1], it[2]]
            }

    tumour_for_match = run_BuildBamIndex_Picard_tumour.out.indexed_out
        .map{ it ->
            [it[2].getFileName(), it[0], it[1], it[2]]
            }
        .groupTuple()

    normal_for_match
        .join(tumour_for_match, by: 0)
        .multiMap { it ->
            normal_bam_reheadered: it[1,2,3]
            tumour_bam_reheadered: it[4]
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
