nextflow.enable.dsl=2

include { run_GatherBamFiles_GATK; run_BuildBamIndex_Picard; run_view_SAMtools } from './bam-processing.nf' addParams(
    is_output_bam: true
    )
include { remove_intermediate_files } from './intermediate-cleanup.nf' addParams(
    save_intermediate_files: false
    )

workflow gather_bams {
    take:
    bams
    intervals
    bam_type
    identifiers

    main:
    // Add index to BAM and intervals channel for merging
    counter_bams = 0
    bams_with_index = bams
        .map{ it ->
            [counter_bams = counter_bams + 1, it]
            }

    counter_intervals = 0
    intervals_with_index = intervals
        .map{ it ->
            [counter_intervals = counter_intervals + 1, it.getFileName()]
            }
    
    // Merge and sort channels by interval and extract only the BAMs
    sorted_bams = bams_with_index
        .join(intervals_with_index, by: 0)
        .map{ it ->
            [it[1], it[2]]
            }
        .toSortedList( { a, b -> a[1] <=> b[1] } )
        .flatten()
        .filter{ it.toString().endsWith('.bam') }
    
    run_GatherBamFiles_GATK(
        sorted_bams.collect(),
        bam_type,
        identifiers
        )
    
    run_view_SAMtools(
        run_GatherBamFiles_GATK.out.merged_bam
        )
    
    // Delete the unconverted, gathered BAM after conversion
    remove_intermediate_files(
        run_GatherBamFiles_GATK.out.merged_bam,
        run_view_SAMtools.out.bam
        )

    // Build the merged BAM index. The intervals parameter is just a placeholder
    run_BuildBamIndex_Picard(
        run_view_SAMtools.out.bam,
        intervals.collect().map{ it -> it[0] }
        )

    emit:
    merged_bam = run_view_SAMtools.out.bam
    merged_bam_index = run_BuildBamIndex_Picard.out.indexed_out
        .map{ it ->
            it[1]
            }
}