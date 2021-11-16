nextflow.enable.dsl=2

include { run_GatherBamFiles_GATK } from './bam-processing.nf'

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
    
    // Merge the channels for sorting by interval
    bam_interval_merged_sorted = bams_with_index
        .join(intervals_with_index, by: 0)
        .map{ it ->
            [it[1], it[2]]
            }
        .toSortedList( { a, b -> a[1] <=> b[1] } )

    
    run_GatherBamFiles_GATK(
        bam_interval_merged_sorted,
        bam_type,
        identifiers
        )

    emit:
    merged_bam = run_GatherBamFiles_GATK.out.merged_bam
    merged_bam_index = run_GatherBamFiles_GATK.out.merged_bam_index
}