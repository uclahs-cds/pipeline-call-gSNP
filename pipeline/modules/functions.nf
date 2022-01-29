// def flatten_samples(groovyx.gpars.dataflow.DataflowBroadcast ich) {
//     ich.map{ it ->
//         outer_tuple = []
//         s = it.size
//         while (!(s instanceof Integer)) {
//             s = s.size
//             }
//         for(i = 0; i < s; i = i + 1) {
//             outer_tuple = outer_tuple + ["$i": it[i]]
//             }
//         outer_tuple
//         }
//         .flatten()
//         .map{ it ->
//             it.values().flatten()
//             }
//     }

workflow flatten_samples {
    take:
    ich

    main:
    ich.map{ it ->
        outer_tuple = []
        for(elem in it) {
            outer_tuple = outer_tuple + ["key": elem]
            }
        outer_tuple
        }
        .flatten()
        .map{ it ->
            it.values().flatten()
            }
        .set{ och }

    emit:
    och = och    
}

// Set workflow to process "element_index" since the process is identical between BAMs and BAIs
workflow gather_vcf_input {
    take:
    normal_ich
    tumour_ich
    element_index

    main:
    ich.map{ it ->
        outer_tuple = []
        for(elem in it) {
            outer_tuple = outer_tuple + ["key": elem]
            }
        outer_tuple
        }
        .flatten()
        .map{ it ->
            it.values().flatten()
            }
        .set{ och }

    emit:
    och = och    
}