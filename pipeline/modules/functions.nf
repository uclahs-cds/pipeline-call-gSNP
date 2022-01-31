

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
