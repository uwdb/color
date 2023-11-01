function build_experiments(experiment_params_list::Vector{ExperimentParams})
    for experiment_params in experiment_params_list
        build_times = [("Dataset", "Partitioner", "NumColors",  "BuildPhase", "BuildTime", "MemoryFootprint")]
        dataset = experiment_params.dataset
        summary_params = experiment_params.summary_params
        data = load_dataset(dataset)
        summary_name = params_to_summary_filename(experiment_params)
        summary_file_location = "Experiments/SerializedSummaries/" * summary_name
        println("Building Color Summary: ", summary_name)
        timing_vec = Float64[]
        results = @timed generate_color_summary(data, summary_params; verbose=1, timing_vec=timing_vec)
        summary_size = Base.summarysize(results.value)
        serialize(summary_file_location, results.value)
        push!(build_times, (string(dataset),
                             string(summary_params.partitioner),
                             string(summary_params.num_colors),
                             "FullTime",
                             string(results.time),
                             string(summary_size)))
        push!(build_times, (string(dataset),
                             string(summary_params.partitioner),
                             string(summary_params.num_colors),
                             "Coloring",
                             string(timing_vec[1]),
                             string(summary_size)))
        push!(build_times, (string(dataset),
                             string(summary_params.partitioner),
                             string(summary_params.num_colors),
                             "CycleCounting",
                             string(timing_vec[2]),
                             string(summary_size)))
        push!(build_times, (string(dataset),
                             string(summary_params.partitioner),
                             string(summary_params.num_colors),
                             "BloomFilter",
                             string(timing_vec[3]),
                             string(summary_size)))
        push!(build_times, (string(dataset),
                             string(summary_params.partitioner),
                             string(summary_params.num_colors),
                             "CardinalityCounting",
                             string(timing_vec[4]),
                             string(summary_size)))
        push!(build_times, (string(dataset),
                             string(summary_params.partitioner),
                             string(summary_params.num_colors),
                             "EdgeStats",
                             string(timing_vec[5]),
                             string(summary_size)))
        results_filename = params_to_results_filename(experiment_params)
        result_file_location = "Experiments/Results/Build_" * results_filename
        writedlm(result_file_location, build_times, ",")
    end
end
