function run_estimation_experiments(experiment_params_list::Vector{ExperimentParams})
    for experiment_params in experiment_params_list
        dataset = experiment_params.dataset
        all_queries = load_querysets([dataset]; require_true_cardinality = true)
        summary_file_location = "Experiments/SerializedSummaries/" * params_to_summary_filename(experiment_params)
        !isfile(summary_file_location) && error("The summary has not been built yet! \n Attempted File Location: $(summary_file_location)")
        summary::ColorSummary = deserialize(summary_file_location)
        experiment_results = []
        push!(experiment_results, ("UpperBound", "Estimate", "LowerBound", "TrueCard", "EstimationTime", "QueryType", "QueryPath"))
        for i in 1:length(all_queries[dataset])
            query::QueryGraph = all_queries[dataset][i].query
            query_path = all_queries[dataset][i].query_path
            exact_size = all_queries[dataset][i].exact_size
            estimate_times = [(@timed get_cardinality_bounds(query, summary;
                                    max_partial_paths = experiment_params.inference_max_paths,
                                    use_partial_sums=experiment_params.use_partial_sums, usingStoredStats=true,
                                    sampling_strategy=experiment_params.sampling_strategy,
                                    only_shortest_path_cycle= experiment_params.only_shortest_path_cycle)).time for _ in 1:3]
            estimate_time = median(estimate_times) # Convert back to seconds from nano seconds

            bounds = get_cardinality_bounds(query, summary;
                                max_partial_paths = experiment_params.inference_max_paths,
                                use_partial_sums=experiment_params.use_partial_sums, usingStoredStats=true,
                                sampling_strategy=experiment_params.sampling_strategy,
                                only_shortest_path_cycle= experiment_params.only_shortest_path_cycle)
            upper_bound = bounds[3]
            estimate = max(1, bounds[2])
            lower_bound = bounds[1]
            query_type = all_queries[dataset][i].query_type
            push!(experiment_results, (upper_bound, estimate, lower_bound, exact_size, estimate_time, query_type, query_path))
        end
        results_file_location = "Experiments/Results/Estimation_"  * params_to_results_filename(experiment_params)
        writedlm(results_file_location, experiment_results, ",")
    end
end
