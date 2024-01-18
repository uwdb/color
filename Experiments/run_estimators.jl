function run_estimation_experiments(experiment_params_list::Vector{ExperimentParams}; timeout::Float64=Inf)
    for experiment_params in experiment_params_list
        dataset = experiment_params.dataset
        replications = experiment_params.n_replications
        all_queries = load_querysets([dataset]; require_true_cardinality = true)
        summary_file_location = "Experiments/SerializedSummaries/" * params_to_summary_filename(experiment_params)
        !isfile(summary_file_location) && error("The summary has not been built yet! \n Attempted File Location: $(summary_file_location)")
        summary::ColorSummary = deserialize(summary_file_location)
        experiment_results = SharedArray{Tuple{Float64, Float64, Float64, String255, String255, Float64, Bool, Int64}}(length(all_queries[dataset]))
        @sync @distributed for i in shuffle(collect(eachindex(experiment_results)))
            query::QueryGraph = all_queries[dataset][i].query
            query_path = all_queries[dataset][i].query_path
            # println(query_path)
            exact_size = all_queries[dataset][i].exact_size
            estimate_results = [(@timed get_cardinality_bounds(query,
                                    summary;
                                    max_partial_paths = experiment_params.inference_max_paths,
                                    use_partial_sums=experiment_params.use_partial_sums,
                                    usingStoredStats=true,
                                    sampling_strategy=experiment_params.sampling_strategy,
                                    only_shortest_path_cycle= experiment_params.only_shortest_path_cycle,
                                    timeout=timeout)) for _ in 1:replications]
            estimate_time = median([x.time for x in  estimate_results]) # Convert back to seconds from nano seconds
            estimate_failure = isnan(estimate_results[1].value) || isinf(estimate_results[1].value) || estimate_results[1].value == 0
            estimate = max(1, estimate_results[1].value)
            if isinf(estimate)
                estimate = 10^35
            end
            if isnan(estimate)
                estimate = 1.0
            end
            query_type = all_queries[dataset][i].query_type
            path_width = get_min_width_node_order(query.graph, return_width=true)
            experiment_results[i] = (estimate, exact_size, estimate_time, query_type, query_path, nv(query.graph), estimate_failure, path_width)
        end
        final_results = [(x[1], x[2], x[3], String(x[4]), String(x[5]), x[6], x[7], x[8]) for x in experiment_results]
        final_results = [("Estimate", "TrueCard", "EstimationTime", "QueryType", "QueryPath", "QuerySize", "Failure", "PathWidth"); final_results]
        results_file_location = "Experiments/Results/Estimation_"  * params_to_results_filename(experiment_params)
        writedlm(results_file_location, final_results, ",")
    end
end
