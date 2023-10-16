
function run_estimation_experiments(experiment_params_list::Vector{ExperimentParams})
    for experiment_params in experiment_params_list
        dataset = experiment_params.dataset
        all_queries = load_querysets([dataset]; require_true_cardinality = true)
        summary_file_location = "Experiments/SerializedSummaries/" * params_to_summary_filename(experiment_params)
        !isfile(summary_file_location) && error("The summary has not been built yet! \n Attempted File Location: $(summary_file_location)")
        summary::ColorSummary = deserialize(summary_file_location)
        experiment_results = []
        push!(experiment_results, ("UpperBound", "Estimate", "LowerBound", "TrueCard", "EstimationTime", "QueryType"))
        for i in 1:length(all_queries[dataset])
            query = all_queries[dataset][i].query
            query_path = all_queries[dataset][i].query_path
            exact_size = all_queries[dataset][i].exact_size
            results = @timed get_cardinality_bounds(query, summary;
                                max_partial_paths = experiment_params.inference_max_paths,
                                use_partial_sums=experiment_params.use_partial_sums, usingStoredStats=true,
                                sampling_strategy=experiment_params.sampling_strategy)
            upper_bound = results.value[3]
            estimate = max(1, results.value[2])
            lower_bound = results.value[1]
            estimate_time = results.time
            query_type = all_queries[dataset][i].query_type
            push!(experiment_results, (upper_bound, estimate, lower_bound, exact_size, estimate_time, query_type))
        end
        results_file_location = "Experiments/Results/Estimation_"  * params_to_results_filename(experiment_params)
        writedlm(results_file_location, experiment_results, ",")
    end
end

# run_custom_queries(experiment_params_list::Vector{ExperimentParams}, custom_queries::Vector{QueryGraph})
#     for experiment_params in experiment_params_list
#         dataset = experiment_params.dataset
#         all_queries = custom_queries
#         summary_file_location = "Experiments/SerializedSummaries/" * params_to_summary_filename(experiment_params)
#         !isfile(summary_file_location) && error("The summary has not been built yet! \n Attempted File Location: $(summary_file_location)")
#         summary::ColorSummary = deserialize(summary_file_location)
#         experiment_results = []
#         push!(experiment_results, ("UpperBound", "Estimate", "LowerBound", "TrueCard", "EstimationTime", "QueryType"))
#         for i in 1:length(all_queries[dataset])
#             query = all_queries[dataset][i].query
#             exact_size = all_queries[dataset][i].exact_size
#             results = @timed get_cardinality_bounds(query, summary;
#                                 max_partial_paths = experiment_params.inference_max_paths,
#                                 use_partial_sums=experiment_params.use_partial_sums, usingStoredStats=true)
#             upper_bound = results.value[3]
#             estimate = max(1, results.value[2])
#             lower_bound = results.value[1]
#             estimate_time = results.time
#             query_type = all_queries[dataset][i].query_type
#             push!(experiment_results, (upper_bound, estimate, lower_bound, exact_size, estimate_time, query_type))
#         end
#         results_file_location = "Experiments/Results/Estimation_"  * params_to_results_filename(experiment_params)
#         writedlm(results_file_location, experiment_results, ",")
#     end
