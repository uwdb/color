include("../Experiments.jl")

#datasets = [human, aids]
datasets = [human, aids, lubm80, yeast, dblp, youtube, eu2005, patents]
queries = load_querysets(datasets)
num_queries = Dict(string(dataset)=>length(queries[dataset]) for dataset in datasets)

methods, comparison_results = comparison_dataset()

failure_counts = Dict()
failure_probabilities = Dict()
for method in methods
    failure_counts[method] = counter(String)
    failure_probabilities[method] = Dict()
    for dataset in datasets
        string_dataset = string(dataset)
        for query in queries[dataset]
            qid = get_query_id(string_dataset, query.query_path)
            comp_key = (string_dataset, method, qid)
            if !haskey(comparison_results, comp_key)
                inc!(failure_counts[method], string_dataset)
            elseif comparison_results[comp_key].Estimate == 0
                inc!(failure_counts[method], string_dataset)
            elseif comparison_results[comp_key].Estimate == Inf
                inc!(failure_counts[method], string_dataset)
            elseif comparison_results[comp_key].Estimate == NaN
                inc!(failure_counts[method], string_dataset)
            end
        end
        failure_probabilities[method][dataset] = failure_counts[method][string_dataset] / num_queries[string_dataset]
    end
end

failure_counts["BSK"] = counter(String)
failure_counts["BSK++"] = counter(String)
failure_counts["AvgQ64"] = counter(String)
failure_probabilities["BSK"] = Dict()
failure_probabilities["BSK++"] = Dict()
failure_probabilities["AvgQ64"] = Dict()
for dataset in datasets
    string_dataset = string(dataset)
    bsk_params = ExperimentParams(deg_stats_type=MaxDegStats,
                                    dataset=dataset,
                                    partitioning_scheme=[(Hash, 64)],
                                    max_cycle_size = -1,
                                    inference_max_paths = 10^30,
                                    use_partial_sums = false,
                                    description = "BSK",
                                    n_replications = 1)
    run_estimation_experiments([bsk_params]; timeout=TIMEOUT_SEC)
    bsk_filename = params_to_results_filename(bsk_params)
    bsk_path = "Experiments/Results/Estimation_" * bsk_filename
    bsk_df = CSV.read(bsk_path, DataFrame; normalizenames=true)
    for i in 1:nrow(bsk_df)
        if bsk_df[i, :Failure]
            inc!(failure_counts["BSK"], string_dataset)
        end
    end
    failure_probabilities["BSK"][string_dataset] = failure_counts["BSK"][string_dataset] / num_queries[string_dataset]


    bsk_agg_params = ExperimentParams(deg_stats_type=MaxDegStats,
                                    dataset=dataset,
                                    partitioning_scheme=[(Hash, 64)],
                                    max_cycle_size = -1,
                                    inference_max_paths = 10^30,
                                    use_partial_sums = true,
                                    description = "BSK++",
                                    n_replications=1)
    run_estimation_experiments([bsk_agg_params]; timeout=TIMEOUT_SEC)
    bsk_agg_filename = params_to_results_filename(bsk_agg_params)
    bsk_agg_path = "Experiments/Results/Estimation_" * bsk_agg_filename
    bsk_agg_df = CSV.read(bsk_agg_path, DataFrame; normalizenames=true)
    for i in 1:nrow(bsk_agg_df)
        if bsk_agg_df[i, :Failure]
            inc!(failure_counts["BSK++"], string_dataset)
        end
    end
    failure_probabilities["BSK++"][string_dataset] = failure_counts["BSK++"][string_dataset] / num_queries[string_dataset]



    avg_params = ExperimentParams(dataset=dataset, n_replications=1)
    run_estimation_experiments([avg_params]; timeout=TIMEOUT_SEC)
    avg_filename = params_to_results_filename(avg_params)
    avg_path = "Experiments/Results/Estimation_" * avg_filename
    avg_df = CSV.read(avg_path, DataFrame; normalizenames=true)
    for i in 1:nrow(avg_df)
        if avg_df[i, :Failure]
            inc!(failure_counts["AvgQ64"], string_dataset)
        end
    end
    failure_probabilities["AvgQ64"][string_dataset] = failure_counts["AvgQ64"][string_dataset] / num_queries[string_dataset]
end
