include("../Experiments.jl")

#datasets = [human, aids]
datasets = [human, aids, lubm80, yeast, dblp, youtube, eu2005, patents]
#datasets = [human, aids, yeast, dblp, youtube, eu2005, patents]
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
        failure_probabilities[method][string_dataset] = failure_counts[method][string_dataset] / num_queries[string_dataset]
    end
end

failure_counts["BSK"] = counter(String)
failure_counts["BSK++"] = counter(String)
failure_counts["AvgMix64"] = counter(String)
failure_probabilities["BSK"] = Dict()
failure_probabilities["BSK++"] = Dict()
failure_probabilities["AvgMix64"] = Dict()
for dataset in datasets
    string_dataset = string(dataset)
    bsk_params = ExperimentParams(deg_stats_type=MaxDegStats,
                                    dataset=dataset,
                                    partitioning_scheme=[(Hash, 64)],
                                    max_cycle_size = -1,
                                    inference_max_paths = 10^30,
                                    summary_max_paths=1000,
                                    use_partial_sums = false,
                                    description = "BSK",
                                    n_replications = 1)
#    run_estimation_experiments([bsk_params]; timeout=TIMEOUT_SEC)
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
                                    summary_max_paths=1000,
                                    use_partial_sums = true,
                                    description = "BSK++",
                                    n_replications=1)
#    run_estimation_experiments([bsk_agg_params]; timeout=TIMEOUT_SEC)
    bsk_agg_filename = params_to_results_filename(bsk_agg_params)
    bsk_agg_path = "Experiments/Results/Estimation_" * bsk_agg_filename
    bsk_agg_df = CSV.read(bsk_agg_path, DataFrame; normalizenames=true)
    for i in 1:nrow(bsk_agg_df)
        if bsk_agg_df[i, :Failure]
            inc!(failure_counts["BSK++"], string_dataset)
        end
    end
    failure_probabilities["BSK++"][string_dataset] = failure_counts["BSK++"][string_dataset] / num_queries[string_dataset]


    mix_scheme = [(QuasiStable, 32), (NeighborNodeLabels, 16), (NodeLabels, 16)]
    avg_params = ExperimentParams(dataset=dataset,
                                    n_replications=2,
                                    partitioning_scheme=mix_scheme)
#    build_experiments([avg_params])
#    run_estimation_experiments([avg_params]; timeout=TIMEOUT_SEC)
    avg_filename = params_to_results_filename(avg_params)
    avg_path = "Experiments/Results/Estimation_" * avg_filename
    avg_df = CSV.read(avg_path, DataFrame; normalizenames=true)
    for i in 1:nrow(avg_df)
        if avg_df[i, :Failure]
            inc!(failure_counts["AvgMix64"], string_dataset)
        end
    end
    failure_probabilities["AvgMix64"][string_dataset] = failure_counts["AvgMix64"][string_dataset] / num_queries[string_dataset]
end

estimators = ["cs", "wj", "jsub", "impr", "cset", "alley", "BSK", "BSK++", "sumrdf", "AvgMix64"]

global latex_table = """
\\begin{table*}[]
\\begin{tabular}{|l|l|l|l|l|l|l|l|l|l|l|}
\\hline
\\textbf{Dataset\\textbackslash{}Method} """
for estimator in estimators
    global latex_table *= """& \\textbf{""" * string(estimator) * """} """
end
global latex_table *= """\\\\
 \\hline"""
for dataset in datasets
    global latex_table *= """\\textbf{""" * string(dataset) * """} """
    for estimator in estimators
        global latex_table *= " & " * @sprintf("%.2f", failure_probabilities[estimator][string(dataset)])
    end
    global latex_table *= """\\\\ \\hline """
end
global latex_table *= """
\\end{tabular}
\\caption{Estimator Failure Rates}
\\label{tbl:estimator-failure}
\\end{table*}
"""
println(latex_table)
