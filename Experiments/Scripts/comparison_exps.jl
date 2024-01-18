
using Profile
include("../Experiments.jl")

datasets = [human, aids, lubm80, yeast, dblp, youtube, eu2005, patents]
#datasets = [human, aids, yeast, dblp, youtube, eu2005, patents]
#datasets = [human, aids]
# datasets = [yeast]

mix_scheme = [(Degree, 8), (QuasiStable, 8), (NeighborNodeLabels, 8), (NodeLabels, 8)]

experiment_params = Vector{ExperimentParams}()
for dataset in datasets
    push!(experiment_params, ExperimentParams(deg_stats_type=AvgDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=mix_scheme,
                                                description = "AvgMix32"))


    push!(experiment_params, ExperimentParams(deg_stats_type=MinDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=mix_scheme,
                                                max_cycle_size = -1,
                                                description = "MinMix32"))

    push!(experiment_params, ExperimentParams(deg_stats_type=MaxDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=mix_scheme,
                                                max_cycle_size = -1,
                                                description = "MaxMix32"))

    push!(experiment_params, ExperimentParams(deg_stats_type=MaxDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(Hash, 64)],
                                                max_cycle_size = -1,
                                                inference_max_paths = 10^30,
                                                summary_max_paths=1000,
                                                use_partial_sums = true,
                                                description = "BSK++"))
    push!(experiment_params, ExperimentParams(deg_stats_type=AvgDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(QuasiStable, 1)],
                                                max_cycle_size = -1,
                                                description = "TradEst"))
end

println("Building...")

build_experiments(experiment_params)

println("Estimating...")

run_estimation_experiments(experiment_params; timeout=TIMEOUT_SEC)

comparison_methods =  ["alley", "wj", "impr", "jsub", "cs", "cset", "sumrdf"]
x_order = [string(data) for data in datasets]
legend_order = [params.description for params in experiment_params][1:Int(length(experiment_params)/ length(datasets))]
legend_order = vcat(legend_order, comparison_methods)

println("Graphing figures 2 and 3...")

graph_grouped_boxplot_with_comparison_methods(experiment_params;
                                                ylims=[10^-5, 10^4],
                                                y_ticks=[10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3, 10^4],
                                                y_type = runtime,
                                                x_type = dataset,
                                                x_order = x_order,
                                                legend_order = legend_order,
                                                grouping=description,
                                                dimensions = (1550, 650),
                                                legend_pos=:topleft,
                                                y_label="Inference Latency log\$_{10}\$ (s)",
                                                filename="overall_runtime")

graph_grouped_boxplot_with_comparison_methods(experiment_params;
                                                ylims=[10^-21, 10^21],
                                                y_ticks=[10^-20, 10^-15, 10^-10, 10^-5, 10^-2, 10^0, 10^2, 10^5, 10^10, 10^15, 10^20],
                                                y_type = estimate_error,
                                                x_type = dataset,
                                                x_order = x_order,
                                                legend_order = legend_order,
                                                grouping=description,
                                                dimensions = (1550, 650),
                                                legend_pos=:bottomleft,
                                                y_label="Relative Error log\$_{10}\$",
                                                filename="overall_error")

println("Graphing figure 4")

graph_grouped_boxplot_with_comparison_methods(experiment_params;
                                                ylims=[10^-5, 10^4],
                                                y_ticks=[10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3, 10^4],
                                                x_type = query_size,
                                                y_type = runtime,
                                                grouping=description,
                                                dimensions = (1450, 550),
                                                legend_pos=:topright,
                                                legend_order = legend_order,
                                                y_label="Inference Latency log\$_{10}\$ (s)",
                                                x_label = "Query Size",
                                                filename="query_size_runtime")


graph_grouped_boxplot_with_comparison_methods(experiment_params;
                                                ylims=[10^-21, 10^21],
                                                x_type = query_size,
                                                y_ticks=[10^-20, 10^-15, 10^-10, 10^-5, 10^-2, 10^0, 10^2, 10^5, 10^10, 10^15, 10^20],
                                                y_type = estimate_error,
                                                grouping=description,
                                                dimensions = (1450, 550),
                                                legend_pos=:bottomleft,
                                                legend_order = legend_order,
                                                y_label="Relative Error log\$_{10}\$",
                                                x_label = "Query Size",
                                                filename="query_size_error")
