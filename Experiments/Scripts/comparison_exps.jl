
using Profile
include("../Experiments.jl")

datasets = [human, aids, lubm80, yeast, dblp, youtube, eu2005, patents]
#datasets = [human, aids, yeast, dblp, youtube, eu2005, patents]
datasets = [human, youtube]

mix_scheme = [(QuasiStable, 32), (NeighborNodeLabels, 16), (NodeLabels, 16)]

experiment_params = Vector{ExperimentParams}()
for dataset in datasets
    push!(experiment_params, ExperimentParams(deg_stats_type=AvgDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=mix_scheme,
                                                description = "AvgMix64"))
#=
    push!(experiment_params, ExperimentParams(deg_stats_type=MinDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=mix_scheme,
                                                max_cycle_size = -1,
                                                description = "MinMix64"))

    push!(experiment_params, ExperimentParams(deg_stats_type=MaxDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=mix_scheme,
                                                max_cycle_size = -1,
                                                description = "MaxMix64"))

    push!(experiment_params, ExperimentParams(deg_stats_type=MaxDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(Hash, 64)],
                                                max_cycle_size = -1,
                                                inference_max_paths = 10^30,
                                                summary_max_paths=1000,
                                                use_partial_sums =false,
                                                description = "BSK++")) =#

    push!(experiment_params, ExperimentParams(deg_stats_type=AvgDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(QuasiStable, 1)],
                                                max_cycle_size = -1,
                                                description = "TradEst"))
end

build_experiments(experiment_params)

run_estimation_experiments(experiment_params; timeout=TIMEOUT_SEC)
comparison_methods =  ["alley", "wj", "impr", "jsub", "cs", "cset", "sumrdf"]
x_order = [string(data) for data in datasets]
legend_order = [params.description for params in experiment_params][1:Int(length(experiment_params)/ length(datasets))]
legend_order = vcat(legend_order, comparison_methods)

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
                                                y_label="Inference Latency 10^ (s)",
                                                filename="overall_runtime1")

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
                                                y_label="Relative Error 10^",
                                                filename="overall_error1")

graph_grouped_bar_plot(experiment_params;
                        grouping=description,
                        y_type=memory_footprint,
                        x_order = x_order,
                        legend_order = legend_order,
                        ylims=[0, 50],
                        y_ticks = [10, 20, 30, 40, 50],
                        legend_pos=:topright,
                        dimensions = (1000, 550),
                        y_label="Memory (MBs)",
                        filename="overall_memory1")

graph_grouped_bar_plot(experiment_params;
                        grouping=description,
                        y_type=build_time,
                        x_order = x_order,
                        legend_order = legend_order,
                        legend_pos=:topright,
                        ylims=[0, 3500],
                        y_ticks = [500, 1000, 1500, 2000, 2500, 3000],
                        dimensions = (1000, 550),
                        y_label="Build Time (s)",
                        filename="overall_build_time1")
