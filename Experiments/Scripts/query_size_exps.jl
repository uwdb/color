
using Profile
include("../Experiments.jl")

#datasets = [human, aids, lubm80, yeast, hprd, dblp, youtube, eu2005, patents, wordnet]
datasets = [human, aids, lubm80, yeast, dblp, youtube, eu2005, patents]

experiment_params = Vector{ExperimentParams}()
for dataset in datasets
    push!(experiment_params, ExperimentParams(deg_stats_type=AvgDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(QuasiStable, 64)],
                                                description = "AvgQ64"))
    push!(experiment_params, ExperimentParams(deg_stats_type=AvgDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(QuasiStable, 32), (NeighborNodeLabels, 32),(QuasiStable, 32), (NeighborNodeLabels, 32)],
                                                description = "AvgQ64N64"))

    push!(experiment_params, ExperimentParams(deg_stats_type=MinDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(QuasiStable, 64)],
                                                description = "MinQ64"))

    push!(experiment_params, ExperimentParams(deg_stats_type=MaxDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(QuasiStable, 64)],
                                                description = "MaxQ64"))

    push!(experiment_params, ExperimentParams(deg_stats_type=MaxDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(Hash, 64)],
                                                max_cycle_size = -1,
                                                description = "BSK"))

    push!(experiment_params, ExperimentParams(deg_stats_type=AvgDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(QuasiStable, 1)],
                                                max_cycle_size = -1,
                                                description = "IndEst"))
end

#build_experiments(experiment_params)

#run_estimation_experiments(experiment_params)

graph_grouped_boxplot_with_comparison_methods(experiment_params;
                                                ylims=[10^-5, 10^4],
                                                y_ticks=[10^-5, 10^-2, 10^0, 10^2, 10^4],
                                                x_type = query_size,
                                                y_type = runtime,
                                                grouping=description,
                                                dimensions = (1450, 550),
                                                legend_pos=:topleft,
                                                y_label="Runtime 10^ (s)",
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
                                                y_label="Relative Error 10^",
                                                x_label = "Query Size",
                                                filename="query_size_error")
