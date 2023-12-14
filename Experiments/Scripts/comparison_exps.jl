
using Profile
include("../Experiments.jl")

#datasets = [human, aids, lubm80, yeast, hprd, dblp, youtube, eu2005, patents, wordnet]
datasets = [human, aids, lubm80, yeast, dblp, youtube, eu2005, patents]

experiment_params = Vector{ExperimentParams}()
for dataset in datasets
    #= push!(experiment_params, ExperimentParams(deg_stats_type=CorrDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(QuasiStable, 64)],
                                                description = "CorrQ64"))
    push!(experiment_params, ExperimentParams(deg_stats_type=CorrDegStats,
                                                dataset=dataset,
                                                partitioning_scheme=[(QuasiStable, 32), (NeighborNodeLabels, 32),(QuasiStable, 32), (NeighborNodeLabels, 32)],
                                                description = "CorrQ64N64")) =#
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
                                                description = "MaxH64"))

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
                                                y_type = runtime,
                                                grouping=description,
                                                dimensions = (1450, 550),
                                                legend_pos=:top,
                                                y_label="Runtime (s)",
                                                filename="runtime")

graph_grouped_boxplot_with_comparison_methods(experiment_params;
                                                ylims=[10^-21, 10^21],
                                                y_ticks=[10^-20, 10^-15, 10^-5, 10^-2, 10^0, 10^2, 10^5, 10^10, 10^15, 10^20],
                                                y_type = estimate_error,
                                                grouping=description,
                                                dimensions = (1450, 550),
                                                legend_pos=:bottomleft,
                                                y_label="Relative Error",
                                                filename="error")

graph_grouped_bar_plot(experiment_params;
                        grouping=description,
                        y_type=memory_footprint,
                        y_lims=[0, 30],
                        y_ticks = [5, 10, 15, 20, 25, 30],
                        dimensions = (775, 550),
                        y_label="Memory (MBs)",
                        filename="memory")

graph_grouped_bar_plot(experiment_params;
                        grouping=description,
                        y_type=build_time,
                        y_lims=[0, 720],
                        y_ticks = [0, 100, 200, 300, 400, 500, 600, 700],
                        dimensions = (775, 550),
                        y_label="Build Time (s)",
                        filename="build_time")
