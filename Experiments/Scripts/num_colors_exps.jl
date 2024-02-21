
using Profile
include("../Experiments.jl")

# The goal of this file is to evaluate the effect of using different numbers of colors during summary-building time, 
# across different datasets.
datasets = [yeast]
partitioning_schemes = [
    [(QuasiStable, 1)],
    [(QuasiStable, 4)],
    [(QuasiStable, 16)],
    [(QuasiStable, 32)],
    [(QuasiStable, 64)],
    [(QuasiStable, 128)],
    [(QuasiStable, 256)],
    [(QuasiStable, 512)],
]
experiment_params = Vector{ExperimentParams}()
for dataset in datasets
    for scheme in partitioning_schemes
            push!(experiment_params, ExperimentParams(dataset=dataset, partitioning_scheme=scheme))
    end
end

#build_experiments(experiment_params)

#run_estimation_experiments(experiment_params)

graph_grouped_box_plot(experiment_params;
                        ylims=[10^-20, 10^15],
                        y_ticks =[10^-20,10^-15, 10^-10, 10^-5, 1, 10^5, 10^10, 10^15],
                        grouping=technique,
                        y_label="Relative Error 10^ (s)",
                        dimensions = (1000, 600),
                        legend_pos = :topright,
                        filename="num_colors_error")

graph_grouped_box_plot(experiment_params;
                        ylims=[10^-4.5, 10^1.5],
                        y_ticks =[10^-4, 10^-3, 10^-2, 10^-1, 1, 10],
                        grouping=technique,
                        y_type=runtime,
                        y_label="Runtime 10^ (s)",
                        dimensions = (1000, 600),
                        legend_pos = :topleft,
                        filename="num_colors_runtime")

graph_grouped_bar_plot(experiment_params,
                        y_type=memory_footprint,
                        ylims=[0, 30],
                        y_ticks =[5 ,10, 15, 20, 25, 30],
                        y_label="Memory Footprint (MB)",
                        dimensions = (600, 400),
                        grouping=technique,
                        legend_pos = :topleft,
                        filename="num_colors_memory")
