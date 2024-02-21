
using Profile
include("../Experiments.jl")

# The goal of this file is to evaluate the memory footprint for using the same technique on different datasets.

#datasets = [aids, yeast, hprd, dblp, youtube, wordnet]
datasets = [wordnet]
num_colors = [4, 8, 16, 32, 64, 128, 256]
experiment_params = Vector{ExperimentParams}()
build_params = Vector{ExperimentParams}()
for dataset in datasets
    for n in num_colors
        push!(build_params, ExperimentParams(dataset=dataset, partitioning_scheme=[(QuasiStable, n)]))
    end
end
build_experiments(build_params)

graph_grouped_bar_plot(build_params; grouping=number_of_colors,
                                          y_type=memory_footprint,
                                          ylims=[0, 16],
                                          filename="memory_size_vs_colors_fp32_int16_2")
