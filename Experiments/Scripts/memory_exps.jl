
using Profile
include("../Experiments.jl")

#datasets = [aids, yeast, hprd, dblp, youtube, wordnet]
datasets = [aids, yeast, hprd, dblp, wordnet]
num_colors = [4, 8, 16, 32, 64, 128]
experiment_params = Vector{ExperimentParams}()
build_params = Vector{ExperimentParams}()
for dataset in datasets
    for n in num_colors
        push!(build_params, ExperimentParams(dataset=dataset, num_colors=n))
    end
end
build_experiments(build_params)

graph_grouped_bar_plot(build_params; grouping=number_of_colors,
                                          y_type=memory_footprint,
                                          y_lims=[0, 16],
                                          filename="memory_size_vs_colors_fp32_int16")
