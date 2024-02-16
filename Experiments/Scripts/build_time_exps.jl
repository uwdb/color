
using Profile
include("../Experiments.jl")

# The goal of this file is to create graphs comparing the build time for different datasets.

datasets = [aids, yeast, hprd, dblp, youtube, wordnet]
#datasets = [wordnet]
experiment_params = Vector{ExperimentParams}()
build_params = Vector{ExperimentParams}()
for dataset in datasets
    push!(build_params, ExperimentParams(dataset=dataset))
end
#build_experiments(build_params)

graph_grouped_bar_plot(build_params; grouping=build_phase,
                                          y_type=build_time,
                                          ylims=[0, 360],
                                          filename="build_time")
