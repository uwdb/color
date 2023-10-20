
using Profile
include("../Experiments.jl")

#datasets = [human, aids, yeast, hprd, dblp, youtube, wordnet]
datasets = [dblp]

experiment_params = Vector{ExperimentParams}()
build_params = Vector{ExperimentParams}()
for dataset in datasets
    push!(build_params, ExperimentParams(dataset=dataset))
    for use_corr in [false, true]
        push!(experiment_params, ExperimentParams(dataset=dataset, use_corr=use_corr))
    end
end

#build_experiments(build_params)

run_estimation_experiments(experiment_params)

graph_grouped_box_plot(experiment_params; grouping=correlation_strategy, filename="corr_comparison_2")
