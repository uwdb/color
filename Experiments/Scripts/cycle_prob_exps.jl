
using Profile
include("../Experiments.jl")

datasets = [aids, yeast, hprd, dblp, youtube, wordnet]

experiment_params = Vector{ExperimentParams}()
build_params = Vector{ExperimentParams}()
for dataset in datasets
    push!(build_params, ExperimentParams(dataset=dataset,
                        num_colors=16,
                        label_refining_rounds=2))
    for only_shortest_path_cycle in [false, true]
        push!(experiment_params, ExperimentParams(dataset=dataset,
                                                    num_colors=16,
                                                    label_refining_rounds=2,
                                                    only_shortest_path_cycle=only_shortest_path_cycle))
    end
end

build_experiments(build_params)

run_estimation_experiments(experiment_params)

graph_grouped_box_plot(experiment_params; grouping=cycle_stats, filename="cycle_stats_exps_w_refining")
