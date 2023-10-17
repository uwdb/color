
using Profile
include("../Experiments.jl")

datasets = [aids, yeast, dblp]

experiment_params = Vector{ExperimentParams}()
build_params = Vector{ExperimentParams}()
for dataset in datasets
    push!(build_params, ExperimentParams(dataset=dataset))
    for only_shortest_path_cycle in [false, true]
        push!(experiment_params, ExperimentParams(dataset=dataset,
                                                  only_shortest_path_cycle=only_shortest_path_cycle))
    end
end

#build_experiments(build_params)

run_estimation_experiments(experiment_params)

graph_grouped_box_plot(experiment_params; grouping=cycle_stats, filename="cycle_stats_exps")
