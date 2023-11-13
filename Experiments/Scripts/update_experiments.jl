using Plots.PlotMeasures
using Graphs
include("../Experiments.jl")

datasets::Vector{DATASET} = [wordnet]
# datasets::Vector{DATASET} = [aids, human, yeast, wordnet, youtube, dblp, patents]
# datasets::Vector{DATASET} = [aids, human, lubm80, yeast, hprd, wordnet, dblp, youtube, eu2005, patents]
max_cycles = 6
proportions_not_updated = [0, 0.2, 0.4, 0.6, 0.8, 1]

experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=current_dataset, partitioner=QuasiStable, max_cycle_size=current_cycle, proportion_not_updated=current_proportion) 
                                                    for current_dataset in datasets for current_cycle in 2:max_cycles for current_proportion in proportions_not_updated]
println("started building")
build_experiments(experiment_params_list)
println("started estimating")
run_estimation_experiments(experiment_params_list)
println("started graphing")
# compare how overall accuracy is affected by summary updates
# graph_grouped_box_plot(experiment_params_list, x_type=dataset, y_type=estimate_error, grouping=proportion_not_updated, filename="overall-accuracy-and-updates")
# compare how cycle stat accuracies are affected by summary updates
graph_grouped_box_plot(experiment_params_list, x_type=proportion_not_updated, y_type=estimate_error, grouping=cycle_size, filename="cycle-stats-and-updates")