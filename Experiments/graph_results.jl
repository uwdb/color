using StatsPlots
using CSV, DataFrames
include("utils.jl")
# swag, think of it like a three column table
# x-value | y-value | group it belongs to
# then repeat for every single dot on the graph

@enum GROUP dataset technique cycle_size summary_paths inference_paths
#todo: query type?

@enum VALUE error runtime

#todo: figure out how to constrain query type
function graph_grouped_box_plot(experiment_params_list::Vector{ExperimentParams}; 
                                        x_type::GROUP=dataset, y_type::VALUE=error,
                                        grouping::GROUP=technique, 
                                        x_label=nothing, y_label=nothing, filename=nothing)
    # for now let's just use the dataset as the x-values and the cycle size as the groups
    x_values = []
    y_values = []
    groups = []
    for experiment_params in experiment_params_list
        # load the results
        results_filename = params_to_results_filename(experiment_params)
        results_path = "Experiments/Results/Estimation_" * results_filename
        # println("results path: ", results_path)
        results_df = CSV.read(results_path, DataFrame; normalizenames=true)

        # get the x_value and grouping (same for all results in this experiment param)
        current_x = get_value_from_param(experiment_params, x_type)
        current_group = get_value_from_param(experiment_params, grouping)
        # keep track of the data points
        for i in 1:nrow(results_df)
            current_y = 0
            if y_type == error
                current_y = results_df[i, :Estimate] / results_df[i, :TrueCard]
            else # y_type == runtime
                current_y = results_df[i, :EstimationTime]
            end
            # push the errors and their groupings into the correct vector
            push!(x_values, current_x)
            push!(y_values, current_y)
            push!(groups, current_group)
        end
    end
    println("starting graphs")
    plot = groupedboxplot(x_values, y_values, group = groups, left_margin = 10mm, bottom_margin = 10mm, yscale =:log10,  ylims=[10^-11, 10^11], yticks=[10^-5, 1, 10^5, 10^10])
    x_label !== nothing && xlabel!(plot, x_label)
    y_label !== nothing && ylabel!(plot, y_label)
    plotname = (isnothing(filename)) ? resultsfilename * ".png" : filename * ".png"
    savefig(plot, "Experiments/Results/Figures/" * plotname)
end

# default to grouping by dataset
function get_value_from_param(experiment_param::ExperimentParams, value_type::GROUP)
    if value_type == dataset
        return experiment_param.dataset
    elseif value_type == cycle_size
        return experiment_param.summary_params.max_cycle_size
    elseif value_type == summary_paths
        return experiment_param.summary_params.max_partial_paths
    elseif value_type == inference_paths
        return experiment_param.inference_max_paths
    else
        # default to grouping by technique
        return experiment_param.summary_params.partitioner
    end
end

# function graph_grouped_box_plots_build_samples(experiment_params_list::Vector{ExperimentParams}; graph_title=nothing)
#     # for now let's just use the dataset as the x-values and the cycle size as the groups
#     x_values = []
#     y_values = []
#     groups = []
#     for experiment_params in experiment_params_list
#         # get dataset and partitioner to group with
#         summary_filename = params_to_results_filename(experiment_params)
#         summary_path = "Experiments/SerializedSummaries/" * summary_filename
#         println("summary path: ", summary_path)
#         summary_df = CSV.read(summary_path, DataFrame; normalizenames=true)
#         dataset = experiment_params.dataset

#         println("dataset: ", dataset)
#         max_partial_paths = experiment_params.summary_params.max_partial_paths;
#         # partitioner = summary_df[1, :Partitioner] # this should be substituted by whatever we group by
#         results_filename = params_to_results_filename(experiment_params)
#         results_path = "Experiments/Results/Estimation_" * results_filename
#         println("results path: ", results_path)
#         results_df = CSV.read(results_path, DataFrame; normalizenames=true)
#         # get the errors
#         for i in 1:nrow(results_df)
#             error = results_df[i, :Estimate] / results_df[i, :TrueCard]
#             # push the errors and their groupings into the correct vector
#             push!(x_values, dataset)
#             push!(y_values, error)
#             push!(groups, max_partial_paths)
#         end
#     end
#     println("starting graphs")
#     plot = groupedboxplot(x_values, y_values, group = groups, yscale =:log10,  ylims=[10^-11, 10^11], yticks=[10^-5, 1, 10^5, 10^10])
#     plotname = (isnothing(graph_title)) ? resultsfilename * ".png" : graph_title * ".png"
#     savefig(plot, "Experiments/Results/Figures/" * plotname)
# end

# function graph_grouped_box_plots_inference_samples(experiment_params_list::Vector{ExperimentParams}; graph_title=nothing)
#     # for now let's just use the dataset as the x-values and the cycle size as the groups
#     x_values = []
#     y_values = []
#     groups = []
#     for experiment_params in experiment_params_list
#         # get dataset and partitioner to group with
#         summary_filename = params_to_results_filename(experiment_params)
#         summary_path = "Experiments/SerializedSummaries/" * summary_filename
#         println("summary path: ", summary_path)
#         summary_df = CSV.read(summary_path, DataFrame; normalizenames=true)
#         dataset = experiment_params.dataset

#         println("dataset: ", dataset)
#         max_partial_paths = experiment_params.inference_max_paths;
#         # partitioner = summary_df[1, :Partitioner] # this should be substituted by whatever we group by
#         results_filename = params_to_results_filename(experiment_params)
#         results_path = "Experiments/Results/Estimation_" * results_filename
#         println("results path: ", results_path)
#         results_df = CSV.read(results_path, DataFrame; normalizenames=true)
#         # get the errors
#         for i in 1:nrow(results_df)
#             error = results_df[i, :Estimate] / results_df[i, :TrueCard]
#             # push the errors and their groupings into the correct vector
#             push!(x_values, dataset)
#             push!(y_values, results_df[i, :EstimationTime])
#             push!(groups, max_partial_paths)
#         end
#     end
#     println("starting graphs")
#     plot = groupedboxplot(x_values, y_values, group = groups, yscale =:log10,  ylims=[10^-11, 10^11], yticks=[10^-5, 1, 10^5, 10^10])
#     plotname = (isnothing(graph_title)) ? resultsfilename * ".png" : graph_title * ".png"
#     savefig(plot, "Experiments/Results/Figures/" * plotname)
# end