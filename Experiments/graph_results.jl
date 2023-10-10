using StatsPlots
using CSV, DataFrames
include("utils.jl")

@enum GROUP dataset technique cycle_size summary_paths inference_paths query_type
#todo: query type

@enum VALUE estimate_error runtime

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
        
        # keep track of the data points
        for i in 1:nrow(results_df)
            current_x = x_type == query_type ? results_df[i, :QueryType] : get_value_from_param(experiment_params, x_type)
            current_group = grouping == query_type ? results_df[i, :QueryType] : get_value_from_param(experiment_params, grouping)
            current_y = 0
            if y_type == estimate_error
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
    gbplot = groupedboxplot(x_values, y_values, group = groups, left_margin = 10mm, bottom_margin = 10mm, yscale =:log10,  ylims=[10^-11, 10^11], yticks=[10^-5, 1, 10^5, 10^10])
    x_label !== nothing && xlabel!(gbplot, x_label)
    y_label !== nothing && ylabel!(gbplot, y_label)
    plotname = (isnothing(filename)) ? resultsfilename * ".png" : filename * ".png"
    savefig(gbplot, "Experiments/Results/Figures/" * plotname)
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
