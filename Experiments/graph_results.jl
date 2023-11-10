@enum GROUP dataset technique cycle_size summary_paths inference_paths query_type sampling_type cycle_stats number_of_colors build_phase
#todo: query type

@enum VALUE estimate_error runtime build_time memory_footprint

function graph_grouped_box_plot(experiment_params_list::Vector{ExperimentParams};
                                        x_type::GROUP=dataset, y_type::VALUE=estimate_error,
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
    results_filename = params_to_results_filename(experiment_params_list[1])
    println("starting graphs")

    # This seems to be necessary for using Plots.jl outside of the ipynb framework.
    # See this: https://discourse.julialang.org/t/deactivate-plot-display-to-avoid-need-for-x-server/19359/15
    ENV["GKSwstype"]="100"
    gbplot = groupedboxplot(x_values, y_values, group = groups, yscale =:log10,
                            ylims=[10^-13, 10^11], yticks=[10^-10, 10^-5, 10^-2, 1, 10^2, 10^5, 10^10],
                            legend = :outertopleft, size = (1000, 600))
    x_label !== nothing && xlabel!(gbplot, x_label)
    y_label !== nothing && ylabel!(gbplot, y_label)
    plotname = (isnothing(filename)) ? results_filename * ".png" : filename * ".png"
    savefig(gbplot, "Experiments/Results/Figures/" * plotname)
end

function comparison_dataset()
    comparison_results = DataFrame(Dataset("Experiments/comparison_results.parquet"))
    comparison_results[!, "QueryType"] .= ""
    for i in 1:nrow(comparison_results)
        dataset = comparison_results[i, :Dataset]
        query_path = comparison_results[i, :Query]
        if dataset == "lubm80"
            comparison_results[i, :QueryType] = match(r".*/lubm80_(.*).txt", query_path).captures[1]
        elseif dataset in ["aids", "human", "yago"]
            comparison_results[i, :QueryType] = match(r"(.*)_.*/.*", query_path).captures[1]
        else
            comparison_results[i, :QueryType] = match(r".*/queryset/.*/query_(.*)_.*", query_path).captures[1]
        end
    end
    return comparison_results
end

function get_query_id(dataset, query_path)
    return if dataset == "lubm80"
        println(query_path)
        match(r".*/queryset/(.*/.*)", query_path).captures[1]
    elseif dataset in ["aids", "human", "yago"]
        match(r".*/queryset/.*/(.*/.*)", query_path).captures[1]
    else
        match(r".*/queryset/.*/(.*)", query_path).captures[1]
    end
end

function graph_grouped_boxplot_with_comparison_methods(experiment_params_list::Vector{ExperimentParams};
                                        y_type::VALUE=estimate_error,
                                        grouping::GROUP=technique,
                                        ylims = [10^-7, 10^7],
                                        x_label=nothing, y_label=nothing, filename=nothing)
    # for now let's just use the dataset as the x-values and the cycle size as the groups
    x_values = []
    y_values = []
    estimators = []
    true_card = Dict()
    for experiment_params in experiment_params_list
        # load the results
        results_filename = params_to_results_filename(experiment_params)
        results_path = "Experiments/Results/Estimation_" * results_filename
        results_df = CSV.read(results_path, DataFrame; normalizenames=true)
        # get the x_value and grouping (same for all results in this experiment param)

        # keep track of the data points
        for i in 1:nrow(results_df)
            current_x = string(get_value_from_param(experiment_params, dataset))

            current_group = string(grouping) * '_' * string(grouping == query_type ? results_df[i, :QueryType] : get_value_from_param(experiment_params, grouping))
            current_y = 0
            if y_type == estimate_error
                current_y = results_df[i, :Estimate] / results_df[i, :TrueCard]
            else # y_type == runtime
                current_y = results_df[i, :EstimationTime]
            end
            true_card[(current_x, get_query_id(string(experiment_params.dataset), results_df[i, :QueryPath]))] = results_df[i, :TrueCard]
            # push the errors and their groupings into the correct vector
            push!(x_values, current_x)
            push!(y_values, current_y)
            push!(estimators, current_group)
        end
    end
    results_filename = params_to_results_filename(experiment_params_list[1])
    comparison_results = comparison_dataset()
    for i in 1:nrow(comparison_results)
        current_x = comparison_results[i, :Dataset]
        estimator = comparison_results[i, :Estimator]
        current_y = 0
        if !haskey(true_card, (current_x, comparison_results[i, :Query]))
            continue
        end
        if y_type == estimate_error
            current_y = max(1, comparison_results[i, :Value]) / true_card[(current_x, comparison_results[i, :Query])]
        else # y_type == runtime
            current_y = comparison_results[i, :Runtime] / 1000.0
        end
        # push the errors and their groupings into the correct vector
        push!(x_values, current_x)
        push!(y_values, current_y)
        push!(estimators, estimator)
    end
    println("starting graphs")

    # This seems to be necessary for using Plots.jl outside of the ipynb framework.
    # See this: https://discourse.julialang.org/t/deactivate-plot-display-to-avoid-need-for-x-server/19359/15
    ENV["GKSwstype"]="100"
    gbplot = groupedboxplot(x_values, y_values, group = estimators, yscale =:log10,
                            ylims=ylims,
                            legend = :outertopleft, size = (1000, 600))
    x_label !== nothing && xlabel!(gbplot, x_label)
    y_label !== nothing && ylabel!(gbplot, y_label)
    plotname = (isnothing(filename)) ? results_filename * ".png" : filename * ".png"
    savefig(gbplot, "Experiments/Results/Figures/" * plotname)
end

function graph_grouped_bar_plot(experiment_params_list::Vector{ExperimentParams};
                                        x_type::GROUP=dataset,
                                        y_type::VALUE=estimate_error,
                                        grouping::GROUP=technique,
                                        x_label=nothing,
                                        y_label=nothing,
                                        y_lims=[0, 10],
                                        filename=nothing)
    # for now let's just use the dataset as the x-values and the cycle size as the groups
    x_values = []
    y_values = Float64[]
    groups = []
    for experiment_params in experiment_params_list
        # load the results
        results_filename = params_to_results_filename(experiment_params)
        prefix = "Experiments/Results/Estimation_"
        if y_type == memory_footprint || y_type == build_time
            prefix = "Experiments/Results/Build_"
        end
        results_path = prefix * results_filename
        results_df = CSV.read(results_path, DataFrame; normalizenames=true)

        # get the x_value and grouping (same for all results in this experiment param)
        println(results_df)
        # keep track of the data points
        for i in 1:nrow(results_df)
            current_x = nothing
            if x_type == query_type
                current_x = results_df[i, :QueryType]
            else
                current_x = get_value_from_param(experiment_params, x_type)
            end
            current_group = nothing
            if grouping == query_type
                current_group = results_df[i, :QueryType]
            elseif grouping == build_phase
                current_group = results_df[i, :BuildPhase]
            else
                current_group = get_value_from_param(experiment_params, grouping)
            end
            current_y = 0
            if y_type == estimate_error
                current_y = results_df[i, :Estimate] / results_df[i, :TrueCard]
            elseif y_type == memory_footprint
                current_y = results_df[i, :MemoryFootprint]/(10^6)
            elseif y_type == build_time
                current_y = results_df[i, :BuildTime]
            else
                     # y_type == runtime
                current_y = results_df[i, :EstimationTime]
            end
            # push the errors and their groupings into the correct vector
            push!(x_values, current_x)
            push!(y_values, current_y)
            push!(groups, current_group)
        end
    end
    results_filename = params_to_results_filename(experiment_params_list[1])
    println("starting graphs")

    # This seems to be necessary for using Plots.jl outside of the ipynb framework.
    # See this: https://discourse.julialang.org/t/deactivate-plot-display-to-avoid-need-for-x-server/19359/15
    ENV["GKSwstype"]="100"
    println(x_values)
    println(y_values)
    println(groups)
    gbplot = StatsPlots.groupedbar(x_values,
                            y_values,
                            group = groups,
                            ylims=y_lims,
#                            yticks = [1, 10^.5, 10, 10^2, 10^2.5],
                            legend = :outertopleft,
                            size = (1000, 600))
    x_label !== nothing && xlabel!(gbplot, x_label)
    y_label !== nothing && ylabel!(gbplot, y_label)
    plotname = (isnothing(filename)) ? results_filename * ".png" : filename * ".png"
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
    elseif value_type == sampling_type
        return experiment_param.sampling_strategy
    elseif value_type == cycle_stats
        return experiment_param.only_shortest_path_cycle
    elseif value_type == number_of_colors
        return experiment_param.summary_params.num_colors
    else
        # default to grouping by technique
        return (experiment_param.summary_params.partitioner, experiment_param.summary_params.label_refining_rounds)
    end
end
