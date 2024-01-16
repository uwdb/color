@enum GROUP dataset technique cycle_size summary_paths inference_paths query_type sampling_type cycle_stats number_of_colors build_phase proportion_updated proportion_deleted deg_stat_type description query_size path_width

@enum VALUE estimate_error runtime build_time memory_footprint update_time total_time

function graph_box_plot(experiment_params_list::Vector{ExperimentParams};
                                x_type::GROUP=dataset,
                                y_type::VALUE=estimate_error,
                                grouping::GROUP=technique,
                                x_order = nothing,
                                legend_pos = :outertopleft,
                                x_label=nothing,
                                y_label=nothing,
                                filename=nothing,
                                ylims=[0, 10^2.5],
                                y_ticks=[10^-10, 10^-5, 10^-2, 1, 10^2, 10^5, 10^10],
                                dimensions = (1000, 800))
    # for now let's just use the dataset as the x-values and the cycle size as the groups
    x_values = []
    y_values = []
    for experiment_params in experiment_params_list
        # load the results
        results_filename = params_to_results_filename(experiment_params)
        results_path = "Experiments/Results/Estimation_" * results_filename
        # println("results path: ", results_path)
        results_df = CSV.read(results_path, DataFrame; normalizenames=true)

        # keep track of the data points
        for i in 1:nrow(results_df)
            current_x = nothing
            if x_type == query_type
                current_x = results_df[i, :QueryType]
            elseif x_type == path_width
                current_x = results_df[i, :PathWidth]
            else
                current_x = get_value_from_param(experiment_params, x_type)
            end
            # current_x = x_type == query_type ? results_df[i, :QueryType] : get_value_from_param(experiment_params, x_type)
            current_y = 0
            if y_type == estimate_error
                current_y = results_df[i, :Estimate] / results_df[i, :TrueCard]
            else # y_type == runtime
                current_y = results_df[i, :EstimationTime]
            end
            # push the errors and their groupings into the correct vector
            push!(x_values, current_x)
            push!(y_values, current_y)
        end
    end

    if isnothing(x_order)
        x_order = sort(unique(x_values))
    end
    x_ticks = ([x for x in 1:length(x_order)], x_order)
    x_order_hash = [hash(x) for x in x_order]
    x_values = [only(indexin(hash(x), x_order_hash)) for x in x_values]
    results_filename = params_to_results_filename(experiment_params_list[1])
    println("starting graphs")

    # This seems to be necessary for using Plots.jl outside of the ipynb framework.
    # See this: https://discourse.julialang.org/t/deactivate-plot-display-to-avoid-need-for-x-server/19359/15
    ENV["GKSwstype"]="100"
    gbplot = boxplot(x_values,
                            [log10(y)  for y in y_values],
                            x_ticks = x_ticks,
                            xlims = [0, length(x_order) + .5],
                            ylims =  (log10(ylims[1]),log10(ylims[2])),
                            y_ticks = [log10(y) for y in y_ticks],
                            legend = false,
                            size = dimensions,
                            bottom_margin = 20px,
                            top_margin = 20px,
                            left_margin = 10mm,
                            titlefont = (12, :black),
                            tickfont = (12, :black),
                            guidefont = (15, :black),
                            whisker_range=2)
    x_label !== nothing && xlabel!(gbplot, x_label)
    y_label !== nothing && ylabel!(gbplot, y_label)
    plotname = (isnothing(filename)) ? results_filename * ".png" : filename * ".png"
    savefig(gbplot, "Experiments/Results/Figures/" * plotname)
end

function graph_grouped_box_plot(experiment_params_list::Vector{ExperimentParams};
                                        x_type::GROUP=dataset,
                                        y_type::VALUE=estimate_error,
                                        grouping::GROUP=technique,
                                        x_order = nothing,
                                        legend_order = nothing,
                                        legend_pos = :outertopleft,
                                        x_label=nothing,
                                        y_label=nothing,
                                        filename=nothing,
                                        ylims=[0, 10^2.5],
                                        y_ticks=[10^-10, 10^-5, 10^-2, 1, 10^2, 10^5, 10^10],
                                        dimensions = (1000, 800))
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
            push!(x_values, string(current_x))
            push!(y_values, current_y)
            push!(groups, current_group)
        end
    end

    if isnothing(x_order)
        x_order = sort(unique(x_values))
    end
    if isnothing(legend_order)
        legend_order = Vector{String}(collect(reverse(sort(unique(groups)))))
    end
    x_ticks = ([x for x in 1:length(x_order)], x_order)
    x_order_hash = [hash(x) for x in x_order]
    x_values = [only(indexin(hash(x), x_order_hash)) for x in x_values]
    sorted_vals = sort(collect(zip(x_values, y_values, groups)), by=(x)->x[1])
    x_values = [x[1] for x in sorted_vals]
    y_values = [x[2] for x in sorted_vals]
    groups = [x[3] for x in sorted_vals]
    group_order_hash = [hash(x) for x in legend_order]
    groups = [only(indexin(hash(x), group_order_hash)) for x in groups]
    results_filename = params_to_results_filename(experiment_params_list[1])
    println("starting graphs")

    # This seems to be necessary for using Plots.jl outside of the ipynb framework.
    # See this: https://discourse.julialang.org/t/deactivate-plot-display-to-avoid-need-for-x-server/19359/15
    ENV["GKSwstype"]="100"
    gbplot = groupedboxplot(x_values,
                            [log10(y)  for y in y_values],
                            group = groups,
                            x_ticks = x_ticks,
                            xlims = [0.5, length(x_order)+.5],
                            ylims =  (log10(ylims[1]),log10(ylims[2])),
                            y_ticks = [log10(y) for y in y_ticks],
                            legend = legend_pos,
                            labels = reshape(legend_order, 1, length(legend_order)),
                            size = dimensions,
                            bottom_margin = 20px,
                            top_margin = 20px,
                            left_margin = 10mm,
                            legend_column = 2,
                            titlefont = (12, :black),
                            legendfont = (11, :black),
                            tickfont = (12, :black),
                            guidefont = (15, :black),
                            outliers=false,
                            whisker_range=2)
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
            comparison_results[i, :QueryType] = match(r".*/query_(.*)_.*", query_path).captures[1]
        end
    end
    results_dict = Dict()
    for i in 1:nrow(comparison_results)
        dataset = comparison_results[i, :Dataset]
        estimator = comparison_results[i, :Estimator]
        query_path = comparison_results[i, :Query]
        results_dict[(dataset, estimator, query_path)] = (Estimate=comparison_results[i, :Value],
                                                            Runtime=comparison_results[i, :Runtime],
                                                            QueryType=comparison_results[i,:QueryType])
    end
    estimators = unique(comparison_results[:, :Estimator])
    return estimators, results_dict
end

function get_query_id(dataset, query_path)
    return if dataset == "lubm80"
        match(r".*/queryset/(.*/.*)", query_path).captures[1]
    elseif dataset in ["aids", "human", "yago"]
        match(r".*/queryset/.*/(.*/.*)", query_path).captures[1]
    else
        match(r".*/(.*/.*).graph", query_path).captures[1]*".txt"
    end
end

function query_size_category(s)
    categories = [3, 4, 6, 9, 12, 16, 24, 32]
    for cat in categories
        if s <= cat
            return cat
        end
    end
end

function graph_grouped_boxplot_with_comparison_methods(experiment_params_list::Vector{ExperimentParams};
                                        x_order = nothing,
                                        legend_order = nothing,
                                        x_type::GROUP=dataset,
                                        y_type::VALUE=estimate_error,
                                        grouping::GROUP=technique,
                                        ylims = [10^-7, 10^7],
                                        y_ticks = [10^-5, 10^-2, 10^0, 10^2, 10^5],
                                        x_label=nothing,
                                        y_label=nothing,
                                        filename=nothing,
                                        legend_pos = :outertopleft,
                                        dimensions = (1000, 600))
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
            data = string(get_value_from_param(experiment_params, dataset))
            current_x = if x_type == dataset
                data
            elseif x_type == query_size
                @sprintf "%02i" query_size_category(results_df[i, :QuerySize])
            end
            current_group = string(grouping == query_type ? results_df[i, :QueryType] : get_value_from_param(experiment_params, grouping))
            current_y = if y_type == estimate_error
                min(10^30, max(1, results_df[i, :Estimate])) / results_df[i, :TrueCard]
            else # y_type == runtime
                results_df[i, :EstimationTime]
            end
            true_card[(data, get_query_id(string(experiment_params.dataset), results_df[i, :QueryPath]))] = (results_df[i, :TrueCard], current_x)
            # push the errors and their groupings into the correct vector
            push!(x_values, string(current_x))
            push!(y_values, current_y)
            push!(estimators, current_group)
        end
    end
    results_filename = params_to_results_filename(experiment_params_list[1])
    estimator_types, comparison_results = comparison_dataset()
    estimator_dataset_missing = Set()
    for (query_key, query_card_and_size) in true_card
        data = query_key[1]
        query_path = query_key[2]
        card = query_card_and_size[1]
        size = query_card_and_size[2]
        for estimator in estimator_types
            comp_key = (data, estimator, query_path)
            (estimate, runtime) = 1, 10 # TODO: We shouldn't use an arbitrary number for runtime here
            if haskey(comparison_results, comp_key)
                result = comparison_results[comp_key]
                estimate = result.Estimate
                runtime = result.Runtime
            else
                push!(estimator_dataset_missing, (estimator, data))
                continue
            end

            current_x = if x_type == dataset
                data
            elseif x_type == query_size
                size
            end

            current_y = if y_type == estimate_error
                min(10^30, max(1, estimate)) / card
            else # y_type == runtime
                runtime / 1000.0
            end
            # push the errors and their groupings into the correct vector
            push!(x_values, string(current_x))
            push!(y_values, current_y)
            push!(estimators, estimator)
        end
    end

    if isnothing(x_order)
        x_order = sort(unique(x_values))
    end
    x_ticks = ([x for x in 1:length(x_order)], x_order)
    x_order_hash = [hash(x) for x in x_order]
    x_values = [only(indexin(hash(x), x_order_hash)) for x in x_values]
    sorted_vals = sort(collect(zip(x_values, y_values, estimators)), by=(x)->x[1])
    x_values = [x[1] for x in sorted_vals]

    if isnothing(legend_order)
        legend_order = Vector{String}(collect(reverse(sort(unique(groups)))))
    end
    groups = [x[3] for x in sorted_vals]
    group_order_hash = [hash(x) for x in legend_order]
    groups = [only(indexin(hash(x), group_order_hash)) for x in groups]

    y_values = [x[2] for x in sorted_vals]
    println("starting graphs")
    # This seems to be necessary for using Plots.jl outside of the ipynb framework.
    # See this: https://discourse.julialang.org/t/deactivate-plot-display-to-avoid-need-for-x-server/19359/15
    ENV["GKSwstype"]="100"
    gbplot =  groupedboxplot(x_values,
                        [log10(y)  for y in y_values],
                        group = groups,
                        x_ticks = x_ticks,
                        xlims = [0.5, length(x_order)+.5],
                        ylims =  (log10(ylims[1]),log10(ylims[2])),
                        y_ticks = [log10(y) for y in y_ticks],
                        legend = legend_pos,
                        labels = reshape(legend_order, 1, length(legend_order)),
                        size = dimensions,
                        bottom_margin = 40px,
                        top_margin = 20px,
                        left_margin = 10mm,
                        legend_column = 2,
                        titlefont = (12, :black),
                        legendfont = (11, :black),
                        tickfont = (12, :black),
                        guidefont = (15, :black),
                        whisker_range=2,
                        outliers = false)
    x_label !== nothing && xlabel!(gbplot, x_label)
    y_label !== nothing && ylabel!(gbplot, y_label)
    y_type == estimate_error && hline!([0], label="exact", linestyle=:solid, lw=2, alpha=.4, color="red")
    plotname = (isnothing(filename)) ? results_filename * ".png" : filename * ".png"
    savefig(gbplot, "Experiments/Results/Figures/" * plotname)
end

function graph_grouped_bar_plot(experiment_params_list::Vector{ExperimentParams};
                                        x_type::GROUP=dataset,
                                        y_type::VALUE=estimate_error,
                                        grouping::GROUP=technique,
                                        x_order = nothing,
                                        legend_order = nothing,
                                        x_label=nothing,
                                        y_label=nothing,
                                        ylims=[0, 10^2.5],
                                        y_ticks = [1, 10^.5, 10, 10^2, 10^2.5],
                                        dimensions = (800, 300),
                                        legend_pos = :topleft,
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
                if grouping != build_phase && results_df[i, :BuildPhase] != "FullTime"
                    continue
                end
                current_y = results_df[i, :BuildTime]
            else
                     # y_type == runtime
                current_y = results_df[i, :EstimationTime]
            end
            # push the errors and their groupings into the correct vector
            push!(x_values, string(current_x))
            push!(y_values, current_y)
            push!(groups, current_group)
        end
    end
    if isnothing(x_order)
        x_order = collect(sort(unique(x_values)))
    end
    x_ticks = ([x for x in 1:length(x_order)], x_order)
    x_order_hash = [hash(x) for x in x_order]
    x_values = [only(indexin(hash(x), x_order_hash)) for x in x_values]

    if isnothing(legend_order)
        legend_order = Vector{String}(collect(reverse(sort(unique(groups)))))
    end
    group_order_hash = [hash(x) for x in legend_order]
    groups = [only(indexin(hash(x), group_order_hash)) for x in groups]

    results_filename = params_to_results_filename(experiment_params_list[1])
    println("starting graphs")

    # This seems to be necessary for using Plots.jl outside of the ipynb framework.
    # See this: https://discourse.julialang.org/t/deactivate-plot-display-to-avoid-need-for-x-server/19359/15
    ENV["GKSwstype"]="100"
    gbplot = StatsPlots.groupedbar(x_values,
                            y_values,
                            group = groups,
                            x_ticks = x_ticks,
                            ylims=ylims,
                            xlims = [0.5, length(x_order)+.5],
                            y_ticks = y_ticks,
                            legend = legend_pos,
                            labels = reshape(legend_order, 1, length(legend_order)),
                            size = dimensions,
                            thickness_scaling=1.25,
                            titlefont = (12, :black),
                            legendfont = (11, :black),
                            tickfont = (12, :black),
                            guidefont = (15, :black),
                            left_margin = 10mm,
                            legend_column = 2,)
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
    elseif value_type == proportion_updated
        return experiment_param.summary_params.proportion_updated
    elseif value_type == proportion_deleted
        return experiment_param.summary_params.proportion_deleted
    elseif value_type == deg_stat_type
        return experiment_param.summary_params.deg_stats_type
    elseif value_type == description
        return experiment_param.description
    else
        # default to grouping by technique
        return experiment_param.summary_params.partitioning_scheme
    end
end
