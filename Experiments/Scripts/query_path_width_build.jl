using Plots.PlotMeasures
include("../Experiments.jl")

current_dataset = youtube
max_paths = -1

# The goal of this file is to demonstrate the significance of the partial sum optimization.
# We use the same datasets and summaries but we try estimating without partial sums, with partial
# sums, and with partial sums and sampling, then record how the inference time changes.

experiment_params_list::Vector{ExperimentParams} = [ExperimentParams(dataset=current_dataset, use_partial_sums=true, inference_max_paths=max_paths),
                                                    ExperimentParams(dataset=current_dataset, use_partial_sums=false, inference_max_paths=max_paths),
                                                    ExperimentParams(dataset=current_dataset, use_partial_sums=true)]
println("started building")
build_experiments(experiment_params_list)
println("started estimating")
run_estimation_experiments(experiment_params_list, timeout=TIMEOUT_SEC)
println("started graphing")

x_values = []
y_values = []
groups = []
for experiment_params in experiment_params_list
    # load the results
    results_filename = params_to_results_filename(experiment_params)
    results_path = "Experiments/Results/Estimation_" * results_filename
    results_df = CSV.read(results_path, DataFrame; normalizenames=true)
    # keep track of the data points
    for i in 1:nrow(results_df)
        failure = results_df[i, :Failure]
        if !failure
            current_x = results_df[i, :PathWidth]
            current_group = "No Partial Agg"
            if experiment_params.use_partial_sums
                current_group = "Partial Agg"
                if experiment_params.inference_max_paths != max_paths
                    current_group *= " + Sampling"
                end
            end
            current_y = results_df[i, :EstimationTime]
            push!(x_values, current_x)
            push!(y_values, current_y)
            push!(groups, current_group)
        end
    end
end
x_order = sort(unique(x_values))
x_ticks = ([x for x in 1:length(x_order)], x_order)
x_order_hash = [hash(x) for x in x_order]
x_values = [only(indexin(hash(x), x_order_hash)) for x in x_values]
results_filename = params_to_results_filename(experiment_params_list[1])

# This seems to be necessary for using Plots.jl outside of the ipynb framework.
# See this: https://discourse.julialang.org/t/deactivate-plot-display-to-avoid-need-for-x-server/19359/15
ENV["GKSwstype"]="100"
ylims=[10^-3, 10^4]
y_ticks=[10^-2, 10^-1, 1, 10, 10^2, 10^2, 10^3]
dimensions = (600, 400)
gbplot = groupedboxplot(x_values,
                        [log10(y) for y in y_values],
                        group = groups,
                        x_ticks = x_ticks,
                        xlims = [0, length(x_order) + .5],
                        ylims =  (log10(ylims[1]),log10(ylims[2])),
                        y_ticks = [log10(y) for y in y_ticks],
                        legend = :topleft,
                        size = dimensions,
                        bottom_margin = 20px,
                        top_margin = 20px,
                        left_margin = 10mm,
                        legend_column = 1,
                        titlefont = (12, :black),
                        legendfont = (11, :black),
                        tickfont = (12, :black),
                        guidefont = (15, :black),
                        whisker_range=2)
xlabel!(gbplot, "Query Path Width")
ylabel!(gbplot, "Inference Latency log\$_{10}\$ (s)")
plotname = "fig_13.png"
savefig(gbplot, "Experiments/Results/Figures/" * plotname)
