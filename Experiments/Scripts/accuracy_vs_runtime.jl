# The goal of this file is to produce graphs which illustrate the runtime/accuracy
# tradeoff for different configurations. For this, each configuration will be a point on
# the 2d graph where one axis is total runtime and one axis is median/75th percentile q-error

# The first graph demonstrates this tradeoff when we vary the number of number colors used
# to build the summary. The tricky aspect is how to handle the number of partial paths
# to keep during inference. If we keep this number constant, then the runtime won't change
# significantly and the error would just increase. Instead, we choose 2 * num_colors to
# scale it up proportionally.

# The second graph demonstrates this tradeoff when we vary the number of partial paths alone
# during inference.

using Profile
include("../Experiments.jl")

function generate_num_colors_graph(dataset::DATASET)
    num_colors = [4, 8, 16, 32, 64, 128]
    p50_q_errors = []
    p50_runtimes = []
    p95_q_errors = []
    p95_runtimes = []
    p99_q_errors = []
    p99_runtimes = []

    for colors in num_colors
        experiment_params = ExperimentParams(dataset=dataset; partitioning_scheme=[(QuasiStable, colors)])
        build_experiments([experiment_params])
        println("Num Colors: ", colors)
        run_estimation_experiments([experiment_params])
        results_filename = params_to_results_filename(experiment_params)
        results_path = "Experiments/Results/Estimation_" * results_filename
        results_df = CSV.read(results_path, DataFrame; normalizenames=true)
        p50_runtime = quantile([runtime for runtime in results_df[:, :EstimationTime]], .50)
        p95_runtime = quantile([runtime for runtime in results_df[:, :EstimationTime]], .95)
        p99_runtime = quantile([runtime for runtime in results_df[:, :EstimationTime]], .99)
        push!(p50_runtimes, p50_runtime)
        push!(p95_runtimes, p95_runtime)
        push!(p99_runtimes, p99_runtime)
        estimates = [estimate for estimate in results_df[:, :Estimate]]
        true_cards = [true_card for true_card in results_df[:, :TrueCard]]
        p50_q_error = quantile([10.0^abs(log10(estimates[i]/float(true_cards[i]))) for i in eachindex(estimates)], .50)
        p95_q_error = quantile([10.0^abs(log10(estimates[i]/float(true_cards[i]))) for i in eachindex(estimates)], .95)
        p99_q_error = quantile([10.0^abs(log10(estimates[i]/float(true_cards[i]))) for i in eachindex(estimates)], .99)
        push!(p50_q_errors, p50_q_error)
        push!(p95_q_errors, p95_q_error)
        push!(p99_q_errors, p99_q_error)
    end

    # This seems to be necessary for using Plots.jl outside of the ipynb framework.
    # See this: https://discourse.julialang.org/t/deactivate-plot-display-to-avoid-need-for-x-server/19359/15
    ENV["GKSwstype"]="100"
    fig = Plot(p50_runtimes, p50_q_errors, seriestype=:scatter, yscale=:log10, xscale=:log10,
                xlabel="P50 Runtime (sec)", ylabel="P50 Q-Error", title="Num Colors "*string(dataset))
    savefig(fig, "Experiments/Results/Figures/NumColors_Accuracy_vs_Runtime_P50_" * string(dataset) * ".png")
    fig = Plot(p95_runtimes, p95_q_errors, seriestype=:scatter, yscale=:log10, xscale=:log10,
                xlabel="P95 Runtime (sec)", ylabel="P95 Q-Error", title="Num Colors "*string(dataset))
    savefig(fig, "Experiments/Results/Figures/NumColors_Accuracy_vs_Runtime_P95_" * string(dataset) * ".png")
    fig = Plot(p99_runtimes, p99_q_errors, seriestype=:scatter, yscale=:log10, xscale=:log10,
                xlabel="P99 Runtime (sec)", ylabel="P99 Q-Error", title="Num Colors "*string(dataset))
    savefig(fig, "Experiments/Results/Figures/NumColors_Accuracy_vs_Runtime_P99_" * string(dataset) * ".png")
end

function generate_partial_paths_graph(dataset::DATASET; num_colors = 64)
    partial_paths = [2, 4, 8, 16, 32, 64, 128, 256, 512]
    p50_q_errors = []
    p50_runtimes = []
    p95_q_errors = []
    p95_runtimes = []
    p99_q_errors = []
    p99_runtimes = []
    build_params = ExperimentParams(dataset=dataset, partitioning_scheme=[(QuasiStable, num_colors)])
    build_experiments([build_params])

    for pp in partial_paths
        println("Partial Paths: ", pp)
        inference_params = ExperimentParams(dataset=dataset, partitioning_scheme=[(QuasiStable, num_colors)], inference_max_paths=pp)
        run_estimation_experiments([inference_params])
        results_filename = params_to_results_filename(inference_params)
        results_path = "Experiments/Results/Estimation_" * results_filename
        results_df = CSV.read(results_path, DataFrame; normalizenames=true)
        p50_runtime = quantile([runtime for runtime in results_df[:, :EstimationTime]], .50)
        p95_runtime = quantile([runtime for runtime in results_df[:, :EstimationTime]], .95)
        p99_runtime = quantile([runtime for runtime in results_df[:, :EstimationTime]], .99)
        push!(p50_runtimes, p50_runtime)
        push!(p95_runtimes, p95_runtime)
        push!(p99_runtimes, p99_runtime)
        estimates = [estimate for estimate in results_df[:, :Estimate]]
        true_cards = [true_card for true_card in results_df[:, :TrueCard]]
        p50_q_error = quantile([10.0^abs(log10(estimates[i]/float(true_cards[i]))) for i in eachindex(estimates)], .50)
        p95_q_error = quantile([10.0^abs(log10(estimates[i]/float(true_cards[i]))) for i in eachindex(estimates)], .95)
        p99_q_error = quantile([10.0^abs(log10(estimates[i]/float(true_cards[i]))) for i in eachindex(estimates)], .99)
        push!(p50_q_errors, p50_q_error)
        push!(p95_q_errors, p95_q_error)
        push!(p99_q_errors, p99_q_error)
    end

    # This seems to be necessary for using Plots.jl outside of the ipynb framework.
    # See this: https://discourse.julialang.org/t/deactivate-plot-display-to-avoid-need-for-x-server/19359/15
    ENV["GKSwstype"]="100"
    println(p95_runtimes)
    println(p95_q_errors)

    fig = Plot(p50_runtimes, p50_q_errors, seriestype=:scatter, yscale=:log10, xscale=:log10,
    xlabel="P50 Runtime (sec)", ylabel="P50 Q-Error", title="Partial Paths "*string(dataset))
    savefig(fig, "Experiments/Results/Figures/PP_Accuracy_vs_Runtime_P50_" * string(dataset) * ".png")
    fig = Plot(p95_runtimes, p95_q_errors, seriestype=:scatter, yscale=:log10, xscale=:log10,
        xlabel="P95 Runtime (sec)", ylabel="P95 Q-Error", title="Partial Paths "*string(dataset))
    savefig(fig, "Experiments/Results/Figures/PP_Accuracy_vs_Runtime_P95_" * string(dataset) * ".png")
    fig = Plot(p99_runtimes, p99_q_errors, seriestype=:scatter, yscale=:log10, xscale=:log10,
        xlabel="P99 Runtime (sec)", ylabel="P99 Q-Error", title="Partial Paths "*string(dataset))
    savefig(fig, "Experiments/Results/Figures/PP_Accuracy_vs_Runtime_P99_" * string(dataset) * ".png")
end

#generate_num_colors_graph(human)
#generate_num_colors_graph(aids)
#generate_num_colors_graph(hprd)
#generate_num_colors_graph(yeast)
generate_partial_paths_graph(human)
generate_partial_paths_graph(aids)
generate_partial_paths_graph(hprd)
generate_partial_paths_graph(yeast)
