# Overall Experiments Harness Include File
using Serialization: serialize, deserialize
using Plots
using Plots.PlotMeasures
using StatsPlots
using CSV, DataFrames
using Parquet2: Dataset
using DelimitedFiles: writedlm
using BenchmarkTools

include("../Source/CardinalityWithColors.jl")
include("utils.jl")
include("load_datasets.jl")
include("load_querysets.jl")
include("build_color_summaries.jl")
include("get_true_cardinalities.jl")
include("run_estimators.jl")
include("graph_results.jl")
