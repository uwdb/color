# Overall Experiments Harness Include File
using BenchmarkTools
using Plots
using Plots.PlotMeasures
using Printf
using Random
using StatsPlots
using Distributed
@everywhere using CSV
@everywhere using DataFrames
@everywhere using DelimitedFiles: writedlm
@everywhere using Parquet2: Dataset
@everywhere using Random
@everywhere using Serialization: serialize, deserialize
@everywhere using SharedArrays
@everywhere using WeakRefStrings


@everywhere include("../src/CardinalityWithColors.jl")
@everywhere include("utils.jl")
@everywhere include("load_datasets.jl")
include("load_querysets.jl")
@everywhere include("build_color_summaries.jl")
include("get_true_cardinalities.jl")
@everywhere include("run_estimators.jl")
include("graph_results.jl")

const TIMEOUT_SEC::Float64 = 60.0
