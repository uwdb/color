using Probably: BloomFilter, constrain
using DataStructures: counter, Dict, Set, Vector, inc!, Queue
using AutoHashEquals
using QuasiStableColors
using StatsBase

QSC = QuasiStableColors
using Graphs: SimpleDiGraphFromIterator, Edge, DiGraph, edges, nv, ne, add_edge!,
                add_vertex!, vertices, all_neighbors, src, dst, outneighbors, inneighbors

BoolPath = Vector{Bool}
NodeId = Int
Color = Int16
StartEndColorPair = Tuple{Color, Color}
@auto_hash_equals mutable struct CyclePathAndColors
    path::BoolPath
    colors::StartEndColorPair
end

@enum PARTITIONER QuasiStable Hash Degree NeighborNodeLabels

function partitioner_to_string(x::PARTITIONER)
    return if x == QuasiStable
        "QS"
    elseif x == Hash
        "H"
    elseif x == Degree
        "D"
    elseif x == NeighborNodeLabels
        "NNL"
    end
end

PartitioningScheme = Vector{Tuple{PARTITIONER, Int}}

function Base.show(io::IO, x::Vector{Tuple{PARTITIONER, Int}})
    output = "["
    prefix = ""
    for (p, n) in x
        output*= prefix * partitioner_to_string(p) * ":" * string(n)
        prefix = ";"
    end
    output *= "]"
    show(io, output)
end

struct ColorSummaryParams
    deg_stats_type::Type
    num_colors::Int
    max_cycle_size::Int
    max_partial_paths::Int
    partitioning_scheme::PartitioningScheme
    weighting::Bool
    proportion_not_updated::Float16
    proportion_deleted::Float16

    function ColorSummaryParams(;deg_stats_type = AvgDegStats, max_cycle_size=4, max_partial_paths=1000,
        partitioning_scheme::PartitioningScheme = [(QuasiStable, 64)], weighting=true, proportion_not_updated = 1.0, proportion_deleted=0.0)
        num_colors = sum([x[2] for x in partitioning_scheme])
        return new(deg_stats_type, num_colors, max_cycle_size, max_partial_paths, partitioning_scheme, weighting, proportion_not_updated, proportion_deleted)
    end
end

function params_to_string(params::ColorSummaryParams)
    summary_name = "ColorSummary_" * string(params.deg_stats_type) * "_"
    summary_name *= string(params.partitioning_scheme) * "_"
    summary_name *= string(params.max_cycle_size) * "_"
    summary_name *= string(params.max_partial_paths)* "_"
    summary_name *= string(params.proportion_not_updated) * "_"
    summary_name *= string(params.proportion_deleted)
    return summary_name
end

include("PropertyGraph.jl")
include("datasets.jl")
include("utils.jl")
include("ExactSizeCalculator.jl")
include("ColoringMethods.jl")
include("DegreeStats.jl")
include("ColorSummary.jl")
include("QuasiStableCardinalityEstimator.jl")
