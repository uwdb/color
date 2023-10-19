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
abstract type Comparable end
import Base .==
function ==(a::T, b::T) where T <: Comparable
    (a.path == b.path) && (a.colors == b.colors)
end
@auto_hash_equals mutable struct CyclePathAndColors
    path::BoolPath
    colors::StartEndColorPair
end

@enum PARTITIONER QuasiStable Hash Degree DirectedDegree SimpleLabel InOut LabelInOut NeighborEdges MostNeighbors


struct ColorSummaryParams
    num_colors::Int
    max_cycle_size::Int
    max_partial_paths::Int
    partitioner::PARTITIONER
    weighting::Bool
    label_refining_rounds::Int

    function ColorSummaryParams(;num_colors::Int=64, max_cycle_size=4, max_partial_paths=1000,
            partitioner::PARTITIONER = QuasiStable, weighting=true, label_refining_rounds = 0)
        return new(num_colors, max_cycle_size, max_partial_paths, partitioner, weighting, label_refining_rounds)
    end
end

function params_to_string(params::ColorSummaryParams)
    summary_name = "ColorSummary_" * string(params.partitioner) * "_"
    summary_name *= string(params.num_colors) * "_"
    summary_name *= string(params.max_cycle_size) * "_"
    summary_name *= string(params.max_partial_paths)* "_"
    summary_name *= string(params.label_refining_rounds)
    return summary_name
end



include("PropertyGraph.jl")
include("datasets.jl")
include("utils.jl")
include("ExactSizeCalculator.jl")
include("ColoringMethods.jl")
include("ColorSummary.jl")
include("QuasiStableCardinalityEstimator.jl")
