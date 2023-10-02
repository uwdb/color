using Probably: BloomFilter, constrain
using DataStructures: counter, Dict, Set, Vector, inc!, Queue
using AutoHashEquals
using QuasiStableColors
using StatsBase

QSC = QuasiStableColors
using Graphs: SimpleDiGraphFromIterator, Edge, DiGraph, edges, nv, ne, add_edge!,
                add_vertex!, vertices, all_neighbors, src, dst, outneighbors, inneighbors





include(pwd() * "/Source/PropertyGraph.jl")
include(pwd() * "/Source/datasets.jl")
include(pwd() * "/Source/utils.jl")
include(pwd() * "/Source/ExactSizeCalculator.jl")
include(pwd() * "/Source/ColorSummary.jl")
include(pwd() * "/Source/QuasiStableCardinalityEstimator.jl")
