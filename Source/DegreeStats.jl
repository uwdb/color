# This file defines the interface and several instantiations of degree statistics. These are
# the core statistic that is used in the estimation procedure.


abstract type DegreeStats end

# A degree statistic is derived from a set of edges. We include the data graph as well in
# order to allow access to things like the broader degree of the nodes and the labels.
function DegreeStats(g::DataGraph, edges::Vector{Tuple{NodeId}})
    throw(ErrorException("DegreeStats is an abstract type, you probably meant to call a particular instance."))
end

# In some places in the code, we need to access an estimate of the degree. The nature of this
# estimate will depend on the kind of statistic (e.g. it might be upper bound, avg, or
# lower bound).
function get_out_deg_estimate(d::DegreeStats)
    throw(ErrorException("DegreeStats is an abstract type, you probably meant to call a particular instance."))
end

function scale_out_deg(d::DegreeStats, s::Float64)
    throw(ErrorException("DegreeStats is an abstract type, you probably meant to call a particular instance."))
end

function scale_in_deg(d::DegreeStats, s::Float64)
    throw(ErrorException("DegreeStats is an abstract type, you probably meant to call a particular instance."))
end


abstract type StatAccumulator end

# Every accumulator needs a way to be initialized to represent a single starting color based
# on the cardinality of that color.
function StatAccumulator(c::Float64)
    throw(ErrorException("StatAccumulator is an abstract type, you probably meant to call a particular instance."))
end

# This is the core operation for extending partial colorings to an additional node.
# 'out_edge' denotes whether the edge is going from w->d (true) or d->w (false).
function extend_coloring(w::StatAccumulator, d::DegreeStats, out_edge::Bool)
    throw(ErrorException("DegreeStats & StatAccumulator are an abstract type, you probably meant to call a particular instance."))
end

# This summation is necessary for partial aggregation.
function sum_colorings(w1::StatAccumulator, w2::StatAccumulator)
    throw(ErrorException("StatAccumulator is an abstract type, you probably meant to call a particular instance."))
end

# Multiplication of this kind happens during sampling and handle_extra_edges to scale the
# weight of a particular path up or down.
function scale_coloring(w::StatAccumulator, s::Float64)
    throw(ErrorException("StatAccumulator is an abstract type, you probably meant to call a particular instance."))
end

# This function is the final output of the accumulator and denotes the number of paths in the
# data graph which are represented by this partial coloring.
function get_count(w::StatAccumulator)
    throw(ErrorException("StatAccumulator is an abstract type, you probably meant to call a particular instance."))
end

function stat_type_to_accumulator(D::Type)
    return if D == MinDegStats
        MinAccumulator
    elseif D == AvgDegStats
        AvgAccumulator
    elseif D == MaxDegStats
        MaxAccumulator
    else
        throw(ErrorException(string(D)* " doesn't have an associated accumulator!"))
    end
end

############################ MinDegStats ################################################

mutable struct MinDegStats <:DegreeStats
    min_in::Float64
    min_out::Float64
end

function MinDegStats(g::DataGraph, edges::Vector{Tuple{NodeId, NodeId, Bool}}, color_size::Int)
    length(edges) == 0 && return MinDegStats(0,0)
    in_counter = counter(NodeId)
    out_counter = counter(NodeId)
    for edge in edges
        if edge[3]
            inc!(out_counter, edge[1])
        else
            inc!(in_counter, edge[1])
        end
    end
    return MinDegStats(minimum(values(in_counter); init=0), minimum(values(out_counter); init=0))
end
get_in_deg_estimate(d::MinDegStats) = d.min_in
get_out_deg_estimate(d::MinDegStats) = d.min_out


# In case p is negative, we to reduce the min degree. Otherwise, keep it.
function add_in_deg(d::MinDegStats, p::Float64, c1::Float64, c2::Float64)
    return MinDegStats(max(0, min(d.min_in + ceil(p), d.min_in)), d.min_out)
end
function add_out_deg(d::MinDegStats, p::Float64, c1::Float64, c2::Float64)
    return MinDegStats(d.min_in, max(0, min(d.min_out + ceil(p), d.min_out)))
end


struct MinAccumulator <:StatAccumulator
    weight::Float64
    MinAccumulator(c::Union{Float64, Int64}) = new(c)
end
get_count(w::MinAccumulator) = w.weight
sum_colorings(w1::MinAccumulator, w2::MinAccumulator) = MinAccumulator(w1.weight + w2.weight)



# Because the minimum degree estimator aims to produce a lower bound, we generally can't
# scale up weights during sampling and we have to treat cycle closure as 0 probability.
function scale_coloring(w::MinAccumulator, s::Float64)
    if s >= 1.0
        return w
    else
        return MinAccumulator(0)
    end
end

function extend_coloring(w::MinAccumulator, d::MinDegStats, out_edge::Bool)
    if out_edge
        return MinAccumulator(w.weight*d.min_out)
    else
        return MinAccumulator(w.weight*d.min_in)
    end
end


############################ AvgDegStats ################################################

mutable struct AvgDegStats <:DegreeStats
    avg_in::Float64
    avg_out::Float64
end

function AvgDegStats(g::DataGraph, edges::Vector{Tuple{NodeId, NodeId, Bool}}, color_size::Int)
    length(edges) == 0 && return AvgDegStats(0,0)
    in_counter = counter(NodeId)
    out_counter = counter(NodeId)
    for edge in edges
        if edge[3]
            inc!(out_counter, edge[1])
        else
            inc!(in_counter, edge[1])
        end
    end
    avg_in = sum([x for x in values(in_counter)]; init=0)/color_size
    avg_out = sum([x for x in values(out_counter)]; init=0)/color_size
    return AvgDegStats(avg_in, avg_out)
end
get_in_deg_estimate(d::AvgDegStats) = d.avg_in
get_out_deg_estimate(d::AvgDegStats) = d.avg_out

function add_in_deg(d::AvgDegStats, p::Float64, c1::Float64, c2::Float64)
    return AvgDegStats(max(0, min(d.avg_in + p/c1, c2)), d.avg_out)
end

function add_out_deg(d::AvgDegStats, p::Float64, c1::Float64, c2::Float64)
    return AvgDegStats(d.avg_in, max(0, min(d.avg_out + p/c1, c2)))
end

struct AvgAccumulator <:StatAccumulator
    weight::Float64
    AvgAccumulator(c::Union{Float64, Int64}) = new(c)
end

get_count(w::AvgAccumulator) = w.weight
sum_colorings(w1::AvgAccumulator, w2::AvgAccumulator) = AvgAccumulator(w1.weight + w2.weight)
scale_coloring(w::AvgAccumulator, s::Float64) = AvgAccumulator(w.weight * s)

function extend_coloring(w::AvgAccumulator, d::AvgDegStats, out_edge::Bool)
    if out_edge
        return AvgAccumulator(w.weight*d.avg_out)
    else
        return AvgAccumulator(w.weight*d.avg_in)
    end
end

############################ MaxDegStats ################################################

mutable struct MaxDegStats <:DegreeStats
    max_in::Float64
    max_out::Float64
end

function MaxDegStats(g::DataGraph, edges::Vector{Tuple{NodeId, NodeId, Bool}}, color_size::Int)
    length(edges) == 0 && return MaxDegStats(0,0)
    in_counter = counter(NodeId)
    out_counter = counter(NodeId)
    for edge in edges
        if edge[3]
            inc!(out_counter, edge[1])
        else
            inc!(in_counter, edge[1])
        end
    end
    return MaxDegStats(maximum(values(in_counter); init=0), maximum(values(out_counter); init=0))
end
get_in_deg_estimate(d::MaxDegStats) = d.max_in
get_out_deg_estimate(d::MaxDegStats) = d.max_out

# In case p is negative, we make sure to never reduce the max degree
function add_in_deg(d::MaxDegStats, p::Float64, c1::Float64, c2::Float64)
    return MaxDegStats(min(max(d.max_in, d.max_in + ceil(p)), c2), d.max_out)
end

function add_out_deg(d::MaxDegStats, p::Float64, c1::Float64, c2::Float64)
    return MaxDegStats(d.max_in, min(max(d.max_out, d.max_out) + ceil(p), c2))
end

struct MaxAccumulator <:StatAccumulator
    weight::Float64
    MaxAccumulator(c::Union{Float64, Int64}) = new(c)
end
get_count(w::MaxAccumulator) = w.weight
sum_colorings(w1::MaxAccumulator, w2::MaxAccumulator) = MaxAccumulator(w1.weight + w2.weight)

# Because the max degree estimator aims to produce an upper bound, we have to treat cycle
# closure as 1.0 probability.
function scale_coloring(w::MaxAccumulator, s::Float64)
    if s >= 1.0
        return MaxAccumulator(w.weight * s)
    else
        return w
    end
end

function extend_coloring(w::MaxAccumulator, d::MaxDegStats, out_edge::Bool)
    if out_edge
        return MaxAccumulator(w.weight*d.max_out)
    else
        return MaxAccumulator(w.weight*d.max_in)
    end
end

############################ VarianceDegreeStats ################################################
# TODO: Add the variance stats here
