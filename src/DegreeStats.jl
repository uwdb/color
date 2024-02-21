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

function get_in_deg_estimate(d::DegreeStats)
    throw(ErrorException("DegreeStats is an abstract type, you probably meant to call a particular instance."))
end

abstract type StatAccumulator end

# Every accumulator needs a way to be initialized to represent a single starting color based
# on the cardinality of that color.
function StatAccumulator(c)
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
function scale_coloring(w::StatAccumulator, s)
    throw(ErrorException("StatAccumulator is an abstract type, you probably meant to call a particular instance."))
end

# This function is the final output of the accumulator and denotes the number of paths in the
# data graph which are represented by this partial coloring.
function get_count(w::StatAccumulator)
    throw(ErrorException("StatAccumulator is an abstract type, you probably meant to call a particular instance."))
end

function get_default_count(D::Type)
    return if D == MinDegStats
        0
    elseif D == AvgDegStats
        0
    elseif D == MaxDegStats
        Inf
    elseif D == CorrDegStats
        0
    else
        throw(ErrorException(string(D)* " doesn't have an associated default count"))
    end
end

function stat_type_to_accumulator(D::Type)
    return if D == MinDegStats
        MinAccumulator
    elseif D == AvgDegStats
        AvgAccumulator
    elseif D == MaxDegStats
        MaxAccumulator
    elseif D == CorrDegStats
        CorrAccumulator
    else
        throw(ErrorException(string(D)* " doesn't have an associated accumulator!"))
    end
end

############################ MinDegStats ################################################

struct MinDegStats <:DegreeStats
    min_in::Float32
    min_out::Float32
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

struct MinAccumulator <:StatAccumulator
    weight::Float32
end
get_count(w::MinAccumulator) = w.weight
sum_colorings(w1::MinAccumulator, w2::MinAccumulator) = MinAccumulator(w1.weight + w2.weight)

# Because the minimum degree estimator aims to produce a lower bound, we generally can't
# scale up weights during sampling and we have to treat cycle closure as 0 probability.
function scale_coloring(w::MinAccumulator, s)
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

struct AvgDegStats <:DegreeStats
    avg_in::Float32
    avg_out::Float32
end

function AvgDegStats(g::DataGraph, edges::Vector{Tuple{NodeId, NodeId, Bool}}, color_size::Int)
    in_edges = 0
    out_edges = 0
    for edge in edges
        if edge[3]
            out_edges += 1
        else
            in_edges += 1
        end
    end
    avg_in = in_edges/color_size
    avg_out = out_edges/color_size
    return AvgDegStats(avg_in, avg_out)
end
get_in_deg_estimate(d::AvgDegStats) = d.avg_in
get_out_deg_estimate(d::AvgDegStats) = d.avg_out

struct AvgAccumulator <:StatAccumulator
    weight::Float32
end

get_count(w::AvgAccumulator) = w.weight
sum_colorings(w1::AvgAccumulator, w2::AvgAccumulator) = AvgAccumulator(w1.weight + w2.weight)
scale_coloring(w::AvgAccumulator, s) = AvgAccumulator(w.weight * s)

function extend_coloring(w::AvgAccumulator, d::AvgDegStats, out_edge::Bool)
    if out_edge
        return scale_coloring(w, d.avg_out)
    else
        return scale_coloring(w, d.avg_in)
    end
end

############################ MaxDegStats ################################################

struct MaxDegStats <:DegreeStats
    max_in::Float32
    max_out::Float32
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

struct MaxAccumulator <:StatAccumulator
    weight::Float32
end
get_count(w::MaxAccumulator) = w.weight
sum_colorings(w1::MaxAccumulator, w2::MaxAccumulator) = MaxAccumulator(w1.weight + w2.weight)

# Because the max degree estimator aims to produce an upper bound, we have to treat cycle
# closure as 1.0 probability.
function scale_coloring(w::MaxAccumulator, s)
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



############################ CorrDegStats ################################################

struct CorrDegStats <:DegreeStats
    avg_in::Float32
    var_in::Float32
    corr_in::Float32
    avg_out::Float32
    var_out::Float32
    corr_out::Float32
end

function CorrDegStats(g::DataGraph, edges::Vector{Tuple{NodeId, NodeId, Bool}}, color_size::Int)
    length(edges) == 0 && return CorrDegStats(0,0,0,0,0)

    l_1st = [0.0, 0.0]
    l_2nd = [0.0, 0.0]
    r_1st = [0.0, 0.0]
    r_2nd = [0.0, 0.0]
    e_1st = [0.0, 0.0]
    e_count = [0.0, 0.0]
    for (x, y, out_edge) in edges
        if out_edge
            l_1st[1] += degree(g.graph, x)
            l_2nd[1] += degree(g.graph, x)^2
            r_1st[1] += degree(g.graph, y)
            r_2nd[1] += degree(g.graph, y)^2
            e_1st[1] += degree(g.graph, x) * degree(g.graph, y)
            e_count[1] += 1
        else
            l_1st[2] += degree(g.graph, x)
            l_2nd[2] += degree(g.graph, x)^2
            r_1st[2] += degree(g.graph, y)
            r_2nd[2] += degree(g.graph, y)^2
            e_1st[2] += degree(g.graph, x) * degree(g.graph, y)
            e_count[2] += 1
        end
    end
    corrs = [0.0, 0.0]
    for i in 1:2
        l_1st[i] /= e_count[i]
        l_2nd[i] /= e_count[i]
        r_1st[i] /= e_count[i]
        r_2nd[i] /= e_count[i]
        e_1st[i] /= e_count[i]
        r_denom = sqrt(max(0.0, r_2nd[i]-r_1st[i]^2))
        l_denom = sqrt(max(0.0, l_2nd[i]-l_1st[i]^2))
        if r_denom > 0 && l_denom > 0
            corrs[i] = (e_1st[i] - l_1st[i]*r_1st[i])/(l_denom * r_denom)
        else
            corrs[i] = 0
        end
    end

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
    var_in = sum([x^2 for x in values(in_counter)]) / color_size - avg_in^2
    avg_out = sum([x for x in values(out_counter)]; init=0)/color_size
    var_out = sum([x^2 for x in values(out_counter)]) / color_size - avg_out^2
    return CorrDegStats(avg_in, var_in, corrs[1], avg_out, var_out, corrs[2])
end
get_in_deg_estimate(d::CorrDegStats) = d.avg_in
get_out_deg_estimate(d::CorrDegStats) = d.avg_out

struct CorrAccumulator <:StatAccumulator
    weight::Float64
    var::Float64
    CorrAccumulator(w) = new(w, 0)
    CorrAccumulator(w, v) = new(w, v)
end

get_count(w::CorrAccumulator) = w.weight
sum_colorings(w1::CorrAccumulator, w2::CorrAccumulator) = CorrAccumulator(w1.weight + w2.weight, w1.var + w2.var)
scale_coloring(w::CorrAccumulator, s) = CorrAccumulator(w.weight * s, w.var * s^2)

@inline function extend_coloring(w::CorrAccumulator, d::CorrDegStats, out_edge::Bool)
    w_std = sqrt(w.var)
    w_2nd = w.var + w.weight^2
    if out_edge
        d_std = sqrt(d.var_out)
        d_2nd = d.var_out + d.avg_out^2
        new_w = w.weight * d.avg_out + w_std * d_std * d.corr_out
        new_2nd = w_2nd * d_2nd + 4 * w.weight * d.avg_out * w_std * d_std *d.corr_out + 2 * d.corr_out^2 * w.var * d.var_out
        new_var = max(0.0, new_2nd - new_w^2) # Needed to handle FP errors & epsilon negatives
        return CorrAccumulator(new_w, new_var)
    else
        d_std = sqrt(d.var_in)
        d_2nd = d.var_in + d.avg_in^2
        new_w = w.weight * d.avg_in + w_std * d_std * d.corr_in
        new_2nd = w_2nd * d_2nd + 4 * w.weight * d.avg_in * w_std * d_std * d.corr_in + 2 * d.corr_in^2 * w.var * d.var_in
        new_var = max(0.0, new_2nd - new_w^2) # Needed to handle FP errors & epsilon negatives
        return CorrAccumulator(new_w, new_var)
    end
end
