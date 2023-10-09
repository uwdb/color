include("../Experiments.jl")


function get_width_distribution(dataset::DATASET; max_exact_width = 10)
    queries = load_querysets([dataset])[dataset]
    treewidths = counter(Int)
    for query in queries
        treewidth = get_min_width_node_order(query.query.graph; return_width=true, max_exact_width=max_exact_width)
        inc!(treewidths, treewidth)
    end
    return treewidths
end
