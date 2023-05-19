using Graphs: SimpleDiGraphFromIterator, Edge

function load_dataset(path; subgraph_matching_data=false)
    n = 0
    edges::Array{Tuple{Tuple{Int, Int}, Int}} = []
    vertices::Array{Tuple{Int, Array{Int}}} = []
    for line in eachline(path)
        if length(line) == 0
            continue
        end
        if line[1] == 'v'
            parts = split(line)
            labels = []
            if (subgraph_matching_data)
                push!(labels, parse(Int, parts[3]))
                push!(vertices, (parse(Int, parts[2])+1, labels))
            else
                for x in parts[3:length(parts)]
                    push!(labels, parse(Int, x))
                end
                push!(vertices, (parse(Int,  parts[2]) + 1, labels))
            end
            n += 1
        elseif line[1] == 'e'
            parts = split(line)
            if (subgraph_matching_data)
                e1, e2 = parse(Int, parts[2])+1, parse(Int, parts[3])+1
                push!(edges, ((e1, e2), 0))
            else
                e1, e2, l1 = parse(Int, parts[2])+1, parse(Int, parts[3])+1, parse(Int, parts[4])
                push!(edges, ((e1, e2), l1))
            end
        end
    end
    g = DataGraph(n)
    for vertex_and_labels in vertices
        update_node_labels!(g, vertex_and_labels[1], vertex_and_labels[2])
    end
    for edge_and_label in edges
        add_labeled_edge!(g, edge_and_label[1], edge_and_label[2])
    end
    return g
end


function load_query(path; subgraph_matching_data=false)
    n = 0
    queryID = 0
    edges::Array{Tuple{Tuple{Int, Int}, Int}} = []
    vertices::Array{Tuple{Int, Int, Int}} = []
    for line in eachline(path)
        if length(line) == 0
            continue
        end
        if line[1] == 't'  
            continue
        elseif line[1] == 'v'
            parts = split(line)
            if (subgraph_matching_data)
                data_label = -1
                label = parse(Int, parts[3])
                push!(vertices, (parse(Int, parts[2]) + 1, label, data_label))
            else
                data_label = parse(Int, parts[4])
                label = parse(Int, parts[3])
            push!(vertices, (parse(Int,  parts[2]) + 1, label, data_label))
            end
            n += 1
        elseif line[1] == 'e'
            parts = split(line)
            if (subgraph_matching_data)
                e1, e2 = parse(Int, parts[2])+1, parse(Int, parts[3])+1
                push!(edges, ((e1, e2), 0))
            else
                e1, e2, l1 = parse(Int, parts[2]) + 1, parse(Int, parts[3]) + 1, parse(Int, parts[4])
                push!(edges, ((e1, e2), l1))
            end
        end
    end
    g = QueryGraph(n)
    for vertex_and_labels in vertices
        update_node_labels!(g, vertex_and_labels[1], vertex_and_labels[2])
        update_data_labels!(g, vertex_and_labels[1], vertex_and_labels[3])
    end
    for edge_and_label in edges
        add_labeled_edge!(g, edge_and_label[1], edge_and_label[2])
    end
    return (queryID, g)
end

function load_true_cardinality(path)
    true_cardinality = 0
    for line in eachline(path)
        true_cardinality = parse(Int, split(line)[1])
    end
    return true_cardinality
end
