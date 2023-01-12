using Graphs: SimpleGraphFromIterator, Edge

function load_dataset(path)
    n = 0
    edges = Vector{Edge{Int}}()
    for line in eachline(path)
        if length(line) == 0
            continue
        end
        if line[1] == 'v'
            n += 1
        elseif line[1] == 'e'
            parts = split(line)
            e1, e2 = parse(Int, parts[2]), parse(Int, parts[3])
            push!(edges, Edge(e1, e2))
        end
    end
    return SimpleGraphFromIterator(edges)
end
