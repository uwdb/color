function load_querysets(datasets::Vector{DATASET}=Vector{DATASET}();require_true_cardinality=true)
    if length(datasets) == 0
        datasets = instances(DATASET)
    end

    query_directories = Dict()
    query_directories[aids] = ["/queryset/aids/Chain_3/",
        "/queryset/aids/Chain_6/",
        "/queryset/aids/Chain_9/",
        "/queryset/aids/Chain_12/",
        "/queryset/aids/Cycle_3/",
        "/queryset/aids/Cycle_6/",
        "/queryset/aids/Flower_6/",
        "/queryset/aids/Flower_9/",
        "/queryset/aids/Flower_12/",
        "/queryset/aids/Graph_3/",
        "/queryset/aids/Graph_6/",
        "/queryset/aids/Graph_9/",
        "/queryset/aids/Graph_12/",
        "/queryset/aids/Petal_6/",
        "/queryset/aids/Petal_9/",
        "/queryset/aids/Petal_12/",
        "/queryset/aids/Star_3/",
        "/queryset/aids/Star_6/",
        "/queryset/aids/Star_9/",
        "/queryset/aids/Tree_3/",
        "/queryset/aids/Tree_6/",
        "/queryset/aids/Tree_9/",
        "/queryset/aids/Tree_12/"]
    query_directories[human] = ["/queryset/human/Chain_3/",
    "/queryset/human/Graph_3/",
    "/queryset/human/Star_3/",
    "/queryset/human/Tree_3/"]
    query_directories[lubm80] = ["/queryset/lubm80"]
    query_directories[yago] = ["/queryset/yago/Chain_3",
    "/queryset/yago/Chain_6",
    "/queryset/yago/Chain_9",
    "/queryset/yago/Chain_12",
    "/queryset/yago/Clique_6",
    "/queryset/yago/Clique_10",
    "/queryset/yago/Cycle_3",
    "/queryset/yago/Cycle_6",
    "/queryset/yago/Cycle_9",
    "/queryset/yago/Flower_6",
    "/queryset/yago/Graph_3",
    "/queryset/yago/Graph_6",
    "/queryset/yago/Graph_9",
    "/queryset/yago/Graph_12",
    "/queryset/yago/Petal_6",
    "/queryset/yago/Petal_9",
    "/queryset/yago/Petal_12",
    "/queryset/yago/Star_3",
    "/queryset/yago/Star_6",
    "/queryset/yago/Star_9",
    "/queryset/yago/Star_12",
    "/queryset/yago/Tree_3",
    "/queryset/yago/Tree_6",
    "/queryset/yago/Tree_9",
    "/queryset/yago/Tree_12",
    ]

    query_directories[yeast] = ["/queryset/yeast"]
    query_directories[hprd] = ["/queryset/hprd"]
    query_directories[wordnet] = ["/queryset/wordnet"]
    query_directories[youtube] = ["/queryset/youtube"]
    query_directories[patents] = ["/queryset/patents"]
    query_directories[eu2005] = ["/queryset/eu2005"]
    query_directories[dblp] = ["/queryset/dblp"]

    query_paths = Dict()
    for dataset in instances(DATASET)
        query_paths[dataset] = [readdir(pwd() * dir, join=true) for dir in query_directories[dataset]]
        query_paths[dataset] = [(query_paths[dataset]...)...]
    end

    all_queries = Dict()
    for dataset in datasets
        all_queries[dataset] = []
        println("Loading Queries For: ", dataset)
        for query_path in query_paths[dataset]
            exact_size = -1
            if dataset == lubm80
                query_type = match(r".*/queryset/.*/lubm80_(.*).txt", query_path).captures[1]
                exact_size = load_true_cardinality(replace(query_path, "queryset"=>"TrueCardinalities"))
                query = load_query(query_path)
            elseif IS_GCARE_DATASET[dataset]
                query_type = match(r".*/queryset/.*/(.*)_.*/.*", query_path).captures[1]
                exact_size = load_true_cardinality(replace(query_path, "queryset"=>"TrueCardinalities"))
                query = load_query(query_path)
            else
                query_type = match(r".*/queryset/.*/query_(.*)_.*", query_path).captures[1]
                if isfile(replace(query_path, "queryset"=>"TrueCardinalities"))
                    exact_size = load_true_cardinality(replace(query_path, "queryset"=>"TrueCardinalities"))
                end
                query = load_query(query_path, subgraph_matching_data=true)
            end
            exact_size == -1 && require_true_cardinality && continue
            push!(all_queries[dataset], (query=query, exact_size=exact_size,
                                            query_type=query_type, query_path=query_path))
        end
    end
    return all_queries
end
