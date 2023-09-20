include("load_querysets.jl")


function calculate_true_cardinalities(datasets::Vector{DATASET}; timeout=120, query_types_to_exclude=[])
    queries = load_querysets(datasets; require_true_cardinality=false)
    for dataset in datasets
        data = load_dataset(dataset)
        for i in 1:length(queries[dataset])
            query = queries[dataset][i].query
            query_type =  queries[dataset][i].query_type
            query_path = queries[dataset][i].query_path
            true_card_path = replace(query_path, "queryset"=>"TrueCardinalities")

            skip_query = isfile(true_card_path)
            for exclusion in query_types_to_exclude
                if occursin(exclusion, query_type)
                    skip_query = true
                end
            end
            skip_query && continue

            println("Query: ", query_path)
            exact_size = get_exact_size(query, data, timeout=timeout)
            if exact_size < 0
                println("Timed Out!")
                continue
            elseif exact_size == 0
                println("Zero Matches!")
                continue
            end
            true_card_file = open(true_card_path, "w")
            show(true_card_file , exact_size)
            close(true_card_file)
        end
    end
end
