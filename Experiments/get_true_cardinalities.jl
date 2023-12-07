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

function verify_true_cardinalities(datasets::Vector{DATASET}; timeout=120)
    queries = load_querysets(datasets; require_true_cardinality=true)
    for dataset in datasets
        if !IS_GCARE_DATASET[dataset]
            error("Can't verify a dataset that's not G-Care!")
        end
        data = load_dataset(dataset)
        for i in 1:length(queries[dataset])
            query = queries[dataset][i].query
            query_path = queries[dataset][i].query_path
            existing_size = queries[dataset][i].exact_size

            println("Query: ", query_path)
            new_size = get_exact_size(query, data, timeout=timeout)
            if new_size < 0
                println("Timed Out!")
            elseif new_size!= existing_size
                println("Exact Size Calc: DISAGREES")
            else
                println("Exact Size Calc: AGREES")
            end
        end
    end
end
