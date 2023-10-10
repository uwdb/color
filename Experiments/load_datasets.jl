function load_dataset(dataset::DATASET)
    if dataset == aids
        aids_data_file_path = "dataset/aids/aids.txt"
        return load_dataset(aids_data_file_path)
    elseif dataset == human
        human_data_file_path = "dataset/human/human.txt"
        return load_dataset(human_data_file_path)
    elseif dataset == lubm80
        lubm80_data_file_path = "dataset/lubm80/lubm80.txt"
        return load_dataset(lubm80_data_file_path)
    elseif dataset == yago
        yago_data_file_path = "dataset/yago/yago.txt"
        return load_dataset(yago_data_file_path)
    elseif dataset == yeast
        yeast_data_file_path = "dataset/yeast/yeast.graph"
        return load_dataset(yeast_data_file_path, subgraph_matching_data=true)
    elseif dataset == hprd
        hprd_data_file_path = "dataset/hprd/hprd.graph"
        return load_dataset(hprd_data_file_path, subgraph_matching_data=true)
    elseif dataset == wordnet
        wordnet_data_file_path = "dataset/wordnet/wordnet.graph"
        return load_dataset(wordnet_data_file_path, subgraph_matching_data=true)
    elseif dataset == dblp
        dblp_data_file_path = "dataset/dblp/dblp.graph"
        return load_dataset(dblp_data_file_path, subgraph_matching_data=true)
    elseif dataset == youtube
        youtube_data_file_path = "dataset/youtube/youtube.graph"
        return load_dataset(youtube_data_file_path, subgraph_matching_data=true)
    elseif dataset == patents
        patents_data_file_path = "dataset/patents/patents.graph"
        return load_dataset(patents_data_file_path, subgraph_matching_data=true)
    elseif dataset == eu2005
        eu2005_data_file_path = "dataset/eu2005/eu2005.graph"
        return load_dataset(eu2005_data_file_path, subgraph_matching_data=true)
    end
end
