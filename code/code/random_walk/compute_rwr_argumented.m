function walks = compute_rwr_argumented(network_files, ngene, gene_clusters, options)
    
    addpath code/random_walk;
    addpath code/constraints;
    
    restart_prob = options.walk.restart_prob;
    
    % weights for mustlinks
    mustlink_penalty  = options.embedding.mustlink_penalty;
    % weights for cannot links
    cannotlink_penalty = options.embedding.cannotlink_penalty; 

    walks = unsupervised_rwr_argumented(network_files, ngene, gene_clusters, ...
        restart_prob, mustlink_penalty, cannotlink_penalty);
    
end

