function x = compute_embedding(walks, gene_clusters, options)
    addpath code/embed
    
    embedding = options.embedding;
    svd_approx = embedding.svd_approx;
    use_clusters = ~(options.num_clusters == -1);
    
    
    if svd_approx
        if ~use_clusters
            fprintf('[Base svd embedding]');
            x = svd_embed(walks, ndim);
        else
            fprintf('[Cluster svd embedding]');
            x = svd_cluster_embed(walks, gene_clusters, options);
        end
    else
        frprintf('[Base KL embedding]')
        x = unsupervised_embed(walks, ndim, 1000); % 3rd arg is max iterations
    end
end

