function walks = unsupervised_rwr_argumented(network_files, ngene, gene_clusters, restart_prob, mustlink_penalty, cannotlink_penalty)
  
    % number of files
    nfiles = length(network_files);
    
    % number of clusters
    nclusters = size(gene_clusters,1);
    
    % initilize space
    walks = zeros(nfiles, ngene+nclusters, ngene+nclusters);
    
    % main loop
    for i=1:nfiles
        
        % load the original network
        A = load_network_sparse(network_files{i}, ngene);
        
        % argument the network with coarse nodes
        %A = [A, mustlink_penalty*gene_clusters'; mustlink_penalty*gene_clusters, sparse(nclusters, nclusters)];
        
        Acc = (-cannotlink_penalty)*ones(nclusters,nclusters);
        Acc = Acc - diag(diag(Acc));
        Acc = sparse(Acc);
        A = [A, mustlink_penalty*gene_clusters'; mustlink_penalty*gene_clusters, Acc];
        
        % form markov transition matrix
        P = markov_mat_sparse(A);
        
        % form random walk matrix (diffusion states)
        if (ngene >= 10000)
            W = rwr_sparse(P, restart_prob);
            W = full(W);
        else 
            P = full(P);
            W = rwr(P, restart_prob);
        end
        
        walks(i,:,:) = W;
        
    end
    
end
