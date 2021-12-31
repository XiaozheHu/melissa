function walks = semi_supervised_rwr_argumented(network_files, ngene, gene_clusters, restart_prob, mustlink_penalty, cannotlink_penalty)
  
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
        
        % get size
        N = size(A,1);
        
        % make sure the diagonal entries of A are zeros
        A = A - spdiags(diag(A),0,N,N);
        
        % compute weighted degree matrix
        d = sum(abs(A),2);
        idx = (d==0);
        d(idx) = 1;
        D = spdiags(d, 0, N, N);
        
        % graph Laplacian (with restart)
        L = D - (1-restart_prob)*A;
        
        % compute diffusion states
        if (ngene >= 10000)
            W = restart_prob*(D/L);
            W = full(W);
        else
            L = full(L);
            W = restart_prob*(D/L);
        end
             
        walks(i,:,:) = W;
        
    end
    
end
