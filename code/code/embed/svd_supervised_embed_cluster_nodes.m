function x = svd_supervised_embed_cluster_nodes(walks, ndim, num_clusters, labels, ... 
    ml_penalty)

    [nnetworks, ngene, ~] = size(walks);
    
    % create the must-link and cannot-link matrices (0s and 1s)
    % this is used as "cannot link" between dummy nodes
    ml = mustlink(labels(:,(end-num_clusters+1):end));
     % reweight them using 
    n_ml = sum(sum(ml))/2;
     % normally we don't need this check but for debugging I have it
    if n_ml > 0
        ml = ml_penalty * ml;
    else
        ml = zeros(ngene);
    end
     
    
    % creates the combined constraint matrix
    S = -ml;
    nz = nnz(S);
    % Checking how dense supervised matrix is
    fprintf('  Supervised penalty has %d non-zero (%f percent)\n', nz, 100 * nz / (ngene * ngene));
    LS = diag(sum(S))-S;
    
    RR_sum = zeros(ngene);
    for i = 1:nnetworks
        W = squeeze(walks(i,:,:));
        R = log(W + 1/ngene); % smoothing
        RR_sum = RR_sum + R * R';
    end
    clear R
    [V, d] = eigs(RR_sum + LS, ndim);
    x = diag(sqrt(sqrt(diag(d)))) * V';
end
