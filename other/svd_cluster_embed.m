function x = svd_cluster_embed(walks, gene_clusters, options)

    % we may want to adjust what we do with the augmenting nodes.

    embedding = options.embedding;
    ndim = embedding.ndim;
    cl_penalty = embedding.cannotlink_penalty;
    ml_penalty = embedding.mustlink_penalty;
    use_unsupervised = options.embedding.use_unsupervised;
    num_clusters = options.num_clusters;
    

    [nnetworks, ngene, ~] = size(walks);
    
    RR_sum = zeros(ngene);
    % add the networks together without dummy nodes
    for i = 1:nnetworks
        W = squeeze(walks(i,:,:));
        R = log(W + 1/ngene); % smoothing
        RR_sum = RR_sum + R * R';
    end
    clear R
    
    % normalize to have mean 0
    % RR_Sum = RR_sum - mean(RR_sum(:));
    % normalize to have std 1
    % RR_sum = RR_sum / std(RR_sum(:));
    
    
    nnodes = ngene + num_clusters;
    augmented = zeros(nnodes);
    augmented(1:ngene, 1:ngene) = RR_sum;
R
    % treat the dummy nodes as their own things
    % the augmenting nodes should be comparable in length to the data
    max_length = max(vecnorm(RR_sum));
    augmented(ngene+1:nnodes, ngene+1:nnodes) = eye(num_clusters) * max_length;
    % create a constraint Laplacian
    A = zeros(nnodes);

    % optional ?  Add the unsupervised links as well
    if use_unsupervised
        A(1:ngene,1:ngene) = 1 / (ngene * ngene);
    end
    
    % add the ML and CL constraints
    A(ngene+1:nnodes,ngene+1:nnodes) = cl_penalty - (cl_penalty * eye(num_clusters));
    A(ngene+1:nnodes,1:ngene) = -ml_penalty * gene_clusters;
    A(1:ngene,ngene+1:nnodes) = -ml_penalty * gene_clusters';
    
    D = diag(sum(A));
    L = D - A;

    [V, d] = eigs(augmented + L, ndim);
    y = diag(sqrt(sqrt(diag(d)))) * V';
    x = y(:, 1:ngene);
end

