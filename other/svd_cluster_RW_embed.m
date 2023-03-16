function x = svd_cluster_RW_embed(walks, gene_clusters, CL_step, options)


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
      
    nnodes = ngene + num_clusters*CL_step;
    augmented = zeros(nnodes);
    augmented(1:ngene, 1:ngene) = RR_sum;

    % treat the dummy nodes as their own things
    % the augmenting nodes should be comparable in length to the data
    max_length = max(vecnorm(RR_sum));
    augmented(ngene+1:nnodes, ngene+1:nnodes) = eye(num_clusters*CL_step) * max_length;
    % create a constraint Laplacian
    A = zeros(nnodes);

    % optional ?  Add the unsupervised links as well
    if use_unsupervised
        A(1:ngene,1:ngene) = 1 / (ngene * ngene);
    end
    
    % add the ML and CL constraints
    
    
    %A(ngene+1:nnodes,ngene+1:nnodes) = cl_penalty - (cl_penalty * eye(num_clusters));
    
    A((ngene+CL_step):CL_step:nnodes,(ngene+CL_step):CL_step:nnodes) = cl_penalty - (cl_penalty * eye(num_clusters));
    A((ngene+1):CL_step:nnodes,1:ngene) = -ml_penalty * gene_clusters;
    A(1:ngene,(ngene+1):CL_step:nnodes) = -ml_penalty * gene_clusters';
    
    toeplitz_vec = zeros(1,CL_step); % creates a matrix of size CL_step x CL_step with pattern
    toeplitz_vec(2) = -ml_penalty;   % 0   -ml   0    0   0
    block = toeplitz(toeplitz_vec);  % -ml  0    -ml  0   0
                                     % 0    -ml  0    -ml 0
                                     % 0    0    -ml  0   -ml
                                     % 0    0    0    -ml 0
    for i=1:num_clusters
        start_idx = ngene+1+(i-1)*CL_step;
        end_idx = start_idx+CL_step-1;
       A(start_idx:end_idx,start_idx:end_idx) = block; 
    end
    
                                     
    D = diag(sum(A));
    L = D - A;

    [V, d] = eigs(augmented + L, ndim);
    y = diag(sqrt(sqrt(diag(d)))) * V';
    x = y(:, 1:ngene);
end

