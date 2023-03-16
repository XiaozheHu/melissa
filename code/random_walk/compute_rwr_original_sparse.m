% compute random walk on original networks
% use sparse matrix
% Modified by Xiaozhe

function walks = compute_rwr_original_sparse(network_files, ngene, options)

    addpath code/random_walk;
    
    restart_prob = options.walk.restart_prob;

    walks = unsupervised_rwr_original_sparse(network_files, ngene, restart_prob);
    
end

