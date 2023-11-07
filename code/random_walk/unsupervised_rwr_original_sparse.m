% compute unserpervised random walk on original networks
% use sparse matrix
% Modified by Xiaozhe

function walks = unsupervised_rwr_original_sparse(network_files, ngene, restart_prob)

  nfiles = length(network_files);
  %walks = zeros(nfiles, ngene, ngene);
  walks = cell(nfiles,1);
  
  for i=1:nfiles
      
    A = load_network_sparse(network_files{i}, ngene);
    
    P = markov_mat_sparse(A);
    
    W = rwr_sparse(P, restart_prob);
    
    walks{i} = W;
    
  end
  
end
