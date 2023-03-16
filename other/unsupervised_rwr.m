function walks = unsupervised_rwr(network_files, ngene, restart_prob)

  nfiles = length(network_files);
  walks = zeros(nfiles, ngene, ngene);
  
  for i=1:nfiles
    A = load_network(network_files{i}, ngene);
    P = markov_mat(A);
    W = rwr(P, restart_prob);
    walks(i,:,:) = W;
  end
end
