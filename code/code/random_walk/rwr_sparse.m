function W = rwr_sparse(P, restart_prob)
  
  % get size
  N = size(P, 1);
  
  % compute random walk matrix (diffusion state)
  Identity = spdiags(ones(N),0,N,N);
  
  W = (Identity - (1 - restart_prob) * P) \ (restart_prob * Identity);
  
end