function P = markov_mat_sparse(A)

  % get size 
  N = size(A,1);

  % make diagonal zero
  A = A - spdiags(diag(A),0,N,N);
  
  %A = A + diag(sum(A) == 0); %(Again, I do not think this is correct -- Xiaozhe)
  
  % compute the transition matrix P = AD^{-1}
  d = sum(A);
  idx = (d==0);
  
  %A = A + spdiags(idx', 0, N, N);
  
  d(idx) = 1;
  
  P = bsxfun(@rdivide, A, d);
end
