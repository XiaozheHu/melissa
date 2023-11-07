function P = markov_mat_sparse(A)

  % compute the transition matrix P = AD^{-1}
  d = sum(A);
  idx = (d==0);  
  d(idx) = 1;
  
  P = bsxfun(@rdivide, A, d);
end
