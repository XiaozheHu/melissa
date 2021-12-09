function P = markov_mat(A)
  A = A - diag(diag(A));
  A = A + diag(sum(A) == 0);
  P = bsxfun(@rdivide, A, sum(A));
end
