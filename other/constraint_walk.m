function c_walk = constraint_walk(A,restart_prob)
    if ~isequal(A, A') % symmetrize
      A = A + A';
    end
    A = A + diag(sum(A, 2) == 0);
    P = markov_mat(A);
    c_walk = rwr(P, restart_prob); 
end
