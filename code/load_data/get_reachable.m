function R = get_reachable(A)
  nterm = size(A, 1);
  R = sparse(nterm, nterm);
  visit = speye(nterm);
  %what if endless loop?
  while nnz(visit) > 0
    % if both elements are zero, result will be zero
    R = R | visit;
    %logical array
    visit = double(visit) * A > 0;
  end
end
