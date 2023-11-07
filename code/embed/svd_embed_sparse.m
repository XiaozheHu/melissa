function [x, d] = svd_embed_sparse(walks, ngene, ndim)
% For sparse format
% Created by Xiaozhe

  % get number of networks
  nnetworks  = size(walks, 1);
  
  % initilize 
  RR_sum = sparse(ngene, ngene);

  for i = 1:nnetworks
    % get i-th random walk
    W = walks{i};

    %R = log(W + 1/ngene); % smoothing
    [row_idx, col_idx, W_val] = find(W);
    idx = (W_val > 0); % random walk matrix should not have negetive entries
    R = sparse(row_idx(idx), col_idx(idx), log(W_val(idx)), ngene, ngene);
  
    %RR_sum = RR_sum + R * R';
    RR_sum = RR_sum + R' * R;

    clear W R row_idx col_idx W_val idx;
  end

  %clear R Q A;
  [V, d] = eigs(RR_sum, ndim);
  x = spdiags(sqrt(sqrt(diag(d))),0,ndim,ndim) * V';

end