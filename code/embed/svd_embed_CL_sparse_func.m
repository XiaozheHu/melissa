function [x, d] = svd_embed_CL_sparse_func(walks, N, ndim)
% when we have cannot link (CL), the entry of the random walks might be
% negative. In this case, smoothing step should be done carefully. First,
% shift the matrix so that the minimal entry is 0.  Then do the smoothing
%
% For sparse format
% Created by Xiaozhe

  % get number of networks
  nnetworks  = size(walks, 1);

  % initilize 
  %RR_sum = sparse(ngene, ngene);

  for i = 1:nnetworks
      % get i-th random walk
      W = walks{i};
    
      % find the minimal entry, if it is non-negative, set it to zero. 
      min_entry = min(min(W)); 
      if min_entry > 0
        min_entry = 0.0;
      end
      
      % smooth
      %R = log(W - min_entry + 1/ngene); % smoothing
      [row_idx, col_idx, W_val] = find(W);
      W_val = (W_val - min_entry);
      idx = (W_val > 0); % random walk matrix should not have negetive entries
      %R = sparse(row_idx(idx), col_idx(idx), log(W_val(idx)), N, N);
      walks{i} = sparse(row_idx(idx), col_idx(idx), log(W_val(idx)), N, N);

      %RR_sum = RR_sum + R * R';
      %RR_sum = RR_sum + R' * R;
      
      clear W row_idx col_idx W_val idx;
  end
  
  %clear R Q A;
  %[V, d] = eigs(RR_sum, ndim);
  [V, d] = eigs(@(xx)RRsum_x(xx, walks, N), N, ndim);
  x = spdiags(sqrt(sqrt(diag(d))),0,ndim,ndim) * V';
  
end

function y = RRsum_x(xx, Rall, N)
    
    % get number of networks
    nnetworks  = size(Rall, 1);

    % initilize 
    y = zeros(N, 1);

    for i = 1:nnetworks

        R = Rall{i};

        y = y + (R'*(R*xx));

    end

end
