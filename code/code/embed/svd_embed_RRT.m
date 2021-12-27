function [x, d] = svd_embed_RRT(walks, ndim)
  [nnetworks, ngene, ~]  = size(walks);
  RR_sum = zeros(ngene);
  for i = 1:nnetworks
    W = squeeze(walks(i,:,:));
    R = log(W + 1/ngene); % smoothing
    %R = full(W);
    RR_sum = RR_sum + R * R';
    clear W R;
  end
  %clear R Q A;
  [V, d] = eigs(RR_sum, ndim);
  x = diag(sqrt(sqrt(diag(d)))) * V';
end
