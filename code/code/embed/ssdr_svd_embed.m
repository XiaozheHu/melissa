function [x, d] = ssdr_svd_embed(walks, ndim, L)

  [nnetworks, ngene, ~]  = size(walks);
  
  %RR_sum = zeros(ngene);
  epsilon = 0.2;  % accuracy 
  k = ceil((6/(epsilon^2/2 - epsilon^3/3))*log10(ngene));  % low dimension
  k = max(k, ceil(ndim/nnetworks));
  
  R_all = [];
 
  for i = 1:nnetworks
    W = squeeze(walks(i,:,:));
    
    %Q = generate_rand_proj( k, ngene, 'dense');
    
    %R_all = [R_all; Q*(log(W + 1/ngene))]; % smoothing + random projection
    R_all = [R_all; (log(W + 1/ngene))];
    clear W;
    
    %R = log(W + 1/ngene); % smoothing
    %R = full(W);
    %RR_sum = RR_sum + R * R';
    %RR_sum = RR_sum + R' * L * R;
    %clear W R;
  end
  %clear R Q A;
  [V, d] = eigs(R_all*L*R_all', ndim);
  %x = diag(sqrt(sqrt(diag(d)))) * V';
  x = V'*R_all;
end
