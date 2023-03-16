function [x, d] = svd_embed_CL(walks, ndim)
% when we have cannot link (CL), the entry of the random walks might be
% negative. In this case, smoothing step should be done carefully. First,
% shift the matrix so that the minimal entry is 0.  Then do the smoothing
%

  [nnetworks, ngene, ~]  = size(walks);
  RR_sum = zeros(ngene);
  for i = 1:nnetworks
      
      W = squeeze(walks(i,:,:));
    
      % find the minimal entry, if it is non-negative, set it to zero. 
      min_entry = min(min(W)); 
      if min_entry > 0
        min_entry = 0.0;
      end
      
      % smooth
      R = log(W - min_entry + 1/ngene); % smoothing
    
      %RR_sum = RR_sum + R * R';
      RR_sum = RR_sum + R' * R;
      
      clear W R;
  end
  
  %clear R Q A;
  [V, d] = eigs(RR_sum, ndim);
  x = diag(sqrt(sqrt(diag(d)))) * V';
  
end
