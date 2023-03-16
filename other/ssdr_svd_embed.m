function [x, d] = ssdr_svd_embed(walks, ndim, L, flag_JL)

  [nnetworks, ngene, ~]  = size(walks);
  
  %RR_sum = zeros(ngene);
  if (flag_JL == 1)
    epsilon = 0.5;  % accuracy 
    k = ceil((6/(epsilon^2/2 - epsilon^3/3))*log10(ngene));  % low dimension
    k = max(k, ceil(ndim/nnetworks));
    %k = ndim;
    R_all = zeros(nnetworks*k, ngene);
  else
    R_all = zeros(nnetworks*ngene, ngene);
  end
  
  for i = 1:nnetworks

    W = squeeze(walks(i,:,:));
    
    if (flag_JL == 1)
        Q = generate_rand_proj( k, ngene, 'dense');
        %R_all = [R_all; Q*(log(W + 1/ngene))]; % smoothing + random projection
        R_all((i-1)*k+1:i*k, :) = Q*(log(W + 1/ngene));
    else
        %R_all = [R_all; (log(W + 1/ngene))];
        R_all((i-1)*ngene+1:i*ngene, :) = (log(W + 1/ngene));
    end

    clear W;
  end
  
  %Q = generate_rand_proj( ndim, ngene*nnetworks, 'dense');
  %R_all = Q*R_all;

  temp_mat = R_all*L*R_all';
  [V, d] = eigs(temp_mat, ndim);
  x = V'*R_all;

end
