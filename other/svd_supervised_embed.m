function x = svd_supervised_embed(walks, ndim, labels, ... 
    cl_penalty, ml_penalty)

    [nnetworks, ngene, ~] = size(walks);
    
    % create the must-link and cannot-link matrices (0s and 1s)
    ml = mustlink(labels);
    cl = cannotlink(ml, labels);
    % reweight them using 
    n_ml = sum(sum(ml))/2;
    n_cl = sum(sum(cl))/2;
    % normally we don't need this check but for debugging I have it
    if n_ml > 0
        ml = ml_penalty * ml;
    else
        ml = zeros(ngene);
    end
    if n_cl > 0
        cl = cl_penalty * cl;
    else
        cl = zeros(ngene);
    end
    
    % creates the combined constraint matrix
    S = ml - cl;
    nz = nnz(S);
    % Checking how dense supervised matrix is
    fprintf('  Supervised penalty has %d non-zero (%f percent)\n', nz, 100 * nz / (ngene * ngene));
    LS = diag(sum(S))-S;
    
    RR_sum = zeros(ngene);
    for i = 1:nnetworks
        W = squeeze(walks(i,:,:));
        R = log(W + 1/ngene); % smoothing
        RR_sum = RR_sum + R * R';
    end
    clear R
    [V, d] = eigs(RR_sum + LS, ndim);
    x = diag(sqrt(sqrt(diag(d)))) * V';
end
