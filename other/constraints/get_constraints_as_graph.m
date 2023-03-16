function mladj = get_constraints_as_graph(walks,labels, ml_penalty)

    [~, ngene, ~] = size(walks);
    
    [l1, l2] = size(labels);
    rf = rand(l1,l2) > 0.95; % use 1-p fraction of the labels
    filtered_labels = (rf .* labels) > 0;

    % create the must-link and cannot-link matrices (0s and 1s)
    ml = mustlink(filtered_labels);
 
    ml = ml_penalty * ml;
     
    mladj = [];
    
    for i = 1:size(ml,1)
        for j = (i+1):size(ml,2)
            if ml(i,j)~=0
                mladj = [mladj; [i,j,ml(i,j)]];
            end    
        end
    end
    nz = nnz(ml);
    % Checking how dense supervised matrix is
    fprintf('  Supervised penalty has %d non-zero (%f percent)\n', nz, 100 * nz / (ngene * ngene));
end
