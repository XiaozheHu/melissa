function ml = mustlink(labels)
    [c, n] = size(labels);
    % for each class create a clique adj. mat.
    ml = zeros(n,n);
    for i = 1:c
        % get the indicator vector for that class
        indicator = labels(i,:) == 1;
        % fills in the ml matrix for that class
        ml(indicator, indicator) = 1;
    end
    % remove diagonal
    ml = ml - diag(diag(ml));
end