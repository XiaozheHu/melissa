function ml = mustlink_weighted(labels, levels)
    [c, n] = size(labels);
    % for each class create a clique adj. mat.
    ml = zeros(n,n);
    weights = zeros(c,1);
    weights(levels <= 5) = 0;
    weights(levels == 6) = 0.25;
    weights(levels == 7) = 0.5;
    weights(levels == 8) = 0.75;
    weights(levels >= 9) = 1;
    for i = 1:c
        % get the indicator vector for that class
        indicator = labels(i,:) == 1;
        % fills in the ml matrix for that class
        ml(indicator, indicator) = weights(i);
    end
    % remove diagonal
    ml = ml - diag(diag(ml));
end