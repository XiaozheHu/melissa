function cl = cannotlink(ml, labels)
    [n,~] = size(ml);
    cl = ones(n,n);
    % indicator vector for if a vertex has a label
    unlabelled = sum(labels) == 0;
    % remove the ml constraints
    cl = cl - ml;
    % remove the unlabeled constraints
    cl(unlabelled, :) = 0;
    cl(:, unlabelled) = 0;
    % remove diagonal
    cl = cl - diag(diag(cl));
end
