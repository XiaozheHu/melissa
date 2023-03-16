function [link_mat] = go_link_matrix(path, ngene, filter, fraction)
    [i1, i2, levels] = textread(path, '%d%d%d');
    
    % write the rules here for how much weight at each level
    weights = zeros(1, length(levels));
    
%     weights(weights <= 5) = 0;
%     weights(weights == 6) = 0.25;
%     weights(wieghts == 7) = 0.5;
%     weights(weights == 8) = 0.75;
%     weights(weights >= 9) = 1;
    
    weights(levels <= 5) = 0;
    weights(levels == 6) = 0.25;
    weights(levels == 7) = 0.5;
    weights(levels == 8) = 0.75;
    weights(levels >= 9) = 1;
    
    inuse = weights > 0;
    
    % save the work and prune the unused weights
    i1 = i1(inuse);
    i2 = i2(inuse);
    weights = weights(inuse);

    if fraction < 1
        random_filter = rand(length(i1), 1) < fraction;
        i1 = i1(random_filter);
        i2 = i2(random_filter);
        weights = weights(random_filter);
    end

    % create the sparse matrix of weights
    %link_mat = sparse(i1,i2,levels,ngene,ngene);
    
    link_mat = sparse(i1,i2,weights,ngene,ngene);
    
    % filter out the links involving data in the test set
    link_mat(~filter, :) = 0;
    link_mat(:, ~filter) = 0;
end

