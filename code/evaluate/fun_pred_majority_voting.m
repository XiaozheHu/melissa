function [acc, f1, auc, class_score] = fun_pred_majority_voting(anno, test_filt,train_filt, knn, dist_mat, weighted)
%----------------------------------------------------
% function prediction with majority voting by kNN
% adapted from matrix_majority_voting.m
% if any node has no voting neighbors among the first k NN, ask its nearest
% labelled neighbor to vote.
%----------------------------------------------------

    addpath code/evaluate;
    
    test_ind = find(test_filt);
    test_anno = anno(:, test_filt);

    labelled = sum(anno, 1) > 0;
    
    % only look at training indices with labels
    % this may be redundant since this filter happens earlier as well
    train_ind = intersect(find(train_filt), find(labelled));
    
    ntest = length(test_ind);
    
    % creates a matrix for scoring the test points
    class_score = zeros(size(anno, 1), ntest);
    
    for p = 1:ntest
        i = test_ind(p);    % gene id
        
        % indices of nearest neighbors
        iknn = knn(:,i);
        % only care about neighbors with labels
        voting_knn = intersect(iknn, train_ind);
        
        if isempty(voting_knn)
            % find the nearest labelled node
            [dist,index]=sort(dist_mat(i,:),'descend');
        
            iknn=index(find(dist,size(dist_mat, 1),'last'));

            voting_knn = intersect(iknn, train_ind, 'stable');
            %voting_knn = voting_knn((end-k2+1):end);
            voting_knn = voting_knn(end);
        end
        
        votes = anno(:, voting_knn);
        if weighted 
            % weight is inversely proportional to distance
            weights = 1 ./ dist_mat(i, voting_knn);
            class_score(:, p) = votes * weights.';
        else
            class_score(:, p) = sum(votes,2);
        end
    end
    
%     % filter out any nodes that had no voting neighbors.
%     had_votes = sum(class_score, 1) > 0;
%     had_votes_anno = test_anno(:, had_votes);
%     
%     % remove nodes with no voters
%     class_score = class_score(:, had_votes);
    
    % make all columns sum to one
    class_score = class_score ./ sum(class_score);
    
    %fprintf('%f out of %f had voting neighbors\n', sum(had_votes), length(had_votes));
    
    [acc, f1, auc] = evaluate_performance(class_score.', test_anno.');
end
