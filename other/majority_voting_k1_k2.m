function [accuracy, predictions] = majority_voting_k1_k2(dist_mat,k1,k2, anno, test_filt,train_filt)
%----------------------------------------------------
% function prediction with majority voting by kNN
%----------------------------------------------------
    %dist_mat = squareform(pdist(X'));

    labelled = sum(anno, 1) > 0;
    test_ind = find(test_filt);
    % only look at training indices with labels
    % this may be redundant since this filter happens earlier as well
    train_ind = intersect(find(train_filt), find(labelled));
    ntest = length(test_ind);
    acc = zeros(ntest,1);
    predictions = zeros(ntest,1);

    parfor p = 1:ntest
        i = test_ind(p);    % gene id
        
        % collect all function labels from kNN neighbors
        fun = [];
        fun_weight = [];
        
        
        %---------------------------
        % to change
        %---------------------------
        
        % have k1 nearest neighbor to vote
        
        [dist,index]=sort(dist_mat(i,:),'descend');
        % use of find makes it so a node doesn't use itself as a neighbor.
        k1nn=index(find(dist,k1,'last'));
        voting_knn = intersect(k1nn, train_ind);
        
        if length(voting_knn)==0
            [dist,index]=sort(dist_mat(i,:),'descend');
            iknn=index(find(dist,size(dist_mat, 1),'last'));
            voting_knn = intersect(iknn, train_ind, 'stable');
            voting_knn = voting_knn((end-k2+1):end);    
        end
   
        for j = 1:length(voting_knn)
            voter = voting_knn(j);
            voter_fun = find(anno(:, voter));
            voter_num_fun = length(voter_fun);
            
            if voter_num_fun == 0 % voter has no label
                continue;
            end
            
            dist = dist_mat(i,voter);
            voter_weight = 1/dist*ones(1,voter_num_fun);
            fun = [fun;voter_fun]; %array of all funcitons from knn, including repeat terms.
            fun_weight = [fun_weight,voter_weight];
        end
        
        % count frequency of all functions voted by kNN neighbors 
        fun_list=unique(fun);  
        fun_freq=zeros(size(fun_list)); % array of weight/frequency of functions from neighbors
        
        for j = 1:length(fun_list)
            fun_freq(j) = sum(fun_weight(find(fun==fun_list(j))));
        end
        
        [~,max_freq_ind] = max(fun_freq);
        pred_fun = fun_list(max_freq_ind);
        predictions(p) = pred_fun; % store the prediction
        true_fun = find(anno(:,i));
        if ~isempty(intersect(pred_fun, true_fun)) 
            acc(p)=1;
        else
            acc(p)=0;
        end 
    end
    %acc = acc(acc ~= -1); % filter out the nodes with no labelled neighbors
    accuracy = mean(acc);
%     mean(acc(indicator==1))
%     mean(acc(indicator==0))
%     sum(indicator==1)
%     sum(indicator==0)
end
