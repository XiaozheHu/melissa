function [best_score] = resnik(anno, full_anno,ontgraph, used_ind, class_score, test_filt)
    nterms = size(ontgraph,1);
    % uses the calculation above to find all the parents
    parents = get_reachable(ontgraph);
    
    term_counts = sum(full_anno, 2);
    % fills in the Pr(v | P(v))
    pr_given_parents = zeros(nterms,1);
    for i = 1:nterms
        
        i_parents = ontgraph(:,i);
        if sum(i_parents) == 0 % handles the root
            pr_given_parents(i) = 1;
            continue;
        end
        
        % the divisor just counts the number of nodes with a label from
        % any one of the parents, and doesn't double count.
        
        anno_by_parents = sum( max( full_anno( i_parents, :), [], 1) );
        
        % this happens when the child and parent are both not used
        % and gives a Nan
        if anno_by_parents == 0 
            pr_given_parents(i) = 1;
            continue;
        end
        
        pr_given_parents(i) = term_counts(i) / anno_by_parents;
    end
    
    % we never care about zero use terms but we don't want log 0 later
    % this flips them to 1s
    pr_given_parents = pr_given_parents + (pr_given_parents == 0); 
    
    % computes i(l) for each term
    info_content = zeros(nterms,1); 
    for i = 1:nterms
        i_parents = ontgraph(:,i);
        
        if sum(i_parents) == 0 % handles the root
            info_content(i) = 0;
            continue;
        end
        
        info_content(i) = -sum(log(pr_given_parents(i_parents)));
    end
    
    
    nused_terms = max(size(used_ind)); % robust to row vecs & col vecs
    % for each of the terms we're using finds the least common ancestor
    lcas = zeros(nused_terms);
    % parents of the used terms
    rel_parents = parents(:,used_ind); % rel for relevant
    
    % BFS down the tracking least common ancestor
    [~,root_index] = min(sum(ontgraph,1));
    queue = [root_index];
    visited = zeros(nterms,1);
    visited(root_index) = 1;
    
    while isempty(queue) == false
        current = queue(1);
        queue = queue(2:end); % pops the first thing off the queue
        
        children = find(ontgraph(current,:));
        desc = rel_parents(current, :);
        
        % this works because we're doing BFS
        lcas((desc.' * desc) > 0) = current;
        
        for child = children
            if visited(child) == 1
                continue
            end
            visited(child) = 1;
            queue = [queue, child];
        end
    end
    
    % use the lcas and info content to get pairwise resnik
    resnik_tt = zeros(nused_terms);
    for i = 1:nused_terms
        for j = 1:nused_terms
            resnik_tt(i,j) = info_content(lcas(i,j));
        end
    end
    
    % actually computes the resnik score
    test_genes = find(test_filt);
    ntest_genes = max(size(test_genes));
    best_score = -inf;
    for t = 1:20 % sweep from 0 to 1 for threshold. (Exact is too slow)
        thresh = t / 20;
        
        score = 0;
        for i = 1:ntest_genes
            gene = test_genes(i);
            true_labels = anno(:,gene) > 0; 
            pred_labels = class_score(:,i) > thresh;
            
            % used in the denominator
            ntrue = sum(true_labels); 
            npred = sum(pred_labels);
            
            if npred == 0
                % this happens when we fail to predict,
                % and avoids weird empty arrays later.
                continue;
            end
            
            % fast way of getting the two sums
            res_submat = resnik_tt(true_labels, pred_labels);
            first_sum = sum(max(res_submat, [], 1));
            second_sum = sum(max(res_submat, [], 2));
            
            score = score + (first_sum + second_sum) / (ntrue + npred);
        end
        best_score = max(best_score, score);
    end
end

