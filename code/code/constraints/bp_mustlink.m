function [ml_graph] = bp_mustlink(anno)
    [nterms, nnodes] = size(anno);
    ml_graph = zeros(nnodes);
    
    % nodes with no labels
    unlabelled = sum(anno,1) == 0;
    
    nlabelled = sum(sum(anno,1) > 0);
    nunlabelled = nnodes - nlabelled;
    
    total_labels = sum(sum(anno));
    % iterate over terms adding rank 2 matrices to the graph
    
    ratio = 1 + nunlabelled/nlabelled;
    
    for i = 1:nterms
        i_labels = anno(i,:) > 0;
        
        nuses = sum(i_labels);
        
        %if nuses == 0
        %    continue;
        %end
        
        % both have the label
        term_graph = i_labels.' * i_labels / (nuses * ratio);
        
        % column has label, row does not
        term_graph(unlabelled, i_labels) = 1 / nnodes;
        
        % row has label, column does not
        term_graph(i_labels, unlabelled) = 1 / (total_labels * ratio);
        
        % neither have label
        term_graph(unlabelled, unlabelled) = nuses / (total_labels * nnodes);
       
        ml_graph = ml_graph + term_graph;
    end
    
    scaling = diag(1 ./ (sum(anno,1) + (sum(anno,1) == 0))); % check this works
    
    ml_graph = ml_graph * scaling;
end

