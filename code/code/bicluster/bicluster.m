function [gene_clusters, label_clusters] = bicluster(anno, train_filt, options)

    num_clusters = options.num_clusters;
    if num_clusters == -1
        gene_clusters = -1;
        label_clusters = -1;
        return
    end

    %clusters_data = spectralCoClustering(anno(:,train_filt), num_clusters, 0);
    clusters_data = ITL(anno(:,train_filt), num_clusters, 0);
    
    gc = clusters_data.NumxCol;
    
    % gene_clusters is a (num_clusters)x(num_labels) binary matrix, 
    % each row is an indicator vector for a cluster (I don't think we use
    % this at all)
    label_clusters = clusters_data.RowxNum;
    
    % converts the gene clusters back into binary
    gene_clusters = sparse(num_clusters, length(train_filt));
    
    % gene_clusters is a (num_clusters)x(num_gene) binary matrix, 
    % each row is an indicator vector for a cluster
    gene_clusters(:,train_filt) = gc;
end

