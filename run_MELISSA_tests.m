clear;clc;
addpath(genpath('code'));
addpath(genpath('data'));

%------------------------
% Example parameters
%------------------------

% use human or yeast data
options.org = 'yeast';

% which type of annotations to use
% options: {bp, mf, cc} for human and yeast GO labels
options.onttype = 'mf'; 

% filter terms into a specific size range 
% examples: [11 30], [31 100], [101 300]
options.ontsize = [31 100];
     
% number of kNN voters
k=10;

% Number of bi-clusters to create (-1 as not bi-cluster)
options.num_clusters = 4; 

% use SVD approximation for Mashup
% recommended: true for human, false for yeast
options.embedding.svd_approx = true;

% whether to stack the matrics really tall for svd or use the
% log transofrm and sum - needs implementation of dummy nodes
% options.embedding.svd_full = false;

% number of dimensions
% this is used as the maximal embeded dimension, tests will be run on
% different dimensions 
options.embedding.ndim =1000; 

% the weight of the edges connecting dummy nodes to true nodes
% choose value between 0 and 1
options.embedding.mustlink_penalty = 1; 

% the weight of the edges connecting dummy nodes to dummy nodes
% choose value between 0 and 1
options.embedding.cannotlink_penalty = 0; 

% chance that the random walk restarts itself
options.walk.restart_prob = 0.5;

% number of folds for k-fold cross validation
% if set to 1 or less then a single experiment is run
options.kfolds = 5;

% if options.kfolds is set to 1 or less then this is
% the portion of the labelled vertices to go in the 
% testing set. (1 - test_fraction) is used to train
options.test_fraction = 0.2;
                    
%Logs the options so we can see the parameters used in log file later
log_options(options);

%------------------------ 
% Construct network file paths
%------------------------
string_nets = {'neighborhood', 'fusion', 'cooccurence', 'coexpression', ...
               'experimental', 'database'};
network_files = cell(1, length(string_nets));
for i = 1:length(string_nets)
     network_files{i} = sprintf('data/networks/%s/%s_string_%s_adjacency.txt', ...
         options.org, options.org, string_nets{i}); 
end

%------------------------
% Load gene list
%------------------------
fprintf('[Loading annotations]\n');
[genes, ngene, anno] = load_anno(options);

fprintf('Number of functional labels: %d\n', size(anno, 1));

if options.kfolds <= 1
    % Generate training and testing sets
    fprintf('Acquiring test filter using %d testing fraction\n', options.test_fraction);
    folds = create_kfolds(anno, options);
else
    folds = create_kfolds(anno, options);
end

% compute the network walks for both mu and f.
fprintf('Compute random walk on original networks \n');

if (ngene >= 10000) % yeast use full matrix and human use sparse matrix
   walks = compute_rwr_original_sparse(network_files, ngene, options);
else 
   walks = compute_rwr_original(network_files, ngene, options);
end

% Mashup embedding is shared across all folds so we take it before cross-val
fprintf('Compute Mashup embedding \n');

[x_mu, d_mu] = svd_embed(walks, options.embedding.ndim);

clear walks;

%-----------------------------------------------
% Perform function prediction by majority voting
%-----------------------------------------------
% step size of embeded dimension 
dim_step = 25;
n_dim_test = floor(options.embedding.ndim / dim_step);

% Mashup results
acc_mu = zeros(length(folds), n_dim_test);
f1_mu  = zeros(length(folds), n_dim_test);
auc_mu = zeros(length(folds), n_dim_test);

% MELISSA results
acc_melissa = zeros(length(folds), n_dim_test);
f1_melissa  = zeros(length(folds), n_dim_test);
auc_melissa = zeros(length(folds), n_dim_test);

weighted = true;
fprintf('weighted: true \n');

%------------------------
% main loop
%------------------------
for i = 1:length(folds)
    
    fprintf('Fold %d / %d \n', i, options.kfolds);
    fprintf('Using k = %f \n', k);
    
    % prepare train and test sets
    train_filt = folds(i).train_filt;
    test_filt = folds(i).test_filt;
    training_labels = anno.*(train_filt.');
    
    % performe biclustering
    [gene_clusters, label_clusters] = bicluster(training_labels, train_filt, options);

    % argument the network and run random walk
    fprintf('Compute argumented random walk \n');
    tic;
    walks_argumented = compute_rwr_argumented(network_files, ngene, gene_clusters, options);
    toc;
    
    % compute melissa embedding
    fprintf('Compute Melissa embedding \n');
    tic;
    if (options.embedding.cannotlink_penalty == 0)
        [x_melissa, d_melissa] = svd_embed(walks_argumented, options.embedding.ndim);
    else
        [x_melissa, d_melissa] = svd_embed_CL(walks_argumented, options.embedding.ndim);
    end
    toc;
    
    clear walks_argumented;
    
    %-------------------------------
    % Perform prediction
    %-------------------------------
    for j = 1:n_dim_test
        
        embed_dim = j*dim_step;
        fprintf('Using embeded dimension = %f \n', embed_dim);
         
        %-------------------------------
        % Perform Mashup for comparison
        %-------------------------------
        fprintf('[Performing Mashup]\n');
        tic;
        [dist_mat,knn] = compute_knn_labelled(x_mu(1:embed_dim,:), k, train_filt);
        [acc, f1, auc, class_score_mashup] = fun_pred_majority_voting(anno, test_filt,train_filt, knn, dist_mat, weighted);
        acc_mu(i,j) = acc;
        f1_mu(i,j)  = f1;
        auc_mu(i,j) = auc;
        toc;
    
        %-------------------------------
        % Perform Melissa 
        %-------------------------------
        fprintf('[Perfoming Melissa]\n');
        tic;
        [dist_mat,knn] = compute_knn_labelled(x_melissa(1:embed_dim,1:ngene), k, train_filt);
        [acc, f1, auc] = fun_pred_majority_voting_mashup_backup(anno, test_filt,train_filt, knn, dist_mat, class_score_mashup, weighted);

        acc_melissa(i,j) = acc;
        f1_melissa(i,j)  = f1;
        auc_melissa(i,j) = auc;
        toc;
    
    end

end 

%-------------------------------
% Summary Performance
%-------------------------------
acc_mu_ave = mean(acc_mu,1);
[val, ind] = max(acc_mu_ave);
fprintf('Mashup Accuracy = %0.4f at embedding dimension = %d \n', val, ind*25);

f1_mu_ave = mean(f1_mu,1);
[val, ind] = max(f1_mu_ave);
fprintf('Mashup F1 = %0.4f at embedding dimension = %d \n', val, ind*25);

auc_mu_ave = mean(auc_mu,1);
[val, ind] = max(auc_mu_ave);
fprintf('Mashup AUROC = %0.4f at embedding dimension = %d \n', val, ind*25);

acc_melissa_ave = mean(acc_melissa,1);
[val, ind] = max(acc_melissa_ave);
fprintf('MELISSA Accuracy = %0.4f at embedding dimension = %d \n', val, ind*25);

f1_melissa_ave = mean(f1_melissa,1);
[val, ind] = max(f1_melissa_ave);
fprintf('MELISSA F1 = %0.4f at embedding dimension = %d \n', val, ind*25);

auc_melissa_ave = mean(auc_melissa,1);
[val, ind] = max(auc_melissa_ave);
fprintf('MELISSA AUROC = %0.4f at embedding dimension = %d \n', val, ind*25);


