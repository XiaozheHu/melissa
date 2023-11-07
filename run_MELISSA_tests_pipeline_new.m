clear;clc;
addpath(genpath('code'));
addpath(genpath('data'));
addpath(genpath('mtba'));

%------------------------
% Example parameters
%------------------------

% random seed
rand_list = {2021,2023,0,272,275};

% use human or yeast data
options.org = 'yeast';

% which type of annotations to use
% options: {bp, mf, cc} for human and yeast GO labels
% options.onttype = 'mf'; 
type_list = {"bp", "mf", "cc"};

% filter terms into a specific size range 
% examples: [11 30], [31 100], [101 300]
% options.ontsize = [31 100];
ont_size_list = {[11 30], [31 100], [101 300]};
     
% number of kNN voters
k=10;

% Number of bi-clusters to create (-1 as not bi-cluster)
options.num_clusters = 2; 

% use SVD approximation for Mashup
% recommended: true for human, false for yeast
options.embedding.svd_approx = false;

% whether to stack the matrics really tall for svd or use the
% log transofrm and sum - needs implementation of dummy nodes
% options.embedding.svd_full = false;

% number of dimensions
% this is used as the maximal embeded dimension, tests will be run on
% different dimensions 
options.embedding.ndim =400; 

% the weight of the edges connecting dummy nodes to true nodes
% choose value between 0 and 1
options.embedding.mustlink_penalty = 1; 

% the weight of the edges connecting dummy nodes to dummy nodes
% choose value between 0 and 1
options.embedding.cannotlink_penalty = 1; 

% chance that the random walk restarts itself
options.walk.restart_prob = 0.5;

% number of folds for k-fold cross validation
% if set to 1 or less then a single experiment is run
options.kfolds = 5;

% if options.kfolds is set to 1 or less then this is
% the portion of the labelled vertices to go in the 
% testing set. (1 - test_fraction) is used to train
options.test_fraction = 0.2;

fileID = fopen('run_result_log_human_wofilter_2023.txt','w');
                    
%Logs the options so we can see the parameters used in log file later
for k3 = 1:length(rand_list)
    rand_num = rand_list{k3};
    fprintf(fileID,"#################################### \n");
    fprintf(fileID,"rand num: %d \n",rand_num);
    fprintf(fileID,"#################################### \n");
    for k1 = 1:length(type_list)
        for k2 = 1:length(ont_size_list)
            options.onttype = type_list{k1};
            options.ontsize = ont_size_list{k2};
            log_options(options);
            fprintf(fileID,"ont type: %s \n",options.onttype);
            fprintf(fileID,"ont size: %d - %d  \n",options.ontsize(1),options.ontsize(2));

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
            % [genes, ngene, anno] = load_anno(options);
            [genes, ngene, anno] = load_anno_wo_filter(options);

            fprintf('Number of functional labels: %d\n', size(anno, 1));

            if options.kfolds <= 1
                % Generate training and testing sets
                fprintf('Acquiring test filter using %d testing fraction\n', options.test_fraction);
                folds = create_kfolds(anno, options,rand_num);
            else
                folds = create_kfolds(anno, options,rand_num);
            end

            % save("variables");

            % compute the network walks for both mu and f.
            fprintf('Compute random walk on original networks \n');
            walks = compute_rwr_original_sparse(network_files, ngene, options);
                               
            % % Mashup embedding is shared across all folds so we take it before cross-val
            fprintf('Compute Mashup embedding \n');

            [x_mu, d_mu] = svd_embed_sparse_func(walks, ngene, options.embedding.ndim);

            clear walks;

            %-----------------------------------------------
            % Perform function prediction by majority voting
            %-----------------------------------------------
            % step size of embeded dimension 
            % dim_step = 25;
            % n_dim_test = floor(options.embedding.ndim / dim_step);
            n_dim_test = 1;
            % 
            % Mashup results
            acc_mu = zeros(length(folds), n_dim_test);
            f1_mu  = zeros(length(folds), n_dim_test);
            auc_mu = zeros(length(folds), n_dim_test);

            % MELISSA results
            acc_melissa = zeros(length(folds), n_dim_test);
            f1_melissa  = zeros(length(folds), n_dim_test);
            auc_melissa = zeros(length(folds), n_dim_test);

            acc_deepNF = zeros(length(folds), n_dim_test);
            f1_deepNF  = zeros(length(folds), n_dim_test);
            auc_deepNF = zeros(length(folds), n_dim_test);

            weighted = true;
            fprintf('weighted: true \n');

            %------------------------
            % main loop
            %------------------------
            %file_name = strcat('deepNF_MDA_arch_15000-9000-1200-9000-15000_features.csv'); %human
            file_name = strcat('deepNF_MDA_arch_12000-600-12000_features.csv'); %yeast
            x_deepNF = readtable(file_name, "FileType","text",'Delimiter', ',');
            x_deepNF = table2array(x_deepNF);
            x_deepNF = transpose(x_deepNF);

            for i = 1:length(folds)

                fprintf('Fold %d / %d \n', i, options.kfolds);
                fprintf('Using k = %f \n', k);

                % prepare train and test sets
                train_filt = folds(i).train_filt;
                test_filt = folds(i).test_filt;
                training_labels = anno.*(train_filt.');

                % perform biclustering
                [gene_clusters, label_clusters] = bicluster(training_labels, train_filt, options);
                % [gene_clusters, label_clusters] = bicluster_spc(training_labels, train_filt, options);
                
                % augment the network and run random walk
                fprintf('Compute argumented random walk \n');
                tic;
                walks_argumented = compute_rwr_argumented_sparse(network_files, ngene, gene_clusters, options);
                toc;

                % compute melissa embedding
                fprintf('Compute Melissa embedding \n');
                tic;
                if (options.embedding.cannotlink_penalty == 0)
                    [x_melissa, d_melissa] = svd_embed_sparse_func(walks_argumented, ngene+options.num_clusters, options.embedding.ndim);
                else
                    [x_melissa, d_melissa] = svd_embed_CL_sparse_func(walks_argumented, ngene+options.num_clusters, options.embedding.ndim);
                end
                toc;
                
                clear walks_argumented;

                
                %-------------------------------
                % Perform prediction
                %-------------------------------
                for j = 1:n_dim_test
                
                    % embed_dim = j*dim_step;
                    embed_dim = options.embedding.ndim;
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

                    fprintf('[Perfoming deepNF embedding]\n');
                    tic;
                    [dist_mat,knn] = compute_knn_labelled(x_deepNF(:,1:ngene), k, train_filt);
                    [acc, f1, auc,class_score_deepNF] = fun_pred_majority_voting(anno, test_filt,train_filt, knn, dist_mat, weighted);
                    % [acc, f1, auc] = fun_pred_majority_voting_mashup_backup(anno, test_filt,train_filt, knn, dist_mat, class_score_mashup, weighted);
                
                    acc_deepNF(i,j) = acc;
                    f1_deepNF(i,j)  = f1;
                    auc_deepNF(i,j) = auc;
                    toc;
                
                end

            end 

            % -------------------------------
            % Summary Performance
            % -------------------------------
            acc_mu_ave = mean(acc_mu,1);
            [val, ind] = max(acc_mu_ave);
            fprintf('Mashup Accuracy = %0.4f at embedding dimension = 400 \n', val);
            fprintf(fileID, 'Mashup Accuracy = %0.4f at embedding dimension = 400 \n', val);

            f1_mu_ave = mean(f1_mu,1);
            [val, ind] = max(f1_mu_ave);
            fprintf('Mashup F1 = %0.4f at embedding dimension = 400 \n', val);
            fprintf(fileID, 'Mashup F1 = %0.4f at embedding dimension = 400 \n', val);

            auc_mu_ave = mean(auc_mu,1);
            [val, ind] = max(auc_mu_ave);
            fprintf('Mashup AUROC = %0.4f at embedding dimension = 400 \n', val);
            fprintf(fileID, 'Mashup AUROC = %0.4f at embedding dimension = 400 \n', val);

            acc_melissa_ave = mean(acc_melissa,1);
            [val, ind] = max(acc_melissa_ave);
            fprintf('MELISSA Accuracy = %0.4f at embedding dimension = 400 \n', val);
            fprintf(fileID, 'MELISSA Accuracy = %0.4f at embedding dimension = 400 \n', val);

            f1_melissa_ave = mean(f1_melissa,1);
            [val, ind] = max(f1_melissa_ave);
            fprintf('MELISSA F1 = %0.4f at embedding dimension = 400 \n', val);
            fprintf(fileID, 'MELISSA F1 = %0.4f at embedding dimension = 400 \n', val);

            auc_melissa_ave = mean(auc_melissa,1);
            [val, ind] = max(auc_melissa_ave);
            fprintf('MELISSA AUROC = %0.4f at embedding dimension = 400 \n', val);
            fprintf(fileID, 'MELISSA AUROC = %0.4f at embedding dimension = 400 \n', val);

            acc_deepNF_ave = mean(acc_deepNF,1);
            [val, ind] = max(acc_deepNF_ave);
            fprintf('original deepNF Accuracy = %0.4f  \n', val);
            fprintf(fileID,'original deepNF Accuracy = %0.4f  \n', val);

            f1_deepNF_ave = mean(f1_deepNF,1);
            [val, ind] = max(f1_deepNF_ave);
            fprintf('original deepNF F1 = %0.4f \n', val);
            fprintf(fileID,'original deepNF F1 = %0.4f  \n', val);

            auc_deepNF_ave = mean(auc_deepNF,1);
            [val, ind] = max(auc_deepNF_ave);
            fprintf('original deepNF AUROC = %0.4f  \n', val);
            fprintf(fileID,'original deepNF AUROC = %0.4f  \n', val);
            
        end 
    end 
end
fclose(fileID);
