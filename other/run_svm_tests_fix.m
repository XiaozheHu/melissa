%clear;clc;
addpath(genpath('code'));
addpath(genpath('data'));
addpath(genpath('mtba'));
addpath(genpath('libsvm'));

%------------------------
% Example parameters
%------------------------
% use human or yeast data
options.org = 'yeast';

% which type of annotations to use
% options: {bp, mf, cc} for human GO,
%          {level1, level2, level3} for yeast MIPS
options.onttype = 'bp'; 

% consider terms in a specific size range (GO only)
% examples: [11 30], [31 100], [101 300]
options.ontsize = [11 30];
     
% number of kNN
k=10;

% Number of bi-clusters to create (-1 to not bi-cluster)
options.num_clusters = 2; 

% use SVD approximation for Mashup
% recommended: true for human, false for yeast
options.embedding.svd_approx = true;

% whether to stack the matrics really tall for svd or use the
% log transofrm and sum - needs implementation of dummy nodes
% options.embedding.svd_full = false;

% number of dimensions
% this is used as the maximal embeded dimension, tests will be run on
% different dimensions 
options.embedding.ndim = 1000; 

% the weight of the edges connecting dummy nodes to true nodes
options.embedding.mustlink_penalty = 1; 

% the weight of the edges connecting dummy nodes to dummy nodes
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

% fix parameters for SVM
gmax = 3; 
cmax = 4;

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

% for z1 = 1:3
%     
% if z1 == 1
%     options.onttype = 'bp';
% elseif z1 == 2
%     options.onttype = 'mf';
% else
%     options.onttype = 'cc';
% end
% 
% fprintf(strcat("STARTING ONT: ", options.onttype, "\n"));
% 
% for z2 = 1:3
% 
% if z2 == 1
%     options.ontsize = [11 30];
% elseif z2 == 2
%     options.ontsize = [31 100];
% else
%     options.ontsize = [101 300];
% end
% 
% fprintf(strcat("STARTING RANGE: ", num2str(z2), "\n"));
%------------------------
% Load gene list
%------------------------
%fprintf('[Loading annotations]\n');
%[genes, ngene, anno] = load_anno(options);
[genes, ngene, anno, ~, ~, filt] = load_anno_full(options);

%fprintf('Number of functional labels: %d\n', size(anno, 1));

if options.kfolds <= 1
    % Generate training and testing sets
    %fprintf('Acquiring test filter using %d testing fraction\n', options.test_fraction);
    folds = create_kfolds(anno, options);
else
    folds = create_kfolds(anno, options);
end

% compute the network walks for both mu and f.
fprintf('Compute random walk on original networks \n');
%tic;
%if (ngene >= 10000) % yeast use full matrix and human use sparse matrix
%    walks = compute_rwr_original_sparse(network_files, ngene, options);
%else 
%    walks = compute_rwr_original(network_files, ngene, options);
%end
%toc;
load('walks.mat')
%load('human_walks.mat');

% mu embedding is shared across all folds so we take it before cross-val
fprintf('Compute Mashup embedding \n');
tic;
[x_mu, d_mu] = svd_embed(walks, options.embedding.ndim);
toc;
%load('human_mashup_embed.mat');

clear walks;

%-----------------------------------------------
% Perform function prediction by majority voting
%-----------------------------------------------
% step size of embeded dimension 
dim_step = 25;
n_dim_test = floor(options.embedding.ndim / dim_step);

% mashup results
acc_mu = zeros(length(folds), n_dim_test);
f1_mu  = zeros(length(folds), n_dim_test);
auc_mu = zeros(length(folds), n_dim_test);

% melissa results
acc_melissa = zeros(length(folds), n_dim_test);
f1_melissa  = zeros(length(folds), n_dim_test);
auc_melissa = zeros(length(folds), n_dim_test);
x_melissa_all = cell(length(folds),1);

weighted = true;
fprintf('weighted: true \n');

%------------------------
% main loop
%------------------------
for i = 1:length(folds)
    
    fprintf('Fold %d / %d \n', i, options.kfolds);
    %fprintf('Using k = %f \n', k);
    
    % prepare train and test sets
    train_filt = folds(i).train_filt;
    test_filt = folds(i).test_filt;
     
    training_labels = anno.*(train_filt.');
    
    % performe biclustering
    [gene_clusters, label_clusters] = bicluster(training_labels, train_filt, options);

    % argument the network and run random walk
    fprintf('Compute argumented random walk \n')
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
    x_melissa_all{i} = x_melissa;
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
        [acc, f1, aupr] = run_svm_fixed_hp(x_mu(1:embed_dim,1:ngene), anno, test_filt, gmax, cmax);
        acc_mu(i,j) = acc;
        f1_mu(i,j)  = f1;
        auc_mu(i,j) = aupr;
        toc;
         
        %-------------------------------
        % Perform Melissa 
        %-------------------------------
        fprintf('[Perfoming Melissa]\n');
        tic;
        [acc, f1, aupr] =  run_svm_fixed_hp(x_melissa(1:embed_dim,1:ngene), anno, test_filt, gmax, cmax);
        acc_melissa(i,j) = acc;
        f1_melissa(i,j)  = f1;
        auc_melissa(i,j) = aupr;
        toc;
        
     end    
    
end 
%end
%end

dim = dim_step:dim_step:options.embedding.ndim;

figure(1);
plot(dim, mean(acc_mu));
hold on;
plot(dim, mean(acc_melissa), '--');

figure(2);
plot(dim, mean(f1_mu));
hold on;
plot(dim, mean(f1_melissa), '--');

figure(3);
plot(dim, mean(auc_mu));
hold on;
plot(dim, mean(auc_melissa), '--');

