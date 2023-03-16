%clear;clc;
addpath(genpath('code'));
addpath(genpath('data'));
addpath(genpath('mtba'));

%% Example parameters

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

% Coefficient of strength of CL constraint
% In the paper, ML coefficient = 1, CL coefficient need to be tuned
CL = 8;

% Dimension of SSDR embedding
SSDR_embedding = 50;

% use 1-p percent of the training labels to form extra graph
extra_graph_filt = 0.95;

% Number of bi-clusters to create (-1 to not bi-cluster)
options.num_clusters = 4;

% use SVD approximation for Mashup
% recommended: true for human, false for yeast
options.embedding.svd_approx = true;

% whether to stack the matrics really tall for svd or use the
% log transofrm and sum - needs implementation of dummy nodes
% options.embedding.svd_full = false;

% number of dimensions
% recommended: 800 for human, 500 for yeast
options.embedding.ndim = 500;

% the weight of the edges connecting dummy nodes to true nodes
options.embedding.mustlink_penalty = 1;

% the weight of the edges connecting dummy nodes to dummy nodes
options.embedding.cannotlink_penalty = 64;

% whether to add 1/ngene^2 strength constraint between all genes
options.embedding.use_unsupervised = false;

% when using go, whether or not to append the extra link matrix
% generated from the labels
options.walk.use_go_link = false;

% when using go and the link matrix, what fraction of links to use
options.walk.go_link_fraction = 1.0;

% chance that the random walk restarts itself
options.walk.restart_prob = 0.5;

% number of folds for k-fold cross validation
% if set to 1 or less then a single experiment is run
options.kfolds = 5;

% if options.kfolds is set to 1 or less then this is
% the portion of the labelled vertices to go in the 
% testing set. (1 - test_fraction) is used to train
options.test_fraction = 0.2;

%% Logs the options so we can see the parameters used in log file later
log_options(options);

%% Construct network file paths
string_nets = {'neighborhood', 'fusion', 'cooccurence', 'coexpression', ...
               'experimental', 'database'};
network_files = cell(1, length(string_nets));
for i = 1:length(string_nets)
     network_files{i} = sprintf('data/networks/%s/%s_string_%s_adjacency.txt', ...
         options.org, options.org, string_nets{i});
end

%% Load gene list
fprintf('[Loading annotations]\n');
[genes, ngene, anno] = load_anno(options);

fprintf('Number of functional labels: %d\n', size(anno, 1));

if options.kfolds <= 1
    %% Generate training and testing sets
    fprintf('Acquiring test filter using %d testing fraction\n', options.test_fraction);
    folds = create_kfolds(anno, options);
else
    folds = create_kfolds(anno, options);
end

%% Performs the specified variant of RWR
fprintf('[Performing RWR step]\n');

% compute the network walks for both mu and f.
%walks = compute_rwr(network_files, ngene, -1, options);
load('walks.mat')

% mu embedding is shared across all folds so we take it before cross-val
%x_mu = svd_embed(walks, options.embedding.ndim);
x_full = svd_embed(walks,ngene);
x_mu = x_full(1:options.embedding.ndim, :);

%% Perform function prediction by majority voting
% mashup results
acc_mu = zeros(length(folds), 1);
f1_mu = zeros(length(folds), 1);
auc_mu = zeros(length(folds), 1);

% SSDR results
acc_SSDR = zeros(length(folds), 1);
f1_SSDR = zeros(length(folds), 1);
auc_SSDR = zeros(length(folds), 1);

weighted = true;
fprintf('weighted: true \n');

for i = 1:length(folds)
    fprintf('Fold %d / %d \n', i, options.kfolds);
    fprintf('Using k = %f \n', k);

    train_filt = folds(i).train_filt;
    test_filt = folds(i).test_filt;

    training_labels = anno.*(train_filt.');
     
    % add extra links
    S = ones(ngene)/(ngene^2);
    % add cannotlink if no label overlaps at all
    % add mustlink if share m_min or more labels
    m_min = 1;
 
    % count number of CL and ML constraint
    Nm = 0; % number of ML constraints
    Nc = 0; % number of CL constraints

    for p = 1:size(anno,2)
        annoi = find(training_labels(:,p)); % label of gene i
        if ~isempty(annoi)
            for q = p:size(anno,2)
                annoj = find(training_labels(:,q)); % label of gene i
                if ~isempty(annoj)
                    shared = intersect(annoi,annoj);
                    % add cannotlink  
                    if isempty(shared)
                        S(p,q) = -2;
                        Nc = Nc + 1;
%                         if rand>0.8  % only add 20% of CL constraints
%                             S(p,q) = -2;
%                             Nc = Nc + 1;
%                        end
                       % add mustlink                      
                     elseif length(shared)>=m_min
                         S(p,q) = 1;
                         Nm = Nm + 1;
          
                    end
                end
            end
        end
    end
    
    S(S == -2) =1/(ngene)^2 - CL/Nc;
    S(S == 1) = 1/(ngene)^2 + 1/Nm;
    
    if ~isequal(S, S') % symmetrize
      S = (S + S')/2;
    end

    S = S - diag(diag(S));
    L = sum(S,2)-S;
    
    [V, d] = eigs(x_mu*L*x_mu', SSDR_embedding);
    x_SSDR = V'*x_mu;


    %% Performs Mashup for comparison
    fprintf('[Performing mu version]\n');
    
    [dist_mat,knn] = compute_knn_labelled(x_mu, k, train_filt);
    
%     % also, compare to mashup embedding at the same dimension
%     x_mu1 = x_full(1:SSDR_embedding, :);
%     [dist_mat,knn] = compute_knn_labelled(x_mu1, k, train_filt);
    
    [acc, f1, auc, class_score_mashup] = fun_pred_majority_voting(anno, test_filt,train_filt, knn, dist_mat, weighted);
    acc_mu(i) = acc;
    f1_mu(i) = f1;
    auc_mu(i) = auc;

%% Perform SSDR version
    fprintf('[Perfoming f version]\n');
    [dist_mat,knn] = compute_knn_labelled(x_SSDR, k, train_filt);
    %[acc, f1, auc] = fun_pred_majority_voting_mashup_backup(anno, test_filt,train_filt, knn, dist_mat, class_score_mashup, weighted);

    [acc, f1, auc, ~] = fun_pred_majority_voting(anno, test_filt,train_filt, knn, dist_mat, weighted);

   acc_SSDR(i) = acc;
   f1_SSDR(i) = f1;
   auc_SSDR(i) = auc;
end


fprintf('[mu accuracy mean = %f ]\n', mean(acc_mu));
fprintf('[mu accuracy std = %f ]\n', std(acc_mu));
fprintf('[mu f1 mean = %f ]\n', mean(f1_mu));
fprintf('[mu f1 std = %f ]\n', std(f1_mu));
fprintf('[mu auc mean = %f ]\n', mean(auc_mu));
fprintf('[mu auc std = %f ]\n', std(auc_mu));

fprintf('[SSDR accuracy mean = %f ]\n', mean(acc_SSDR));
fprintf('[SSDR accuracy std = %f ]\n', std(acc_SSDR));
fprintf('[SSDR f1 mean = %f ]\n', mean(f1_SSDR));
fprintf('[SSDR f1 std = %f ]\n', std(f1_SSDR));
fprintf('[SSDR auc mean = %f ]\n', mean(auc_SSDR));
fprintf('[SSDR auc std = %f ]\n', std(auc_SSDR));


