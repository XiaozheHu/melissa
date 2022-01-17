% get k dimensional Mashup embedding and then to SSDR to update the
% embedding (same dimension)

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
options.onttype = 'mf';

% consider terms in a specific size range (GO only)
% examples: [11 30], [31 100], [101 300]
options.ontsize = [31 100];

% number of kNN
%options.kNN_neighbor=10;

% Coefficient of strength of CL constraint
% In the paper, ML coefficient = 1, CL coefficient need to be tuned
% CL = 8;

% Dimension of SSDR embedding
options.SSDR_embedding = 50;

% use 1-p percent of the training labels to form extra graph
options.extra_graph_filt = 0.95;

% Number of bi-clusters to create (-1 to not bi-cluster)
%options.num_clusters = 4;

% use SVD approximation for Mashup
% recommended: true for human, false for yeast
%options.embedding.svd_approx = true;

% whether to stack the matrics really tall for svd or use the
% log transofrm and sum - needs implementation of dummy nodes
% options.embedding.svd_full = false;

% number of dimensions
% recommended: 800 for human, 500 for yeast
options.embedding.ndim = 1000;

% the weight of the edges connecting dummy nodes to true nodes
options.embedding.mustlink_penalty = 1;

% the weight of the edges connecting dummy nodes to dummy nodes
options.embedding.cannotlink_penalty = 8;

% whether to add 1/ngene^2 strength constraint between all genes
%options.embedding.use_unsupervised = false;

% when using go, whether or not to append the extra link matrix
% generated from the labels
%options.walk.use_go_link = false;

% when using go and the link matrix, what fraction of links to use
%options.walk.go_link_fraction = 1.0;

% chance that the random walk restarts itself
options.walk.restart_prob = 0.5;

% number of folds for k-fold cross validation
% if set to 1 or less then a single experiment is run
options.kfolds = 5;

% if options.kfolds is set to 1 or less then this is
% the portion of the labelled vertices to go in the 
% testing set. (1 - test_fraction) is used to train
options.test_fraction = 0.2;

% fix parameters for SVM
options.svm.gmax = 3; 
options.svm.cmax = 4;

%-----------------------------------------------
% Logs the options so we can see the parameters used in log file later
%-----------------------------------------------
log_options(options);

%-----------------------------------------------
% Construct network file paths
%-----------------------------------------------
string_nets = {'neighborhood', 'fusion', 'cooccurence', 'coexpression', ...
               'experimental', 'database'};
network_files = cell(1, length(string_nets));
for i = 1:length(string_nets)
     network_files{i} = sprintf('data/networks/%s/%s_string_%s_adjacency.txt', ...
         options.org, options.org, string_nets{i});
end

%-----------------------------------------------
% Load gene list
%-----------------------------------------------
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

%-----------------------------------------------
% Performs the specified variant of RWR
%-----------------------------------------------
fprintf('[Performing RWR step]\n');

% compute the network walks for both mu and f.
fprintf('Compute random walk on original networks \n');
%walks = compute_rwr(network_files, ngene, -1, options);
load('walks.mat')

% mu embedding is shared across all folds so we take it before cross-val
fprintf('Compute Mashup embedding \n');
tic;
x_mu = svd_embed(walks, options.embedding.ndim);
toc;
%x_full = svd_embed(walks,ngene);
%x_mu = x_full(1:options.embedding.ndim, :);

%-----------------------------------------------
% Perform function prediction by majority voting
%-----------------------------------------------
% step size of embeded dimension 
dim_step = 25;
n_dim_test = floor(options.embedding.ndim / dim_step);

% mashup results
acc_mu = zeros(length(folds), n_dim_test);
f1_mu = zeros(length(folds), n_dim_test);
auc_mu = zeros(length(folds), n_dim_test);

% SSDR results
acc_SSDR = zeros(length(folds), n_dim_test);
f1_SSDR = zeros(length(folds), n_dim_test);
auc_SSDR = zeros(length(folds), n_dim_test);

%weighted = true;
%fprintf('weighted: true \n');

for i = 1:length(folds)
    fprintf('Fold %d / %d \n', i, options.kfolds);
    
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
    
    %S(S == -2) = 1/(ngene)^2 + 1/Nc;
    %S(S == 1) = 1/(ngene)^2 - CL/Nm;
    
    S(S == -2) = 1/(ngene)^2 - options.embedding.cannotlink_penalty/Nc;
    S(S == 1) = 1/(ngene)^2 + options.embedding.mustlink_penalty/Nm;
    
    if ~isequal(S, S') % symmetrize
      S = (S + S')/2;
    end

    S = S - spdiags(diag(S),0,ngene,ngene);
    L = sum(S,2)-S;
    
    %-------------------------------
    % Perform prediction
    %-------------------------------
    for j = 1:n_dim_test
        
        embed_dim = j*dim_step;
        fprintf('Using embeded dimension = %f \n', embed_dim);
        
        %-------------------------------
        % Perform Mashup for comparison
        %-------------------------------
        fprintf('[Performing Mashup version]\n');   
        tic;
        %[acc, f1, aupr, gmax, cmax] = run_svm(x_mu(1:embed_dim,:), anno, test_filt);  
        [acc, f1, aupr] =  run_svm_fixed_hp(x_mu(1:embed_dim,:),...
            anno, test_filt, options.svm.gmax, options.svm.cmax);
        acc_mu(i,j) = acc;
        f1_mu(i,j) = f1;
        auc_mu(i,j) = aupr;
        toc;

        %-------------------------------
        % Compute SSDR
        %-------------------------------
        fprintf('Compute SSDR embedding \n');
        tic;
        X = x_mu(1:embed_dim,:)*L*(x_mu(1:embed_dim,:))';
        [V, d] = eigs(X, embed_dim);
        x_SSDR = V'*x_mu(1:embed_dim,:);
        toc;
        
        %-------------------------------
        % Perform SSDR
        %-------------------------------
        fprintf('[Perfoming SSDR version]\n');
        tic;
        %[acc, f1, aupr, gmax, cmax] = run_svm(x_SSDR(1:embeded_dim,:), anno, test_filt);
        [acc, f1, aupr] =  run_svm_fixed_hp(x_SSDR(1:embed_dim,:), ...
            anno, test_filt, options.svm.gmax, options.svm.cmax);
        acc_SSDR(i,j) = acc;
        f1_SSDR(i,j) = f1;
        auc_SSDR(i,j) = aupr;
        toc;
    
    end
    
end

dim = dim_step:dim_step:options.embedding.ndim;

figure(1);
plot(dim, mean(acc_mu));
hold on;
plot(dim, mean(acc_SSDR), '--');

figure(2);
plot(dim, mean(f1_mu));
hold on;
plot(dim, mean(f1_SSDR), '--');

figure(3);
plot(dim, mean(auc_mu));
hold on;
plot(dim, mean(auc_SSDR), '--');
