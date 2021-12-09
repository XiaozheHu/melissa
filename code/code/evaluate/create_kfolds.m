function folds = create_kfolds(anno, options)

addpath code/evaluate;

kfolds = options.kfolds;

ng = size(anno, 2);
% filters proteins with no labels
label_filt = (sum(anno) > 0).'; 

ntest = floor(ng / kfolds);
rng(2021)
perm = randperm(ng);


for i = kfolds:-1:1 % handle the last fold to avoid bad indexing
    last = i * ntest;
    if i == kfolds % handles ng not being a multiple of kfold
        last = ng;
    end

    test_ind = perm((i - 1) * ntest + 1 : last);
    test_filt = false(ng,1);
    test_filt(test_ind) = true;

    train_filt = (~test_filt) & label_filt;
    test_filt = test_filt & label_filt;
    
    folds(i).train_filt   = train_filt;
    folds(i).test_filt    = test_filt;
    folds(i).ntrain       = sum(test_filt);
    folds(i).ntest        = sum(train_filt);
    folds(i).train_labels = anno.*(train_filt.');
    folds(i).test_labels  = anno.*(test_filt.');
end
