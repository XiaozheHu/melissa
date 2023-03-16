function folds = create_train_test(anno, options)
    
    addpath code/evaluate;
    
    test_fraction = options.test_fraction;
    
    [~, test_filt] = cv_partition(anno, test_fraction); 
    % filters proteins with no labels
    label_filt = (sum(anno) > 0).'; 
    % removes the unlabeled proteins from the train and test sets
    folds(1).train_filt =(~test_filt) & label_filt;
    folds(1).test_filt = test_filt & label_filt;

    folds(1).ntest = sum(test_filt);
    folds(1).ntrain = sum(train_filt);

    % applies the filters as masks to the annotations
    folds(1).train_labels = anno.*(train_filt.');
    folds(1).test_labels = anno.*(test_filt.');
end

