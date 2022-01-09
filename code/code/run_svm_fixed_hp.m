function [acc, f1, aupr] = run_svm(x, anno, test_filt, gmax, cmax)
    addpath code/evaluate;
    % Parameters
    gvec = -3:1:0;
    cvec = -2:1:2;
  
    % Scale features
    maxval = max(x, [], 2);
    minval = min(x, [], 2);
    x = bsxfun(@times, bsxfun(@minus, x, minval), 1 ./ (maxval - minval));
  
    % Filter genes with no annotations
    filt = sum(anno) > 0;
    anno = anno(:,filt);
    x = x(:,filt);
  
    [nclass, ngene] = size(anno);
      
    test_filt = test_filt(filt);
    ntest = sum(test_filt);
    ntrain = ngene - ntest;
  
    K = rbf_kernel(x', x', 10^gvec(gmax));
  
    % use the best parameters we computed above
    Ktrain = K(~test_filt,~test_filt);
    Ktest = K(test_filt,~test_filt);
  
    class_score = zeros(ntest, nclass);
    parfor s = 1:nclass
        Ytrain = full(double(anno(s,~test_filt)') * 2 - 1);
        Ytest = full(double(anno(s,test_filt)') * 2 - 1);
  
        model = svmtrain(Ytrain, [(1:ntrain)', Ktrain], ['-t 4 -b 1 -q -c ', num2str(10^cvec(cmax))]);
        posind = find(model.Label > 0);
        if ~isempty(posind)
            [~, ~, dec] = svmpredict(Ytest, [(1:ntest)', Ktest], model, '-q');
            class_score(:,s) = dec(:,posind);
        end
    end
  

    [acc, f1, aupr] = evaluate_performance(class_score, anno(:,test_filt)');
    fprintf('acc: %f, f1: %f, aupr: %f\n', acc, f1, aupr);
end
  
