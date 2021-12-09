function [acc, f1, aupr] = run_svm(x, anno, test_filt)
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
  
    rbfK = cell(length(gvec), 1);
    for i = 1:length(gvec)
        rbfK{i} = rbf_kernel(x', x', 10^gvec(i));
    end
  
      
    retmax = -inf;
    gmax = 1;
    cmax = 1;
    for gi = 1:length(gvec)
        for ci = 1:length(cvec)
            tt = tic;

            Ktrain = rbfK{gi}(~test_filt,~test_filt);
            Ktest = rbfK{gi}(test_filt,~test_filt);
  
            class_score = zeros(ntest, nclass);
            parfor s = 1:nclass
                Ytrain = full(double(anno(s,~test_filt)') * 2 - 1);
                Ytest = full(double(anno(s,test_filt)') * 2 - 1);
  
                model = svmtrain(Ytrain, [(1:size(Ktrain,1))', Ktrain], ['-t 4 -b 1 -q -c ', num2str(10^cvec(ci))]);
                posind = find(model.Label > 0);
                if ~isempty(posind)
                    [~, ~, dec] = svmpredict(Ytest, [(1:size(Ktest,1))', Ktest], model, '-q');
                    class_score(:,s) = dec(:,posind);
                end
            end
            % TODO make it just the passed indices
            [~, ~, cv_result] = evaluate_performance(class_score, anno(:,test_filt)');

            if retmax < cv_result
                retmax = cv_result;
                gmax = gi;
                cmax = ci;
            end
  
            fprintf('gi:%d, ci:%d, cv_result:%f, ', gi, ci, cv_result); toc(tt)
        end
    end
  

    % use the best parameters we computed above
    Ktrain = rbfK{gmax}(~test_filt,~test_filt);
    Ktest = rbfK{gmax}(test_filt,~test_filt);
  
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
  
