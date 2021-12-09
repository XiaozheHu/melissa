function [acc, f1, auprc] = evaluate_performance(class_score, label)
%-------------------------------------------------------
% compute the accuracy, f1, and auprc scores for the predictions
% Input:    class_score: matrix of label confidences. Size is 
%                        [nelements, nclasses]
%           label: matrix of true labels Size is 
%                  [nelements, nclasses]
% Output:   acc: percent of top guesses that hit a correct label
%           f1: f1 score of the predictions
%           auprc: area under the precision recall curve
%-------------------------------------------------------
  alpha = 3;

  label = label > 0;

  [ncase, nclass] = size(class_score);

  [~,o] = sort(class_score,2,'descend');
  p = sub2ind(size(label),(1:ncase)',o(:,1));
  acc = mean(label(p));

  a = repmat((1:ncase)',1,alpha);
  pred = sparse(a, o(:,1:alpha), 1, size(label,1),size(label,2));

  % creates a 2x2 matrix with (1,1) = False negs, (1,2) = True negs
  %                           (2,1) = False pos,  (2,2) = True pos
  tab = crosstab([0;1;pred(:)],[0;1;label(:)]) - eye(2);
  
  f1 = 2*tab(2,2) / (2*tab(2,2)+tab(1,2)+tab(2,1));

  [~, auprc] = auc(label(:), class_score(:));
end
