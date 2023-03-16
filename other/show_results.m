[val,idx] = max(mean(acc_mu));
fprintf('Mashup:  acc = %f, dim = %d\n', val, dim(idx));

%[val,idx] = max(mean(acc_melissa));
%fprintf('Melissa: acc = %f, dim = %d\n', val, dim(idx));

[val,idx] = max(mean(acc_SSDR));
fprintf('SSDR: acc = %f, dim = %d\n', val, dim(idx));

[val,idx] = max(mean(f1_mu));
fprintf('Mashup:  f1 = %f, dim = %d\n', val, dim(idx));

%[val,idx] = max(mean(f1_melissa));
%fprintf('Melissa: f1 = %f, dim = %d\n', val, dim(idx));

[val,idx] = max(mean(f1_SSDR));
fprintf('SSDR: f1 = %f, dim = %d\n', val, dim(idx));

[val,idx] = max(mean(auc_mu));
fprintf('Mashup:  auc = %f, dim = %d\n', val, dim(idx));

%[val,idx] = max(mean(auc_melissa));
%fprintf('Melissa: auc = %f, dim = %d\n', val, dim(idx));

[val,idx] = max(mean(auc_SSDR));
fprintf('SSDR: auc = %f, dim = %d\n', val, dim(idx));
