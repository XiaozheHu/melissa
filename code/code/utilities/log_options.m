function success = log_options(options)
  %% Logs the options so we can see the parameters used in log file later
  tf = {'false','true'};
  fprintf('\n\n\n---------------PARAMETERS--------------\n');
  fprintf('options.org: %s\n', options.org);
  fprintf('options.onttype: %s\n', options.onttype);
  fprintf('options.ontsize: [%d, %d]\n', options.ontsize(1), options.ontsize(2));
  fprintf('options.num_clusters: %d\n', options.num_clusters);
  fprintf('options.embedding.ndim: %d\n', options.embedding.ndim);
  fprintf('options.embedding.svd_approx: %s\n', tf{options.embedding.svd_approx + 1});
  fprintf('options.embedding.mustlink_penalty: %d\n', options.embedding.mustlink_penalty);
  fprintf('options.embedding.cannotlink_penalty: %d\n', options.embedding.cannotlink_penalty);
  %fprintf('options.embedding.use_unsupervised: %s\n', tf{options.embedding.use_unsupervised + 1});
  %fprintf('options.walk.use_go_link: %s\n', tf{options.walk.use_go_link + 1});
  %fprintf('options.walk.go_link_fraction: %f\n', options.walk.go_link_fraction);
  fprintf('options.walk.restart_prob: %f\n', options.walk.restart_prob);
  fprintf('options.kfolds: %f\n', options.kfolds);
  fprintf('options.test_fraction: %f\n', options.test_fraction);
  fprintf('\n\n');

  success = true;
end

