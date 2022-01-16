function success = log_options(options)
  %% Logs the options so we can see the parameters used in log file later
  tf = {'false','true'};
  fprintf('\n\n\n---------------PARAMETERS--------------\n');
  if isfield(options, 'org') 
      fprintf('options.org: %s\n', options.org);
  end
  if isfield(options, 'onttype') 
      fprintf('options.onttype: %s\n', options.onttype);
  end
  if isfield(options, 'ontsize') 
      fprintf('options.ontsize: [%d, %d]\n', options.ontsize(1), options.ontsize(2));
  end
  if isfield(options, 'kNN_neighbor')
      fprintf('options.kNN_neighbor: %d\n', options.kNN_neighbor);
  end
  if isfield(options, 'num_clusters')
      fprintf('options.num_clusters: %d\n', options.num_clusters);
  end
  
  if isfield(options, 'SSDR_embedding')
      fprintf('options.SSDR_embedding: %d\n', options.SSDR_embedding);
  end
  if isfield(options, 'extra_graph_filt')
      fprintf('options.extra_graph_filt: %d\n', options.extra_graph_filt);
  end
  if isfield(options, 'kfolds')
      fprintf('options.kfolds: %d\n', options.kfolds);
  end
  if isfield(options, 'test_fraction')
      fprintf('options.test_fraction: %f\n', options.test_fraction);
  end
  if isfield(options, 'embedding')
      if isfield(options.embedding, 'ndim')
          fprintf('options.embedding.ndim: %d\n', options.embedding.ndim);
      end
      if isfield(options.embedding, 'mustlink_penalty')
          fprintf('options.embedding.mustlink_penalty: %d\n', options.embedding.mustlink_penalty);
      end
      if isfield(options.embedding, 'cannotlink_penalty')
          fprintf('options.embedding.cannotlink_penalty: %d\n', options.embedding.cannotlink_penalty);
      end
      %fprintf('options.embedding.svd_approx: %s\n', tf{options.embedding.svd_approx + 1});
      %fprintf('options.embedding.use_unsupervised: %s\n', tf{options.embedding.use_unsupervised + 1});
  end
  if isfield(options, 'walk')
      %fprintf('options.walk.use_go_link: %s\n', tf{options.walk.use_go_link + 1});
      %fprintf('options.walk.go_link_fraction: %f\n', options.walk.go_link_fraction);
      if isfield(options.walk, 'restart_prob')
          fprintf('options.walk.restart_prob: %f\n', options.walk.restart_prob);
      end
  end
  if isfield(options, 'svm')
      if isfield(options.svm, 'gmax')
          fprintf('options.svm.gmax: %f\n', options.svm.gmax);
      end
      if isfield(options.svm, 'cmax')
          fprintf('options.svm.cmax: %f\n', options.svm.gmax);
      end 
  end
  
  fprintf('\n\n');

  success = true;
end

