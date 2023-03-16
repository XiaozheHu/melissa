function walks = semisupervised_rwr(network_files, ngene, restart_prob,...
      teleport_prob, labels)
  nfiles = length(network_files);
  walks = zeros(nfiles, ngene, ngene);
  for i=1:nfiles
    A = load_network(network_files{i}, ngene);
    P = markov_mat(A);
    P = semisupervise_markov(P, labels, teleport_prob);
    W = rwr(P, restart_prob);
    walks(i,:,:) = W;
  end
  clear A P
end
  

function Ps = semisupervise_markov(P, labels, teleport_prob)
  Ps = P;
  n = size(P,1);
  nlabels = size(labels,1);

  % indicator vector for vertices in the training set
  train_filt = (sum(labels) > 0).';

  % number of classes for each vertex
  vertex_label_counts = sum(labels, 1);
  % number of vertices in each class
  class_label_counts = sum(labels, 2);

  % start by cutting out alpha percent of the labeled transitions
  for i = 1:n
    if ~train_filt(i)
      % ignore data not in the training set
      continue
    end
    Ps(i,:) = (1 - teleport_prob) * Ps(i,:);
  end

  % iterate over each class and vertex and adds weak edges between vertices
  % that share a label. The edge weight from u to v is defined to be
  % sum_{l labels u and v} teleport_prob / ((num labels of u)(num vertices with label l))
  % roughly, this uniformly picks a class to teleport using, then picks from the vertices
  % in that class randomly  
  for i = 1:nlabels
    % grabs just the labels for class i
    % and normalizes the vector to sum to alpha
    class_labels = labels(i,:) * teleport_prob / class_label_counts(i);
    for j = 1:n
      if class_labels(j) == 0
        % ignore the non-class members
        continue
      end
      vlc = vertex_label_counts(j);
      % adds the class label weighted by the 
      Ps(j,:) = Ps(j,:) + class_labels / vlc;
    end
  end
end