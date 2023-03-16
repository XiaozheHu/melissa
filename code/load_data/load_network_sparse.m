function A = load_network_sparse(filename, ngene)
    
    % load data
    M = dlmread(filename);
    
    % generate the adjacent matrix in sparse format
    A = sparse(M(:,1), M(:,2), M(:,3), ngene, ngene);
    
    % symmetrize A
    A = (A + A')/2;
    
    % add ones to the diagonal 
    % A = A + diag(sum(A, 2) == 0);
  end