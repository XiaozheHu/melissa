function [dist_mat,knn]=compute_knn(X,k)
%-------------------------------------------------------
% compute the pairwise distances and find kNN
% Input:    X: data representation matrix, ndim x ngenes
%           k: first k nearest neighbors
% Output:   dist_mat: distance matrix, dist_mat(i,j) is the distance between 
%                     node i and j, ngenes x ngenes
%           knn: k nearest neighbor matrix, k x ngenes
%-------------------------------------------------------

dist_mat = squareform(pdist(X'));
 
n=size(X,2);
knn=zeros(k,n);

% find knn neighbors for each of the nodes
for i=1:n
    [dist,index]=sort(dist_mat(i,:),'descend');
    % use of find makes it so a node doesn't use itself as a neighbor.
    nn=index(find(dist,k,'last'));
    knn(:,i)=nn;
end
end