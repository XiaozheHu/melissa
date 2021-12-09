% load network data
%
% [Input]
% weight: 0 represents the all the edges have same weight. 1 represents all
% the deges have different weight
% nnode: number of nodes
% bp: back propogation probability
%
% [Output]
% A : the sparse matrix of network.
%
function [A] = load_network(net_prefix,weight,nnode,bp)


if weight
    [g1,g2,s] = textread(net_prefix, '%d%d%f');
    n = nnode;
    A = sparse(g1,g2,s,n,n);
else
   [g1,g2] = textread(net_prefix, '%d%d');
    n = max( max(g1),max(g2));
    A = sparse(g1,g2,true,n,n);
end

is_sym = isequal(A,A');

if ~is_sym
  A = A + bp*A';
end

for i=1:n
  A(i,i) = 0;
end

  fprintf('Size of network: %d\n', n);

end
