function link_mat = select_p_percent(link_mat, p, rng_seed)
[row,col,nn] = find(link_mat);
n = length(row);
nsample = floor(n*p);
rng = rng_seed;
ind = randsample(n,nsample);
ind = sort(ind);
n1 = size(link_mat,1);
link_mat = sparse(row(ind), col(ind), nn(ind), n1,n1);
end