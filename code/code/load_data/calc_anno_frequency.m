function [least_anno, most_anno]=calc_anno_frequency(anno)
%-----------------------------------------------------------------
% calculate frequency of function labels and assign multilabeled 
% genes with the labels that annotate the most/least genes.
% E.G. if a gene X has label A and B, A labels more genes than B does.
%       in least_anno X has label B and in most_anno X has label B.
%  Input:   anno: nlabel x ngenes matrix 
%  Output:  least_anno; most_anno
%------------------------------------------------------------------
ngene = size(anno,2);
least_anno = zeros(ngene,1);
most_anno = zeros(ngene,1);
for i = 1:ngene
    a = anno(:,i);
    if a(a~=0)~=0
       nnz_ind = find(a);
       least_anno(i) = nnz_ind(1);
       most_anno(i)=nnz_ind(end);
    else
        least_anno(i) = 0;
        most_anno(i) = 0;
    end
end
end