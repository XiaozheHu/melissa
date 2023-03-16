function [anno, golevels, filt] = filter_anno(anno, golevels, thres) % add extra argument and return.
  anno = double(anno > 0);
  termsize = sum(anno, 2);
  
  [~, ord] = sort(termsize, 'ascend');
  anno = anno(ord,:);
  golevels = golevels(ord);
  % apply the same ordering

  in = anno * anno';
  un = size(anno, 2) - (1 - anno) * (1 - anno)';
  jacc = in ./ un;

  max_jacc = max(triu(jacc, 1), [], 2);
  anno = anno(max_jacc <= thres,:);
  filt = max_jacc <= thres;
  golevels = golevels(max_jacc <= thres);
end
