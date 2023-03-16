%% gotype: 'bp', 'mf', or 'cc'
%% genes: cell array of query gene names
%% ontsize (optional): [min_size max_size], only return terms within a size range
%% no_overlap (optional): remove redundant terms
%%
function [anno, golevels, anno_full, filt] = load_go_full(org, gotype, genes, ontsize, no_overlap)
  prefix = sprintf('go_%s_ref', org);
  
  go_path = sprintf('data/annotations/%s', org);
  gogene = textread(sprintf('%s/%s_genes.txt', go_path, prefix), '%s');
  filt = ismember(genes, gogene);
  
  goterm = textread(sprintf('%s/%s_%s_terms.txt', go_path, prefix, gotype), '%s');
  [g, t] = textread(sprintf('%s/%s_%s_adjacency.txt', go_path, prefix, gotype), '%d%d');
  goanno = sparse(t, g, true, length(goterm), length(gogene));
  
  anno = zeros(length(goterm), length(genes));
  
  genemap = containers.Map(gogene, 1:length(gogene));
  s2goind = cell2mat(values(genemap, genes(filt)));
  anno(:,filt) = goanno(:,s2goind);
  
  %% Use ontology graph to propagation annotations
  termfile = sprintf('%s/graph/go_%s.terms', go_path, gotype);
  nterm = length(textread(termfile, '%s'));
  
  % each entry is the go level of the term with same index in termfile.
  golevels = textread(sprintf('%s/graph/go_%s.levels', go_path, gotype), '%d');
  
  mapfile = sprintf('%s/graph/go_%s.map', go_path, gotype);
  [t, i] = textread(mapfile, '%s\t%d');
  m = containers.Map(t, i);
  
  linkfile = sprintf('%s/graph/go_%s.links', go_path, gotype);
  M = dlmread(linkfile);
  ontgraph = sparse(M(:,1), M(:,2), true, nterm, nterm);
  reachable = get_reachable(ontgraph);
  
  
  goind = cell2mat(values(m, goterm(isKey(m, goterm))));
  anno = (reachable(:,goind) * anno(isKey(m, goterm), :)) > 0; % propagate
  anno_full = anno;

  % filter out any unused go terms.
  filt1 = sum(anno, 2) > 0;
  anno = anno(filt1,:);
  golevels = golevels(filt1);
  
%   if exist('no_overlap', 'var') && no_overlap
%     [anno_full, ~]  = filter_anno(anno_full, golevels, 0.1);
%   end
  
  
%% Select terms in the given size range
  if exist('ontsize', 'var')
    term_size = sum(anno, 2);
    filt2 = ontsize(1) <= term_size & term_size <= ontsize(2);
    anno = anno(filt2,:);
    golevels = golevels(filt2);
  end
  %% Remove redundant terms
  if exist('no_overlap', 'var') && no_overlap
   [anno, golevels, filt3]  = filter_anno(anno, golevels, 0.1);
  end 
 
  filt = {filt1, filt2, filt3};
  %filt = {filt1, filt2};
end
