function [genes, ngene, anno, levels] = load_anno_wo_filter(options)
    addpath code/load_data
    
    org = options.org;
    onttype = options.onttype;
    ontsize = options.ontsize;

    gene_file = sprintf('data/networks/%s/%s_string_genes.txt', org, org);
    genes = textread(gene_file, '%s');
    ngene = length(genes);
    
    if contains(onttype, 'level') % check if using go or mips
        anno = load_mips(onttype, genes);
        levels = 0;
    else
        % [anno, levels] = load_go(org, onttype, genes, ontsize, true); 
        [anno, levels] = load_go(org, onttype, genes, ontsize, false); 
    end
end