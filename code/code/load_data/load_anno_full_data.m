function [genes, ngene, anno, levels, full_anno, ontgraph, used_ind] = load_anno(options)
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
        [anno, levels, full_anno, ontgraph, used_ind] = load_go_full_data(org, onttype, genes, ontsize, true); 
    end
end

