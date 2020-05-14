#vangl2_genes.R
source('rnaseq_functions.R')

#check how genes differentially expressed by vangl2 genotype respond to other treatments



# UPLOAD DATASETS -------------------------------------------------

#vangl2 response
vdat = read.table('deseq/genotype_results.tsv')
vdat$gene = rownames(vdat)
sig = vdat %>% 
  filter(padj < 0.1)
rownames(sig)=sig$gene
sig = merge_gene_names(sig)
sig_genes = sig$gene


read_compare = function(filePath){
  cdat = read.table(filePath)
  cdat$gene = rownames(cdat)
  return(cdat)
}


#ethanol
edat = read_compare('deseq/ethanol_results.tsv')

#temperature
tdat = read_compare('deseq/temperature_results.tsv')


#function to plot agreement between two deseq results sets
plot_agreement = function(cdat, sig, xlab, ylab, subtitle){
  mdat = cdat %>% 
    dplyr::rename(log2Compare = log2FoldChange) %>% 
    filter(gene %in% sig_genes) %>% 
    left_join(sig, by = 'gene')
  
  mdat %>% 
    mutate(vpos = log2FoldChange > 0,
           cpos = log2Compare > 0,
           agree = vpos==cpos) %>% 
    ggplot(aes(x=log2FoldChange, y=log2Compare, color=agree)) +
    geom_point(size=3) +
    geom_hline(yintercept = 0, lty=2, color='grey') +
    geom_vline(xintercept = 0, lty=2, color='grey') +
    scale_color_manual(values=c('red', 'black')) +
    labs(x=xlab,
         y=ylab,
         subtitle=subtitle)
  
}

#agreement for ethanol
eplt = plot_agreement(edat,
               sig,
               xlab='vangl2 genotype response',
               ylab='ethanol treatment response',
               subtitle='')

#agreement for temp
tplt = plot_agreement(tdat,
               sig,
               xlab='vangl2 genotype response',
               ylab='temperature treatment response',
               subtitle='')

plot_grid(eplt, tplt)


# BUILD A LARGE MERGED TABLE ----------------------------------------------

sub_sig = sig %>% 
  dplyr::select(log2FoldChange,
                padj) %>% 
  rename(vangl2.log2FoldChange = log2FoldChange,
         vangl2.padj = padj) %>% 
  mutate(gene = rownames(sig))

#format the other two
head(edat)
sub_edat = edat %>% 
  dplyr::select(log2FoldChange,
                padj) %>% 
  rename(etoh.log2FoldChange = log2FoldChange,
         etoh.padj = padj) %>% 
  mutate(gene=rownames(edat)) %>% 
  filter(gene %in% sig_genes)

sub_tdat = tdat %>% 
  dplyr::select(log2FoldChange,
                padj) %>% 
  rename(temp.log2FoldChange = log2FoldChange,
         temp.padj = padj) %>% 
  mutate(gene = rownames(tdat)) %>% 
  filter(gene %in% sig_genes)

#upload the wgcna results
wdat = read.csv('wgcna/geneModuleMembership.csv') %>% 
  rename(gene=X) %>% 
  filter(gene %in% sig_genes)


#merg up
to_merge = list(sub_sig, sub_edat, sub_tdat, wdat)
fmdat = purrr::reduce(to_merge, left_join, by = 'gene')
fmdat = data.frame(fmdat)
rownames(fmdat) = fmdat$gene
fmdat2 = merge_gene_names(fmdat)
write_csv(fmdat2, path='deseq/merged_sig_genotype.csv')




