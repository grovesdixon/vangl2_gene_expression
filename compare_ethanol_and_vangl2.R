#compare_ethanol_and_vangle2.R

setwd('~/gitreps/vangl2_gene_expression/')
source('rnaseq_functions.R')

edat=read.table('deseq/ethanol_results.tsv')
gdat=read.table('deseq/genotype_results.tsv')

mdat = merge(edat, gdat, by = 0)
mdat$ensemble=mdat$Row.names
mdat$Row.names<-NULL
mdat$sig.g = mdat$padj.y < 0.1
sdat = mdat %>% 
  filter(sig.g==TRUE)
rownames(sdat)=sdat$ensemble
sdat2 = merge_gene_names(sdat, sort.column='pvalue.y')


lm1=lm(mdat$log2FoldChange.x ~mdat$log2FoldChange.y)
r2 = round(summary(lm1)$r.squared, digits=2)
p=
r=round(cor(x=mdat$log2FoldChange.x, y=mdat$log2FoldChange.y), digits=2)

mdat %>% 
  ggplot(aes(x=log2FoldChange.x, y=log2FoldChange.y)) +
  geom_hline(yintercept = 0, lty=2, lwd=0.5) +
  geom_vline(xintercept = 0, lty=2, lwd=0.5) +
  geom_point(color='black') +
  geom_smooth(method='lm', se=FALSE) +
  geom_point(data=sdat2, aes(x=log2FoldChange.x, y=log2FoldChange.y), fill='red', color='black', pch=21, size=3) +
  geom_text(data=sdat2, aes(label=external_gene_name),hjust=0, vjust=0, color='red') +
  labs(x=bquote(ethanol~log[2]~foldChange), y=bquote(vangl[2]~log[2]~foldChange)) +
  annotate("text", x = -2.5, y = 2,
           label = paste('italic(r) ==', r), parse=TRUE, color='black',
           hjust=0)




mean(abs(mdat$log2FoldChange.x))
mean(abs(mdat$log2FoldChange.y))
