#look_at_genes_within_GO_terms.R
setwd('./wgcna/go_mwu')

#------------ LOOK AT GENES WITHIN GO TERMS ------------#
library(cowplot)
library(tidyverse)

#upload the results for each genes output by the functions above
#re-pick input if you need to
input = 'brown_moduleInput.csv'
goDivision = 'BP'

#upload GO enrichment results
resName = paste(paste('MWU', goDivision, sep = "_"), input, sep = "_")
go.res=read.table(resName, header = T)
sig = go.res[go.res$p.adj < 0.1,]
head(sig, n=20)


#upload the gene-GO associations
geneResName=paste(goDivision, input, sep='_')
gene.res=read.table(geneResName, header=T)
head(gene.res)


#search the GO results for a particular term
searchString = 'development'
sig[grep(searchString, sig$name),] 


#select the GO term from your results file (sigGo above)
go="GO:0004984" #olfactory receptor activity
go='GO:0004930' #G-protein coupled receptor activity
go='GO:0007409;GO:0048812' #neuron projection morphogenesis
go='GO:0043009;GO:0009792;GO:0009790' #embryo development

#subset for that GO term
go.genes = gene.res[gene.res$term == go, 'seq']
length(go.genes)

#get gene names
geneSet = as.character(go.genes)

#gather the names from ensemble using biomarRt
library("biomaRt")
embl = useMart("ensembl", dataset="drerio_gene_ensembl")
names = getBM(attributes = c('ensembl_gene_id', 'external_gene_name','description'), filters = c('ensembl_gene_id'), values = geneSet, mart = embl)
head(names)

#named genes are here:
names

#merge with deseq results
ll=load('../../deseq/ethanol_results.Rdata')
mdat = res.eth %>% 
  data.frame() %>% 
  mutate(ensembl_gene_id = rownames(res.eth)) %>% 
  right_join(names, by = 'ensembl_gene_id')
head(mdat)
geneSet = mdat$ensembl_gene_id

#check the log2 values match expectation from heatmap
mdat %>% 
  ggplot(aes(y=log2FoldChange)) +
  geom_boxplot() +labs(title='Brown module')



mdat %>% 
  write_tsv(path='~/Desktop/brownModule_BP_embryDevelopment_genes.tsv')

