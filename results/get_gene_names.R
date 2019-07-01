#get_gene_names.R

library(tidyverse)
source('rnaseq_functions.R')

#read in dataframe with ensemble IDs to get gene names for
inputFile = 'deseq/genotype_results.tsv'
dat = read.table(inputFile)

#merge in gene names
namedDat = merge_gene_names(dat, col=F, sort.column='pvalue')
head(namedDat)


#write out the named dataframe
write.table(namedDat, file='~/Desktop/named_genotype_results.tsv', sep='\t', row.names=F, quote=F)
