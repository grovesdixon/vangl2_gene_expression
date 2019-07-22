

library(cowplot)
library(tidyverse)


ll=load("initialize_counts/rld.Rdata")
ll=load("deseq/deseqInput.Rdata")
source('rnaseq_functions.R')

coldata = coldata %>% 
  mutate(Temperature = as.character(Temperature),
         Temperature = if_else(Temperature=='roomTemp',
                               'room temp.',
                               Temperature),
         Genotype = as.character(Genotype),
         Genotype = if_else(Genotype=='het_wt',
                            'wt/het',
                            Genotype))

#plot pca
NTOP=10000
epca <- mod.plotPCA.df(rld.df, coldata, intgroup = 'Ethanol', returnData=F, ntop=NTOP) + labs(colour='', subtitle='Ethanol treatment')
tpca <- mod.plotPCA.df(rld.df, coldata, intgroup = 'Temperature', returnData=F, ntop=NTOP) + labs(colour='', subtitle='Temperature treatment')
gpca <- mod.plotPCA.df(rld.df, coldata, intgroup = 'Genotype', returnData=F, ntop=NTOP) + labs(colour='', subtitle='Genotype')
bpca <- mod.plotPCA.df(rld.df, coldata, intgroup = 'Batch', returnData=F, ntop=NTOP) + labs(colour='', subtitle='Batch')
plot_grid(epca, tpca, gpca, bpca, nrow=2)
