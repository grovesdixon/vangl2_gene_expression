#mut_vs_wthet.R
#look at the ethanol effect in each of the genotype groups individually


library(tidyverse)
library(cowplot)

#load full dataset
ll=load('deseq/ethanol_results.Rdata')
fdat = data.frame(res.eth)

#load wt/het data
ll=load('deseq/wtHetOnly_ethanol_results.Rdata')
wdat = data.frame(res.eth)

#load mutant data
ll=load('deseq/mutantOnly_ethanol_results.Rdata')
mdat = data.frame(res.eth)


#merge and get difference
cdat = merge(wdat, mdat, by = 0) %>% 
  mutate(fcw = log2FoldChange.x,
         fcm = log2FoldChange.y,
         diff = fcm - fcw)
head(cdat)
cdat %>% 
  ggplot(aes(x=diff)) +
  geom_density()
cutoff = 1
cdat = cdat %>% 
  mutate(different = abs(diff)>=cutoff)
sum(cdat$different)


#plot wt vs mutant
cdat %>% 
  ggplot(aes(x=log2FoldChange.x, y=log2FoldChange.y, color=different)) +
  geom_hline(yintercept = 0, lty=2, size=0.5) +
  geom_vline(xintercept = 0, lty=2, size=0.5) +
  geom_point() +
  geom_abline(slope=1, intercept=0,
              na.rm = FALSE, show.legend = NA) +
  labs(x='wt/het log2 diff.', y='mutant log2 diff.') +
  scale_color_manual(values=c('black', 'red'))


#get the gene names
source('rnaseq_functions.R')
rownames(cdat)=cdat$Row.names
ndat = merge_gene_names(cdat)
head(ndat)


#write out
ndat %>% 
  write_tsv(path='deseq/mutWt.tsv')







  