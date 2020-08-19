#study_comparison.R
source('rnaseq_functions.R')
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

#function to prepare deseq results file for plotting
prep_res = function(res){
  x=data.frame(res)
  x$gene=rownames(res)
  return(x)
}


#LOAD ETHANOL RESULTS FOR THIS STUDY (vangl2 rnaseq)
ll=load('deseq/ethanol_results.Rdata')
ll
new=prep_res(res.eth)


#LOAD RESULTS FROM PREVIOUS STUDY (wildtype rnaseq)
#for full dataset (all timepoints together)
ll=load('~/gitreps/Drerio_early_ethanol_RNAseq/deseq/ethanol_full_LRT_results.Rdata')
ll
#8 hour
ll=load('~/gitreps/Drerio_early_ethanol_RNAseq/deseq/ethanol_8hr_LRT_results.Rdata')
ll
#10 hour
ll=load('~/gitreps/Drerio_early_ethanol_RNAseq/deseq/ethanol_10hr_LRT_results.Rdata')
ll
#14 hour
ll=load('~/gitreps/Drerio_early_ethanol_RNAseq/deseq/ethanol_14hr_LRT_results.Rdata')
ll



#organize subsets of previous study and labels into lists
oldRes = list(res.eth, res.e, res.t, res.f)
oldf = lapply(oldRes, function(x) prep_res(x) %>% left_join(new, by = 'gene'))
names = c('Full', '8hr', '10hr', '14hr')
ylabs = c('Old Full Ethanol',
          'Old 8hr Ethanol',
          'Old 10hr Ethanol',
          'Old 14 hr Ethanol')
xlab = 'New Full Ethanol'

#function to build scatterplot
do_scatter = function(x, xlab, ylab, main){
  lm1=lm(x$log2FoldChange.y~x$log2FoldChange.x)
  print(summary(lm1))
  r2=round(summary(lm1)$r.squared, digits=2)
  x %>% 
    ggplot(aes(x=log2FoldChange.y, y=log2FoldChange.x)) +
    geom_point() +
    geom_smooth(method='lm') +
    annotate("text", x = 1.5, y = -2,
             label = paste('italic(R) ^ 2 ==', r2), parse=TRUE, color='blue') +
    labs(x=xlab, y=ylab, subtitle=main)
    
}


#plot them
pltList = list()
for (i in 1:length(oldf)){
  main=names[i]
  df=oldf[[i]]
  plt=do_scatter(df, xlab='New', ylab='Old', main=main) + lims(x=c(-3,3),y=c(-3,3))
  pltList[[names[i]]]=plt
}

plot_grid(plotlist = pltList)



#PLOT OLD 10 HR VS NEW 10HR

head(res.t)
head(new)

#make calls for significance accross datasets
CUT=0.1
x=prep_res(res.t) %>% 
  data.frame() %>% 
  left_join(data.frame(new), by = 'gene') %>% 
  mutate(significant = if_else(padj.x<CUT,
                               'old',
                               'none'),
         significant = if_else(padj.y<CUT & padj.x < CUT,
                               'both',
                               significant),
         significant = if_else(padj.y<CUT & padj.x >=CUT,
                               'new',
                               significant),
         significant = factor(significant, levels=c('none', 'old', 'new', 'both')))

#pull just the significant ones for getting names
sx = x %>% 
  filter(significant!='none') %>% 
  rename
sx2 = merge_gene_names(sx, col = 'gene')
head(sx2)
sx2 %>% 
  write_csv(path='~/Desktop/oldVnew10hrEthanol.csv')

#build plot with signiciance combinations color coded
x %>% 
  filter(!is.na(significant)) %>% 
  ggplot(aes(x=log2FoldChange.x, y=log2FoldChange.y,
             fill=significant)) +
  geom_point(alpha=0.5,
             color='black',
             pch=21) +
  scale_fill_manual(values=c('black', 'dodgerblue', 'forestgreen', 'firebrick')) +
  labs(x='old', y='new', subtitle='10 hpf') +
  geom_point(data=sx2[sx2$significant %in% c('old'),],
             size=2,
             fill='dodgerblue',
             color='black',
             pch=21) +
  geom_point(data=sx2[sx2$significant %in% c('new'),],
             size=2,
             fill='forestgreen',
             color='black',
             pch=21) +
  geom_point(data=sx2[sx2$significant %in% c('both'),],
             size=3,
             fill='firebrick',
             color='black',
             pch=21)


#build plot highlighting genes signficant in both
x %>% 
  filter(!is.na(significant)) %>% 
  mutate(`Both significant`=if_else(significant=='both',
                                    TRUE,
                                    FALSE)) %>% 
  ggplot(aes(x=log2FoldChange.x, y=log2FoldChange.y,
             fill=`Both significant`)) +
  geom_point(alpha=0.5,
             pch=21) +
  scale_fill_manual(values=c('black', 'red')) +
  labs(x='old', y='new', subtitle='10 hpf') +
  geom_point(data=sx2[sx2$significant %in% c('both'),],
             size=2,
             fill='firebrick',
             color='black',
             pch=21) +
  labs(fill='Both significant')



# compare with WGCNA ------------------------------------------------------

#isolate significance calls
sdat = x %>% 
  dplyr::select(gene, significant)

#upload wgcna module membership values and merge with significance calls
wdat = read.csv('wgcna/geneModuleMembership.csv') %>% 
  rename(gene = X) %>% 
  dplyr::select(gene, assignment) %>% 
  tibble() %>% 
  inner_join(sdat, by = 'gene')

#plot amounts splits by module
wdat %>% 
  filter(!is.na(significant)) %>% 
  ggplot(aes(x=assignment, fill=significant)) +
  geom_bar()


#plot proportion 'both'
pdat = wdat %>% 
  filter(!is.na(significant)) %>% 
  group_by(assignment) %>% 
  summarize(nboth = sum(significant == 'both'),
            N = n()) %>% 
  mutate(prop = nboth / N) %>% 
  arrange(prop) 
lvl_cols = unique(pdat$assignment)
pdat$assignment = factor(pdat$assignment, levels = lvl_cols)
pdat %>% 
  ggplot(aes(x=assignment, y=prop, fill=assignment)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = as.character(lvl_cols)) +
  labs(y='proportion of module',
       x='module',
       fill='module') +
  coord_flip()

