#deseq_treatment_effects.R 

#SET UP THE DATA TO RUN DESEQ
library(DESeq2)
library(cowplot)
library(tidyverse)
source("rnaseq_functions.R")


#load the data
ll = load('deseq/deseqInput.Rdata')
ll
coldata
head(counts)
sum(colnames(counts) == coldata$sample)==ncol(counts)


#get variance stabilized counts for timepoint 2
dds<-DESeqDataSetFromMatrix(counts,
	colData = coldata, 
	design = formula(~Ethanol +
	                   Temperature +
	                   Genotype +
	                   Temperature:Genotype +
	                   Ethanol:Genotype)
	)

INDEPDENDENT_FILTERING = FALSE
dds <- DESeq(dds)
resultsNames(dds)
e.res = results(dds, contrast = c('Ethanol', 'ethanol', 'control'), independentFiltering=INDEPDENDENT_FILTERING)
t.res = results(dds, contrast = c('Temperature', 'roomTemp', 'incubator'), independentFiltering=INDEPDENDENT_FILTERING)
m.res = results(dds, contrast = c('Genotype', 'mutant', 'het_wt'), independentFiltering=INDEPDENDENT_FILTERING)
i.res = results(dds, name='Ethanolethanol.Genotypemutatnt', independentFiltering=INDEPDENDENT_FILTERING)


#VOLCANO PLOTS
res = e.res;TITLE='Ethanol effect'
res=t.res;TITLE='Temperature effect'
res=m.res;TITLE='Genotype effect'
res=i.res;TITLE='Ethanol:Genotype interaction'

sdf = data.frame(res) %>% 
  arrange(pvalue) %>% 
  mutate(sig = !is.na(padj) & padj < 0.1) %>% 
  as_tibble()
sdf


g=ggplot(data=sdf, aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
	geom_point(alpha=0.1) +
	scale_color_manual(values=c('black', 'red')) + 
	labs(subtitle=TITLE,
	     y=bquote(-log[10]~pvalue),
	     x=bquote(log[2]~fold~difference)) +
  lims(y=c(0,35))
plot(g)




res.eth = e.res
save(dds, res.eth, file='deseq/ethanol_results.Rdata')
write.table(e.res, file='deseq/ethanol_results.tsv', quote=F, sep='\t')
write.table(t.res, file='deseq/temperature_results.tsv', quote=F, sep='\t')
write.table(m.res, file='deseq/genotype_results.tsv', quote=F, sep='\t')


#WRITE OUT FOR GO ENRICHMENT
write_out_go(e.res, './go_mwu/ethanolForGO.csv')
write_out_go(t.res, './go_mwu/tempForGO.csv')




# repeat for mutants only -------------------------------------------------

mcoldata = coldata %>% 
  filter(Genotype=='mutant')
mcounts = counts[, mcoldata$sample]
dim(mcoldata)
dim(mcounts)


#get variance stabilized counts for timepoint 2
dds<-DESeqDataSetFromMatrix(mcounts,
                            colData = mcoldata, 
                            design = formula(~Ethanol +
                                               Temperature)
)

INDEPDENDENT_FILTERING = FALSE
dds <- DESeq(dds)
resultsNames(dds)
e.res = results(dds, contrast = c('Ethanol', 'ethanol', 'control'), independentFiltering=INDEPDENDENT_FILTERING)
t.res = results(dds, contrast = c('Temperature', 'roomTemp', 'incubator'), independentFiltering=INDEPDENDENT_FILTERING)


#VOLCANO PLOTS
res = e.res;TITLE='Ethanol effect'
res=t.res;TITLE='Temperature effect'


sdf = data.frame(res) %>% 
  arrange(pvalue) %>% 
  mutate(sig = !is.na(padj) & padj < 0.1) %>% 
  as_tibble()
sdf


g=ggplot(data=sdf, aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
  geom_point(alpha=0.1) +
  scale_color_manual(values=c('black', 'red')) + 
  labs(subtitle=TITLE,
       y=bquote(-log[10]~pvalue),
       x=bquote(log[2]~fold~difference)) +
  lims(y=c(0,35))
plot(g)




res.eth = e.res
save(dds, res.eth, file='deseq/mutantOnly_ethanol_results.Rdata')
write.table(e.res, file='deseq/mutantOnly_ethanol_results.tsv', quote=F, sep='\t')
write.table(t.res, file='deseq/mutantOnly_temperature_results.tsv', quote=F, sep='\t')




# repeat for wild type/het only -------------------------------------------


mcoldata = coldata %>% 
  filter(Genotype=='het_wt')
mcounts = counts[, mcoldata$sample]
dim(mcoldata)
dim(mcounts)


#get variance stabilized counts for timepoint 2
dds<-DESeqDataSetFromMatrix(mcounts,
                            colData = mcoldata, 
                            design = formula(~Ethanol +
                                               Temperature)
)

INDEPDENDENT_FILTERING = FALSE
dds <- DESeq(dds)
resultsNames(dds)
e.res = results(dds, contrast = c('Ethanol', 'ethanol', 'control'), independentFiltering=INDEPDENDENT_FILTERING)
t.res = results(dds, contrast = c('Temperature', 'roomTemp', 'incubator'), independentFiltering=INDEPDENDENT_FILTERING)


#VOLCANO PLOTS
res = e.res;TITLE='Ethanol effect'
res=t.res;TITLE='Temperature effect'


sdf = data.frame(res) %>% 
  arrange(pvalue) %>% 
  mutate(sig = !is.na(padj) & padj < 0.1) %>% 
  as_tibble()
sdf


g=ggplot(data=sdf, aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
  geom_point(alpha=0.1) +
  scale_color_manual(values=c('black', 'red')) + 
  labs(subtitle=TITLE,
       y=bquote(-log[10]~pvalue),
       x=bquote(log[2]~fold~difference)) +
  lims(y=c(0,35))
plot(g)




res.eth = e.res
save(dds, res.eth, file='deseq/wtHetOnly_ethanol_results.Rdata')
write.table(e.res, file='deseq/wtHetOnly_ethanol_results.tsv', quote=F, sep='\t')
write.table(t.res, file='deseq/wtHetOnly_temperature_results.tsv', quote=F, sep='\t')




