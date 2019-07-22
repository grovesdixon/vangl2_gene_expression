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
i.res = results(dds, name="Ethanolethanol.Genotypemutant", independentFiltering=INDEPDENDENT_FILTERING)


#VOLCANO PLOTS
res = e.res;TITLE='Ethanol effect'
res=t.res;TITLE='Temperature effect'
res=m.res;TITLE='Genotype effect'
res=i.res;TITLE='Ethanol:Genotype interaction'


ev=volcano_from_deseq_res(e.res, 'Ethanol effect')
tv=volcano_from_deseq_res(t.res, 'Temperature effect')
gv=volcano_from_deseq_res(m.res, 'Genotype effect')
iv=volcano_from_deseq_res(i.res, 'Ethanol:Genotype interaction')
legend <- cowplot::get_legend(ev+theme(legend.position='top'))
lplt = plot_grid(legend)

pltList=list(ev,tv,gv,iv)
modVolcs = function(x){
  x= x + theme(
    axis.title.x=element_blank(),
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    # axis.text.y=element_blank(),
    # axis.ticks.y=element_blank(),
    legend.position='none')
  return(x)
}
mList = lapply(pltList, function(x) modVolcs(x))


#build multipanel
ylab = ggdraw() + draw_label(bquote(-log[10]~pvalue), angle=90)
plts=plot_grid(plotlist = mList)
ylabPlts = plot_grid(ylab, plts, nrow=1, rel_widths=c(0.05, 1))
top = plot_grid(lplt, ylabPlts, nrow=2, rel_heights=c(0.1, 1))
xlab = plot_grid(ggdraw() + draw_label(bquote(log[2]~fold~difference)))
plot_grid(top, xlab, nrow=2, rel_heights=c(1, 0.05))




res.eth = e.res
save(dds, res.eth, file='deseq/ethanol_results.Rdata')
save(e.res, t.res, m.res, i.res, file='deseq/all_res.Rdata')
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




