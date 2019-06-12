#deseq_treatment_effects.R 

#SET UP THE DATA TO RUN DESEQ
library(DESeq2)
library(cowplot)
library(tidyverse)

#SAVE YOUR DIRECTORY NAME AND SET WD
directory<-"/Users/grovesdixon/lab_files/projects/vangle_zebrafish"


setwd(directory)
source("rnaseq_functions.R")


#load the data
ll = load('deseqInput.Rdata')
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
	                   Ethanol:Genotype)
	)

INDEPDENDENT_FILTERING = FALSE
dds <- DESeq(dds)
resultsNames(dds)
e.res = results(dds, contrast = c('Ethanol', 'ethanol', 'control'), independentFiltering=INDEPDENDENT_FILTERING)
t.res = results(dds, contrast = c('Temperature', 'roomTemp', 'incubator'), independentFiltering=INDEPDENDENT_FILTERING)
m.res = results(dds, contrast = c('Genotype', 'mutatnt', 'het_wt'), independentFiltering=INDEPDENDENT_FILTERING)
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
save(dds, res.eth, file='ethanol_results.Rdata')
write.table(e.res, file='ethanol_results.tsv', quote=F, sep='\t')
write.table(t.res, file='temperature_results.tsv', quote=F, sep='\t')
write.table(m.res, file='genotype_results.tsv', quote=F, sep='\t')





