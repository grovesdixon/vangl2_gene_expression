#genes_for_clueICV.R
library("biomaRt")
hsembl = useMart("ensembl", dataset="hsapiens_gene_ensembl") #get the human dataset
setwd("/Users/grovesdixon/lab_files/projects/vangle_zebrafish")
source('~/gitreps/zebrafish_early_ethanol_RNASeq/scripts/zebrafish_RNAseq_functions.R')



#read in ethanol DESEq results
lnames=load('ethanol_results.Rdata')
head(res.eth)
res.eth=res.eth[!is.na(res.eth$padj),]
sdat=data.frame(res.eth[res.eth$padj<0.1,])
dim(sdat)


#upregulated by ethanol
up=sdat[sdat$log2FoldChange > 0,]
sub=merge_gene_names(up, sort.column='pvalue')
geneSet=sub$external_gene_name
geneSet=sub(" (1 of many)", "", geneSet, fixed=T)
h.upGeneSet=get_human_ids(geneSet)



#downregulated by ethanol
down=sdat[sdat$log2FoldChange < 0,]
sub=merge_gene_names(down, sort.column='pvalue')
geneSet=sub$external_gene_name
geneSet=sub(" (1 of many)", "", geneSet, fixed=T)
h.downGeneSet=get_human_ids(geneSet)

#check no genes that show up in both
h.upGeneSet[h.upGeneSet %in% h.downGeneSet]

#write them out
write.table(h.upGeneSet, "ethanol_significant_upregulated.txt", quote=F, row.names=F, col.names=F)
write.table(h.downGeneSet, "ethanol_significant_downregulated.txt", quote=F, row.names=F, col.names=F)







#gather genes based on GO terms
godat = read.table("datasets/zebrafish_embl_go.tsv", header = T)
godat = read.table("goMWU/BP_ethanol_goMWU_input.csv", header = T) #output from GO_MWU_ethanol.R--connects genes to merged GO terms
godat = read.table("goMWU/MF_ethanol_goMWU_input.csv", header = T) #output from GO_MWU_ethanol.R--connects genes to merged GO terms

head(godat)
go='GO:0009653' #anatomical structure morphogenesis
go='GO:0048598' #embryonic morphogenesis
go='GO:0048513' #organ development (animal organ development)
go='GO:0009887' #organ morphogenesis (animal organ morphogenesis)

termSet=c()



sub = godat[grep('organ morphogenesis', godat$name),]
sub = godat[grep(go, godat$term),]
head(sub)
nrow(sub)

sgo = sdat[rownames(sdat) %in% sub$seq,]
head(sgo)
nrow(sgo)

















