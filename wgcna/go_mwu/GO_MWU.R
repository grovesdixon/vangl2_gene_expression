
# GO_MWU uses continuous measure of significance (such as fold-change or -log(p-value) ) to identify GO categories that are significantly enriches with either up- or down-regulated genes. The advantage - no need to impose arbitrary significance cutoff.

# If the measure is binary (0 or 1) the script will perform a typical "GO enrichment" analysis based Fisher's exact test: it will show GO categories over-represented among the genes that have 1 as their measure. 

# On the plot, different fonts are used to indicate significance and color indicates enrichment with either up (red) or down (blue) regulated genes. No colors are shown for binary measure analysis.

# The tree on the plot is hierarchical clustering of GO categories based on shared genes. Categories with no branch length between them are subsets of each other.

# The fraction next to GO category name indicates the fracton of "good" genes in it; "good" genes being the ones exceeding the arbitrary absValue cutoff (option in gomwuPlot). For Fisher's based test, specify absValue=0.5. This value does not affect statistics and is used for plotting only.

# Stretch the plot manually to match tree to text

# Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu

################################################################
# First, press command-D on mac or ctrl-shift-H in Rstudio and navigate to the directory containing scripts and input files. Then edit, mark and execute the following bits of code, one after another.
setwd("./wgcna/go_mwu")

#set variables to run for each module
ll=load('moduleInputFiles.Rdata')
ll
divisions = c('CC', 'MF', 'BP')
goAnnotations="zebrafish_embl_go.tsv" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
source("gomwu.functions.R")



# run go_mwu for each -----------------------------------------------------


for (goDivision in divisions){
  print('==============')
  print('==============')
  print('==============')
  print(goDivision)
  for (input in inputFiles){
    print('--------------')
    print(paste(input, '...', sep='.'))
    # Calculating stats. It might take ~3 min for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
    gomwuStats(input, goDatabase, goAnnotations, goDivision,
               perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
               largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
               smallest=5,   # a GO category should contain at least this many genes to be considered
               clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
               Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
               # Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
               #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
    )
  }
}




# record how many significant and save thsoe ---------------------------------------------
library(tidyverse)
divRec = c()
inRec = c()
sigRec = c()
for (goDivision in divisions){
  for (input in inputFiles){
    resName = paste(paste('MWU', goDivision, sep = "_"), input, sep = "_")
    go.res=read.table(resName, header = T)
    totSig = sum(go.res$p.adj < 0.1)
    divRec = append(divRec, goDivision)
    inRec = append(inRec, input)
    sigRec = append(sigRec, totSig)
    sig=go.res[go.res$p.adj < 0.1,]
    sigOut=paste( c('./resultsFiles/', goDivision, input, '_sigGos.tsv'), collapse='')
    if (nrow(sig)>0){
      sig %>% 
        write_tsv(path=sigOut)
    }
  }
}
res = tibble('goDivision'=divRec,
                 'input'=inRec,
                 'nSig'=sigRec)
res %>% 
  write_tsv(path='./gomwu_results_summary.tsv')



# plot results for each module with any significant enrichment -------------------
sig = res %>% 
  filter(nSig>2)

for (i in 1:nrow(sig)){
  row=sig[i,]
  goDivision=row['goDivision']
  input = row['input']
  figFileName = paste('./resultsFiles/', sep='', paste(paste(goDivision, input, sep='_'), 'tree.pdf', sep='_'))
  pdf(figFileName)
  gomwuPlot(input,goAnnotations,goDivision,
            absValue=0.00001,  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
            level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
            level2=0.05, # FDR cutoff to print in regular (not italic) font.
            level3=0.01, # FDR cutoff to print in large bold font.
            txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
            treeHeight=0.5, # height of the hierarchical clustering tree
            #	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
  )
  dev.off()
}
    

#Code below is included in separate script look_at_genes_within_GO_terms.R

#------------ LOOK AT GENES WITHIN GO TERMS ------------#
library(cowplot)
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
go='GO:0043009;GO:0009792;GO:0009790'

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
  geom_boxplot()


#------------ BUILD HEATMAP FOR A GIVEN GO TERM ------------#
library(pheatmap)


#you'll need the variance stabilized counts
ll=load('../../initialize_counts/rld.Rdata')
ll
head(rld.df)


#subset the variance stabilized counts to include only the GO of interest
g=rld.df[geneSet,]
head(g)
dim(g)

#gather gene names from names above
# names$description[]

# colnames(labs)=c('ensembl_gene_id')
# labs=merge(labs, names, by = 'ensembl_gene_id', all.x=T)
# labs$ensembl_gene_id == rownames(g)

#plot the heatmap
pheatmap(g,cluster_cols=T,border_color=NA,clustering_distance_rows="correlation")









