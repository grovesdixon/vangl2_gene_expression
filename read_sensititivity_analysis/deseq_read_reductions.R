#deseq_read_reductions.R

rm(list=ls())
library(tidyverse)
library('cowplot')
theme_set(theme_cowplot())

# load the data -----------------------------------------------------------

#load the coldata
ll = load('deseq/deseqInput.Rdata')
ll
head(counts)
head(coldata)

# BUILD REDUCED COUNT SETS ------------------------------------------------

# #build a step-wise reduced counts matrices
reductionMultipliers0 = sapply(4:1, function(x) return(2^-x))
reductionMultipliers = append(reductionMultipliers0, seq(0.6, 1, by=0.05))
reductionMultipliers = seq(0.05, 1, by=0.05)
resampledList = list()

for (redM in reductionMultipliers){
  pctLab = paste('pct', redM*100, sep='.')
  print('---------------')
  print(paste('reduction =', redM))
  resampled = round(counts*redM, digits=0)
  print('mean resampled proportion:')
  print(sum(resampled) / sum(counts))
  resampledList[[pctLab]] = resampled
}
names(resampledList)
resampledList[['original']]=counts #load the original as the set you started with (no resampling)

# double-check things worked -----------------------------------------------------

#check totals
reductionMultipliers
map(resampledList, function(x) return(sum(x)/sum(counts)))

#check lineup
samples = coldata$sample
check_lineup = function(x){
  sum(samples==colnames(x)) == ncol(x)
}
map(resampledList, head)
check = unlist(map(resampledList, check_lineup))
sum(check)==length(resampledList)

# run deseq on each -------------------------------------------------------
library(DESeq2)
#function to run deseq as in main script
run_deseq = function(counts){
  dds<-DESeqDataSetFromMatrix(counts,
                              colData = coldata, 
                              design = formula(~Ethanol +
                                                 Temperature +
                                                 Genotype +
                                                 Temperature:Genotype +
                                                 Ethanol:Genotype))
  dds <- DESeq(dds)
  return(dds)
}

#loop through each resampled counts matrix and run

total_count_list = list()
for (n in names(resampledList)){
  print('---------')
  print(n)
  counts = resampledList[[n]]
  total_count_list[[n]] = sum(as.matrix(counts), na.rm=TRUE)
}


dds_list = list()
names(resampledList)
for (n in names(resampledList)){
  print('---------')
  print(n)
  counts = resampledList[[n]]
  ndds = run_deseq(counts)
  dds_list[[n]] = ndds
}


# pull results ------------------------------------------------------------

#choose results parameters
INDEPDENDENT_FILTERING = FALSE
contrast = c('Ethanol', 'ethanol', 'control')

#loop through deseq objects to pull results
res_list = list()
for (n in names(dds_list)){
  print('----------')
  print(n)
  ndds = dds_list[[n]]
  r = results(ndds, contrast = contrast, independentFiltering=INDEPDENDENT_FILTERING)
  res_list[[n]] = r
}



# save --------------------------------------------------------------------

save(res_list, total_count_list, file='read_sensititivity_analysis/ethanol_res_list.Rdata')


