#plot_read_reductions.R

CUT = 0.1

# load deseq results ------------------------------------------------------

#load the coldata
ll = load('deseq/deseqInput.Rdata')
ll
head(counts)
head(coldata)

ll = load('read_sensititivity_analysis/ethanol_res_list.Rdata')
ll


# get significance counts -------------------------------------------------

get_nsig = function(res){
  nsig = sum(res$padj < CUT, na.rm=TRUE)
  return(nsig)
}

get_measured = function(res){
  m = sum(!is.na(res$log2FoldChange))
  return(m)
}


sig_counts = unlist(map(res_list, get_nsig))
measured_counts = unlist(map(res_list, get_measured))
ori_med = sum(counts)/nrow(coldata)

rdat = data.frame('pct' = names(sig_counts),
                  'nMeasured' = measured_counts,
                  'nSig' = sig_counts) %>% 
  mutate(prop = sub('pct.', '', pct, fixed=TRUE),
         prop = as.numeric(prop)/100,
         prop = if_else(pct=='original',
                        1,
                        prop),
         count = prop*ori_med,
         million = count/1e6)

splt = rdat %>% 
  ggplot(aes(x=million, y=nSig)) +
  geom_point() +
  geom_line() +
  labs(x='median reads per sample (million)',
       y='N significant genes detected')


# look at correlation drop ------------------------------------------------

original = res_list[['original']]
ol2 = original$log2FoldChange
get_cor = function(res){
  lm1 = lm(res$log2FoldChange ~ ol2)
  r2=summary(lm1)$r.squared
  return(r2)
}


cors = unlist(map(res_list, get_cor))
rdat$cors = cors


cplt = rdat %>% 
  ggplot(aes(x=million, y=cors)) +
  geom_point() +
  geom_line() +
  labs(x='median reads per sample (million)',
       y=bquote(R^2~'with full dataset'))

plot_grid(cplt, splt, nrow=1)



apply(counts, 2, sum)/1e6
