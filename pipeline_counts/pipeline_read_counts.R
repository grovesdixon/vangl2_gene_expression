

library(tidyverse)
library(cowplot)


input = 'pipeline_counts/all_pipeline_counts.txt'


sdat0 = read.table(input, header=F, col.names=c('Run', 'value', 'stat'), stringsAsFactors=F) %>% 
  as_tibble()

laneSummed = sdat0 %>% 
  filter(stat %in% c('rawCounts', 'trimmedCounts')) %>% 
  group_by(Run, stat) %>% 
  summarize(value=sum(value)) 


statLevels = c("rawCounts", "trimmedCounts", "predupPropPaired", "dedupPropPair", "geneCounted")

sdat = sdat0 %>% 
  filter(!stat %in% c('rawCounts', 'trimmedCounts')) %>% 
  bind_rows(laneSummed) %>% 
  filter(stat %in% statLevels) %>% 
  mutate(stat=factor(stat, levels=statLevels),
         temp = substr(Run, 1,1),
         ethanol=grepl('E', Run),
         mutant=grepl('M', Run))


#order the stats in pipline order
order = sdat %>% 
  filter(stat=='rawCounts') %>% 
  arrange(by=value) %>% 
  pull(Run) %>% 
  rev() 
sdat$Run = factor(sdat$Run, levels=order)


#get table of results
sdat %>% 
  group_by(stat) %>% 
  summarize(`Mean(million)`=mean(value)/1e6) %>% 
  write_tsv('pipeline_counts/summary_table.tsv')

#plot full barplot
bp<-sdat %>% 
  ggplot(aes(x=stat, y=value, color=ethanol, fill=ethanol, group=Run)) +
    geom_bar(stat='identity', position='dodge') +
    labs(y='Read count', x='Pipeline step') +
    theme(axis.text.x=element_text(angle=20, vjust=0.75))
bp


#plot gene counted only
gc = sdat %>% 
  filter(stat=='geneCounted') %>% 
  arrange(by=value)
gc$Run = factor(gc$Run, levels = rev(pull(gc, Run)))

histMeans = gc %>% 
  ggplot(aes(x=value)) +
  geom_histogram(bins=12) +
  labs(x='Mean gene counted reads accross projects')

rankedMeans = gc %>% 
  ggplot(aes(x=stat, y=value, color=ethanol, fill=ethanol, group=Run)) +
  geom_bar(stat='identity', position='dodge') +
  labs(y='Read count', x='Pipeline step', title='Mean Counts')
rankedMeans

plot_grid(histMeans, rankedMeans, nrow=1, rel_widths=c(.5,1))

#plot scatter abs
quartz()
lp<-sdat %>% 
  ggplot(aes(x=stat, y=value, color=ethanol, group=Run)) +
  geom_point() +
  geom_line(aes(group=Run)) +
  theme(legend.position='none', axis.text.x=element_text(angle=20, vjust=0.75)) +
  labs(y='Read count', x='Pipeline step', subtitle='read counts')
  

#plot scatter prop raw
pp<-sdat %>% 
  group_by(Run) %>% 
  mutate(prop=value/max(value)) %>% 
  ggplot(aes(x=stat, y=prop, color=ethanol, group=Run)) +
    geom_point() +
    geom_line(aes(group=Run)) +
    labs(y='Proportion raw reads', x='Pipeline step', subtitle='read proportions') +
    theme(axis.text.x=element_text(angle=20, vjust=0.75))

plot_grid(lp, pp, histMeans, rankedMeans, nrow=2, rel_widths=c(1,1.3))


#PLOT ALL 4 TOGETHER
plot_grid(lp, pp, )






quartz()
lp = sdat %>% 
  filter(stat %in% c("rawCounts", "trimmedCounts", "predupMapped", "dedupMapped", "geneCounted")) %>% 
  mutate(stat=factor(stat, levels=c("rawCounts", "trimmedCounts", "predupMapped", "dedupMapped", "geneCounted"))) %>% 
  ggplot(aes(x=stat, y=fixValue, color=my_title)) +
  geom_point() +
  geom_line(aes(group=Run)) +
  theme(legend.position='none', axis.text.x=element_text(angle=20, vjust=0.75)) +
  labs(y='Read count', x='Pipeline step', subtitle='read counts')


pp<-sdat %>% 
  filter(stat %in% c("rawCounts", "trimmedCounts", "predupMapped", "dedupMapped", "geneCounted")) %>% 
  mutate(stat=factor(stat, levels=c("rawCounts", "trimmedCounts", "predupMapped", "dedupMapped", "geneCounted"))) %>% 
  group_by(Run) %>% 
  mutate(prop = fixValue/max(fixValue) ) %>% 
  ggplot(aes(x=stat, y=prop, color=my_title)) +
    geom_point() +
    geom_line(aes(group=Run)) +
    labs(y='Proportion raw reads', x='Pipeline step', subtitle='read proportions') +
    theme(legend.position='none', axis.text.x=element_text(angle=20, vjust=0.75))
plot_grid(lp, pp)
