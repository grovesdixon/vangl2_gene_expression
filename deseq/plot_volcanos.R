#plot_volcanos.R

ll=load('deseq/all_res.Rdata')
source('rnaseq_functions.R')

summary(e.res)
summary(t.res)
summary(m.res)


ev=volcano_from_deseq_res(e.res, 'Ethanol effect')
tv=volcano_from_deseq_res(t.res, 'Temperature effect')
gv=volcano_from_deseq_res(m.res, 'Genotype effect')
iv=volcano_from_deseq_res(i.res, 'Etoh:Genotype interaction')
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