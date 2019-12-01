#zebrafish_RNAseq_functions.R
library(genefilter)
library(ggplot2)
write_out_go = function(df, outPath){
  go = data.frame(df) %>% 
    mutate(gene=rownames(df),
           upregulated = log2FoldChange>0,
           logp= if_else(upregulated,
                         -log(pvalue, 10),
                         -log(pvalue, 10)*-1)) %>%
    dplyr::select(gene,logp)
  go %>% 
    write_csv(outPath)
}


volcano_from_deseq_res = function(res, TITLE, YLIM=c(0,35), XLIM=c(-2.5, 2.5), ALPHA=0.5){
  sdf = data.frame(res) %>% 
    arrange(pvalue) %>% 
    mutate(sig = !is.na(padj) & padj < 0.1) %>% 
    as_tibble()
  g=ggplot(data=sdf, aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
    geom_point(alpha=ALPHA) +
    scale_color_manual(values=c('black', 'red')) + 
    labs(subtitle=TITLE,
         y=bquote(-log[10]~pvalue),
         x=bquote(log[2]~fold~difference)) +
    lims(y=YLIM, x=XLIM) +
    guides(color = guide_legend(title='FDR < 0.1'))
  return(g)
}



#modified version of the function provided with the DESeq package
#allows for output of more principal components
mod.plotPCA <- function (object, intgroup = "condition", ntop = 25000, returnData = F, pcs = 10, pc1 = 1, pc2 = 2, main = "\n") 
{
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
        length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
        drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = " : "))
    }
    else {
        colData(object)[[intgroup]]
    }
    d <- data.frame(pca$x[,1:pcs], group = group, 
        intgroup.df, name = colnames(object))
    attr(d, "percentVar") <- percentVar[1:2]
    g = ggplot(data = d, aes_string(x = paste('PC', pc1, sep = ''),
    y = paste('PC', pc2, sep = ''), color = "group")) + 
        geom_point(size = 3) + xlab(paste0(paste0(paste0("PC", pc1), ": "), round(percentVar[pc1] * 
        100), "% variance")) + ylab(paste0(paste0(paste0("PC", pc2), ": "), round(percentVar[pc2] * 
        100), "% variance")) + coord_fixed()
   g = g + ggtitle(main)
   g = g + theme_bw()
   print(g)
   return(d)
}



printPCA = function(dat, pc1, pc2){
	ggplot(data = dat, aes_string(x = paste('PC', pc1, sep = ''), y = paste('PC', pc2, sep = ''), color = "group")) + 
        geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 
        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 
        100), "% variance")) + coord_fixed()
}




#modified version of the function provided in the DESeq package that can be run from a dataframe
#instead of DESeqTransform object output from
mod.plotPCA.df <- function (df, coldat, intgroup = "condition", ntop = 25000, returnData = F, pcs = 10, pc1 = 1, pc2 = 2, main = "\n", SIZE = 3) 
{
    rv <- rowVars(df)
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
        length(rv)))]
    pca <- prcomp(t(df[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    intgroup.df <- as.data.frame(coldat[, intgroup, 
        drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = " : "))
    }
    else {
        coldat[[intgroup]]
    }
    d <- data.frame(pca$x[,1:pcs], group = group, 
        intgroup.df, name = colnames(df))
    attr(d, "percentVar") <- percentVar[1:2]
    g = ggplot(data = d, aes_string(x = paste('PC', pc1, sep = ''),
    y = paste('PC', pc2, sep = ''), color = "group")) + 
        geom_point(size = SIZE) + xlab(paste0(paste0(paste0("PC", pc1), ": "), round(percentVar[pc1] * 
        100), "% variance")) + ylab(paste0(paste0(paste0("PC", pc2), ": "), round(percentVar[pc2] * 
        100), "% variance")) + coord_fixed()
   g = g + ggtitle(main) +
     theme(legend.position='bottom')
   # g = g + theme_bw()
   print(g)
   if (returnData == T){
   return(d)
   } else {
     return(g)
   }
}




#Function to write an input file for GO MWU analysis
#dat = a DESeq results object (eg: results(dds, contrast=c('treatment', 'control', 'treated')))
#out.name = the desired name for the output file
write.gomwu.input = function(dat, out.name){
	sign = dat[,'log2FoldChange'] > 0
	sign[sign == TRUE] <- 1
	sign[sign == FALSE] <- -1
	# logp = (-log(dat[,'pvalue'], 10) ) * sign
	stat = abs(dat$stat)*sign
	out = data.frame(rownames(dat), stat)
	colnames(out) = c('gene', 'stat')
	head(out)
	write.csv(out, out.name, row.names = F, quote = F)
}




merge_gene_names = function(dat, col=F, sort.column=F){
	library("biomaRt")
	if(col!=FALSE){
		geneSet=dat[,col]
	}
	else{
		geneSet=rownames(dat)
		}
	embl = useMart("ensembl", dataset="drerio_gene_ensembl")
	names = getBM(attributes = c('ensembl_gene_id', 'external_gene_name','description', 'chromosome_name', 'start_position', 'end_position'), filters = c('ensembl_gene_id'), values = geneSet, mart = embl)
	dat2=dat
	dat2$ensembl_gene_id<-geneSet
	dat3=merge(dat2, names, by = 'ensembl_gene_id', all.x=T)
	if(sort.column!=F){
		dat3=dat3[order(dat3[,sort.column]),]
	}
	return(dat3)
}


get_human_ids = function(geneSet){
	hnames = getBM(attributes = c('ensembl_gene_id', 'external_gene_name','description'), filters = c('external_gene_name'), values = geneSet, mart = hsembl)
	hgeneSet=hnames$external_gene_name
	prop=round(length(hgeneSet) / length(geneSet), digits=2)*100
	pct=paste(prop, "%", sep='')
	print(paste("Percentage of gene set with human IDs =", pct))
	return(hgeneSet)
}


get_gene_names = function(geneSet){
	library("biomaRt")
	embl = useMart("ensembl", dataset="drerio_gene_ensembl")
	names = getBM(attributes = c('ensembl_gene_id', 'external_gene_name','description'), filters = c('ensembl_gene_id'), values = geneSet, mart = embl)
	return(names)
}


output_module_for_gomwu = function(dat){
	
}


