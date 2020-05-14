#README
#This folder contains R scripts to testing the effect of randomly reducing the fold coverage

1. reduce the reads by randomly resampling proportionally to original dataset
	read_reductions.R
2. re-run DESeq on the reduced datasets
	deseq_read_reductions.R
3. plot results for comparison
	plot_read_reductions.R