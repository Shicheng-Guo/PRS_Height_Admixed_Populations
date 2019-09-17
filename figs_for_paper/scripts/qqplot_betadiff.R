#!/usr/bin/env Rscript
library(data.table)
library(snpStats)
library(ggplot2)
library(dplyr)
library(parallel)
library(mgcViz)
#########################
cat('checkpoint number 1\n')

pdf('~/height_prediction/figs_for_paper/figs/qqplot_test.pdf')
qq.chisq(combo$Beta_Diff_Chisq, df=2,  main="QQ plot ALL SNPs",  xlab="Expected", ylab="Observed")
qq.chisq(combo_prs$Beta_Diff_Chisq, df=2,  main="QQ plot PRS",  xlab="Expected", ylab="Observed")
dev.off()

