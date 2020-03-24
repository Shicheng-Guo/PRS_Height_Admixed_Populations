#!/usr/bin/env Rscript
######################################
# Format summary statistics file to  #
# use with ldpred                    #  
######################################
library(data.table)
library(dplyr)

input<-fread("zcat ~/height_prediction/gwas/input/50_raw_filtered.txt.gz")
input[,c("CHR", "POS","A2","A1") := tstrsplit(variant, ":", fixed=TRUE)][,variant:=NULL]  #fix columns. In the UKB, variants are listed as CHR:POS: REF:ALT, where ALT is the effect allele. So, in order to be compatilble with my scripts (where allele 1 is effect all
input[CHR %in% seq(1:22)]-> input
gc()
input[, N:=n_complete_samples][, AC:=NULL][, BETA:=beta][, PVAL:=pval]
input[,n_complete_samples:=NULL][, beta:=NULL][, SE:=se]
input[,.(CHR, POS, A1,A2,PVAL, BETA, SE,N)]-> input
gc()
input<-unique(input, by=c('CHR', 'POS'))
gc()
as.integer(input$CHR)-> input$CHR
as.integer(input$POS)-> input$POS
arrange(input, CHR, POS) %>% as.data.table -> input
input[,SNP:=paste0(CHR, ":", POS)]
fwrite(input, file="~/height_prediction/ldpred/output/Height.QC.txt", sep="\t") #sumstats file needs to be \t separated for ldpred
system('bgzip ~/height_prediction/ldpred/output/Height.QC.txt')
system('mv ~/height_prediction/ldpred/output/Height.QC.txt.gz ~/height_prediction/ldpred/output/Height.QC.gz')


