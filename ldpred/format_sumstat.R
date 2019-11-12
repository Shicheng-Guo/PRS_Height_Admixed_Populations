#!/usr/bin/env Rscript
######################################
# Format summary statistics file to  #
# use with ldpred                    #  
######################################
library(data.table)
library(dplyr)

input<-fread("zcat ~/height_prediction/gwas/input/50.assoc.tsv.gz")
input[,c("CHR", "POS","A2","A1") := tstrsplit(variant, ":", fixed=TRUE)][,variant:=NULL]  #fix columns. In the UKB, variants are listed as CHR:POS: REF:ALT, where ALT is the effect allele. So, in order to be compatilble with my scripts (where allele 1 is effect all
input[, SNP_ID:=rsid][, N:=nCompleteSamples][, AC:=NULL][, BETA:=beta][, PVAL:=pval]
input[,rsid:=NULL][,nCompleteSamples:=NULL][, beta:=NULL][, SE:=se]
input[,.(CHR, POS, SNP_ID, A1,A2,PVAL BETA, SE,N)]-> input
#ainput[,CHR:=paste0('chr', CHR)]
arrange(input, CHR, POS) %>% as.data.table -> input
fwrite(input, file="~/height_prediction/ldpred/output/Height.QC.txt", sep="\t") #sumstats file needs to be \t separated for ldpred
system('bgzip ~/height_prediction/ldpred/output/Height.QC.txt')
system('mv ~/height_prediction/ldpred/output/Height.QC.txt.gz ~/height_prediction/ldpred/output/Height.QC.gz')

