#!/usr/bin/env Rscript
############################
library(data.table)
library(dplyr)


##add P from gwas to sib file
###
#UKB
fread('zcat ~/height_prediction/gwas/input/50.assoc.tsv.gz')-> ukb_height #read in GWAS summary statistics for height from the Uk Biobank
ukb_height[,c("CHR", "POS","Allele2","Allele1") := tstrsplit(variant, ":", fixed=TRUE)][,variant:=NULL]
ukb_height[, MarkerName:=rsid][, N:=nCompleteSamples][, AC:=NULL][, b_u:=beta][,p_u:=pval][, SE:=se]
ukb_height[,.(MarkerName,Allele1,Allele2, SE, p_u, N, CHR, POS)]-> ukb_height
ukb_height[,CHR:=as.integer(CHR)]
ukb_height[,POS:=as.integer(POS)]
ukb_height[order(CHR,POS)]-> ukb_height
setkey(ukb_height,CHR,POS, MarkerName)
#sib
gc()
fread("~/height_prediction/sib_betas/input/sibestimates_50.tsv")-> ukb_height_sib
ukb_height_sib[,SE:= -abs(beta)/qnorm(pval/2)]
ukb_height_sib[, MarkerName:=SNP]
ukb_height_sib[,CHR:=as.integer(CHR)]
ukb_height_sib[,POS:=as.integer(POS)]
ukb_height_sib[order(CHR,POS)]-> ukb_height_sib
setkey(ukb_height_sib,CHR,POS, MarkerName)
ukb_height[ukb_height_sib, nomatch=0]-> res1
res1[,SNP:=MarkerName][,BETA:=beta][,P:=p_u]
res1[,.(SNP, CHR, POS,A1,A2, n,BETA,P)]-> res1

#
fwrite('~/height_prediction/sib_betas/input/sib_beta_gwas_p.txt', sep="\t")
