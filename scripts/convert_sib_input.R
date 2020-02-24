#!/usr/bin/env Rscript
############################
library(data.table)
library(dplyr)


##add P from gwas to sib file
###
#UKB
fread('zcat ~/height_prediction/gwas/input/50_raw_filtered.txt.gz')-> ukb_height #read in GWAS summary statistics for height from the Uk Biobank
ukb_height[,c("CHR", "POS","REF","ALT") := tstrsplit(variant, ":", fixed=TRUE)][,variant:=NULL]
gc()
ukb_height[, N:=n_complete_samples][, AC:=NULL][, b_u:=beta][,p_u:=pval][, SE:=se]
ukb_height[,.(REF,ALT, SE, p_u, N, CHR, POS)]-> ukb_height
gc()
ukb_height[CHR %in% 1:22]-> ukb_height
ukb_height[,CHR:=as.integer(CHR)]
ukb_height[,POS:=as.integer(POS)]
gc()
ukb_height[order(CHR,POS)]-> ukb_height
setkey(ukb_height,CHR,POS)
#sib
gc()
fread("~/height_prediction/sib_betas/input/sibestimates_50.tsv")-> ukb_height_sib
ukb_height_sib[,SE:= -abs(beta)/qnorm(pval/2)]
ukb_height_sib[, MarkerName:=SNP]
ukb_height_sib[,CHR:=as.integer(CHR)]
ukb_height_sib[,POS:=as.integer(POS)]
ukb_height_sib[order(CHR,POS)]-> ukb_height_sib
setkey(ukb_height_sib,CHR,POS)
ukb_height[ukb_height_sib, nomatch=0]-> res1
res1[,SNP:=MarkerName][,BETA:=beta][,P:=p_u]
res1[,.(SNP, CHR, POS,A1,A2, n,BETA,P)]-> res1

#
fwrite(res1, '~/height_prediction/sib_betas/input/sib_beta_gwas_p.txt', sep="\t")


JHS<-fread('~/height_prediction/input/JHS/JHS_b37_strand.bim')[,c(1,2,4)]
colnames(JHS)<-c('CHR', 'SNP', 'POS')
setkey(JHS, CHR, POS)
res1[JHS,nomatch=0][, SNP:=i.SNP][,i.SNP:=NULL][,N:=n][, n:=NULL]-> JHS1
JHS1[,.(SNP,A1,A2, BETA, P, N, CHR, POS)]-> JHS1
fwrite(JHS1,'~/height_prediction/sib_betas/JHS/plink2/output/betas_for_plink.txt', sep = "\t")
remove(JHS, JHS1)
#WHI
WHI<-fread('~/height_prediction/input/WHI/WHI_b37_strand_include.bim')[,c(1,2,4)]
colnames(WHI)<-c('CHR', 'SNP', 'POS')
setkey(WHI, CHR, POS)
res1[WHI,nomatch=0][, SNP:=i.SNP][,i.SNP:=NULL][,N:=n][, n:=NULL]-> WHI1
WHI1[,.(SNP,A1,A2, BETA, P, N, CHR, POS)]-> WHI1
fwrite(WHI1,'~/height_prediction/sib_betas/WHI/plink2/output/betas_for_plink.txt', sep = "\t")
remove(WHI, WHI1)
#HRS_afr
HRS_afr<-fread('~/height_prediction/input/HRS_afr/HRS_AFR_b37_strand_include.bim')[,c(1,2,4)]
colnames(HRS_afr)<-c('CHR', 'SNP', 'POS')
setkey(HRS_afr, CHR, POS)
res1[HRS_afr,nomatch=0][, SNP:=i.SNP][,i.SNP:=NULL][,N:=n][, n:=NULL]-> HRS_afr1
HRS_afr1[,.(SNP,A1,A2, BETA, P, N, CHR, POS)]-> HRS_afr1
fwrite(HRS_afr1,'~/height_prediction/sib_betas/HRS_afr/plink2/output/betas_for_plink.txt', sep = "\t")
remove(HRS_afr, HRS_afr1)
#
HRS_eur<-fread('~/height_prediction/input/HRS_eur/HRS_EUR_b37_strand_include.bim')[,c(1,2,4)]
colnames(HRS_eur)<-c('CHR', 'SNP', 'POS')
setkey(HRS_eur, CHR, POS)
res1[HRS_eur,nomatch=0][, SNP:=i.SNP][,i.SNP:=NULL][,N:=n][, n:=NULL]-> HRS_eur1
HRS_eur1[,.(SNP,A1,A2, BETA, P, N, CHR, POS)]-> HRS_eur1
fwrite(HRS_eur1,'~/height_prediction/sib_betas/HRS_eur/plink2/output/betas_for_plink.txt', sep = "\t")
remove(HRS_eur, HRS_eur1)
#
UKB_eur<-fread('~/height_prediction/input/ukb_eur/UKB_EUR.bim')[,c(1,2,4)]
colnames(UKB_eur)<-c('CHR', 'SNP', 'POS')
setkey(UKB_eur, CHR, POS)
res1[UKB_eur,nomatch=0][, SNP:=i.SNP][,i.SNP:=NULL][,N:=n][, n:=NULL]-> UKB_eur1
UKB_eur1[,.(SNP,A1,A2, BETA, P, N, CHR, POS)]-> UKB_eur1
fwrite(UKB_eur1,'~/height_prediction/sib_betas/ukb_eur/plink2/output/betas_for_plink.txt', sep = "\t")
remove(UKB_eur, UKB_eur1)
#
UKB_afr<-fread('~/height_prediction/input/ukb_afr/UKB_AFR.bim')[,c(1,2,4)]
colnames(UKB_afr)<-c('CHR', 'SNP', 'POS')
setkey(UKB_afr, CHR, POS)
res1[UKB_afr,nomatch=0][, SNP:=i.SNP][,i.SNP:=NULL][,N:=n][, n:=NULL]-> UKB_afr1
UKB_afr1[,.(SNP,A1,A2, BETA, P, N, CHR, POS)]-> UKB_afr1
fwrite(UKB_afr1,'~/height_prediction/sib_betas/ukb_afr/plink2/output/betas_for_plink.txt', sep = "\t")
remove(UKB_afr, UKB_afr1)

