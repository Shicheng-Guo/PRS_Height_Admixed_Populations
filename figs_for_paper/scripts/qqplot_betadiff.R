#!/usr/bin/env Rscript
library(readr)
library(data.table)
library(reshape2)
library(asbio)
library(ggplot2)
library(dplyr)
library(parallel)
library(mgcViz)
#########################
source('~/height_prediction/strat_prs/scripts/Rsq_R2.R')
source('~/height_prediction/gwas/WHI/scripts/PolygenicScore_v3.R')
source('~/height_prediction/gwas/WHI/scripts/short_fun.R')
########################
cat('checkpoint number 1\n')

plink<-fread('~/height_prediction/runSmartpCA-master/UKB_AFR/association_v3.Res.Height.glm.linear.adjusted', header=T, fill=T)
colnames(plink)[2]<-'MarkerName'
colnames(plink)[1]<-'CHR'
N<-8816*2
#
plink2<-fread('~/height_prediction/runSmartpCA-master/UKB_AFR/test3.txt', fill=T)
colnames(plink2)<-c("CHR","POS", "MarkerName","REF","ALT","A1","TEST","OBS_CT","PLINK", "SE","T_STAT", "UNADJ") #plink is BETA
setkey(plink, MarkerName, CHR, UNADJ)
setkey(plink2, MarkerName, CHR, UNADJ)
plink[plink2, nomatch=0]-> final_plink
final_plink$CHR<-as.numeric(final_plink$CHR)
arrange(final_plink, CHR,POS) %>% as.data.table -> final_plink
setkey(final_plink, MarkerName, CHR, POS)

ukb_height<-fread('zcat ~/height_prediction/gwas/input/50.assoc.tsv.gz')
ukb_height[,c("CHR", "POS","Allele2","Allele1") := tstrsplit(variant, ":", fixed=TRUE)][,variant:=NULL]
ukb_height[, MarkerName:=rsid][, N:=nCompleteSamples][, AC:=NULL][, b:=beta][,p:=pval]
ukb_height[,rsid:=NULL][,nCompleteSamples:=NULL][, beta:=NULL][, pval:=NULL][, SE:=se]
ukb_height[,.(MarkerName,Allele1,Allele2, b, SE, p, N, CHR, POS)]-> ukb_height
ukb_height$CHR<-as.numeric(ukb_height$CHR)
ukb_height$POS<-as.numeric(ukb_height$POS)
setkey(ukb_height, MarkerName, CHR, POS)

#combine
final_plink[ukb_height, nomatch=0]-> combo
combo[,TEST:=NULL]
combo[,Beta_Diff:=abs(PLINK-b)]
combo[, Beta_Diff_Chisq:=((PLINK-b)/(sqrt(((SE)^2+(i.SE)^2))))^2]

fread('/project/mathilab/bbita/gwas_admix/new_height/ukb_afr_betas_100000_0.0005.txt')-> prs_snps

combo_prs<-combo[MarkerName %in% prs_snps$MarkerName]


#qqplot


pdf('~/height_prediction/figs_for_paper/figs/qqplot_test.pdf')
qq.chisq(combo$Beta_Diff_Chisq, df=2,  main="QQ plot ALL SNPs",  xlab="Expected", ylab="Observed")
qq.chisq(combo_prs$Beta_Diff_Chisq, df=2,  main="QQ plot PRS",  xlab="Expected", ylab="Observed")
dev.off()
