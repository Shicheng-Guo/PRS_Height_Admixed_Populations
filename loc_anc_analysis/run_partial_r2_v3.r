#!/usr/bin/env Rscript

############################
old <- Sys.time()
#args<-c('phys_100000_0.0005', 21)
## Load libraries
source('~/height_prediction/strat_prs/scripts/Rsq_R2.R')
library(optparse)
library(data.table)
library(dplyr)
library(readr)
library(tidyr)
library(parallel)
library(reshape2)
library(asbio)
library(ggplot2)
#options(scipen=999)
source('~/height_prediction/scripts/mclapply2.R')
##Load other PRSs

prs<-vector('list', 22)
for(chr in 1:22){
	prs[[chr]]<-readRDS(paste0('~/height_prediction/loc_anc_analysis/output/chr', chr,'_0_prs_EuR.Rds'))
}
samp_names<-names(prs[[1]])
dt<-data.table(SUBJID=samp_names, PRS_0=NA)

prs2<-vector('list', length(samp_names))
names(prs2)<-samp_names


for (S in samp_names){
	prs2[[S]]<-sum(prs[[1]][[S]],prs[[2]][[S]],prs[[3]][[S]],prs[[4]][[S]],prs[[5]][[S]],prs[[6]][[S]], prs[[7]][[S]], prs[[8]][[S]], prs[[9]][[S]], prs[[10]][[S]],prs[[11]][[S]], prs[[12]][[S]], prs[[13]][[S]],prs[[14]][[S]], prs[[15]][[S]], prs[[16]][[S]], prs[[17]][[S]], prs[[18]][[S]], prs[[19]][[S]], prs[[20]][[S]], prs[[21]][[S]],prs[[22]][[S]], na.rm=T)
}

dt[, PRS_0:=scale(unlist(prs2))]

#phenotype
fread('~/height_prediction/input/WHI/WHI_phenotypes.txt')-> Pheno_WHI

Pheno_WHI[, SUBJID:=paste0("0_", as.character(Pheno_WHI[, SUBJID]))]
setkey(Pheno_WHI, SUBJID)
#add ancestry
ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/WHI/rfmix_anc_chr', X, '.txt'))))
anc_WHI<-ancestry %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual
setkey(anc_WHI, SUBJID)
##

#PRS_EUR<-data.table(PRS_EUR=unlist(readRDS('~/height_prediction/gwas/WHI/output/PGS_WHI_phys_100000_0.0005.Rds')),SUBJID=names(unlist(readRDS('~/height_prediction/gwas/WHI/output/PGS_WHI_phys_100000_0.0005.Rds'))))
#a<-merge(dt, PRS_EUR, by="SUBJID")
setkey(dt, SUBJID)
setkey(Pheno_WHI, SUBJID)
dt[Pheno_WHI][anc_WHI]-> final
final2<-final[,c("SUBJID","PRS_0","AGE", "EUR_ANC", "HEIGHTX")]

final2[,AGE2:=AGE^2]

#
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0, data=final2))*100 # 1%, compared to 3.87% for when using all SNPs and alfa=0
#
saveRDS(final2, '~/height_prediction/loc_anc_analysis/output/all_PRS_WHI_Eur.Rds')
remove(list=ls())
###################
##JHS

####
prs<-vector('list', 22)
for(chr in 1:22){
        prs[[chr]]<-readRDS(paste0('~/height_prediction/loc_anc_analysis/output/chr', chr,'_0_prs_JHS_EuR.Rds'))
}
samp_names<-names(prs[[1]])
dt<-data.table(SUBJID=samp_names, PRS_0=NA)

prs2<-vector('list', length(samp_names))
names(prs2)<-samp_names


for (S in samp_names){
        prs2[[S]]<-sum(prs[[1]][[S]],prs[[2]][[S]],prs[[3]][[S]],prs[[4]][[S]],prs[[5]][[S]],prs[[6]][[S]], prs[[7]][[S]], prs[[8]][[S]], prs[[9]][[S]], prs[[10]][[S]],prs[[11]][[S]], prs[[12]][[S]], prs[[13]][[S]],prs[[14]][[S]], prs[[15]][[S]], prs[[16]][[S]], prs[[17]][[S]], prs[[18]][[S]], prs[[19]][[S]], prs[[20]][[S]], prs[[21]][[S]],prs[[22]][[S]], na.rm=T)
}

dt[, PRS_0:=scale(unlist(prs2))]
dt$SUBJID<-substr(dt$SUBJID, 1,7))
#phenotype
fread('~/height_prediction/input/JHS/JHS_phenotypes.txt')-> Pheno_JHS
Pheno_JHS[,AGE2:=age_baseline^2]
setkey(Pheno_JHS, SUBJID)
#add ancestry
anc_JHS<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/JHS/rfmix_anc_chr', X, '.txt'))))
anc_JHS$SUBJID<-substr(anc_JHS[,SUBJID],3,9)
anc_JHS<-anc_JHS %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual
setkey(anc_JHS, SUBJID)
##
setkey(dt, SUBJID)
setkey(Pheno_JHS, SUBJID)
dt[Pheno_JHS][anc_JHS]-> final
final[, AGE:=age_baseline]
final[, HEIGHTX:=height_baseline]
final2<-final[,c("SUBJID","PRS_0","EUR_ANC", "SEX", "HEIGHTX", "AGE2", "AGE")]

partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_0, data=final2))*100 # 0.75%
#
saveRDS(final2, '~/height_prediction/loc_anc_analysis/output/all_PRS_JHS_Eur.Rds')
remove(list=ls())
