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

readRDS("~/height_prediction/imputed/output/all_prs.Rds")-> res_all
names(res_all)<-c("PLINK", "PLINK_TSTAT_1")
data.table(SUBJID=names(res_all[[1]][[1]]), PRS_afr=(unlist(res_all[[1]][[1]])+unlist(res_all[[2]][[1]])+unlist(res_all[[3]][[1]])+unlist(res_all[[4]][[1]])+unlist(res_all[[5]][[1]])+unlist(res_all[[6]][[1]])+unlist(res_all[[7]][[1]])+unlist(res_all[[8]][[1]])+unlist(res_all[[9]][[1]])+unlist(res_all[[10]][[1]])+unlist(res_all[[11]][[1]])+unlist(res_all[[12]][[1]])+unlist(res_all[[13]][[1]])+unlist(res_all[[14]][[1]])+unlist(res_all[[15]][[1]])+unlist(res_all[[16]][[1]])+unlist(res_all[[17]][[1]])+ unlist(res_all[[18]][[1]])+unlist(res_all[[19]][[1]])+unlist(res_all[[20]][[1]])+unlist(res_all[[21]][[1]])+unlist(res_all[[22]][[1]])),PRS_EUR=unlist(readRDS('~/height_prediction/gwas/WHI/output/PGS_WHI_phys_100000_0.0005.Rds')))-> a

a[, PRS_afr:=scale(PRS_afr)]

#read in the prss and add them
prs<-vector('list', 29)
names(prs)<-c("0", "0.05","0.1", "0.15", "0.2","0.21", "0.22", "0.23", "0.24","0.25","0.26", "0.27", "0.28", "0.29", "0.3","0.35","0.4","0.45","0.5","0.55", "0.6","0.65", "0.7", "0.75","0.8", "0.85","0.9","0.95","1")
for (A in names(prs)){
	prs[[A]]<-vector('list',22)
		for(chr in 1:22){
		prs[[A]][[chr]]<-readRDS(paste0('~/height_prediction/loc_anc_analysis/output/chr', chr,'_',A, '_prs.Rds'))
	}
}
samp_names<-names(prs[[1]][[1]])
dt<-data.table(SUBJID=samp_names, PRS_0=NA, PRS_0.05=NA, PRS_0.1=NA, PRS_0.15=NA, PRS_0.2=NA, PRS_0.21=NA, PRS_0.22=NA, PRS_0.23=NA, PRS_0.24=NA,  PRS_0.25=NA, PRS_0.26=NA, PRS_0.27=NA, PRS_0.28=NA, PRS_0.29=NA, PRS_0.3=NA, PRS_0.35=NA, PRS_0.4=NA, PRS_0.45=NA, PRS_0.5=NA, PRS_0.55=NA, PRS_0.6=NA,PRS_0.65=NA, PRS_0.7=NA, PRS_0.75=NA,PRS_0.8=NA, PRS_0.85=NA, PRS_0.9=NA, PRS_0.95=NA, PRS_1=NA)

prs2<-vector('list', 29)
names(prs2)<-c("0", "0.05","0.1", "0.15", "0.2","0.21", "0.22", "0.23", "0.24","0.25","0.26", "0.27", "0.28", "0.29","0.3","0.35","0.4","0.45","0.5","0.55", "0.6","0.65", "0.7", "0.75","0.8", "0.85","0.9","0.95","1")
for (A in names(prs2)){
	prs2[[A]]<-vector('list', length(prs[[1]][[1]]))
	names(prs2[[A]])<-gsub("_A", "",samp_names)
}

for (A in names(prs2)){
	for (S in names(prs2[[A]])){
	prs2[[A]][[S]]<-sum(prs[[A]][[1]][[S]],prs[[A]][[2]][[S]],prs[[A]][[3]][[S]],prs[[A]][[4]][[S]],prs[[A]][[5]][[S]],prs[[A]][[6]][[S]], prs[[A]][[7]][[S]], prs[[A]][[8]][[S]], prs[[A]][[9]][[S]], prs[[A]][[10]][[S]],prs[[A]][[11]][[S]], prs[[A]][[12]][[S]], prs[[A]][[13]][[S]],prs[[A]][[14]][[S]], prs[[A]][[15]][[S]], prs[[A]][[16]][[S]], prs[[A]][[17]][[S]], prs[[A]][[18]][[S]], prs[[A]][[19]][[S]], prs[[A]][[20]][[S]], prs[[A]][[21]][[S]],prs[[A]][[22]][[S]], na.rm=T)
	}
}

dt[, PRS_0:=scale(unlist(prs2[[1]]))]
dt[, PRS_0.05:=scale(unlist(prs2[[2]]))]
dt[, PRS_0.1:=scale(unlist(prs2[[3]]))]
dt[, PRS_0.15:=scale(unlist(prs2[[4]]))]
dt[, PRS_0.2:=scale(unlist(prs2[[5]]))]
dt[, PRS_0.21:=scale(unlist(prs2[[6]]))]
dt[, PRS_0.22:=scale(unlist(prs2[[7]]))]
dt[, PRS_0.23:=scale(unlist(prs2[[8]]))]
dt[, PRS_0.24:=scale(unlist(prs2[[9]]))]
dt[, PRS_0.25:=scale(unlist(prs2[[10]]))]
dt[, PRS_0.26:=scale(unlist(prs2[[11]]))]
dt[, PRS_0.27:=scale(unlist(prs2[[12]]))]
dt[, PRS_0.28:=scale(unlist(prs2[[13]]))]
dt[, PRS_0.29:=scale(unlist(prs2[[14]]))]
dt[, PRS_0.3:=scale(unlist(prs2[[15]]))]
dt[, PRS_0.35:=scale(unlist(prs2[[16]]))]
dt[, PRS_0.4:=scale(unlist(prs2[[17]]))]
dt[, PRS_0.45:=scale(unlist(prs2[[18]]))]
dt[, PRS_0.5:=scale(unlist(prs2[[19]]))]
dt[, PRS_0.55:=scale(unlist(prs2[[20]]))]
dt[, PRS_0.6:=scale(unlist(prs2[[21]]))]
dt[, PRS_0.65:=scale(unlist(prs2[[22]]))]
dt[, PRS_0.7:=scale(unlist(prs2[[23]]))]
dt[, PRS_0.75:=scale(unlist(prs2[[24]]))]
dt[, PRS_0.8:=scale(unlist(prs2[[25]]))]
dt[, PRS_0.85:=scale(unlist(prs2[[26]]))]
dt[, PRS_0.9:=scale(unlist(prs2[[27]]))]
dt[, PRS_0.95:=scale(unlist(prs2[[28]]))]
dt[, PRS_1:=scale(unlist(prs2[[29]]))]

merge(dt,a, by="SUBJID")-> dt
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
final2<-final[,c("SUBJID","PRS_0", "PRS_0.05","PRS_0.1", "PRS_0.15", "PRS_0.2","PRS_0.21", "PRS_0.22","PRS_0.23", "PRS_0.24", "PRS_0.25", "PRS_0.26","PRS_0.27","PRS_0.28","PRS_0.29","PRS_0.3","PRS_0.35","PRS_0.4","PRS_0.45","PRS_0.5","PRS_0.55", "PRS_0.6","PRS_0.65", "PRS_0.7", "PRS_0.75","PRS_0.8", "PRS_0.85","PRS_0.9","PRS_0.95","PRS_1", "AGE", "EUR_ANC", "HEIGHTX", "PRS_EUR", "PRS_afr")]

final2[,AGE2:=AGE^2]

#
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2))*100 #4.1%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0, data=final2))*100 #3.872254
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.05, data=final2))*100 #3.909988
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.1, data=final2))*100  # 3.962055
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.15, data=final2))*100 #4.000619
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.2, data=final2))*100 #4.022925
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.21, data=final2))*100 #
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.22, data=final2))*100 #
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.23, data=final2))*100 #
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.24, data=final2))*100 #
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.25, data=final2))*100 #4.026249
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.3, data=final2))*100 #4.008073
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.35, data=final2))*100 #3.966297
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.4, data=final2))*100 #3.899456
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.45, data=final2))*100 #3.806926
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.5, data=final2))*100 #3.704721
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.55, data=final2))*100 #3.547346
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.6, data=final2))*100 #3.384211
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.65, data=final2))*100 # 3.203051
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.7, data=final2))*100  #3.007917
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.75, data=final2))*100 #2.803241
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.8, data=final2))*100 #2.593527
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.85, data=final2))*100 #2.383062
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.9, data=final2))*100  #2.17569
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.95, data=final2))*100 #1.974655
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_1, data=final2))*100 #1.782522
#
saveRDS(final2, '~/height_prediction/loc_anc_analysis/output/all_PRS_WHI.Rds')

part_R2<-data.table(Alfa=c(0,0.05,0.1,0.15,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1),part_R2=c(
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0, data=final2)), partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.05, data=final2)), partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.1, data=final2)),
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.15, data=final2)), partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.2, data=final2)), partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.21, data=final2)),
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.22, data=final2)), partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.23, data=final2)), partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.24, data=final2)),
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.25, data=final2)), partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.26, data=final2)),
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.27, data=final2)), partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.28, data=final2)), partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.29, data=final2)),
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.3, data=final2)),partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.35, data=final2)),partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.4, data=final2)),
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.45, data=final2)),partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.5, data=final2)), partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.55, data=final2)),
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.6, data=final2)),partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.65, data=final2)), partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.7, data=final2)),
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.75, data=final2)),partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.8, data=final2)),partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.85, data=final2)),
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.9, data=final2)),partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0.95, data=final2)),partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_1, data=final2))))

saveRDS(part_R2, '~/height_prediction/loc_anc_analysis/output/part_R2_WHI.Rds')
###################
####
