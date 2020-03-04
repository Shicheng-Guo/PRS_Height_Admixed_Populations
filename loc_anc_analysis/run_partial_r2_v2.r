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
options(scipen=999)
source('~/height_prediction/scripts/mclapply2.R')
library(SOAR)
library(TeachingDemos)
##Load other PRSs
txtStart("~/height_prediction/loc_anc_analysis/README.txt")
readRDS("~/height_prediction/imputed/output/all_prs.Rds")-> res_all
names(res_all)<-c("PLINK", "PLINK_TSTAT_1")
data.table(SUBJID=names(res_all[[1]][[1]]), PRS_afr=(unlist(res_all[[1]][[1]])+unlist(res_all[[2]][[1]])+unlist(res_all[[3]][[1]])+unlist(res_all[[4]][[1]])+unlist(res_all[[5]][[1]])+unlist(res_all[[6]][[1]])+unlist(res_all[[7]][[1]])+unlist(res_all[[8]][[1]])+unlist(res_all[[9]][[1]])+unlist(res_all[[10]][[1]])+unlist(res_all[[11]][[1]])+unlist(res_all[[12]][[1]])+unlist(res_all[[13]][[1]])+unlist(res_all[[14]][[1]])+unlist(res_all[[15]][[1]])+unlist(res_all[[16]][[1]])+unlist(res_all[[17]][[1]])+ unlist(res_all[[18]][[1]])+unlist(res_all[[19]][[1]])+unlist(res_all[[20]][[1]])+unlist(res_all[[21]][[1]])+unlist(res_all[[22]][[1]])),PRS_EUR=unlist(readRDS('~/height_prediction/gwas/WHI/output/PGS_WHI_phys_100000_0.0005.Rds')))-> a

a[, PRS_afr:=scale(PRS_afr)]

#read in the prss and add them
prs<-vector('list', 101)
#names(prs)<-c("0", "0.05","0.1", "0.15", "0.2","0.21", "0.22", "0.23", "0.24","0.25","0.26", "0.27", "0.28", "0.29", "0.3","0.35","0.4","0.45","0.5","0.55", "0.6","0.65", "0.7", "0.75","0.8", "0.85","0.9","0.95","1")
names(prs)<-seq(from=0, to=1, by=0.01)
names(prs)[1]<-"0.00"
names(prs)[11]<-"0.10"
names(prs)[21]<-"0.20"
names(prs)[31]<-"0.30"
names(prs)[41]<-"0.40"
names(prs)[51]<-"0.50"
names(prs)[61]<-"0.60"
names(prs)[71]<-"0.70"
names(prs)[81]<-"0.80"
names(prs)[91]<-"0.90"
names(prs)[101]<-"1.00"
for (A in names(prs)){
	prs[[A]]<-vector('list',22)
		for(chr in 1:22){
		prs[[A]][[chr]]<-readRDS(paste0('~/height_prediction/loc_anc_analysis/output/chr', chr,'_',A, '_prs_WHI.Rds'))
	}
}
samp_names<-names(prs[[1]][[1]])
cols_names<-paste0("PRS_", names(prs))

dt<-data.table(SUBJID=samp_names)
#dt[, (cols_names):=NA]

prs2<-vector('list', 101)
names(prs2)<-names(prs)
for (A in names(prs2)){
	prs2[[A]]<-vector('list', length(prs[[1]][[1]]))
	names(prs2[[A]])<-gsub("_A", "",samp_names)
}

for (A in names(prs2)){
	for (S in names(prs2[[A]])){
	prs2[[A]][[S]]<-sum(prs[[A]][[1]][[S]],prs[[A]][[2]][[S]],prs[[A]][[3]][[S]],prs[[A]][[4]][[S]],prs[[A]][[5]][[S]],prs[[A]][[6]][[S]], prs[[A]][[7]][[S]], prs[[A]][[8]][[S]], prs[[A]][[9]][[S]], prs[[A]][[10]][[S]],prs[[A]][[11]][[S]], prs[[A]][[12]][[S]], 
	prs[[A]][[13]][[S]],prs[[A]][[14]][[S]], prs[[A]][[15]][[S]], prs[[A]][[16]][[S]], prs[[A]][[17]][[S]], prs[[A]][[18]][[S]], prs[[A]][[19]][[S]], prs[[A]][[20]][[S]], prs[[A]][[21]][[S]],prs[[A]][[22]][[S]], na.rm=T)
	cat(S, '\r')
	}
	cat(A, '\n')
        gc()
}

prs3<-lapply(prs2, function(X) unlist(X))
names(prs3)<-cols_names
dt2<-dt #backup
for(I in cols_names){	
	dt[,(I):=scale(prs3[[I]])]
}


merge(dt,a, by="SUBJID")-> dt
#phenotype
fread('~/height_prediction/input/WHI/WHI_phenotypes.txt')-> Pheno_WHI

Pheno_WHI[, SUBJID:=paste0("0_", as.character(Pheno_WHI[, SUBJID]))]
setkey(Pheno_WHI, SUBJID)
#add ancestry
ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/WHI/rfmix_anc_chr', X, '.txt'))))
anc_WHI<-ancestry %>% dplyr::group_by(SUBJID) %>% dplyr::summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual
setkey(anc_WHI, SUBJID)
##

#PRS_EUR<-data.table(PRS_EUR=unlist(readRDS('~/height_prediction/gwas/WHI/output/PGS_WHI_phys_100000_0.0005.Rds')),SUBJID=names(unlist(readRDS('~/height_prediction/gwas/WHI/output/PGS_WHI_phys_100000_0.0005.Rds'))))
#a<-merge(dt, PRS_EUR, by="SUBJID")
setkey(dt, SUBJID)
setkey(Pheno_WHI, SUBJID)
dt[Pheno_WHI][anc_WHI]-> final
final2<-select(final, cols_names, PRS_EUR, PRS_afr,AGE, EUR_ANC,HEIGHTX, SUBJID) %>% as.data.table
final2[,AGE2:=AGE^2]

#
colnames(final2)[1:103]-> my_cols
part_r2<-vector('list', 103)
names(part_r2)<- my_cols
for(I in my_cols){
	part_r2[[I]]<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+get(I), data=final2))*100
}
	
saveRDS(final2, '~/height_prediction/loc_anc_analysis/output/all_PRS_WHI.Rds')

saveRDS(part_r2, '~/height_prediction/loc_anc_analysis/output/part_R2_WHI.Rds')
txtStop()
###################
####
