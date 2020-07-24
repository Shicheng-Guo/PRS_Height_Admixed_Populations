#!/usr/bin/env Rscript
library(readr)
library(data.table)
library(reshape2)
library(asbio)
library(ggplot2)
library(dplyr)
library(parallel)
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

cat('checkpoint number 2\n')

res_all<-vector('list', 22)
for(I in 22:1){
	res_all[[I]]<-short_fun(args=I)
	saveRDS(res_all[[I]], file=paste0("~/height_prediction/gwas/WHI/output/chr", I, "_prs.Rds"))
}
saveRDS(res_all, file="~/height_prediction/gwas/WHI/output/all_prs.Rds")


res_all_HRS<-vector('list', 22)
for(I in 1:22){
        res_all_HRS[[I]]<-short_fun_v2(args=I)
        saveRDS(res_all_HRS[[I]], file=paste0("~/height_prediction/gwas/WHI/output/chr", I, "_prs_HRS.Rds"))
        cat('Chr ', I, ' done\n')
 }
saveRDS(res_all_HRS, file="~/height_prediction/gwas/WHI/output/all_prs_HRS.Rds")

res_all_JHS<-vector('list', 22)
for(I in 1:22){
        short_fun_v3(args=I)-> res_all_JHS[[I]]
        saveRDS(res_all_JHS[[I]], file=paste0("~/height_prediction/gwas/WHI/output/chr", I, "_prs_JHS.Rds"))
        cat('Chr ', I, ' done\n')
 }
saveRDS(res_all_JHS, file="~/height_prediction/gwas/WHI/output/all_prs_JHS.Rds")

readRDS("~/height_prediction/gwas/WHI/output/all_prs.Rds")-> res_all
names(res_all)<-c("POP1","POP2","ALL","ALL_Tstat1", "PLINK", "PLINK_Tstat1")
readRDS("~/height_prediction/gwas/WHI/output/all_prs_HRS.Rds")-> res_all_HRS
names(res_all_HRS)<-c("POP1","POP2","ALL","ALL_Tstat1", "PLINK", "PLINK_Tstat1")
readRDS("~/height_prediction/gwas/WHI/output/all_prs_JHS.Rds")-> res_all_JHS
names(res_all_JHS)<-c("POP1","POP2","ALL","ALL_Tstat1", "PLINK", "PLINK_Tstat1")

cat('checkpoint number 3\n')
data.table(SUBJID=names(res_all[[1]][[1]]), 
PRS_POP1=(unlist(res_all[[1]][[1]])+unlist(res_all[[2]][[1]])+unlist(res_all[[3]][[1]])+unlist(res_all[[4]][[1]])+unlist(res_all[[5]][[1]])+unlist(res_all[[6]][[1]])+unlist(res_all[[7]][[1]])+unlist(res_all[[8]][[1]])+unlist(res_all[[9]][[1]])+unlist(res_all[[10]][[1]])+unlist(res_all[[11]][[1]])+unlist(res_all[[12]][[1]])+unlist(res_all[[13]][[1]])+unlist(res_all[[14]][[1]])+unlist(res_all[[15]][[1]])+unlist(res_all[[16]][[1]])+unlist(res_all[[17]][[1]])+ unlist(res_all[[18]][[1]])+unlist(res_all[[19]][[1]])+unlist(res_all[[20]][[1]])+unlist(res_all[[21]][[1]])+unlist(res_all[[22]][[1]])), 
PRS_POP2=(unlist(res_all[[1]][[2]])+unlist(res_all[[2]][[2]])+unlist(res_all[[3]][[2]])+unlist(res_all[[4]][[2]])+unlist(res_all[[5]][[2]])+unlist(res_all[[6]][[2]])+unlist(res_all[[7]][[2]])+unlist(res_all[[8]][[2]])+unlist(res_all[[9]][[2]])+unlist(res_all[[10]][[2]])+unlist(res_all[[11]][[2]])+unlist(res_all[[12]][[2]])+unlist(res_all[[13]][[2]])+unlist(res_all[[14]][[2]])+unlist(res_all[[15]][[2]])+unlist(res_all[[16]][[2]])+unlist(res_all[[17]][[2]])+ unlist(res_all[[18]][[2]])+unlist(res_all[[19]][[2]])+unlist(res_all[[20]][[2]])+unlist(res_all[[21]][[2]])+unlist(res_all[[22]][[2]])), 
PRS_all=(unlist(res_all[[1]][[3]])+unlist(res_all[[2]][[3]])+unlist(res_all[[3]][[3]])+unlist(res_all[[4]][[3]])+unlist(res_all[[5]][[3]])+unlist(res_all[[6]][[3]])+unlist(res_all[[7]][[3]])+unlist(res_all[[8]][[3]])+unlist(res_all[[9]][[3]])+unlist(res_all[[10]][[3]])+unlist(res_all[[11]][[3]])+unlist(res_all[[12]][[3]])+unlist(res_all[[13]][[3]])+unlist(res_all[[14]][[3]])+unlist(res_all[[15]][[3]])+unlist(res_all[[16]][[3]])+unlist(res_all[[17]][[3]])+ unlist(res_all[[18]][[3]])+unlist(res_all[[19]][[3]])+unlist(res_all[[20]][[3]])+unlist(res_all[[21]][[3]])+unlist(res_all[[22]][[3]])), 
PRS_all_tstat_1=(unlist(res_all[[1]][[4]])+unlist(res_all[[2]][[4]])+unlist(res_all[[3]][[4]])+unlist(res_all[[4]][[4]])+unlist(res_all[[5]][[4]])+unlist(res_all[[6]][[4]])+unlist(res_all[[7]][[4]])+unlist(res_all[[8]][[4]])+unlist(res_all[[9]][[4]])+unlist(res_all[[10]][[4]])+unlist(res_all[[11]][[4]])+unlist(res_all[[12]][[4]])+unlist(res_all[[13]][[4]])+unlist(res_all[[14]][[4]])+unlist(res_all[[15]][[4]])+unlist(res_all[[16]][[4]])+unlist(res_all[[17]][[4]])+ unlist(res_all[[18]][[4]])+unlist(res_all[[19]][[4]])+unlist(res_all[[20]][[4]])+unlist(res_all[[21]][[4]])+unlist(res_all[[22]][[4]])),
PRS_plink=(unlist(res_all[[1]][[5]])+unlist(res_all[[2]][[5]])+unlist(res_all[[3]][[5]])+unlist(res_all[[4]][[5]])+unlist(res_all[[5]][[5]])+unlist(res_all[[6]][[5]])+unlist(res_all[[7]][[5]])+unlist(res_all[[8]][[5]])+unlist(res_all[[9]][[5]])+unlist(res_all[[10]][[5]])+unlist(res_all[[11]][[5]])+unlist(res_all[[12]][[5]])+unlist(res_all[[13]][[5]])+unlist(res_all[[14]][[5]])+unlist(res_all[[15]][[5]])+unlist(res_all[[16]][[5]])+unlist(res_all[[17]][[5]])+ unlist(res_all[[18]][[5]])+unlist(res_all[[19]][[5]])+unlist(res_all[[20]][[5]])+unlist(res_all[[21]][[5]])+unlist(res_all[[22]][[5]])),
PRS_plink_tstat_1=(unlist(res_all[[1]][[6]])+unlist(res_all[[2]][[6]])+unlist(res_all[[3]][[6]])+unlist(res_all[[4]][[6]])+unlist(res_all[[5]][[6]])+unlist(res_all[[6]][[6]])+unlist(res_all[[7]][[6]])+unlist(res_all[[8]][[6]])+unlist(res_all[[9]][[6]])+unlist(res_all[[10]][[6]])+unlist(res_all[[11]][[6]])+unlist(res_all[[12]][[6]])+unlist(res_all[[13]][[6]])+unlist(res_all[[14]][[6]])+unlist(res_all[[15]][[6]])+unlist(res_all[[16]][[6]])+unlist(res_all[[17]][[6]])+ unlist(res_all[[18]][[6]])+unlist(res_all[[19]][[6]])+unlist(res_all[[20]][[6]])+unlist(res_all[[21]][[6]])+unlist(res_all[[22]][[6]])),
PRS_EUR=unlist(readRDS('~/height_prediction/gwas/WHI/output/PGS_WHI_phys_100000_0.0005.Rds')))-> a

a[, PRS_POP1:=scale(PRS_POP1)]
a[, PRS_POP2:=scale(PRS_POP2)]
a[, PRS_all:=scale(PRS_all)]
a[, PRS_EUR:=scale(PRS_EUR)]
a[, PRS_plink:=scale(PRS_plink)]
a[, PRS_plink_tstat_1:=scale(PRS_plink_tstat_1)]
a[, PRS_all_tstat_1:=scale(PRS_all_tstat_1)]

remove(res_all)
gc()
data.table(SUBJID=names(res_all_HRS[[1]][[1]]),
PRS_POP1=(unlist(res_all_HRS[[1]][[1]])+unlist(res_all_HRS[[2]][[1]])+unlist(res_all_HRS[[3]][[1]])+unlist(res_all_HRS[[4]][[1]])+unlist(res_all_HRS[[5]][[1]])+unlist(res_all_HRS[[6]][[1]])+unlist(res_all_HRS[[7]][[1]])+unlist(res_all_HRS[[8]][[1]])+unlist(res_all_HRS[[9]][[1]])+unlist(res_all_HRS[[10]][[1]])+unlist(res_all_HRS[[11]][[1]])+unlist(res_all_HRS[[12]][[1]])+unlist(res_all_HRS[[13]][[1]])+unlist(res_all_HRS[[14]][[1]])+unlist(res_all_HRS[[15]][[1]])+unlist(res_all_HRS[[16]][[1]])+unlist(res_all_HRS[[17]][[1]])+ unlist(res_all_HRS[[18]][[1]])+unlist(res_all_HRS[[19]][[1]])+unlist(res_all_HRS[[20]][[1]])+unlist(res_all_HRS[[21]][[1]])+unlist(res_all_HRS[[22]][[1]])),
PRS_POP2=(unlist(res_all_HRS[[1]][[2]])+unlist(res_all_HRS[[2]][[2]])+unlist(res_all_HRS[[3]][[2]])+unlist(res_all_HRS[[4]][[2]])+unlist(res_all_HRS[[5]][[2]])+unlist(res_all_HRS[[6]][[2]])+unlist(res_all_HRS[[7]][[2]])+unlist(res_all_HRS[[8]][[2]])+unlist(res_all_HRS[[9]][[2]])+unlist(res_all_HRS[[10]][[2]])+unlist(res_all_HRS[[11]][[2]])+unlist(res_all_HRS[[12]][[2]])+unlist(res_all_HRS[[13]][[2]])+unlist(res_all_HRS[[14]][[2]])+unlist(res_all_HRS[[15]][[2]])+unlist(res_all_HRS[[16]][[2]])+unlist(res_all_HRS[[17]][[2]])+ unlist(res_all_HRS[[18]][[2]])+unlist(res_all_HRS[[19]][[2]])+unlist(res_all_HRS[[20]][[2]])+unlist(res_all_HRS[[21]][[2]])+unlist(res_all_HRS[[22]][[2]])),
PRS_all=(unlist(res_all_HRS[[1]][[3]])+unlist(res_all_HRS[[2]][[3]])+unlist(res_all_HRS[[3]][[3]])+unlist(res_all_HRS[[4]][[3]])+unlist(res_all_HRS[[5]][[3]])+unlist(res_all_HRS[[6]][[3]])+unlist(res_all_HRS[[7]][[3]])+unlist(res_all_HRS[[8]][[3]])+unlist(res_all_HRS[[9]][[3]])+unlist(res_all_HRS[[10]][[3]])+unlist(res_all_HRS[[11]][[3]])+unlist(res_all_HRS[[12]][[3]])+unlist(res_all_HRS[[13]][[3]])+unlist(res_all_HRS[[14]][[3]])+unlist(res_all_HRS[[15]][[3]])+unlist(res_all_HRS[[16]][[3]])+unlist(res_all_HRS[[17]][[3]])+ unlist(res_all_HRS[[18]][[3]])+unlist(res_all_HRS[[19]][[3]])+unlist(res_all_HRS[[20]][[3]])+unlist(res_all_HRS[[21]][[3]])+unlist(res_all_HRS[[22]][[3]])),
PRS_all_tstat_1=(unlist(res_all_HRS[[1]][[4]])+unlist(res_all_HRS[[2]][[4]])+unlist(res_all_HRS[[3]][[4]])+unlist(res_all_HRS[[4]][[4]])+unlist(res_all_HRS[[5]][[4]])+unlist(res_all_HRS[[6]][[4]])+unlist(res_all_HRS[[7]][[4]])+unlist(res_all_HRS[[8]][[4]])+unlist(res_all_HRS[[9]][[4]])+unlist(res_all_HRS[[10]][[4]])+unlist(res_all_HRS[[11]][[4]])+unlist(res_all_HRS[[12]][[4]])+unlist(res_all_HRS[[13]][[4]])+unlist(res_all_HRS[[14]][[4]])+unlist(res_all_HRS[[15]][[4]])+unlist(res_all_HRS[[16]][[4]])+unlist(res_all_HRS[[17]][[4]])+ unlist(res_all_HRS[[18]][[4]])+unlist(res_all_HRS[[19]][[4]])+unlist(res_all_HRS[[20]][[4]])+unlist(res_all_HRS[[21]][[4]])+unlist(res_all_HRS[[22]][[4]])),
PRS_plink=(unlist(res_all_HRS[[1]][[5]])+unlist(res_all_HRS[[2]][[5]])+unlist(res_all_HRS[[3]][[5]])+unlist(res_all_HRS[[4]][[5]])+unlist(res_all_HRS[[5]][[5]])+unlist(res_all_HRS[[6]][[5]])+unlist(res_all_HRS[[7]][[5]])+unlist(res_all_HRS[[8]][[5]])+unlist(res_all_HRS[[9]][[5]])+unlist(res_all_HRS[[10]][[5]])+unlist(res_all_HRS[[11]][[5]])+unlist(res_all_HRS[[12]][[5]])+unlist(res_all_HRS[[13]][[5]])+unlist(res_all_HRS[[14]][[5]])+unlist(res_all_HRS[[15]][[5]])+unlist(res_all_HRS[[16]][[5]])+unlist(res_all_HRS[[17]][[5]])+ unlist(res_all_HRS[[18]][[5]])+unlist(res_all_HRS[[19]][[5]])+unlist(res_all_HRS[[20]][[5]])+unlist(res_all_HRS[[21]][[5]])+unlist(res_all_HRS[[22]][[5]])),
PRS_plink_tstat_1=(unlist(res_all_HRS[[1]][[6]])+unlist(res_all_HRS[[2]][[6]])+unlist(res_all_HRS[[3]][[6]])+unlist(res_all_HRS[[4]][[6]])+unlist(res_all_HRS[[5]][[6]])+unlist(res_all_HRS[[6]][[6]])+unlist(res_all_HRS[[7]][[6]])+unlist(res_all_HRS[[8]][[6]])+unlist(res_all_HRS[[9]][[6]])+unlist(res_all_HRS[[10]][[6]])+unlist(res_all_HRS[[11]][[6]])+unlist(res_all_HRS[[12]][[6]])+unlist(res_all_HRS[[13]][[6]])+unlist(res_all_HRS[[14]][[6]])+unlist(res_all_HRS[[15]][[6]])+unlist(res_all_HRS[[16]][[6]])+unlist(res_all_HRS[[17]][[6]])+ unlist(res_all_HRS[[18]][[6]])+unlist(res_all_HRS[[19]][[6]])+unlist(res_all_HRS[[20]][[6]])+unlist(res_all_HRS[[21]][[6]])+unlist(res_all_HRS[[22]][[6]])),PRS_EUR=unlist(readRDS('~/height_prediction/gwas/HRS_afr/output/PGS_HRS_afr_phys_100000_0.0005.Rds')))-> a_HRS

a_HRS[, PRS_POP1:=scale(PRS_POP1)]
a_HRS[, PRS_POP2:=scale(PRS_POP2)]
a_HRS[, PRS_all:=scale(PRS_all)]
a_HRS[, PRS_EUR:=scale(PRS_EUR)]
a_HRS[, PRS_plink:=scale(PRS_plink)]
a_HRS[, PRS_plink_tstat_1:=scale(PRS_plink_tstat_1)]
a_HRS[, PRS_all_tstat_1:=scale(PRS_all_tstat_1)]

remove(res_all_HRS)
gc()
data.table(SUBJID=names(res_all_JHS[[1]][[1]]),
PRS_POP1=(unlist(res_all_JHS[[1]][[1]])+unlist(res_all_JHS[[2]][[1]])+unlist(res_all_JHS[[3]][[1]])+unlist(res_all_JHS[[4]][[1]])+unlist(res_all_JHS[[5]][[1]])+unlist(res_all_JHS[[6]][[1]])+unlist(res_all_JHS[[7]][[1]])+unlist(res_all_JHS[[8]][[1]])+unlist(res_all_JHS[[9]][[1]])+unlist(res_all_JHS[[10]][[1]])+unlist(res_all_JHS[[11]][[1]])+unlist(res_all_JHS[[12]][[1]])+unlist(res_all_JHS[[13]][[1]])+unlist(res_all_JHS[[14]][[1]])+unlist(res_all_JHS[[15]][[1]])+unlist(res_all_JHS[[16]][[1]])+unlist(res_all_JHS[[17]][[1]])+ unlist(res_all_JHS[[18]][[1]])+unlist(res_all_JHS[[19]][[1]])+unlist(res_all_JHS[[20]][[1]])+unlist(res_all_JHS[[21]][[1]])+unlist(res_all_JHS[[22]][[1]])),
PRS_POP2=(unlist(res_all_JHS[[1]][[2]])+unlist(res_all_JHS[[2]][[2]])+unlist(res_all_JHS[[3]][[2]])+unlist(res_all_JHS[[4]][[2]])+unlist(res_all_JHS[[5]][[2]])+unlist(res_all_JHS[[6]][[2]])+unlist(res_all_JHS[[7]][[2]])+unlist(res_all_JHS[[8]][[2]])+unlist(res_all_JHS[[9]][[2]])+unlist(res_all_JHS[[10]][[2]])+unlist(res_all_JHS[[11]][[2]])+unlist(res_all_JHS[[12]][[2]])+unlist(res_all_JHS[[13]][[2]])+unlist(res_all_JHS[[14]][[2]])+unlist(res_all_JHS[[15]][[2]])+unlist(res_all_JHS[[16]][[2]])+unlist(res_all_JHS[[17]][[2]])+ unlist(res_all_JHS[[18]][[2]])+unlist(res_all_JHS[[19]][[2]])+unlist(res_all_JHS[[20]][[2]])+unlist(res_all_JHS[[21]][[2]])+unlist(res_all_JHS[[22]][[2]])),
PRS_all=(unlist(res_all_JHS[[1]][[3]])+unlist(res_all_JHS[[2]][[3]])+unlist(res_all_JHS[[3]][[3]])+unlist(res_all_JHS[[4]][[3]])+unlist(res_all_JHS[[5]][[3]])+unlist(res_all_JHS[[6]][[3]])+unlist(res_all_JHS[[7]][[3]])+unlist(res_all_JHS[[8]][[3]])+unlist(res_all_JHS[[9]][[3]])+unlist(res_all_JHS[[10]][[3]])+unlist(res_all_JHS[[11]][[3]])+unlist(res_all_JHS[[12]][[3]])+unlist(res_all_JHS[[13]][[3]])+unlist(res_all_JHS[[14]][[3]])+unlist(res_all_JHS[[15]][[3]])+unlist(res_all_JHS[[16]][[3]])+unlist(res_all_JHS[[17]][[3]])+ unlist(res_all_JHS[[18]][[3]])+unlist(res_all_JHS[[19]][[3]])+unlist(res_all_JHS[[20]][[3]])+unlist(res_all_JHS[[21]][[3]])+unlist(res_all_JHS[[22]][[3]])),
PRS_all_tstat_1=(unlist(res_all_JHS[[1]][[4]])+unlist(res_all_JHS[[2]][[4]])+unlist(res_all_JHS[[3]][[4]])+unlist(res_all_JHS[[4]][[4]])+unlist(res_all_JHS[[5]][[4]])+unlist(res_all_JHS[[6]][[4]])+unlist(res_all_JHS[[7]][[4]])+unlist(res_all_JHS[[8]][[4]])+unlist(res_all_JHS[[9]][[4]])+unlist(res_all_JHS[[10]][[4]])+unlist(res_all_JHS[[11]][[4]])+unlist(res_all_JHS[[12]][[4]])+unlist(res_all_JHS[[13]][[4]])+unlist(res_all_JHS[[14]][[4]])+unlist(res_all_JHS[[15]][[4]])+unlist(res_all_JHS[[16]][[4]])+unlist(res_all_JHS[[17]][[4]])+ unlist(res_all_JHS[[18]][[4]])+unlist(res_all_JHS[[19]][[4]])+unlist(res_all_JHS[[20]][[4]])+unlist(res_all_JHS[[21]][[4]])+unlist(res_all_JHS[[22]][[4]])),
PRS_plink=(unlist(res_all_JHS[[1]][[5]])+unlist(res_all_JHS[[2]][[5]])+unlist(res_all_JHS[[3]][[5]])+unlist(res_all_JHS[[4]][[5]])+unlist(res_all_JHS[[5]][[5]])+unlist(res_all_JHS[[6]][[5]])+unlist(res_all_JHS[[7]][[5]])+unlist(res_all_JHS[[8]][[5]])+unlist(res_all_JHS[[9]][[5]])+unlist(res_all_JHS[[10]][[5]])+unlist(res_all_JHS[[11]][[5]])+unlist(res_all_JHS[[12]][[5]])+unlist(res_all_JHS[[13]][[5]])+unlist(res_all_JHS[[14]][[5]])+unlist(res_all_JHS[[15]][[5]])+unlist(res_all_JHS[[16]][[5]])+unlist(res_all_JHS[[17]][[5]])+ unlist(res_all_JHS[[18]][[5]])+unlist(res_all_JHS[[19]][[5]])+unlist(res_all_JHS[[20]][[5]])+unlist(res_all_JHS[[21]][[5]])+unlist(res_all_JHS[[22]][[5]])),
PRS_plink_tstat_1=(unlist(res_all_JHS[[1]][[6]])+unlist(res_all_JHS[[2]][[6]])+unlist(res_all_JHS[[3]][[6]])+unlist(res_all_JHS[[4]][[6]])+unlist(res_all_JHS[[5]][[6]])+unlist(res_all_JHS[[6]][[6]])+unlist(res_all_JHS[[7]][[6]])+unlist(res_all_JHS[[8]][[6]])+unlist(res_all_JHS[[9]][[6]])+unlist(res_all_JHS[[10]][[6]])+unlist(res_all_JHS[[11]][[6]])+unlist(res_all_JHS[[12]][[6]])+unlist(res_all_JHS[[13]][[6]])+unlist(res_all_JHS[[14]][[6]])+unlist(res_all_JHS[[15]][[6]])+unlist(res_all_JHS[[16]][[6]])+unlist(res_all_JHS[[17]][[6]])+ unlist(res_all_JHS[[18]][[6]])+unlist(res_all_JHS[[19]][[6]])+unlist(res_all_JHS[[20]][[6]])+unlist(res_all_JHS[[21]][[6]])+unlist(res_all_JHS[[22]][[6]])),
PRS_EUR=unlist(readRDS('~/height_prediction/gwas/JHS/output/PGS_JHS_phys_100000_0.0005.Rds')))-> a_JHS

a_JHS[, PRS_POP1:=scale(PRS_POP1)]
a_JHS[, PRS_POP2:=scale(PRS_POP2)]
a_JHS[, PRS_all:=scale(PRS_all)]
a_JHS[, PRS_EUR:=scale(PRS_EUR)]
a_JHS[, PRS_plink:=scale(PRS_plink)]
a_JHS[, PRS_plink_tstat_1:=scale(PRS_plink_tstat_1)]
a_JHS[, PRS_all_tstat_1:=scale(PRS_all_tstat_1)]


#phenotype
fread('~/height_prediction/input/WHI/WHI_phenotypes.txt')-> Pheno_WHI

Pheno_WHI[, SUBJID:=paste0("0_", as.character(Pheno_WHI[, SUBJID]))]
setkey(Pheno_WHI, SUBJID)
#add ancestry
ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/WHI/rfmix_anc_chr', X, '.txt'))))
anc_WHI<-ancestry %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual
setkey(anc_WHI, SUBJID)
##
setkey(a, SUBJID)
setkey(Pheno_WHI, SUBJID)
a[Pheno_WHI][anc_WHI]-> final

#HRS
fread('~/height_prediction/input/HRS_afr/HRS_AFR_phenotypes.txt', fill=T)[,1:4]-> Pheno_HRS
Pheno_HRS[, SUBJID:=paste0(ID, "_", ID)]
setkey(Pheno_HRS, SUBJID)
#add ancestry
ancestry_HRS<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/HRS_afr/rfmix_anc_chr', X, '.txt'))))
anc_HRS<-ancestry_HRS %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual
anc_HRS[, SUBJID:=paste0(SUBJID, "_", SUBJID)]
setkey(anc_HRS, SUBJID)
##
setkey(a_HRS, SUBJID)
setkey(Pheno_HRS, SUBJID)
a_HRS[Pheno_HRS][anc_HRS]-> final_HRS
final_HRS[,SEX:=ifelse(SEX==1, "Male", "Female")]
final_HRS[, HEIGHTX:=100*HEIGHT]
##JHS

fread('~/height_prediction/input/JHS/JHS_phenotypes.txt', fill=T)-> Pheno_JHS
#Pheno_JHS[, SUBJID:=paste0("0_", SUBJID, "_", FAM)]
setkey(Pheno_JHS, SUBJID)
#add ancestry
ancestry_JHS<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/JHS/rfmix_anc_chr', X, '.txt'))))
anc_JHS<-ancestry_JHS %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual
#anc_JHS[,SUBJID:=gsub(":", "_", SUBJID)]
anc_JHS$SUBJID<-substr(anc_JHS[,SUBJID],3,9)
setkey(anc_JHS, SUBJID)
##
a_JHS[,SUBJID:=substr(a_JHS[,SUBJID],1,7)]
setkey(a_JHS, SUBJID)
setkey(Pheno_JHS, SUBJID)
a_JHS[Pheno_JHS, nomatch=0][anc_JHS, nomatch=0]-> final_JHS
PRS_POP2=(unlist(res_all[[1]][[2]])+unlist(res_all[[2]][[2]])+unlist(res_all[[3]][[2]])+unlist(res_all[[4]][[2]])+unlist(res_all[[5]][[2]])+unlist(res_all[[6]][[2]])+unlist(res_all[[7]][[2]])+unlist(res_all[[8]][[2]])+unlist(res_all[[9]][[2]])+unlist(res_all[[10]][[2]])+unlist(res_all[[11]][[2]])+unlist(res_all[[12]][[2]])+unlist(res_all[[13]][[2]])+unlist(res_all[[14]][[2]])+unlist(res_all[[15]][[2]])+unlist(res_all[[16]][[2]])+unlist(res_all[[17]][[2]])+ unlist(res_all[[18]][[2]])+unlist(res_all[[19]][[2]])+unlist(res_all[[20]][[2]])+unlist(res_all[[21]][[2]])+unlist(res_all[[22]][[2]])), 
PRS_all=(unlist(res_all[[1]][[3]])+unlist(res_all[[2]][[3]])+unlist(res_all[[3]][[3]])+unlist(res_all[[4]][[3]])+unlist(res_all[[5]][[3]])+unlist(res_all[[6]][[3]])+unlist(res_all[[7]][[3]])+unlist(res_all[[8]][[3]])+unlist(res_all[[9]][[3]])+unlist(res_all[[10]][[3]])+unlist(res_all[[11]][[3]])+unlist(res_all[[12]][[3]])+unlist(res_all[[13]][[3]])+unlist(res_all[[14]][[3]])+unlist(res_all[[15]][[3]])+unlist(res_all[[16]][[3]])+unlist(res_all[[17]][[3]])+ unlist(res_all[[18]][[3]])+unlist(res_all[[19]][[3]])+unlist(res_all[[20]][[3]])+unlist(res_all[[21]][[3]])+unlist(res_all[[22]][[3]])), 
PRS_all_tstat_1=(unlist(res_all[[1]][[4]])+unlist(res_all[[2]][[4]])+unlist(res_all[[3]][[4]])+unlist(res_all[[4]][[4]])+unlist(res_all[[5]][[4]])+unlist(res_all[[6]][[4]])+unlist(res_all[[7]][[4]])+unlist(res_all[[8]][[4]])+unlist(res_all[[9]][[4]])+unlist(res_all[[10]][[4]])+unlist(res_all[[11]][[4]])+unlist(res_all[[12]][[4]])+unlist(res_all[[13]][[4]])+unlist(res_all[[14]][[4]])+unlist(res_all[[15]][[4]])+unlist(res_all[[16]][[4]])+unlist(res_all[[17]][[4]])+ unlist(res_all[[18]][[4]])+unlist(res_all[[19]][[4]])+unlist(res_all[[20]][[4]])+unlist(res_all[[21]][[4]])+unlist(res_all[[22]][[4]])),
PRS_plink=(unlist(res_all[[1]][[5]])+unlist(res_all[[2]][[5]])+unlist(res_all[[3]][[5]])+unlist(res_all[[4]][[5]])+unlist(res_all[[5]][[5]])+unlist(res_all[[6]][[5]])+unlist(res_all[[7]][[5]])+unlist(res_all[[8]][[5]])+unlist(res_all[[9]][[5]])+unlist(res_all[[10]][[5]])+unlist(res_all[[11]][[5]])+unlist(res_all[[12]][[5]])+unlist(res_all[[13]][[5]])+unlist(res_all[[14]][[5]])+unlist(res_all[[15]][[5]])+unlist(res_all[[16]][[5]])+unlist(res_all[[17]][[5]])+ unlist(res_all[[18]][[5]])+unlist(res_all[[19]][[5]])+unlist(res_all[[20]][[5]])+unlist(res_all[[21]][[5]])+unlist(res_all[[22]][[5]])),
PRS_plink_tstat_1=(unlist(res_all[[1]][[6]])+unlist(res_all[[2]][[6]])+unlist(res_all[[3]][[6]])+unlist(res_all[[4]][[6]])+unlist(res_all[[5]][[6]])+unlist(res_all[[6]][[6]])+unlist(res_all[[7]][[6]])+unlist(res_all[[8]][[6]])+unlist(res_all[[9]][[6]])+unlist(res_all[[10]][[6]])+unlist(res_all[[11]][[6]])+unlist(res_all[[12]][[6]])+unlist(res_all[[13]][[6]])+unlist(res_all[[14]][[6]])+unlist(res_all[[15]][[6]])+unlist(res_all[[16]][[6]])+unlist(res_all[[17]][[6]])+ unlist(res_all[[18]][[6]])+unlist(res_all[[19]][[6]])+unlist(res_all[[20]][[6]])+unlist(res_all[[21]][[6]])+unlist(res_all[[22]][[6]])),
PRS_EUR=unlist(readRDS('~/height_prediction/gwas/WHI/output/PGS_WHI_phys_100000_0.0005.Rds')))-> a

a[, PRS_POP1:=scale(PRS_POP1)]
a[, PRS_POP2:=scale(PRS_POP2)]
a[, PRS_all:=scale(PRS_all)]
a[, PRS_EUR:=scale(PRS_EUR)]
a[, PRS_plink:=scale(PRS_plink)]
a[, PRS_plink_tstat_1:=scale(PRS_plink_tstat_1)]
a[, PRS_all_tstat_1:=scale(PRS_all_tstat_1)]

remove(res_all)
gc()
data.table(SUBJID=names(res_all_HRS[[1]][[1]]),
PRS_POP1=(unlist(res_all_HRS[[1]][[1]])+unlist(res_all_HRS[[2]][[1]])+unlist(res_all_HRS[[3]][[1]])+unlist(res_all_HRS[[4]][[1]])+unlist(res_all_HRS[[5]][[1]])+unlist(res_all_HRS[[6]][[1]])+unlist(res_all_HRS[[7]][[1]])+unlist(res_all_HRS[[8]][[1]])+unlist(res_all_HRS[[9]][[1]])+unlist(res_all_HRS[[10]][[1]])+unlist(res_all_HRS[[11]][[1]])+unlist(res_all_HRS[[12]][[1]])+unlist(res_all_HRS[[13]][[1]])+unlist(res_all_HRS[[14]][[1]])+unlist(res_all_HRS[[15]][[1]])+unlist(res_all_HRS[[16]][[1]])+unlist(res_all_HRS[[17]][[1]])+ unlist(res_all_HRS[[18]][[1]])+unlist(res_all_HRS[[19]][[1]])+unlist(res_all_HRS[[20]][[1]])+unlist(res_all_HRS[[21]][[1]])+unlist(res_all_HRS[[22]][[1]])),
PRS_POP2=(unlist(res_all_HRS[[1]][[2]])+unlist(res_all_HRS[[2]][[2]])+unlist(res_all_HRS[[3]][[2]])+unlist(res_all_HRS[[4]][[2]])+unlist(res_all_HRS[[5]][[2]])+unlist(res_all_HRS[[6]][[2]])+unlist(res_all_HRS[[7]][[2]])+unlist(res_all_HRS[[8]][[2]])+unlist(res_all_HRS[[9]][[2]])+unlist(res_all_HRS[[10]][[2]])+unlist(res_all_HRS[[11]][[2]])+unlist(res_all_HRS[[12]][[2]])+unlist(res_all_HRS[[13]][[2]])+unlist(res_all_HRS[[14]][[2]])+unlist(res_all_HRS[[15]][[2]])+unlist(res_all_HRS[[16]][[2]])+unlist(res_all_HRS[[17]][[2]])+ unlist(res_all_HRS[[18]][[2]])+unlist(res_all_HRS[[19]][[2]])+unlist(res_all_HRS[[20]][[2]])+unlist(res_all_HRS[[21]][[2]])+unlist(res_all_HRS[[22]][[2]])),
PRS_all=(unlist(res_all_HRS[[1]][[3]])+unlist(res_all_HRS[[2]][[3]])+unlist(res_all_HRS[[3]][[3]])+unlist(res_all_HRS[[4]][[3]])+unlist(res_all_HRS[[5]][[3]])+unlist(res_all_HRS[[6]][[3]])+unlist(res_all_HRS[[7]][[3]])+unlist(res_all_HRS[[8]][[3]])+unlist(res_all_HRS[[9]][[3]])+unlist(res_all_HRS[[10]][[3]])+unlist(res_all_HRS[[11]][[3]])+unlist(res_all_HRS[[12]][[3]])+unlist(res_all_HRS[[13]][[3]])+unlist(res_all_HRS[[14]][[3]])+unlist(res_all_HRS[[15]][[3]])+unlist(res_all_HRS[[16]][[3]])+unlist(res_all_HRS[[17]][[3]])+ unlist(res_all_HRS[[18]][[3]])+unlist(res_all_HRS[[19]][[3]])+unlist(res_all_HRS[[20]][[3]])+unlist(res_all_HRS[[21]][[3]])+unlist(res_all_HRS[[22]][[3]])),
PRS_all_tstat_1=(unlist(res_all_HRS[[1]][[4]])+unlist(res_all_HRS[[2]][[4]])+unlist(res_all_HRS[[3]][[4]])+unlist(res_all_HRS[[4]][[4]])+unlist(res_all_HRS[[5]][[4]])+unlist(res_all_HRS[[6]][[4]])+unlist(res_all_HRS[[7]][[4]])+unlist(res_all_HRS[[8]][[4]])+unlist(res_all_HRS[[9]][[4]])+unlist(res_all_HRS[[10]][[4]])+unlist(res_all_HRS[[11]][[4]])+unlist(res_all_HRS[[12]][[4]])+unlist(res_all_HRS[[13]][[4]])+unlist(res_all_HRS[[14]][[4]])+unlist(res_all_HRS[[15]][[4]])+unlist(res_all_HRS[[16]][[4]])+unlist(res_all_HRS[[17]][[4]])+ unlist(res_all_HRS[[18]][[4]])+unlist(res_all_HRS[[19]][[4]])+unlist(res_all_HRS[[20]][[4]])+unlist(res_all_HRS[[21]][[4]])+unlist(res_all_HRS[[22]][[4]])),
PRS_plink=(unlist(res_all_HRS[[1]][[5]])+unlist(res_all_HRS[[2]][[5]])+unlist(res_all_HRS[[3]][[5]])+unlist(res_all_HRS[[4]][[5]])+unlist(res_all_HRS[[5]][[5]])+unlist(res_all_HRS[[6]][[5]])+unlist(res_all_HRS[[7]][[5]])+unlist(res_all_HRS[[8]][[5]])+unlist(res_all_HRS[[9]][[5]])+unlist(res_all_HRS[[10]][[5]])+unlist(res_all_HRS[[11]][[5]])+unlist(res_all_HRS[[12]][[5]])+unlist(res_all_HRS[[13]][[5]])+unlist(res_all_HRS[[14]][[5]])+unlist(res_all_HRS[[15]][[5]])+unlist(res_all_HRS[[16]][[5]])+unlist(res_all_HRS[[17]][[5]])+ unlist(res_all_HRS[[18]][[5]])+unlist(res_all_HRS[[19]][[5]])+unlist(res_all_HRS[[20]][[5]])+unlist(res_all_HRS[[21]][[5]])+unlist(res_all_HRS[[22]][[5]])),
PRS_plink_tstat_1=(unlist(res_all_HRS[[1]][[6]])+unlist(res_all_HRS[[2]][[6]])+unlist(res_all_HRS[[3]][[6]])+unlist(res_all_HRS[[4]][[6]])+unlist(res_all_HRS[[5]][[6]])+unlist(res_all_HRS[[6]][[6]])+unlist(res_all_HRS[[7]][[6]])+unlist(res_all_HRS[[8]][[6]])+unlist(res_all_HRS[[9]][[6]])+unlist(res_all_HRS[[10]][[6]])+unlist(res_all_HRS[[11]][[6]])+unlist(res_all_HRS[[12]][[6]])+unlist(res_all_HRS[[13]][[6]])+unlist(res_all_HRS[[14]][[6]])+unlist(res_all_HRS[[15]][[6]])+unlist(res_all_HRS[[16]][[6]])+unlist(res_all_HRS[[17]][[6]])+ unlist(res_all_HRS[[18]][[6]])+unlist(res_all_HRS[[19]][[6]])+unlist(res_all_HRS[[20]][[6]])+unlist(res_all_HRS[[21]][[6]])+unlist(res_all_HRS[[22]][[6]])),PRS_EUR=unlist(readRDS('~/height_prediction/gwas/HRS_afr/output/PGS_HRS_afr_phys_100000_0.0005.Rds')))-> a_HRS

a_HRS[, PRS_POP1:=scale(PRS_POP1)]
a_HRS[, PRS_POP2:=scale(PRS_POP2)]
a_HRS[, PRS_all:=scale(PRS_all)]
a_HRS[, PRS_EUR:=scale(PRS_EUR)]
a_HRS[, PRS_plink:=scale(PRS_plink)]
a_HRS[, PRS_plink_tstat_1:=scale(PRS_plink_tstat_1)]
a_HRS[, PRS_all_tstat_1:=scale(PRS_all_tstat_1)]

remove(res_all_HRS)
gc()
data.table(SUBJID=names(res_all_JHS[[1]][[1]]),
PRS_POP1=(unlist(res_all_JHS[[1]][[1]])+unlist(res_all_JHS[[2]][[1]])+unlist(res_all_JHS[[3]][[1]])+unlist(res_all_JHS[[4]][[1]])+unlist(res_all_JHS[[5]][[1]])+unlist(res_all_JHS[[6]][[1]])+unlist(res_all_JHS[[7]][[1]])+unlist(res_all_JHS[[8]][[1]])+unlist(res_all_JHS[[9]][[1]])+unlist(res_all_JHS[[10]][[1]])+unlist(res_all_JHS[[11]][[1]])+unlist(res_all_JHS[[12]][[1]])+unlist(res_all_JHS[[13]][[1]])+unlist(res_all_JHS[[14]][[1]])+unlist(res_all_JHS[[15]][[1]])+unlist(res_all_JHS[[16]][[1]])+unlist(res_all_JHS[[17]][[1]])+ unlist(res_all_JHS[[18]][[1]])+unlist(res_all_JHS[[19]][[1]])+unlist(res_all_JHS[[20]][[1]])+unlist(res_all_JHS[[21]][[1]])+unlist(res_all_JHS[[22]][[1]])),
PRS_POP2=(unlist(res_all_JHS[[1]][[2]])+unlist(res_all_JHS[[2]][[2]])+unlist(res_all_JHS[[3]][[2]])+unlist(res_all_JHS[[4]][[2]])+unlist(res_all_JHS[[5]][[2]])+unlist(res_all_JHS[[6]][[2]])+unlist(res_all_JHS[[7]][[2]])+unlist(res_all_JHS[[8]][[2]])+unlist(res_all_JHS[[9]][[2]])+unlist(res_all_JHS[[10]][[2]])+unlist(res_all_JHS[[11]][[2]])+unlist(res_all_JHS[[12]][[2]])+unlist(res_all_JHS[[13]][[2]])+unlist(res_all_JHS[[14]][[2]])+unlist(res_all_JHS[[15]][[2]])+unlist(res_all_JHS[[16]][[2]])+unlist(res_all_JHS[[17]][[2]])+ unlist(res_all_JHS[[18]][[2]])+unlist(res_all_JHS[[19]][[2]])+unlist(res_all_JHS[[20]][[2]])+unlist(res_all_JHS[[21]][[2]])+unlist(res_all_JHS[[22]][[2]])),
PRS_all=(unlist(res_all_JHS[[1]][[3]])+unlist(res_all_JHS[[2]][[3]])+unlist(res_all_JHS[[3]][[3]])+unlist(res_all_JHS[[4]][[3]])+unlist(res_all_JHS[[5]][[3]])+unlist(res_all_JHS[[6]][[3]])+unlist(res_all_JHS[[7]][[3]])+unlist(res_all_JHS[[8]][[3]])+unlist(res_all_JHS[[9]][[3]])+unlist(res_all_JHS[[10]][[3]])+unlist(res_all_JHS[[11]][[3]])+unlist(res_all_JHS[[12]][[3]])+unlist(res_all_JHS[[13]][[3]])+unlist(res_all_JHS[[14]][[3]])+unlist(res_all_JHS[[15]][[3]])+unlist(res_all_JHS[[16]][[3]])+unlist(res_all_JHS[[17]][[3]])+ unlist(res_all_JHS[[18]][[3]])+unlist(res_all_JHS[[19]][[3]])+unlist(res_all_JHS[[20]][[3]])+unlist(res_all_JHS[[21]][[3]])+unlist(res_all_JHS[[22]][[3]])),
PRS_all_tstat_1=(unlist(res_all_JHS[[1]][[4]])+unlist(res_all_JHS[[2]][[4]])+unlist(res_all_JHS[[3]][[4]])+unlist(res_all_JHS[[4]][[4]])+unlist(res_all_JHS[[5]][[4]])+unlist(res_all_JHS[[6]][[4]])+unlist(res_all_JHS[[7]][[4]])+unlist(res_all_JHS[[8]][[4]])+unlist(res_all_JHS[[9]][[4]])+unlist(res_all_JHS[[10]][[4]])+unlist(res_all_JHS[[11]][[4]])+unlist(res_all_JHS[[12]][[4]])+unlist(res_all_JHS[[13]][[4]])+unlist(res_all_JHS[[14]][[4]])+unlist(res_all_JHS[[15]][[4]])+unlist(res_all_JHS[[16]][[4]])+unlist(res_all_JHS[[17]][[4]])+ unlist(res_all_JHS[[18]][[4]])+unlist(res_all_JHS[[19]][[4]])+unlist(res_all_JHS[[20]][[4]])+unlist(res_all_JHS[[21]][[4]])+unlist(res_all_JHS[[22]][[4]])),
PRS_plink=(unlist(res_all_JHS[[1]][[5]])+unlist(res_all_JHS[[2]][[5]])+unlist(res_all_JHS[[3]][[5]])+unlist(res_all_JHS[[4]][[5]])+unlist(res_all_JHS[[5]][[5]])+unlist(res_all_JHS[[6]][[5]])+unlist(res_all_JHS[[7]][[5]])+unlist(res_all_JHS[[8]][[5]])+unlist(res_all_JHS[[9]][[5]])+unlist(res_all_JHS[[10]][[5]])+unlist(res_all_JHS[[11]][[5]])+unlist(res_all_JHS[[12]][[5]])+unlist(res_all_JHS[[13]][[5]])+unlist(res_all_JHS[[14]][[5]])+unlist(res_all_JHS[[15]][[5]])+unlist(res_all_JHS[[16]][[5]])+unlist(res_all_JHS[[17]][[5]])+ unlist(res_all_JHS[[18]][[5]])+unlist(res_all_JHS[[19]][[5]])+unlist(res_all_JHS[[20]][[5]])+unlist(res_all_JHS[[21]][[5]])+unlist(res_all_JHS[[22]][[5]])),
PRS_plink_tstat_1=(unlist(res_all_JHS[[1]][[6]])+unlist(res_all_JHS[[2]][[6]])+unlist(res_all_JHS[[3]][[6]])+unlist(res_all_JHS[[4]][[6]])+unlist(res_all_JHS[[5]][[6]])+unlist(res_all_JHS[[6]][[6]])+unlist(res_all_JHS[[7]][[6]])+unlist(res_all_JHS[[8]][[6]])+unlist(res_all_JHS[[9]][[6]])+unlist(res_all_JHS[[10]][[6]])+unlist(res_all_JHS[[11]][[6]])+unlist(res_all_JHS[[12]][[6]])+unlist(res_all_JHS[[13]][[6]])+unlist(res_all_JHS[[14]][[6]])+unlist(res_all_JHS[[15]][[6]])+unlist(res_all_JHS[[16]][[6]])+unlist(res_all_JHS[[17]][[6]])+ unlist(res_all_JHS[[18]][[6]])+unlist(res_all_JHS[[19]][[6]])+unlist(res_all_JHS[[20]][[6]])+unlist(res_all_JHS[[21]][[6]])+unlist(res_all_JHS[[22]][[6]])),
PRS_EUR=unlist(readRDS('~/height_prediction/gwas/JHS/output/PGS_JHS_phys_100000_0.0005.Rds')))-> a_JHS

a_JHS[, PRS_POP1:=scale(PRS_POP1)]
a_JHS[, PRS_POP2:=scale(PRS_POP2)]
a_JHS[, PRS_all:=scale(PRS_all)]
a_JHS[, PRS_EUR:=scale(PRS_EUR)]
final2_HRS_JHS<-rbind(final2_HRS[, Dt:='HRS'], final2_JHS[,Dt:='JHS'])
#ok, so for ukb_afr and WHI (both), POP1 is AFR
cat('another checkpoint\n')

final2[,AGE2:=AGE^2][, PRS_eur_afr_local:=(mean(EUR_ANC)*(PRS_EUR))+(mean(AFR_ANC)*PRS_all)]
final2[,AGE2:=AGE^2][, PRS_eur_afr_plink:=(mean(EUR_ANC)*(PRS_EUR))+(mean(AFR_ANC)*PRS_plink)]
final2_HRS_JHS[,AGE2:=AGE^2][, PRS_eur_afr_local:=(mean(EUR_ANC)*(PRS_EUR))+(mean(AFR_ANC)*PRS_all)]
final2_HRS_JHS[,AGE2:=AGE^2][, PRS_eur_afr_plink:=(mean(EUR_ANC)*(PRS_EUR))+(mean(AFR_ANC)*PRS_plink)]

fwrite(final2, file="~/height_prediction/gwas/WHI/output/all_prs_whi.txt", sep="\t")
fwrite(final2_HRS, file="~/height_prediction/gwas/HRS_afr/output/all_prs_hrs.txt", sep="\t")
fwrite(final2_JHS, file="~/height_prediction/gwas/JHS/output/all_prs_jsh.txt", sep="\t")
fwrite(final2_HRS_JHS, file="~/height_prediction/gwas/JHS/output/all_prs_hrs_jsh.txt", sep="\t")

partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2))*100 #4.1%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_all, data=final2))*100 #0.35%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_plink, data=final2))*100 #0.39%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_POP1, data=final2))*100 #0.09 %
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_POP2, data=final2))*100 #0.1896%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_plink_tstat_1, data=final2))*100 #0.34%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_all_tstat_1, data=final2))*100 #0.31%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_eur_afr_local, data=final2))*100 #0.937%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_eur_afr_plink, data=final2))*100 # #0.613%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2),lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR+PRS_eur_afr_plink, data=final2))*100 #0.0715%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2),lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR+PRS_plink, data=final2))*100 #0.071
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2),lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR+PRS_eur_afr_local, data=final2))*100 #0.052
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2),lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR+PRS_all, data=final2))*100 #0.052
#
partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+Dt+EUR_ANC, data=final2_HRS_JHS), lm(HEIGHTX~SEX+AGE+AGE2+Dt+EUR_ANC+PRS_EUR, data=final2_HRS_JHS))*100 #2.93

###################
##################a
my_alpha<-function(alpha=0.9, dt='final2'){
	my_dt<-get(dt)
	my_dt[,PRS_comb:=((1-alpha)*PRS_EUR)+(alpha*PRS_all)]
	res2<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=my_dt), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_comb,data=my_dt ))
	return(res2)
}
my_alpha_conc<-function(alpha=0.9, dt='final2_HRS_JHS'){
	my_dt<-get(dt)
        my_dt[,PRS_comb:=((1-alpha)*PRS_EUR)+(alpha*PRS_all)]
	res2<-partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+Dt+EUR_ANC, data=my_dt), lm(HEIGHTX~SEX+AGE+AGE2+Dt+EUR_ANC+PRS_comb, data=my_dt))
	return(res2)	
}

a_vec<-seq(0,1,0.001)
lapply(a_vec, function(X) my_alpha(alpha=X))-> wei_PRS
lapply(a_vec, function(X) my_alpha_conc(alpha=X))-> wei_PRS_conc
rbind(data.table(part_R2=unlist(wei_PRS), alfa=a_vec, Dataset='WHI'), data.table(part_R2=unlist(wei_PRS_conc), alfa=a_vec, Dataset='JHS+HRS'))-> wanna_plot
#

ggplot(wanna_plot, aes(x=alfa, y=part_R2, colour=Dataset)) + geom_point(size=0.8)+ geom_line() + 
geom_vline(xintercept=optimize(my_alpha, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, col='red', lty=2) + 
geom_vline(xintercept=optimize(my_alpha_conc, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, col='red', lty=2) +
geom_hline(yintercept=optimize(my_alpha, interval=c(0,1), maximum=T, tol = 0.0001)$objective, col='red', lty=2) +
geom_hline(yintercept=optimize(my_alpha_conc, interval=c(0,1), maximum=T, tol = 0.0001)$objective, col='red', lty=2) +
coord_cartesian(ylim = c(0, 0.048), xlim=c(0,0.6))
ggsave('~/height_prediction/gwas/WHI/figs/alfa_all.pdf')
optimize(my_alpha, interval=c(0,1), maximum=T, tol = 0.0001) #one liner for max 0.07690677 for 0.0415066
optimize(my_alpha_conc, interval=c(0,1), maximum=T, tol = 0.0001) #one liner for max  0.1191175 0.02988471
###################
##################
my_alpha_v2<-function(alpha=0.9, dt='final2'){
geom_vline(xintercept=optimize(my_alpha, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, col='red', lty=2) + 
geom_vline(xintercept=optimize(my_alpha_conc, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, col='red', lty=2) +
geom_hline(yintercept=optimize(my_alpha, interval=c(0,1), maximum=T, tol = 0.0001)$objective, col='red', lty=2) +
geom_hline(yintercept=optimize(my_alpha_conc, interval=c(0,1), maximum=T, tol = 0.0001)$objective, col='red', lty=2) +
coord_cartesian(ylim = c(0, 0.048), xlim=c(0,0.6))
ggsave('~/height_prediction/gwas/WHI/figs/alfa_all.pdf')
optimize(my_alpha, interval=c(0,1), maximum=T, tol = 0.0001) #one liner for max 0.07690677 for 0.0415066
optimize(my_alpha_conc, interval=c(0,1), maximum=T, tol = 0.0001) #one liner for max  0.1191175 0.02988471
###################
##################
my_alpha_v2<-function(alpha=0.9, dt='final2'){
        my_dt<-get(dt)
        my_dt[,PRS_comb:=((1-alpha)*PRS_EUR)+(alpha*PRS_plink)]
        res2<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=my_dt), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_comb,data=my_dt ))
        return(res2)
}
my_alpha_conc_v2<-function(alpha=0.9, dt='final2_HRS_JHS'){
        my_dt<-get(dt)
        my_dt[,PRS_comb:=((1-alpha)*PRS_EUR)+(alpha*PRS_plink)]
        res2<-partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+Dt+EUR_ANC, data=my_dt), lm(HEIGHTX~SEX+AGE+AGE2+Dt+EUR_ANC+PRS_comb, data=my_dt))
        return(res2)
}
lapply(a_vec, function(X) my_alpha_v2(alpha=X))-> wei_PRS_v2
lapply(a_vec, function(X) my_alpha_conc_v2(alpha=X))-> wei_PRS_conc_v2
rbind(data.table(part_R2=unlist(wei_PRS_v2), alfa=a_vec, Dataset='WHI'), data.table(part_R2=unlist(wei_PRS_conc_v2), alfa=a_vec, Dataset='JHS+HRS'))-> wanna_plot_v2
#
ggplot(wanna_plot_v2, aes(x=alfa, y=part_R2, colour=Dataset)) + geom_point(size=0.8)+ geom_line() +
geom_vline(xintercept=optimize(my_alpha_v2, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, col='red', lty=2) +
geom_vline(xintercept=optimize(my_alpha_conc_v2, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, col='red', lty=2) +
geom_hline(yintercept=optimize(my_alpha_v2, interval=c(0,1), maximum=T, tol = 0.0001)$objective, col='red', lty=2) +
geom_hline(yintercept=optimize(my_alpha_conc_v2, interval=c(0,1), maximum=T, tol = 0.0001)$objective, col='red', lty=2) +
coord_cartesian(ylim = c(0, 0.048), xlim=c(0,0.6))
ggsave('~/height_prediction/gwas/WHI/figs/alfa_plink.pdf')
optimize(my_alpha_v2, interval=c(0,1), maximum=T, tol = 0.0001) #one liner for max 0.03613721 // 0.04169
optimize(my_alpha_conc_v2, interval=c(0,1), maximum=T, tol = 0.0001) #one liner for max  0.1012799 0.02973544

##################
max_alpha<-function(alpha, dt='final2'){
	my_dt<-get(dt)
	my_dt[,PRS_comb:=(AFR_ANC*alpha*PRS_all)+((1-(alpha*AFR_ANC))*PRS_EUR)]
	res<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=my_dt), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_comb, data=my_dt))
	return(res)
}

max_alpha_conc<-function(alpha, dt='final2_HRS_JHS'){
        my_dt<-get(dt)
        my_dt[,PRS_comb:=(AFR_ANC*alpha*PRS_all)+((1-(alpha*AFR_ANC))*PRS_EUR)]
        res<-partial.R2(lm(HEIGHTX~Dt+SEX+AGE+AGE2+EUR_ANC, data=my_dt), lm(HEIGHTX~Dt+SEX+AGE+AGE2+EUR_ANC+PRS_comb, data=my_dt))
        return(res)
}

a_vec2<-seq(from=0, to=1, by=0.001)
all_prs<-lapply(a_vec, function(X) max_alpha(alpha=X))
all_prs_conc<-lapply(a_vec, function(X) max_alpha_conc(alpha=X))

optimize(max_alpha, interval=c(0,1), maximum=T, tol = 0.0001) # 0.1064444 #0.04157011
optimize(max_alpha_conc, interval=c(0,1), maximum=T, tol = 0.0001) #0.1367484 #0.02981409

temp_dt<-rbind(data.table(part_R2=unlist(all_prs), alfa=a_vec, Dataset='WHI'), data.table(part_R2=unlist(all_prs_conc), alfa=a_vec, Dataset='JHS+HRS'))

saveRDS(temp_dt, file='~/height_prediction/gwas/WHI/output/temp_dt.Rds')
ggplot(temp_dt, aes(x=alfa, y=part_R2, colour=Dataset)) + geom_point(size=0.8)+ geom_line() +
geom_vline(xintercept=optimize(max_alpha, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, col='red', lty=2) +
geom_vline(xintercept=optimize(max_alpha_conc, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, col='red', lty=2) +
geom_hline(yintercept=optimize(max_alpha, interval=c(0,1), maximum=T, tol = 0.0001)$objective, col='red', lty=2) +
geom_hline(yintercept=optimize(max_alpha_conc, interval=c(0,1), maximum=T, tol = 0.0001)$objective, col='red', lty=2) +
coord_cartesian(ylim = c(0, 0.048), xlim=c(0,0.6))
ggsave('~/height_prediction/gwas/WHI/figs/alfa_LA_withAnc.pdf')

###################
max_alpha_v2<-function(alpha, dt='final2'){
        my_dt<-get(dt)
        my_dt[,PRS_comb:=(AFR_ANC*alpha*PRS_plink)+((1-(alpha*AFR_ANC))*PRS_EUR)]
        res<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=my_dt), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_comb, data=my_dt))
        return(res)
}

max_alpha_conc_v2<-function(alpha, dt='final2_HRS_JHS'){
        my_dt<-get(dt)
        my_dt[,PRS_comb:=(AFR_ANC*alpha*PRS_plink)+((1-(alpha*AFR_ANC))*PRS_EUR)]
        res<-partial.R2(lm(HEIGHTX~Dt+SEX+AGE+AGE2+EUR_ANC, data=my_dt), lm(HEIGHTX~Dt+SEX+AGE+AGE2+EUR_ANC+PRS_comb, data=my_dt))
        return(res)
}

all_prs_v2<-lapply(a_vec, function(X) max_alpha_v2(alpha=X))
all_prs_conc_v2<-lapply(a_vec, function(X) max_alpha_conc_v2(alpha=X))
geom_line(size=1.2) + coord_cartesian(xlim=c(0,0.5), ylim=c(0.02, 0.042)) + labs(x="Alpha", y=expression(Partial~R^2)) +
annotate("segment", x = optimize(my_alpha, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, xend = optimize(my_alpha, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, y = 0.04, yend = optimize(my_alpha, interval=c(0,1), maximum=T, tol = 0.0001)$objective,colour='gray', size=1, alpha=0.6, arrow=arrow()) +
annotate("segment", x = optimize(max_alpha, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, xend = optimize(max_alpha, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, y = 0.043, yend = optimize(max_alpha, interval=c(0,1), maximum=T, tol = 0.0001)$objective, colour='gray', size=1, alpha=0.6, arrow=arrow(), lty=2) +
annotate("segment", x = optimize(my_alpha_conc, interval=c(0,1), maximum=T, tol = 0.0001)$maximum,, xend = optimize(my_alpha_conc, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, y = 0.027, yend=optimize(my_alpha_conc, interval=c(0,1), maximum=T, tol = 0.0001)$objective, colour='gray', size=1, alpha=0.6, arrow=arrow()) +
annotate("segment", x = optimize(max_alpha_conc, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, xend = optimize(max_alpha_conc, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, y = 0.03, yend = optimize(max_alpha_conc, interval=c(0,1), maximum=T, tol = 0.0001)$objective, colour = "gray", size=1, alpha=0.6, arrow=arrow(), lty=2) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.key=element_blank(), legend.background=element_blank(),legend.title=element_blank()) + labs(x="Alpha", y=expression(Partial~R^2))

ggsave('~/height_prediction/gwas/WHI/figs/test_this_all.pdf')

ggplot(test5, aes(x=alpha, y=value, colour=variable)) + facet_wrap(~Dataset)+
scale_colour_manual(values=c("#96a8b2", "#101010")) +
geom_line(size=1.2) + coord_cartesian(xlim=c(0,0.5), ylim=c(0.02, 0.042)) + labs(x="Alpha", y=expression(Partial~R^2)) +
#annotate("segment", x = optimize(my_alpha_v2, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, xend = optimize(my_alpha_v2, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, y = 0.04, yend = optimize(my_alpha_v2, interval=c(0,1), maximum=T, tol = 0.0001)$objective,colour='gray', size=1, alpha=0.6, arrow=arrow()) +
#annotate("segment", x = optimize(max_alpha_v2, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, xend = optimize(max_alpha_v2, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, y = 0.043, yend = optimize(max_alpha_v2, interval=c(0,1), maximum=T, tol = 0.0001)$objective, colour='gray', size=1, alpha=0.6, arrow=arrow(), lty=2) +
#annotate("segment", x = optimize(my_alpha_conc_v2, interval=c(0,1), maximum=T, tol = 0.0001)$maximum,, xend = optimize(my_alpha_conc_v2, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, y = 0.027, yend=optimize(my_alpha_conc_v2, interval=c(0,1), maximum=T, tol = 0.0001)$objective, colour='gray', size=1, alpha=0.6, arrow=arrow()) +
#annotate("segment", x = optimize(max_alpha_conc_v2, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, xend = optimize(max_alpha_conc_v2, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, y = 0.03, yend = optimize(max_alpha_conc_v2, interval=c(0,1), maximum=T, tol = 0.0001)$objective, colour = "gray", size=1, alpha=0.6, arrow=arrow(), lty=2) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.key=element_blank(), legend.background=element_blank(),legend.title=element_blank(), legend.text=element_text(size=12)) + labs(x="Alpha", y=expression(Partial~R^2))

ggsave('~/height_prediction/gwas/WHI/figs/test_this_plink.pdf')

