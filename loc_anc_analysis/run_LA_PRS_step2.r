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


#read in the prss and add them

prs<-vector('list',22)
for(chr in 1:22){
	prs[[chr]]<-readRDS(paste0('~/height_prediction/loc_anc_analysis/output/chr', chr, '_prs.Rds'))
}


a<-seq(from=1, to=length(prs[[1]]), by=2)

prs2<-vector('list', 22)

for(chr in 1:22){
	prs2[[chr]]<-vector('list', length(a))
	names(prs2[[chr]])<-a
}

for(chr in 1:22){
	for(i in a){
		prs2[[chr]][[i]]<-prs[[chr]][[i]]+prs[[chr]][[i+1]]
	}
	prs2[[chr]]<-unlist(prs2[[chr]])
}

prs3<-vector('list', length(prs2[[1]]))

for(S in 1:length(prs2[[1]])){
	prs3[[S]]<-sum(prs2[[1]][[S]],prs2[[2]][[S]],prs2[[3]][[S]],prs2[[4]][[S]],prs2[[5]][[S]]+prs2[[6]][[S]], prs2[[7]][[S]], prs2[[8]][[S]], prs2[[9]][[S]], prs2[[10]][[S]],prs2[[11]][[S]], prs2[[12]][[S]], prs2[[13]][[S]],prs2[[14]][[S]], prs2[[15]][[S]], prs2[[16]][[S]], prs2[[17]][[S]], prs2[[18]][[S]], prs2[[19]][[S]], prs2[[20]][[S]], prs2[[21]][[S]],prs2[[22]][[S]], na.rm=T)
}

chr<-20
input_MATRIX<-readRDS(file=paste0('~/height_prediction/loc_anc_analysis/output/chr', chr, '_temp.Rds'))
names(prs3)<-gsub("_A", "", unlist(lapply(input_MATRIX, function(X) colnames(X)[[1]])))



#phenotype
fread('~/height_prediction/input/WHI/WHI_phenotypes.txt')-> Pheno_WHI

Pheno_WHI[, SUBJID:=paste0("0_", as.character(Pheno_WHI[, SUBJID]))]
setkey(Pheno_WHI, SUBJID)
#add ancestry
ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/WHI/rfmix_anc_chr', X, '.txt'))))
anc_WHI<-ancestry %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual
setkey(anc_WHI, SUBJID)
##
a<-data.table(SUBJID=names(prs3), PRS_LocAnc=unlist(prs3), PRS_EUR=unlist(readRDS('~/height_prediction/gwas/WHI/output/PGS_WHI_phys_100000_0.0005.Rds')))
setkey(a, SUBJID)
setkey(Pheno_WHI, SUBJID)
a[Pheno_WHI][anc_WHI]-> final
final[,c("SUBJID", "PRS_LocAnc","PRS_EUR", "SEX", "HEIGHTX", "EUR_ANC","AFR_ANC","AGE")]-> final2
final2[,PRS_LocAnc:=scale(PRS_LocAnc)]
final2[,AGE2:=AGE^2][, PRS_eur_afr_LocAnc:=(mean(EUR_ANC)*(PRS_EUR))+(mean(AFR_ANC)*PRS_LocAnc)]

#


partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2))*100 #4.1%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_LocAnc, data=final2))*100 #0.061 %
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_eur_afr_LocAnc, data=final2))*100 # #0.08896 %
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2),lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR+PRS_eur_afr_LocAnc, data=final2))*100 #0.01268%
#


###################
##################a
a_vec<-seq(0,1,0.001)
##################
my_alpha_v2<-function(alpha=0.9, dt='final2'){
        my_dt<-get(dt)
        my_dt[,PRS_comb:=((1-alpha)*PRS_EUR)+(alpha*PRS_LocAnc)]
        res2<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=my_dt), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_comb,data=my_dt ))
        return(res2)
}
my_alpha_conc_v2<-function(alpha=0.9, dt='final2_HRS_JHS'){
        my_dt<-get(dt)
        my_dt[,PRS_comb:=((1-alpha)*PRS_EUR)+(alpha*PRS_LocAnc)]
        res2<-partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+Dt+EUR_ANC, data=my_dt), lm(HEIGHTX~SEX+AGE+AGE2+Dt+EUR_ANC+PRS_comb, data=my_dt))
        return(res2)
}
lapply(a_vec, function(X) my_alpha_v2(alpha=X))-> wei_PRS_v2

#
optimize(my_alpha_v2, interval=c(0,1), maximum=T, tol = 0.0001) #one liner for max 0.004748 //0.04113066
##################
a_vec2<-seq(from=0, to=1, by=0.001)

###################
max_alpha_v2<-function(alpha, dt='final2'){
        my_dt<-get(dt)
        my_dt[,PRS_comb:=(AFR_ANC*alpha*PRS_LocAnc)+((1-(alpha*AFR_ANC))*PRS_EUR)]
        res<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=my_dt), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_comb, data=my_dt))
        return(res)
}

max_alpha_conc_v2<-function(alpha, dt='final2_HRS_JHS'){
        my_dt<-get(dt)
        my_dt[,PRS_comb:=(AFR_ANC*alpha*PRS_LocAnc)+((1-(alpha*AFR_ANC))*PRS_EUR)]
        res<-partial.R2(lm(HEIGHTX~Dt+SEX+AGE+AGE2+EUR_ANC, data=my_dt), lm(HEIGHTX~Dt+SEX+AGE+AGE2+EUR_ANC+PRS_comb, data=my_dt))
        return(res)
}

all_prs_v2<-lapply(a_vec, function(X) max_alpha_v2(alpha=X))

optimize(max_alpha_v2, interval=c(0,1), maximum=T, tol = 0.0001) # 



