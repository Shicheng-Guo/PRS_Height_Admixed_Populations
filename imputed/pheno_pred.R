#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#**************************************
#*      CALCULATE PARTIAL R2  *********
#**************************************
library("optparse")
library(data.table)
library(dplyr)
library(biomaRt)
library(parallel)
options(scipen=10)
options(digits=10)
library("optparse")
library(ggplot2);
library(reshape2); library(wesanderson)
library(rlist)
library(asbio)
library(GGally)
library(tidyr)
library(hexbin)
library(psychometric)
library(boot)

#read in PGS scores
PGS_HRS_afr<-vector('list', 35)
nam<-paste(c(5000,10000,25000,50000,75000,100000,500000), c(0.0005,0.00005, 0.000005,0.0000005,0.00000005), sep="_")
names(PGS_HRS_afr)<-nam
for(N in nam){
        readRDS(paste0('~/height_prediction/imputed/output/PGS2_', nam, '.Rds'))[[1]]-> PGS_HRS_afr[[N]]
}

#read in phenotype data
fread('~/height_prediction/input/HRS_afr/HRS_AFR_phenotypes.txt', fill=T)[,V5:=NULL]-> Pheno_HRS_afr
#a partial R2 function
source('~/height_prediction/strat_prs/scripts/Rsq_R2.R')


#add PGS to Pheno table in order to be able to make multiple analyses
#Pheno_HRS_afr[,ID:=paste0(ID, "_", ID)]
for (N in nam){
	for(j in 1:length(PGS_HRS_afr)){
        	gsub("[0-9]+_","",names(PGS_HRS_afr[[N]])[[j]])-> names(PGS_HRS_afr[[N]])[[j]]
	}

	as.character(Pheno_HRS_afr[[N]]$ID)-> Pheno_HRS_afr[[N]]$ID
	setkey(Pheno_HRS_afr[[N]], ID)
}
#add ancestry
ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/HRS_afr/rfmix_anc_chr', X, '.txt'))))

anc_HRS_afr<-ancestry %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual

anc_HRS_afr[,ID:=SUBJID][,SUBJID:=NULL]
as.character(anc_HRS_afr$ID)-> anc_HRS_afr$ID
setkey(anc_HRS_afr, ID)

PGS2_HRS_afr<-vector('list', length(PGS_HRS_afr))
names(PGS2_HRS_afr)<-nam

for(N in nam){
	data.table(ID=names(PGS_HRS_afr[[N]]), PGS=unlist(PGS_HRS_afr[[N]]))-> PGS2_HRS_afr[[N]]
	setkey(PGS2_HRS_afr[[N]], ID)
	PGS2_HRS_afr[[N]][Pheno_HRS_afr, nomatch=0]-> PGS2_HRS_afr[[N]]
	PGS2_HRS_afr[[N]][anc_HRS_afr, nomatch=0]-> PGS2_HRS_afr[[N]]
	PGS2_HRS_afr[[N]][,AGE2:=AGE^2]
	PGS2_HRS_afr[[N]][AFR_ANC>=0.05]-> PGS2_HRS_afr[[N]]
	PGS2_HRS_afr[[N]][which(!is.na(PGS2_HRS_afr[,HEIGHT])),]-> PGS2_HRS_afr[[N]]
}
lapply(nam, function(X) lm(HEIGHT~SEX, X))-> lm0_HRS_afr
lapply(nam, function(X) lm(HEIGHT~PGS, X))-> lm1_HRS_afr
lapply(nam, function(X) lm(HEIGHT~AGE, X))-> lm2_HRS_afr
lapply(nam, function(X) lm(HEIGHT~AGE2, X))-> lm3_HRS_afr
lapply(nam, function(X) lm(HEIGHT~EUR_ANC, X))-> lm4_HRS_afr
lapply(nam, function(X) lm(HEIGHT~PGS+AGE, X))-> lm5_HRS_afr
lapply(nam, function(X) lm(HEIGHT~PGS+AGE2,X))-> lm6_HRS_afr
lapply(nam, function(X) lm(HEIGHT~SEX+AGE+AGE2+EUR_ANC, X))-> lm7_HRS_afr
lapply(nam, function(X) lm(HEIGHT~SEX+AGE+AGE2+EUR_ANC+PGS, X))-> lm8_HRS_afr

partial_R2<-lapply(nam, function(X) partial.R2(lm7_HRS_afr[[X]], lm8_HRS_afr[[X]])
names(partial_R2)<-nam


##now HRS_eur
PGS_HRS_afr<-vector('list', 35)
names(PGS_HRS_eur)<-nam
for(N in nam){
	readRDS(paste0('~/height_prediction/imputed/output/PGS2', nam, '.Rds'))[[2]]-> PGS_HRS_eur[[N]]
}
#read in phenotype data
fread('~/height_prediction/input/HRS_eur/HRS_EUR_phenotypes.txt')-> Pheno_HRS_eur
###########
##############

#add PGS to Pheno table in order to be able to make multiple analyses

Pheno_HRS_eur[,ID:=paste0(ID, "_", ID)]
setkey(Pheno_HRS_eur, ID)

PGS2_HRS_eur<-vector('list', 35)
names(PGS2_HRS_eur)<-nam

for(N in nam)){
	data.table(ID=names(PGS_HRS_eur[[N]]), PGS=unlist(PGS_HRS_eur[[N]]))-> PGS2_HRS_eur[[N]]
	setkey(PGS2_HRS_eur[[N]], ID)
	PGS2_HRS_eur[[N]][Pheno_HRS_eur, nomatch=0]-> PGS2_HRS_eur[[N]]
	PGS2_HRS_eur[[N]][,EUR_ANC:=1]
	PGS2_HRS_eur[[N]][,AGE2:=AGE^2]
	PGS2_HRS_eur[[N]][,HEIGHT:=HEIGHT*100]
	PGS2_HRS_eur[[N]]$SEX<-as.factor(PGS2_HRS_eur[[N]]$SEX)
}

lapply(nam, function(X) lm(HEIGHT~SEX, X))-> lm1_HRS_eur
lapply(nam, function(X) lm(HEIGHT~PGS, X))-> lm2_HRS_eur
lapply(nam, function(X) lm(HEIGHT~AGE, X))-> lm3_HRS_eur
lapply(nam, function(X) lm(HEIGHT~AGE2,X))-> lm4_HRS_eur
lapply(nam, function(X) lm(HEIGHT~PGS+AGE, X))-> lm5_HRS_eur
lapply(nam, function(X) lm(HEIGHT~PGS+AGE2,X))-> lm6_HRS_eur
lapply(nam, function(X) lm(HEIGHT~SEX+AGE+AGE2, X))-> lm7_HRS_eur
lapply(nam, function(X) lm(HEIGHT~SEX+AGE+AGE2+PGS, X))-> lm8_HRS_eur

partial_R2<-lapply(nam, function(X) partial.R2(lm7_HRS_eur[[X]],lm8_HRS_eur[[X]])) #
