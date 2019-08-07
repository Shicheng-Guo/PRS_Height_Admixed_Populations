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
PGS_UKB_afr<-vector('list', 35)
nam<-paste(rep(c(5000,10000,25000,50000,75000,100000,500000),5), c(0.0005,0.00005, 0.000005,0.0000005,0.00000005), sep="_")
names(PGS_UKB_afr)<-nam
for(N in nam){
        readRDS(paste0('~/height_prediction/imputed/output/PGS2', N, '_UKB.Rds'))[[1]]-> PGS_UKB_afr[[N]]
}

#read in phenotype data
fread('~/height_prediction/input/ukb_afr/UKB_AFR_pheno.txt', fill=T)[,ANC.PC:=NULL]-> Pheno_UKB_afr
#a partial R2 function
source('~/height_prediction/strat_prs/scripts/Rsq_R2.R')


#add PGS to Pheno table in order to be able to make multiple analyses
#Pheno_UKB_afr[,ID:=paste0(ID, "_", ID)]
for (N in nam){
	for(j in 1:length(PGS_UKB_afr[[N]])){
        	gsub("[0-9]+_","",names(PGS_UKB_afr[[N]])[[j]])-> names(PGS_UKB_afr[[N]])[[j]]
	}
}

as.character(Pheno_UKB_afr$ID)-> Pheno_UKB_afr$ID
setkey(Pheno_UKB_afr, ID)
#add ancestry
ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/ukb_afr/rfmix_anc_chr', X, '.txt'))))

anc_UKB_afr<-ancestry %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual

anc_UKB_afr[,ID:=SUBJID][,SUBJID:=NULL]
as.character(anc_UKB_afr$ID)-> anc_UKB_afr$ID
setkey(anc_UKB_afr, ID)

PGS2_UKB_afr<-vector('list', length(PGS_UKB_afr))
names(PGS2_UKB_afr)<-nam

for(N in nam){
	data.table(ID=names(PGS_UKB_afr[[N]]), PGS=unlist(PGS_UKB_afr[[N]]))-> PGS2_UKB_afr[[N]]
	setkey(PGS2_UKB_afr[[N]], ID)
	PGS2_UKB_afr[[N]][Pheno_UKB_afr, nomatch=0]-> PGS2_UKB_afr[[N]]
	PGS2_UKB_afr[[N]][anc_UKB_afr, nomatch=0]-> PGS2_UKB_afr[[N]]
	PGS2_UKB_afr[[N]][,AGE2:=Age^2]
	PGS2_UKB_afr[[N]][AFR_ANC>=0.05]-> PGS2_UKB_afr[[N]]
	PGS2_UKB_afr[[N]]$Sex<-as.factor(PGS2_UKB_afr[[N]]$Sex)
	PGS2_UKB_afr[[N]][which(!is.na(PGS2_UKB_afr[[N]][,Height])),]-> PGS2_UKB_afr[[N]]
}

lapply(PGS2_UKB_afr, function(X) lm(Height~Sex, X))-> lm0_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~PGS, X))-> lm1_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~Age, X))-> lm2_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~AGE2, X))-> lm3_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~EUR_ANC, X))-> lm4_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~PGS+Age, X))-> lm5_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~PGS+AGE2,X))-> lm6_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~Sex+Age+AGE2+EUR_ANC, X))-> lm7_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~Sex+Age+AGE2+EUR_ANC+PGS, X))-> lm8_UKB_afr

partial_R2<-lapply(nam, function(X) partial.R2(lm7_UKB_afr[[X]], lm8_UKB_afr[[X]]))
names(partial_R2)<-nam


##now UKB_eur
PGS_UKB_eur<-vector('list', 35)
names(PGS_UKB_eur)<-nam
for(N in nam){
	readRDS(paste0('~/height_prediction/imputed/output/PGS2', N, '_UKB.Rds'))[[2]]-> PGS_UKB_eur[[N]]
}
#read in phenotype data
fread('~/height_prediction/input/ukb_eur/UKB_EUR_pheno.txt')-> Pheno_UKB_eur
###########
##############

#add PGS to Pheno table in order to be able to make multiple analyses

Pheno_UKB_eur[,ID:=paste0(ID, "_", ID)]
setkey(Pheno_UKB_eur, ID)

PGS2_UKB_eur<-vector('list', 35)
names(PGS2_UKB_eur)<-nam

for(N in nam){
	data.table(ID=names(PGS_UKB_eur[[N]]), PGS=unlist(PGS_UKB_eur[[N]]))-> PGS2_UKB_eur[[N]]
	setkey(PGS2_UKB_eur[[N]], ID)
	PGS2_UKB_eur[[N]][Pheno_UKB_eur, nomatch=0]-> PGS2_UKB_eur[[N]]
	PGS2_UKB_eur[[N]][,EUR_ANC:=1]
	PGS2_UKB_eur[[N]][,AGE2:=Age^2]
	PGS2_UKB_eur[[N]][,Height:=Height*100]
	PGS2_UKB_eur[[N]]$Sex<-as.factor(PGS2_UKB_eur[[N]]$Sex)
}

lapply(PGS2_UKB_eur, function(X) lm(Height~Sex, X))-> lm1_UKB_eur
lapply(PGS2_UKB_eur, function(X) lm(Height~PGS, X))-> lm2_UKB_eur
lapply(PGS2_UKB_eur, function(X) lm(Height~Age, X))-> lm3_UKB_eur
lapply(PGS2_UKB_eur, function(X) lm(Height~AGE2,X))-> lm4_UKB_eur
lapply(PGS2_UKB_eur, function(X) lm(Height~PGS+Age, X))-> lm5_UKB_eur
lapply(PGS2_UKB_eur, function(X) lm(Height~PGS+AGE2,X))-> lm6_UKB_eur
lapply(PGS2_UKB_eur, function(X) lm(Height~Sex+Age+AGE2, X))-> lm7_UKB_eur
lapply(PGS2_UKB_eur, function(X) lm(Height~Sex+Age+AGE2+PGS, X))-> lm8_UKB_eur

partial_R2_eur<-lapply(nam, function(X) partial.R2(lm7_UKB_eur[[X]],lm8_UKB_eur[[X]])) #
names(partial_R2_eur)<-nam


#combine all in a table

readRDS('~/height_prediction/gwas/ukb_afr/output/Nr_SNPs_UKB_afr.Rds')[Name %in% paste0('phys_', nam)]-> ukb_afr
readRDS('/project/mathilab/bbita/gwas_admix/new_height/Nr_SNPs_UKB.Rds')[Name %in% paste0('phys_', nam)]-> ukb_eur #need to update this path

setkey(ukb_afr, Name)
setkey(ukb_eur, Name)
dt<-data.table(Name=paste0('phys_',nam), UKB_afr_imp=unlist(partial_R2), UKB_eur_imp=unlist(partial_R2_eur), Nr_imp=unlist(lapply(nam, function(X) nrow(do.call(rbind,readRDS(paste0('~/height_prediction/imputed/output/UKB_vec_all_', X,'.Rds')))))))
setkey(dt, Name)
setkey(dt, Name)[ukb_afr][ukb_eur]-> dt2
dt2[, UKB_afr:=Part_R2]
dt2[, UKB_eur:=i.Part_R2]
dt2[,i.Nr:=NULL][,Part_R2:=NULL][,i.Part_R2:=NULL]
dt2[,eur_diff:=UKB_eur_imp-UKB_eur]
dt2[,afr_diff:=UKB_afr_imp-UKB_afr]
dt2[,Nr_diff:=Nr_imp-Nr]
saveRDS(dt2,'~/height_prediction/imputed/output/comparison_ukb.Rds')

