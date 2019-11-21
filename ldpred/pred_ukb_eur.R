#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(ggplot2);library(reshape2); library(wesanderson)
library(rlist)
library(asbio)
library(GGally)
library(tidyr)
library(hexbin)
library(psychometric)
library(boot)

a<-c('_p1.0000e-03', '_p3.0000e-03','_p1.0000e-02', '_p3.0000e-02', '_p1.0000e-01', '_p3.0000e-01','_p1.0000e+00', '-inf')
lapply(a , function(X) fread(paste0('~/height_prediction/ldpred/output/EUR.score_LDpred', X,'.txt')))-> PGS_UKB_eur
names(PGS_UKB_eur)<-a
fread('~/height_prediction/input/ukb_eur/UKB_EUR_pheno.txt')-> Pheno_UKB_eur
#a partial R2 function
source('~/height_prediction/strat_prs/scripts/Rsq_R2.R')


setkey(Pheno_UKB_eur, ID)
for(I in a){
	colnames(PGS_UKB_eur[[I]])[1]<-'ID'
	setkey(PGS_UKB_eur[[I]], ID)
	PGS_UKB_eur[[I]][Pheno_UKB_eur, nomatch=0]-> PGS_UKB_eur[[I]]
	PGS_UKB_eur[[I]][, EUR_ANC:=1]
	PGS_UKB_eur[[I]][,age2:=Age^2]
	PGS_UKB_eur[[I]][-which(is.na(PGS_UKB_eur[[I]][,Height])),]-> PGS_UKB_eur[[I]]
}

lapply(PGS_UKB_eur, function(X) lm(Height~Sex, data=X))-> lm0_UKB_eur
lapply(PGS_UKB_eur, function(X) lm(Height~PRS, data=X))-> lm1_UKB_eur
lapply(PGS_UKB_eur, function(X) lm(Height~Age, data=X))-> lm2_UKB_eur
lapply(PGS_UKB_eur, function(X) lm(Height~age2, data=X))-> lm3_UKB_eur
lapply(PGS_UKB_eur, function(X) lm(Height~EUR_ANC, data=X))-> lm4_UKB_eur
lapply(PGS_UKB_eur, function(X) lm(Height~PRS+Age, data=X))-> lm5_UKB_eur
lapply(PGS_UKB_eur, function(X) lm(Height~PRS+age2, data=X))-> lm6_UKB_eur
lapply(PGS_UKB_eur, function(X) lm(Height~Sex+Age+age2+EUR_ANC, data=X))-> lm7_UKB_eur
lapply(PGS_UKB_eur, function(X) lm(Height~Sex+Age+age2+EUR_ANC+PRS, data=X))-> lm8_UKB_eur


partial_r2_UKB_eur<-lapply(1:length(PGS_UKB_eur), function(X) partial.R2(lm7_UKB_eur[[X]], lm8_UKB_eur[[X]]))
names(partial_r2_UKB_eur)<- names(PGS_UKB_eur)

