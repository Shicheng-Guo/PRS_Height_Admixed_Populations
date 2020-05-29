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

a<-c('_p1.0000e-08', '_p1.0000e-07','_p1.0000e-06', '_p1.0000e-05', '_p3.0000e-05', '_p1.0000e-04','_p3.0000e-04', '_p1.0000e-03','_p3.0000e-03','_p1.0000e-02', '_p3.0000e-02', '_p1.0000e-01', '_p3.0000e-01', '_p1.0000e+00')

a1<-paste0('~/height_prediction/ldpred/output/UKB_EUR.score_P+T_r0.20',a, '.txt')

b<-c('_p1.0000e-03', '_p3.0000e-03','_p1.0000e-02', '_p3.0000e-02', '_p1.0000e-01', '_p3.0000e-01','_p1.0000e+00', '-inf')

b1<-lapply(b, function(X) paste0('~/height_prediction/ldpred/output/UKB_EUR.score_LDpred', X, '.txt'))

d<-c(unlist(a1), unlist(b1))
lapply(a1 , function(X) fread(X))-> PGS_UKB_eur
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

lapply(PGS_UKB_eur, function(X) X[, Quantile:='total'][,Med_Eur_Anc:=1])-> PGS3_UKB_eur

B_UKB_eur<-vector('list', length(PGS3_UKB_eur))
names(B_UKB_eur)<-names(PGS3_UKB_eur)

for (I in names(PGS3_UKB_eur)){
        B_UKB_eur[[I]]<-data.table(Quant="total",R_sq=partial_r2_UKB_eur[[I]],Med_Eur_Anc=1)
        B_UKB_eur[[I]][,N:=nrow(PGS3_UKB_eur[[I]])]
        B_UKB_eur[[I]][,K:=1] #number of predictors. Need to check later if this is correct.
        B_UKB_eur[[I]][, LCL:=CI.Rsq(R_sq, k=K, n=N)[3]]
        B_UKB_eur[[I]][, UCL:=CI.Rsq(R_sq, k=K, n=N)[4]]
}


### add confidence intervals calculated with bootstrap: https://www.statmethods.net/advstats/bootstrapping.html
results.UKB_eur<-vector('list', length(PGS3_UKB_eur))
names(results.UKB_eur)<-names(PGS3_UKB_eur)


for (I in names(PGS3_UKB_eur)){
        results.UKB_eur[[I]] <- boot(data=PGS_UKB_eur[[I]], statistic=rsq.R2, R=1000, formula1=Height~Sex+Age+age2, formula2=Height~Sex+Age+age2+PRS)
        cat(' done\n')
}
saveRDS(PGS3_UKB_eur, file='~/height_prediction/ldpred/output/PGS3_UKB_eur.Rds')
saveRDS(results.UKB_eur, file='~/height_prediction/ldpred/output/results.UKB_eur.Rds')

#confidence intervals

boots.ci.UKB_eur<-lapply(results.UKB_eur, function(X)  boot.ci(X, type = c("norm", 'basic', "perc")))
names(boots.ci.UKB_eur)<-names(results.UKB_eur)

for (I in names(PGS3_UKB_eur)){
        B_UKB_eur[[I]][,HVB_L:=1]
        B_UKB_eur[[I]][,HVB_U:=1]
        B_UKB_eur[[I]][, Dataset:='UKB_EUR']
        B_UKB_eur[[I]][, boots_norm_L:=boots.ci.UKB_eur[[I]]$normal[2]]
        B_UKB_eur[[I]][, boots_norm_U:=boots.ci.UKB_eur[[I]]$normal[3]]
        B_UKB_eur[[I]][, boots_perc_L:=boots.ci.UKB_eur[[I]]$perc[4]]
        B_UKB_eur[[I]][, boots_perc_U:=boots.ci.UKB_eur[[I]]$perc[5]]
        B_UKB_eur[[I]][, boots_basic_L:=boots.ci.UKB_eur[[I]]$basic[4]]
        B_UKB_eur[[I]][, boots_basic_U:=boots.ci.UKB_eur[[I]]$basic[5]]
}

saveRDS(B_UKB_eur, file="~/height_prediction/ldpred/output/B_UKB_eur.Rds")
