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
options(scipen=999)

a<-c('_p1.0000e-08', '_p1.0000e-07','_p1.0000e-06', '_p1.0000e-05', '_p3.0000e-05', '_p1.0000e-04','_p3.0000e-04', '_p1.0000e-03','_p3.0000e-03','_p1.0000e-02', '_p3.0000e-02', '_p1.0000e-01', '_p3.0000e-01', '_p1.0000e+00')

a1<-paste0('~/height_prediction/ldpred/output/HRS_EUR.score_P+T_r0.20',a, '.txt')

b<-c('_p1.0000e-03', '_p3.0000e-03','_p1.0000e-02', '_p3.0000e-02', '_p1.0000e-01', '_p3.0000e-01','_p1.0000e+00', '-inf')
b1<-paste0('~/height_prediction/ldpred/output/HRS_EUR.score_LDpred', b, '.txt')

d<-c(unlist(a1), unlist(b1))
lapply(d , function(X) fread(X))-> PGS_HRS_eur

#read in phenotype data
fread('~/height_prediction/input/HRS_eur/HRS_EUR_phenotypes.txt')-> Pheno_HRS_eur
###########
source('~/height_prediction/strat_prs/scripts/Rsq_R2.R')

my_names<-c(paste0('pt', a), paste0('gibbs', b))
names(PGS_HRS_eur)<-my_names

##############

#add PGS to Pheno table in order to be able to make multiple analyses

setkey(Pheno_HRS_eur, ID)

for (I in my_names){
	colnames(PGS_HRS_eur[[I]])[1]<-'ID'
        setkey(PGS_HRS_eur[[I]], ID)
        PGS_HRS_eur[[I]][Pheno_HRS_eur, nomatch=0]-> PGS_HRS_eur[[I]]
	PGS_HRS_eur[[I]][, EUR_ANC:=1]
        PGS_HRS_eur[[I]][,AGE2:=AGE^2]
	PGS_HRS_eur[[I]][,HEIGHT:=HEIGHT*100]
	PGS_HRS_eur[[I]]$SEX<-as.factor(PGS_HRS_eur[[I]]$SEX)
        dt_f<-PGS_HRS_eur[[I]][SEX==2]
        dt_m<-PGS_HRS_eur[[I]][SEX==1]
        sd1_f<-sd(dt_f$HEIGHT)
        m1_f<-mean(dt_f$HEIGHT)
        sd1_m<-sd(dt_m$HEIGHT)
        m1_m<-mean(dt_m$HEIGHT)
        dt_f<-dt_f[HEIGHT>=m1_f-(2*sd1_f)]
        dt_m<-dt_m[HEIGHT>=m1_m-(2*sd1_m)]
	PGS_HRS_eur[[I]]<-rbind(dt_f, dt_m)
}
lapply(PGS_HRS_eur, function(X) lm(HEIGHT~SEX, X))-> lm1_HRS_eur
lapply(PGS_HRS_eur, function(X) lm(HEIGHT~PRS, X))-> lm2_HRS_eur
lapply(PGS_HRS_eur, function(X) lm(HEIGHT~AGE, X))-> lm3_HRS_eur
lapply(PGS_HRS_eur, function(X) lm(HEIGHT~AGE2, X))-> lm4_HRS_eur
lapply(PGS_HRS_eur, function(X) lm(HEIGHT~PRS+AGE, X))-> lm5_HRS_eur
lapply(PGS_HRS_eur, function(X) lm(HEIGHT~PRS+AGE2, X))-> lm6_HRS_eur
lapply(PGS_HRS_eur, function(X) lm(HEIGHT~SEX+AGE+AGE2, X))-> lm7_HRS_eur
lapply(PGS_HRS_eur, function(X) lm(HEIGHT~SEX+AGE+AGE2+PRS, X))-> lm8_HRS_eur

partial_r2_HRS_eur<-lapply(1:length(PGS_HRS_eur), function(X) partial.R2(lm7_HRS_eur[[X]], lm8_HRS_eur[[X]]))
names(partial_r2_HRS_eur)<-names(PGS_HRS_eur)
#

lapply(PGS_HRS_eur, function(X) X[, Quantile:='total'][,Med_Eur_Anc:=1])-> PGS2_HRS_eur

B_HRS_eur<-vector('list', length(PGS2_HRS_eur))
names(B_HRS_eur)<-names(PGS2_HRS_eur)

for (I in names(PGS2_HRS_eur)){
        B_HRS_eur[[I]]<-data.table(Quant="total",R_sq=partial_r2_HRS_eur[[I]],Med_Eur_Anc=1)
        B_HRS_eur[[I]][,N:=nrow(PGS2_HRS_eur[[I]])]
        B_HRS_eur[[I]][, K:=1] #number of predictors. Need to check later if this is correct.
        B_HRS_eur[[I]][, LCL:=CI.Rsq(R_sq, k=K, n=N)[3]]
        B_HRS_eur[[I]][, UCL:=CI.Rsq(R_sq, k=K, n=N)[4]]
}

### add confidence intervals calculated with bootstrap: https://www.statmethods.net/advstats/bootstrapping.html
results.HRS_eur<-vector('list', length(PGS_HRS_eur))
names(results.HRS_eur)<-names(PGS2_HRS_eur)

for (I in names(PGS2_HRS_eur)){
        results.HRS_eur[[I]]<- boot(data=PGS_HRS_eur[[I]], statistic=rsq.R2, R=1000, formula1=HEIGHT~SEX+AGE+AGE2, formula2=HEIGHT~SEX+AGE+AGE2+PRS)
        #names(results.HRS_eur[[I]])<-c(a1[[I]], "total")
        cat(I, ' done\n')
}

#95% confidence intervals.
saveRDS(results.HRS_eur, file='~/height_prediction/ldpred/output/results.HRS_eur.Rds')
saveRDS(PGS2_HRS_eur, file='~/height_prediction/ldpred/output/PGS3_HRS_eur.Rds')

boots.ci.HRS_eur<-lapply(results.HRS_eur, function(X)  boot.ci(X, type = c("norm", 'basic', "perc")))
names(boots.ci.HRS_eur)<-names(results.HRS_eur)

for (I in names(PGS_HRS_eur)){
        B_HRS_eur[[I]][,HVB_L:=1]
        B_HRS_eur[[I]][,HVB_U:=1]
        B_HRS_eur[[I]][, Dataset:='HRS_eur']
        B_HRS_eur[[I]][, boots_norm_L:=boots.ci.HRS_eur[[I]]$normal[2]]
        B_HRS_eur[[I]][, boots_norm_U:=boots.ci.HRS_eur[[I]]$normal[3]]
        B_HRS_eur[[I]][, boots_perc_L:=boots.ci.HRS_eur[[I]]$perc[4]]
        B_HRS_eur[[I]][, boots_perc_U:=boots.ci.HRS_eur[[I]]$perc[5]]
        B_HRS_eur[[I]][, boots_basic_L:=boots.ci.HRS_eur[[I]]$basic[4]]
        B_HRS_eur[[I]][, boots_basic_U:=boots.ci.HRS_eur[[I]]$basic[5]]
}

saveRDS(B_HRS_eur, file="~/height_prediction/ldpred/output/B_HRS_eur.Rds")


