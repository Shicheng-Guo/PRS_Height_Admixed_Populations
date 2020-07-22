#!/usr/bin/env Rscript
##preable###
library("optparse")
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
library(TeachingDemos)
options(scipen=999)

a<-c('_p1.0000e-08','_p1.0000e-07','_p1.0000e-06','_p1.0000e-05','_p3.0000e-05','_p1.0000e-04','_p3.0000e-04','_p1.0000e-03','_p3.0000e-03','_p1.0000e-02', '_p3.0000e-02', '_p1.0000e-01', '_p3.0000e-01', '_p1.0000e+00')

a1<-paste0('~/height_prediction/ldpred/output/HRS_AFR.score_P+T_r0.20',a, '.txt')

b<-c('_p1.0000e-03', '_p3.0000e-03','_p1.0000e-02', '_p3.0000e-02', '_p1.0000e-01', '_p3.0000e-01','_p1.0000e+00', '-inf')

b1<-paste0('~/height_prediction/ldpred/output/HRS_AFR.score_LDpred', b, '.txt')

d<-c(unlist(a1), unlist(b1))
lapply(d , function(X) fread(X))-> PGS_HRS_afr
#read in PGS scores
#read in phenotype data
fread('~/height_prediction/input/HRS_afr/HRS_AFR_phenotypes.txt', fill=T)[,V5:=NULL]-> Pheno_HRS_afr
#a partial R2 function
source('~/height_prediction/strat_prs/scripts/Rsq_R2.R')

my_names<-c(paste0('pt', a), paste0('gibbs', b))
names(PGS_HRS_afr)<-my_names

##############
#add PGS to Pheno table in order to be able to make multiple analyses
#Pheno_HRS_afr[,ID:=paste0(ID, "_", ID)]
as.character(Pheno_HRS_afr$ID)-> Pheno_HRS_afr$ID
Pheno_HRS_afr[, HEIGHT:=HEIGHT*100]
Pheno_HRS_afr[, ":="(SEX=ifelse(SEX==1, "Male", "Female"))]
setkey(Pheno_HRS_afr, ID)
#add ancestry

ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/HRS_afr/rfmix_anc_chr', X, '.txt'))))

anc_HRS_afr<-ancestry %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual
cat('checkpoint\n')
anc_HRS_afr[,ID:=SUBJID][,SUBJID:=NULL]
as.character(anc_HRS_afr$ID)-> anc_HRS_afr$ID
setkey(anc_HRS_afr, ID)

for (I in my_names){
	colnames(PGS_HRS_afr[[I]])[1]<-'ID'
	PGS_HRS_afr[[I]]$ID<-as.character(PGS_HRS_afr[[I]]$ID)
        setkey(PGS_HRS_afr[[I]], ID)
        PGS_HRS_afr[[I]][Pheno_HRS_afr, nomatch=0]-> PGS_HRS_afr[[I]]
        setkey(PGS_HRS_afr[[I]], ID)
        PGS_HRS_afr[[I]][anc_HRS_afr, nomatch=0]-> PGS_HRS_afr[[I]]
        PGS_HRS_afr[[I]][,AGE2:=AGE^2]
	PGS_HRS_afr[[I]][AFR_ANC>=0.05]-> PGS_HRS_afr[[I]]
        PGS_HRS_afr[[I]][which(!is.na(PGS_HRS_afr[[I]][,HEIGHT])),]-> PGS_HRS_afr[[I]]
	PGS_HRS_afr[[I]]$SEX<-as.factor(PGS_HRS_afr[[I]]$SEX)
	dt_f<-PGS_HRS_afr[[I]][SEX=='Female']
	dt_m<-PGS_HRS_afr[[I]][SEX=='Male']
	sd1_f<-sd(dt_f$HEIGHT)
        m1_f<-mean(dt_f$HEIGHT)
	sd1_m<-sd(dt_m$HEIGHT)
        m1_m<-mean(dt_m$HEIGHT)
	dt_f<-dt_f[HEIGHT>=m1_f-(2*sd1_f)]
	dt_m<-dt_m[HEIGHT>=m1_m-(2*sd1_m)]
	PGS_HRS_afr[[I]]<-rbind(dt_f, dt_m)
	cat(I)
	cat('\n')
}
cat('checkpoint\n')

lapply(PGS_HRS_afr, function(X) lm(HEIGHT~SEX, X))-> lm0_HRS_afr
lapply(PGS_HRS_afr, function(X) lm(HEIGHT~PRS, X))-> lm1_HRS_afr
lapply(PGS_HRS_afr, function(X) lm(HEIGHT~AGE, X))-> lm2_HRS_afr
lapply(PGS_HRS_afr, function(X) lm(HEIGHT~AGE2, X))-> lm3_HRS_afr
lapply(PGS_HRS_afr, function(X) lm(HEIGHT~EUR_ANC, X))-> lm4_HRS_afr
lapply(PGS_HRS_afr, function(X) lm(HEIGHT~PRS+AGE, X))-> lm5_HRS_afr
lapply(PGS_HRS_afr, function(X) lm(HEIGHT~PRS+AGE2, X))-> lm6_HRS_afr
lapply(PGS_HRS_afr, function(X) lm(HEIGHT~SEX+AGE+AGE2+EUR_ANC, X))-> lm7_HRS_afr
lapply(PGS_HRS_afr, function(X) lm(HEIGHT~SEX+AGE+AGE2+EUR_ANC+PRS, X))-> lm8_HRS_afr

partial_r2_HRS_afr<-lapply(1:length(PGS_HRS_afr), function(X) partial.R2(lm7_HRS_afr[[X]], lm8_HRS_afr[[X]])) #
names(partial_r2_HRS_afr)<-names(PGS_HRS_afr)

lapply(PGS_HRS_afr, function(X) X[, Quantile:= cut(EUR_ANC,
                                breaks=quantile(EUR_ANC, probs=seq(0,1, by=1/2), na.rm=TRUE),
                                include.lowest=TRUE)])-> PGS2_HRS_afr

names(PGS2_HRS_afr)<-names(PGS_HRS_afr)

lapply(1:length(PGS2_HRS_afr), function(X) PGS2_HRS_afr[[X]][,Med_Eur_Anc:=median(EUR_ANC),by=Quantile])

lapply(1:length(PGS2_HRS_afr), function(X) as.character(unique((PGS2_HRS_afr[[X]]$Quantile))))-> a

lapply(a, function(X) c(X[2],X[1]))-> a1 #check

names(a1)<-names(PGS2_HRS_afr)

r2_HRS_afr<-vector('list', length(PGS2_HRS_afr))	
names(r2_HRS_afr)<-names(PGS2_HRS_afr)

for(I in names(r2_HRS_afr)){
        r2_HRS_afr[[I]]<-vector('list', length(a1[[I]]))
        names(r2_HRS_afr[[I]])<-a1[[I]]
        for(i in a1[[I]]){
        r2_HRS_afr[[I]][[i]]<-partial.R2(lm(HEIGHT~SEX+AGE+AGE2+EUR_ANC, PGS2_HRS_afr[[I]][Quantile==i]),lm(HEIGHT~SEX+AGE+AGE2+EUR_ANC+PRS, PGS2_HRS_afr[[I]][Quantile==i]))
        }
}

B_HRS_afr<-vector('list', length(r2_HRS_afr))
names(B_HRS_afr)<-names(r2_HRS_afr)

for (I in names(r2_HRS_afr)){
        B_HRS_afr[[I]]<-data.table(Quant=c(a1[[I]], "total"),
        R_sq=c(unlist(r2_HRS_afr[[I]]), partial_r2_HRS_afr[[I]]),
        Med_Eur_Anc=c(unique(PGS2_HRS_afr[[I]][Quantile==a1[[I]][1]][,Med_Eur_Anc]), unique(PGS2_HRS_afr[[I]][Quantile==a1[[I]][2]][,Med_Eur_Anc]), median(PGS2_HRS_afr[[I]][, EUR_ANC])))
        B_HRS_afr[[I]][,N:=c(nrow(PGS2_HRS_afr[[I]][Quantile==a1[[I]][1]]), nrow(PGS2_HRS_afr[[I]][Quantile==a1[[I]][2]]),nrow(PGS2_HRS_afr[[I]]))]
        B_HRS_afr[[I]][,K:=1] #number of predictors. Need to check later if this is correct.
        B_HRS_afr[[I]][, LCL:=CI.Rsq(R_sq, k=K, n=N)[3]]
        B_HRS_afr[[I]][, UCL:=CI.Rsq(R_sq, k=K, n=N)[4]]
}

### add confidence intervals calculated with bootstrap: https://www.statmethods.net/advstats/bootstrapping.html
results.HRS_afr<-vector('list', length(PGS2_HRS_afr))
names(results.HRS_afr)<-names(PGS2_HRS_afr)

for (I in names(PGS2_HRS_afr)){
results.HRS_afr[[I]]<-vector('list', length(a1[[I]])+1)
        lapply(a1[[I]], function(i) boot(data=PGS2_HRS_afr[[I]][Quantile==i], statistic=rsq.R2,
                R=1000, formula1=HEIGHT~SEX+AGE+AGE2+EUR_ANC, formula2=HEIGHT~SEX+AGE+AGE2+EUR_ANC+PRS))-> results.HRS_afr[[I]]
        cat(I)
        cat(' done\n')
}

for (I in names(PGS_HRS_afr)){
        tp <- boot(data=PGS_HRS_afr[[I]], statistic=rsq.R2, R=1000, formula1=HEIGHT~SEX+AGE+AGE2+EUR_ANC, formula2=HEIGHT~SEX+AGE+AGE2+EUR_ANC+PRS)
        list.append(results.HRS_afr[[I]], tp)-> results.HRS_afr[[I]]
        names(results.HRS_afr[[I]])<-c(a1[[I]], "total")
	cat(I)
        cat(' done\n')
}
saveRDS(PGS2_HRS_afr, file='~/height_prediction/ldpred/output/PGS3_HRS_afr.Rds')
saveRDS(results.HRS_afr, file='~/height_prediction/ldpred/output/results.HRS_afr.Rds')
#confidence intervals
boots.ci.HRS_afr<-lapply(results.HRS_afr, function(Y) lapply(Y, function(X) boot.ci(X, type = c("norm", 'basic', "perc"))))
names(boots.ci.HRS_afr)<-names(results.HRS_afr)

for (I in names(PGS2_HRS_afr)){
        B_HRS_afr[[I]][1:2,]-> a
        B_HRS_afr[[I]][3,]-> b
        a[,HVB_L:=sapply(a$Quant, function(X) as.numeric(gsub("\\]","",gsub("\\(","",gsub("\\[","",strsplit(X,",")[[1]])))))[1,]]
        a[,HVB_U:=sapply(a$Quant, function(X) as.numeric(gsub("\\]","",gsub("\\(","",gsub("\\[","",strsplit(X,",")[[1]])))))[2,]]
        b[,HVB_L:=1]
        b[,HVB_U:=1]
        rbind(a,b)->B_HRS_afr[[I]]
        B_HRS_afr[[I]][, Dataset:='HRS_AFR']
        B_HRS_afr[[I]][, boots_norm_L:=sapply(1:3, function(X) boots.ci.HRS_afr[[I]][[X]]$normal[2])]
        B_HRS_afr[[I]][, boots_norm_U:=sapply(1:3, function(X) boots.ci.HRS_afr[[I]][[X]]$normal[3])]
        B_HRS_afr[[I]][, boots_perc_L:=sapply(1:3, function(X) boots.ci.HRS_afr[[I]][[X]]$perc[4])]
        B_HRS_afr[[I]][, boots_perc_U:=sapply(1:3, function(X) boots.ci.HRS_afr[[I]][[X]]$perc[5])]
        B_HRS_afr[[I]][, boots_basic_L:=sapply(1:3, function(X) boots.ci.HRS_afr[[I]][[X]]$basic[4])]
        B_HRS_afr[[I]][, boots_basic_U:=sapply(1:3, function(X) boots.ci.HRS_afr[[I]][[X]]$basic[5])]
}

saveRDS(B_HRS_afr, file="~/height_prediction/ldpred/output/B_HRS_afr.Rds")
