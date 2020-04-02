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

lapply(a , function(X) fread(paste0('~/height_prediction/ldpred/output/UKB_AFR_pt.score_P+T_r0.20', X,'.txt')))-> PGS_UKB_afr
names(PGS_UKB_afr)<-a
fread('~/height_prediction/input/ukb_afr/UKB_AFR_pheno.txt')-> Pheno_UKB_afr
#a partial R2 function
source('~/height_prediction/strat_prs/scripts/Rsq_R2.R')


setkey(Pheno_UKB_afr, ID)

#add ancestry
ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/ukb_afr/rfmix_anc_chr', X, '.txt'))))

anc_ukb_afr<-ancestry %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual

anc_ukb_afr[,ID:=SUBJID][,SUBJID:=NULL]
setkey(anc_ukb_afr, ID)


for (I in a){
	colnames(PGS_UKB_afr[[I]])[1]<-'ID'
        setkey(PGS_UKB_afr[[I]], ID)
        PGS_UKB_afr[[I]][Pheno_UKB_afr, nomatch=0]-> PGS_UKB_afr[[I]]
	setkey(PGS_UKB_afr[[I]], ID)
        PGS_UKB_afr[[I]][anc_ukb_afr, nomatch=0]-> PGS_UKB_afr[[I]]
        PGS_UKB_afr[[I]][,age2:=Age^2]
        PGS_UKB_afr[[I]][AFR_ANC>=0.05]-> PGS_UKB_afr[[I]]
        PGS_UKB_afr[[I]][-which(is.na(PGS_UKB_afr[[I]][,Height])),]-> PGS_UKB_afr[[I]]
	PGS_UKB_afr[[I]][ID!="6007195"]-> PGS_UKB_afr[[I]]
}

lapply(PGS_UKB_afr, function(X) lm(Height~Sex, X))-> lm0_UKB_afr
lapply(PGS_UKB_afr, function(X) lm(Height~PRS, X))-> lm1_UKB_afr
lapply(PGS_UKB_afr, function(X) lm(Height~Age, X))-> lm2_UKB_afr
lapply(PGS_UKB_afr, function(X) lm(Height~age2, X))-> lm3_UKB_afr
lapply(PGS_UKB_afr, function(X) lm(Height~EUR_ANC, X))-> lm4_UKB_afr
lapply(PGS_UKB_afr, function(X) lm(Height~PRS+Age, X))-> lm5_UKB_afr
lapply(PGS_UKB_afr, function(X) lm(Height~PRS+age2, X))-> lm6_UKB_afr
lapply(PGS_UKB_afr, function(X) lm(Height~Sex+Age+age2+EUR_ANC, X))-> lm7_UKB_afr
lapply(PGS_UKB_afr, function(X) lm(Height~Sex+Age+age2+EUR_ANC+PRS, X))-> lm8_UKB_afr

partial_r2_UKB_afr<-lapply(1:length(PGS_UKB_afr), function(X) partial.R2(lm7_UKB_afr[[X]], lm8_UKB_afr[[X]])) 
names(partial_r2_UKB_afr)<-names(PGS_UKB_afr)

lapply(PGS_UKB_afr, function(X) X[, Quantile:= cut(EUR_ANC,
                                breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.25), na.rm=TRUE),
                                include.lowest=TRUE)])-> PGS2_UKB_afr

names(PGS2_UKB_afr)<-names(PGS_UKB_afr)

lapply(1:length(PGS2_UKB_afr), function(X) PGS2_UKB_afr[[X]][,Med_Eur_Anc:=median(EUR_ANC),by=Quantile])

lapply(1:length(PGS2_UKB_afr), function(X) as.character(unique((PGS2_UKB_afr[[X]]$Quantile))))-> a

lapply(a, function(X) c(X[2],X[4], X[3], X[1]))-> a1 #check

names(a1)<-names(PGS2_UKB_afr)

r2_UKB_afr<-vector('list', length(PGS2_UKB_afr))
names(r2_UKB_afr)<-names(PGS2_UKB_afr)

for(I in names(r2_UKB_afr)){
        r2_UKB_afr[[I]]<-vector('list', length(a1[[I]]))
        names(r2_UKB_afr[[I]])<-a1[[I]]
        for(i in a1[[I]]){
        r2_UKB_afr[[I]][[i]]<-partial.R2(lm(Height~Sex+Age+age2+EUR_ANC, PGS2_UKB_afr[[I]][Quantile==i]),lm(Height~Sex+Age+age2+EUR_ANC+PRS, PGS2_UKB_afr[[I]][Quantile==i]))
        }
}

B_UKB_afr<-vector('list', length(r2_UKB_afr))
names(B_UKB_afr)<-names(r2_UKB_afr)

partial_r2_UKB_afr<-lapply(1:length(PGS_UKB_afr), function(X) partial.R2(lm7_UKB_afr[[X]], lm8_UKB_afr[[X]]))
names(partial_r2_UKB_afr)<- names(PGS_UKB_afr)


for (I in names(r2_UKB_afr)){
        B_UKB_afr[[I]]<-data.table(Quant=c(a1[[I]], "total"),
        R_sq=c(unlist(r2_UKB_afr[[I]]), partial_r2_UKB_afr[[I]]),
        Med_Eur_Anc=c(unique(PGS2_UKB_afr[[I]][Quantile==a1[[I]][1]][,Med_Eur_Anc]), unique(PGS2_UKB_afr[[I]][Quantile==a1[[I]][2]][,Med_Eur_Anc]),unique(PGS2_UKB_afr[[I]][Quantile==a1[[I]][3]][,Med_Eur_Anc]),unique(PGS2_UKB_afr[[I]][Quantile==a1[[I]][4]][,Med_Eur_Anc]), median(PGS2_UKB_afr[[I]][, EUR_ANC])))
        B_UKB_afr[[I]][,N:=c(nrow(PGS2_UKB_afr[[I]][Quantile==a1[[I]][1]]), nrow(PGS2_UKB_afr[[I]][Quantile==a1[[I]][2]]), nrow(PGS2_UKB_afr[[I]][Quantile==a1[[I]][3]]), nrow(PGS2_UKB_afr[[I]][Quantile==a1[[I]][4]]),nrow(PGS2_UKB_afr[[I]]))]
        B_UKB_afr[[I]][,K:=1] #number of predictors. Need to check later if this is correct.
        B_UKB_afr[[I]][, LCL:=CI.Rsq(R_sq, k=K, n=N)[3]]
        B_UKB_afr[[I]][, UCL:=CI.Rsq(R_sq, k=K, n=N)[4]]
}

### add confidence intervals calculated with bootstrap: https://www.statmethods.net/advstats/bootstrapping.html
results.UKB_afr<-vector('list', length(PGS2_UKB_afr))
names(results.UKB_afr)<-names(PGS2_UKB_afr)

for (I in names(PGS2_UKB_afr)){
	results.UKB_afr[[I]]<-vector('list', length(a1[[I]])+1)
        lapply(a1[[I]], function(i) boot(data=PGS2_UKB_afr[[I]][Quantile==i], statistic=rsq.R2,
                R=999, formula1=Height~Sex+Age+age2+EUR_ANC, formula2=Height~Sex+Age+age2+EUR_ANC+PRS))-> results.UKB_afr[[I]]
        cat(I)
        cat(' done\n')
}

for (I in names(PGS2_UKB_afr)){
        tp <- boot(data=PGS2_UKB_afr[[I]], statistic=rsq.R2, R=999, formula1=Height~Sex+Age+age2+EUR_ANC, formula2=Height~Sex+Age+age2+EUR_ANC+PRS)
        list.append(results.UKB_afr[[I]], tp)-> results.UKB_afr[[I]]
        names(results.UKB_afr[[I]])<-c(a1[[I]], "total")
        cat(' done\n')
}
saveRDS(PGS2_UKB_afr, file='~/height_prediction/ldpred/output/PGS2_UKB_afr_pt.Rds')
saveRDS(results.UKB_afr, file='~/height_prediction/ldpred/output/results.UKB_afr_pt.Rds')

#confidence intervals

boots.ci.UKB_afr<-lapply(results.UKB_afr, function(Y) lapply(Y, function(X) boot.ci(X, type = c("norm", 'basic', "perc"))))
names(boots.ci.UKB_afr)<-names(results.UKB_afr)

for (I in names(PGS2_UKB_afr)){
        B_UKB_afr[[I]][1:4,]-> a
        B_UKB_afr[[I]][5,]-> b
        a[,HVB_L:=sapply(a$Quant, function(X) as.numeric(gsub("\\]","",gsub("\\(","",gsub("\\[","",strsplit(X,",")[[1]])))))[1,]]
        a[,HVB_U:=sapply(a$Quant, function(X) as.numeric(gsub("\\]","",gsub("\\(","",gsub("\\[","",strsplit(X,",")[[1]])))))[2,]]
        b[,HVB_L:=1]
        b[,HVB_U:=1]
        rbind(a,b)->B_UKB_afr[[I]]
        B_UKB_afr[[I]][, Dataset:='UKB_AFR']
        B_UKB_afr[[I]][, boots_norm_L:=sapply(1:5, function(X) boots.ci.UKB_afr[[I]][[X]]$normal[2])]
        B_UKB_afr[[I]][, boots_norm_U:=sapply(1:5, function(X) boots.ci.UKB_afr[[I]][[X]]$normal[3])]
        B_UKB_afr[[I]][, boots_perc_L:=sapply(1:5, function(X) boots.ci.UKB_afr[[I]][[X]]$perc[4])]
        B_UKB_afr[[I]][, boots_perc_U:=sapply(1:5, function(X) boots.ci.UKB_afr[[I]][[X]]$perc[5])]
        B_UKB_afr[[I]][, boots_basic_L:=sapply(1:5, function(X) boots.ci.UKB_afr[[I]][[X]]$basic[4])]
        B_UKB_afr[[I]][, boots_basic_U:=sapply(1:5, function(X) boots.ci.UKB_afr[[I]][[X]]$basic[5])]
}


saveRDS(B_UKB_afr, file="~/height_prediction/ldpred/output/B_UKB_afr_pt.Rds")

