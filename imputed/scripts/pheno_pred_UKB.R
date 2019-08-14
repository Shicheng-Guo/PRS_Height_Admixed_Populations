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


lapply(PGS2_UKB_afr, function(X) X[, Quantile:= cut(EUR_ANC,
                                breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.25), na.rm=TRUE),
                                include.lowest=TRUE)])-> PGS3_UKB_afr

names(PGS3_UKB_afr)<-names(PGS2_UKB_afr)

lapply(1:length(PGS3_UKB_afr), function(X) PGS3_UKB_afr[[X]][,Med_Eur_Anc:=median(EUR_ANC),by=Quantile])

lapply(1:length(PGS3_UKB_afr), function(X) as.character(unique((PGS3_UKB_afr[[X]]$Quantile))))-> a

lapply(a, function(X) c(X[2],X[4], X[3], X[1]))-> a1 #check

names(a1)<-names(PGS3_UKB_afr)

r2_UKB_afr<-vector('list', length(PGS3_UKB_afr))
names(r2_UKB_afr)<-names(PGS3_UKB_afr)

for(I in names(r2_UKB_afr)){
        r2_UKB_afr[[I]]<-vector('list', length(a1[[I]]))
        names(r2_UKB_afr[[I]])<-a1[[I]]
        for(i in a1[[I]]){
        r2_UKB_afr[[I]][[i]]<-partial.R2(lm(Height~Sex+Age+AGE2+EUR_ANC, PGS3_UKB_afr[[I]][Quantile==i]),lm(Height~Sex+Age+AGE2+EUR_ANC+PGS, PGS3_UKB_afr[[I]][Quantile==i]))
        }
}



B_UKB_afr<-vector('list', length(r2_UKB_afr))
names(B_UKB_afr)<-names(r2_UKB_afr)

for (I in names(r2_UKB_afr)){
        B_UKB_afr[[I]]<-data.table(Quant=c(a1[[I]], "total"),
        R_sq=c(unlist(r2_UKB_afr[[I]]), partial_R2[[I]]),
        Med_Eur_Anc=c(unique(PGS3_UKB_afr[[I]][Quantile==a1[[I]][1]][,Med_Eur_Anc]), unique(PGS3_UKB_afr[[I]][Quantile==a1[[I]][2]][,Med_Eur_Anc]),unique(PGS3_UKB_afr[[I]][Quantile==a1[[I]][3]][,Med_Eur_Anc]),unique(PGS3_UKB_afr[[I]][Quantile==a1[[I]][4]][,Med_Eur_Anc]), median(PGS3_UKB_afr[[I]][, EUR_ANC])))
        B_UKB_afr[[I]][,N:=c(nrow(PGS3_UKB_afr[[I]][Quantile==a1[[I]][1]]), nrow(PGS3_UKB_afr[[I]][Quantile==a1[[I]][2]]), nrow(PGS3_UKB_afr[[I]][Quantile==a1[[I]][3]]), nrow(PGS3_UKB_afr[[I]][Quantile==a1[[I]][4]]),nrow(PGS3_UKB_afr[[I]]))]
        B_UKB_afr[[I]][,K:=1] #number of predictors. Need to check later if this is correct.
        B_UKB_afr[[I]][, LCL:=CI.Rsq(R_sq, k=K, n=N)[3]]
        B_UKB_afr[[I]][, UCL:=CI.Rsq(R_sq, k=K, n=N)[4]]
}
### add confidence intervals calculated with bootstrap: https://www.statmethods.net/advstats/bootstrapping.html
results.UKB_afr<-vector('list', length(PGS3_UKB_afr))
names(results.UKB_afr)<-names(PGS3_UKB_afr)

for (I in names(PGS3_UKB_afr)){
results.UKB_afr[[I]]<-vector('list', length(a1[[I]])+1)
        lapply(a1[[I]], function(i) boot(data=PGS3_UKB_afr[[I]][Quantile==i], statistic=rsq.R2,
                R=999, formula1=Height~Sex+Age+AGE2+EUR_ANC, formula2=Height~Sex+Age+AGE2+EUR_ANC+PGS))-> results.UKB_afr[[I]]
        cat(I)
        cat(' done\n')
}

for (I in names(PGS3_UKB_afr)){
        tp <- boot(data=PGS2_UKB_afr[[I]], statistic=rsq.R2, R=999, formula1=Height~Sex+Age+AGE2+EUR_ANC, formula2=Height~Sex+Age+AGE2+EUR_ANC+PGS)
        list.append(results.UKB_afr[[I]], tp)-> results.UKB_afr[[I]]
        names(results.UKB_afr[[I]])<-c(a1[[I]], "total")
        cat(' done\n')
}


saveRDS(PGS3_UKB_afr, file='~/height_prediction/imputed/output/PGS3_UKB_afr.Rds')
saveRDS(results.UKB_afr, file='~/height_prediction/imputed/output/results.UKB_afr.Rds')

#confidence intervals

boots.ci.UKB_afr<-lapply(results.UKB_afr, function(Y) lapply(Y, function(X) boot.ci(X, type = c("norm", 'basic', "perc"))))
names(boots.ci.UKB_afr)<-names(results.UKB_afr)

for (I in names(PGS3_UKB_afr)){
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

saveRDS(B_UKB_afr, file="~/height_prediction/imputed/output/B_UKB_afr.Rds")
###################################
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

