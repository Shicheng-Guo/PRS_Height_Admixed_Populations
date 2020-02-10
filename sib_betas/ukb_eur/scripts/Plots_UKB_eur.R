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

#read in PGS scores
readRDS('~/height_prediction/sib_betas/ukb_eur/output/PGS_ukb_eur.Rds')-> PGS_UKB_eur
#read in phenotype data
fread('~/height_prediction/input/ukb_eur/UKB_EUR_pheno.txt')-> Pheno_UKB_eur
#a partial R2 function
source('~/height_prediction/strat_prs/scripts/Rsq_R2.R')

##############

#add PGS to Pheno table in order to be able to make multiple analyses
#Pheno_UKB_eur[,ID:=paste0(ID, "_", ID)]
Pheno_UKB_eur[, ID:=paste0(ID, "_", ID)]
setkey(Pheno_UKB_eur, ID)


PGS2_UKB_eur<-vector('list', length(PGS_UKB_eur))
names(PGS2_UKB_eur)<-names(PGS_UKB_eur)


for (I in names(PGS_UKB_eur)){
        data.table(ID=names(PGS_UKB_eur[[I]]), PGS=unlist(PGS_UKB_eur[[I]]))-> PGS2_UKB_eur[[I]]
        setkey(PGS2_UKB_eur[[I]], ID)
        PGS2_UKB_eur[[I]][Pheno_UKB_eur, nomatch=0]-> PGS2_UKB_eur[[I]]
	PGS2_UKB_eur[[I]][, EUR_ANC:=1]
        PGS2_UKB_eur[[I]][,age2:=Age^2]
        PGS2_UKB_eur[[I]][-which(is.na(PGS2_UKB_eur[[I]][,Height])),]-> PGS2_UKB_eur[[I]]
}

lapply(PGS2_UKB_eur, function(X) lm(Height~Sex, X))-> lm0_UKB_eur
lapply(PGS2_UKB_eur, function(X) lm(Height~PGS, X))-> lm1_UKB_eur
lapply(PGS2_UKB_eur, function(X) lm(Height~Age, X))-> lm2_UKB_eur
lapply(PGS2_UKB_eur, function(X) lm(Height~age2, X))-> lm3_UKB_eur
lapply(PGS2_UKB_eur, function(X) lm(Height~EUR_ANC, X))-> lm4_UKB_eur
lapply(PGS2_UKB_eur, function(X) lm(Height~PGS+Age, X))-> lm5_UKB_eur
lapply(PGS2_UKB_eur, function(X) lm(Height~PGS+age2, X))-> lm6_UKB_eur
lapply(PGS2_UKB_eur, function(X) lm(Height~Sex+Age+age2+EUR_ANC, X))-> lm7_UKB_eur
lapply(PGS2_UKB_eur, function(X) lm(Height~Sex+Age+age2+EUR_ANC+PGS, X))-> lm8_UKB_eur

partial.R2(lm7_UKB_eur[[63]],lm8_UKB_eur[[63]]) # 20.43% 2nd highest after LD_block_AFR

partial_r2_UKB_eur<-lapply(1:length(PGS2_UKB_eur), function(X) partial.R2(lm7_UKB_eur[[X]], lm8_UKB_eur[[X]])) #range 8.5 to 21.5%
names(partial_r2_UKB_eur)<-names(PGS2_UKB_eur)


Nr_SNPs<-rep(NA, length(PGS2_UKB_eur))
names(Nr_SNPs)<- names(PGS2_UKB_eur)

for(I in names(Nr_SNPs)){
        readRDS(paste0('~/height_prediction/sib_betas/ukb_eur/output/hei_', I, '.Rds'))-> smtg
        sum(sapply(smtg, function(X) nrow(X)))-> Nr_SNPs[I]
	remove(smtg)
	gc()
        cat(I, ' \n')
}

data.table(Nr=unlist(Nr_SNPs), Name=names(Nr_SNPs), Part_R2=unlist(partial_r2_UKB_eur))-> A_table
saveRDS(A_table, file='~/height_prediction/sib_betas/ukb_eur/output/Nr_SNPs_UKB_eur.Rds')

#
cor.test(unlist(Nr_SNPs), unlist(partial_r2_UKB_eur))# 0.245
summary(lm(Part_R2~Nr,data=A_table))$r.squared #0.06

for(I in 1:length(PGS2_UKB_eur)){
        A<-ggpairs(PGS2_UKB_eur[[I]][,.(Height, Sex, PGS, Age, age2)])
        png(filename=paste0('~/height_prediction/sib_betas/ukb_eur/figs/UKB_eur_ggpairs_', names(PGS2_UKB_eur)[I], '.png'))
        print(A)
        cat(names(PGS2_UKB_eur)[I])
        cat(' done\n')
        dev.off()
}


lapply(PGS2_UKB_eur, function(X) X[, Quantile:='total'][,Med_Eur_Anc:=1])-> PGS3_UKB_eur

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
        results.UKB_eur[[I]] <- boot(data=PGS2_UKB_eur[[I]], statistic=rsq.R2, R=999, formula1=Height~Sex+Age+age2, formula2=Height~Sex+Age+age2+PGS)
#        names(results.UKB_eur[[I]])<-c(a1[[I]], "total")
        cat(' done\n')
}
saveRDS(PGS3_UKB_eur, file='~/height_prediction/sib_betas/ukb_eur/output/PGS3_UKB_eur.Rds')
saveRDS(results.UKB_eur, file='~/height_prediction/sib_betas/ukb_eur/output/results.UKB_eur.Rds')

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

saveRDS(B_UKB_eur, file="~/height_prediction/sib_betas/ukb_eur/output/B_UKB_eur.Rds")
