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

txtStart(paste0("~/height_prediction/gwas/HRS_afr/plots_out.txt"))

#read in PGS scores
readRDS('~/height_prediction/gwas/HRS_afr/output/PGS_HRS_afr.Rds')-> PGS_HRS_afr
#read in phenotype data
fread('~/height_prediction/input/HRS_afr/HRS_AFR_phenotypes.txt', fill=T)[,V5:=NULL]-> Pheno_HRS_afr
#a partial R2 function
source('~/height_prediction/strat_prs/scripts/Rsq_R2.R')

##############

#add PGS to Pheno table in order to be able to make multiple analyses
#Pheno_HRS_afr[,ID:=paste0(ID, "_", ID)]
for(j in 1:length(PGS_HRS_afr)){
	gsub("[0-9]+_","",names(PGS_HRS_afr[[j]]))-> names(PGS_HRS_afr[[j]])
}
as.character(Pheno_HRS_afr$ID)-> Pheno_HRS_afr$ID
Pheno_HRS_afr[, HEIGHT:=HEIGHT*100]
Pheno_HRS_afr[, ":="(SEX=ifelse(SEX==1, "Male", "Female"))]
setkey(Pheno_HRS_afr, ID)
#add ancestry
ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/HRS_afr/rfmix_anc_chr', X, '.txt'))))

anc_HRS_afr<-ancestry %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual

anc_HRS_afr[,ID:=SUBJID][,SUBJID:=NULL]
as.character(anc_HRS_afr$ID)-> anc_HRS_afr$ID
setkey(anc_HRS_afr, ID)

PGS2_HRS_afr<-vector('list', length(PGS_HRS_afr))
names(PGS2_HRS_afr)<-names(PGS_HRS_afr)


for (I in names(PGS_HRS_afr)){
        data.table(ID=names(PGS_HRS_afr[[I]]), PGS=unlist(PGS_HRS_afr[[I]]))-> PGS2_HRS_afr[[I]]
        setkey(PGS2_HRS_afr[[I]], ID)
        PGS2_HRS_afr[[I]][Pheno_HRS_afr, nomatch=0]-> PGS2_HRS_afr[[I]]
        PGS2_HRS_afr[[I]][anc_HRS_afr, nomatch=0]-> PGS2_HRS_afr[[I]]
        PGS2_HRS_afr[[I]][,AGE2:=AGE^2]
	PGS2_HRS_afr[[I]][AFR_ANC>=0.05]-> PGS2_HRS_afr[[I]]
        PGS2_HRS_afr[[I]][which(!is.na(PGS2_HRS_afr[[I]][,HEIGHT])),]-> PGS2_HRS_afr[[I]]
}

lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~SEX, X))-> lm0_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~PGS, X))-> lm1_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~AGE, X))-> lm2_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~AGE2, X))-> lm3_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~EUR_ANC, X))-> lm4_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~PGS+AGE, X))-> lm5_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~PGS+AGE2, X))-> lm6_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~SEX+AGE+AGE2+EUR_ANC, X))-> lm7_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~SEX+AGE+AGE2+EUR_ANC+PGS, X))-> lm8_HRS_afr


partial_r2_HRS_afr<-lapply(1:length(PGS2_HRS_afr), function(X) partial.R2(lm7_HRS_afr[[X]], lm8_HRS_afr[[X]])) #
names(partial_r2_HRS_afr)<-names(PGS2_HRS_afr)

partial.R2(lm7_HRS_afr[[63]],lm8_HRS_afr[[63]]) #0.02376103
partial.R2(lm7_HRS_afr[[67]],lm8_HRS_afr[[67]]) #0.03114019
partial.R2(lm7_HRS_afr[[64]],lm8_HRS_afr[[64]]) #0.02612528
partial.R2(lm7_HRS_afr[[35]],lm8_HRS_afr[[35]]) #0.02863406

Nr_SNPs<-rep(NA, length(PGS2_HRS_afr))
names(Nr_SNPs)<- names(PGS2_HRS_afr)

for(I in names(Nr_SNPs)){
        readRDS(paste0('~/height_prediction/gwas/HRS_afr/output/hei_', I, '.Rds'))-> smtg
        sum(sapply(smtg, function(X) nrow(X)))-> Nr_SNPs[I]
	remove(smtg)
	gc()
        cat(I, ' \n')
}

data.table(Nr=unlist(Nr_SNPs), Name=names(Nr_SNPs), Part_R2=unlist(partial_r2_HRS_afr))-> A_table
saveRDS(A_table, file='~/height_prediction/gwas/HRS_afr/output/Nr_SNPs_HRS_afr.Rds')

#
cor.test(unlist(Nr_SNPs), unlist(partial_r2_HRS_afr))#  0.84
summary(lm(Part_R2~Nr,data=A_table))$r.squared #0.72

for(I in 1:length(PGS2_HRS_afr)){
        A<-ggpairs(PGS2_HRS_afr[[I]][,.(HEIGHT, SEX, PGS, AGE, AGE2,EUR_ANC)])
        png(filename=paste0('~/height_prediction/gwas/HRS_afr/figs/HRS_afr_ggpairs_', names(PGS2_HRS_afr)[I], '.png'))
        print(A)
        cat(names(PGS2_HRS_afr)[I])
        cat(' done\n')
        dev.off()
}


lapply(PGS2_HRS_afr, function(X) X[, Quantile:= cut(EUR_ANC,
                                breaks=quantile(EUR_ANC, probs=seq(0,1, by=1/2), na.rm=TRUE),
                                include.lowest=TRUE)])-> PGS3_HRS_afr

names(PGS3_HRS_afr)<-names(PGS2_HRS_afr)

lapply(1:length(PGS3_HRS_afr), function(X) PGS3_HRS_afr[[X]][,Med_Eur_Anc:=median(EUR_ANC),by=Quantile])

lapply(1:length(PGS3_HRS_afr), function(X) as.character(unique((PGS3_HRS_afr[[X]]$Quantile))))-> a

lapply(a, function(X) c(X[2],X[1]))-> a1 #check

names(a1)<-names(PGS3_HRS_afr)

r2_HRS_afr<-vector('list', length(PGS3_HRS_afr))	
names(r2_HRS_afr)<-names(PGS3_HRS_afr)

for(I in names(r2_HRS_afr)){
        r2_HRS_afr[[I]]<-vector('list', length(a1[[I]]))
        names(r2_HRS_afr[[I]])<-a1[[I]]
        for(i in a1[[I]]){
        r2_HRS_afr[[I]][[i]]<-partial.R2(lm(HEIGHT~SEX+AGE+AGE2+EUR_ANC, PGS3_HRS_afr[[I]][Quantile==i]),lm(HEIGHT~SEX+AGE+AGE2+EUR_ANC+PGS, PGS3_HRS_afr[[I]][Quantile==i]))
        }
}

B_HRS_afr<-vector('list', length(r2_HRS_afr))
names(B_HRS_afr)<-names(r2_HRS_afr)

for (I in names(r2_HRS_afr)){
        B_HRS_afr[[I]]<-data.table(Quant=c(a1[[I]], "total"),
        R_sq=c(unlist(r2_HRS_afr[[I]]), partial_r2_HRS_afr[[I]]),
        Med_Eur_Anc=c(unique(PGS3_HRS_afr[[I]][Quantile==a1[[I]][1]][,Med_Eur_Anc]), unique(PGS3_HRS_afr[[I]][Quantile==a1[[I]][2]][,Med_Eur_Anc]), median(PGS3_HRS_afr[[I]][, EUR_ANC])))
        B_HRS_afr[[I]][,N:=c(nrow(PGS3_HRS_afr[[I]][Quantile==a1[[I]][1]]), nrow(PGS3_HRS_afr[[I]][Quantile==a1[[I]][2]]),nrow(PGS3_HRS_afr[[I]]))]
        B_HRS_afr[[I]][,K:=1] #number of predictors. Need to check later if this is correct.
        B_HRS_afr[[I]][, LCL:=CI.Rsq(R_sq, k=K, n=N)[3]]
        B_HRS_afr[[I]][, UCL:=CI.Rsq(R_sq, k=K, n=N)[4]]
}

### add confidence intervals calculated with bootstrap: https://www.statmethods.net/advstats/bootstrapping.html
results.HRS_afr<-vector('list', length(PGS3_HRS_afr))
names(results.HRS_afr)<-names(PGS3_HRS_afr)

for (I in names(PGS3_HRS_afr)){
results.HRS_afr[[I]]<-vector('list', length(a1[[I]])+1)
        lapply(a1[[I]], function(i) boot(data=PGS3_HRS_afr[[I]][Quantile==i], statistic=rsq.R2,
                R=999, formula1=HEIGHT~SEX+AGE+AGE2+EUR_ANC, formula2=HEIGHT~SEX+AGE+AGE2+EUR_ANC+PGS))-> results.HRS_afr[[I]]
        cat(I)
        cat(' done\n')
}

for (I in names(PGS3_HRS_afr)){
        tp <- boot(data=PGS2_HRS_afr[[I]], statistic=rsq.R2, R=999, formula1=HEIGHT~SEX+AGE+AGE2+EUR_ANC, formula2=HEIGHT~SEX+AGE+AGE2+EUR_ANC+PGS)
        list.append(results.HRS_afr[[I]], tp)-> results.HRS_afr[[I]]
        names(results.HRS_afr[[I]])<-c(a1[[I]], "total")
        cat(' done\n')
}
saveRDS(PGS3_HRS_afr, file='~/height_prediction/gwas/HRS_afr/output/PGS3_HRS_afr.Rds')
saveRDS(results.HRS_afr, file='~/height_prediction/gwas/HRS_afr/output/results.HRS_afr.Rds')

#confidence intervals

boots.ci.HRS_afr<-lapply(results.HRS_afr, function(Y) lapply(Y, function(X) boot.ci(X, type = c("norm", 'basic', "perc"))))
names(boots.ci.HRS_afr)<-names(results.HRS_afr)

for (I in names(PGS3_HRS_afr)){
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


for (I in names(PGS3_HRS_afr)){
       myp<-ggplot(B_HRS_afr[[I]][1:2,], aes(x=Med_Eur_Anc, y=R_sq)) +
        geom_point(size=1.5, shape=21, fill="white") +
        geom_errorbar(aes(x=Med_Eur_Anc, ymin=boots_perc_L, ymax=boots_perc_U), width=0.05, size=0.8, color="cornflowerblue") +
        geom_line(color='lightgray')+
        geom_errorbarh(aes(x=Med_Eur_Anc, xmin=HVB_L, xmax=HVB_U), width=0.05, size=0.5, color="cornflowerblue") +
        labs(title = "HRS_AFR") + ylab("R-squared") + xlab("European Ancestry Proportion")
        print(myp)
        ggsave(paste0('~/height_prediction/gwas/HRS_afr/figs/HRS_afr_error_bars_', I, '.png'))
}

saveRDS(B_HRS_afr, file="~/height_prediction/gwas/HRS_afr/output/B_HRS_afr.Rds")
txtStop()
