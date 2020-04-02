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
nam<-paste(rep(c(5000,10000,25000,50000,75000,100000,500000),5), c(0.0005,0.00005, 0.000005,0.0000005,0.00000005), sep="_")
names(PGS_HRS_afr)<-nam
for(N in nam){
        readRDS(paste0('~/height_prediction/imputed/output/PGS2', N, '.Rds'))[[1]]-> PGS_HRS_afr[[N]]
}

#read in phenotype data
fread('~/height_prediction/input/HRS_afr/HRS_AFR_phenotypes.txt', fill=T)[,V5:=NULL]-> Pheno_HRS_afr
#a partial R2 function
source('~/height_prediction/strat_prs/scripts/Rsq_R2.R')

#add PGS to Pheno table in order to be able to make multiple analyses
#Pheno_HRS_afr[,ID:=paste0(ID, "_", ID)]
for (N in nam){
	for(j in 1:length(PGS_HRS_afr[[N]])){
        	gsub("[0-9]+_","",names(PGS_HRS_afr[[N]])[[j]])-> names(PGS_HRS_afr[[N]])[[j]]
	}
}

as.character(Pheno_HRS_afr$ID)-> Pheno_HRS_afr$ID
setkey(Pheno_HRS_afr, ID)
Pheno_HRS_afr[, HEIGHT:=HEIGHT*100]
Pheno_HRS_afr[, ":="(SEX=ifelse(SEX==1, "Male", "Female"))]
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
	#PGS2_HRS_afr[[N]][which(!is.na(PGS2_HRS_afr[[N]][,HEIGHT])),]-> PGS2_HRS_afr[[N]]
        dt_f<-PGS2_HRS_afr[[N]][SEX=='Female']
        dt_m<-PGS2_HRS_afr[[N]][SEX=='Male']
        sd1_f<-sd(dt_f$HEIGHT)
        m1_f<-mean(dt_f$HEIGHT)
        sd1_m<-sd(dt_m$HEIGHT)
        m1_m<-mean(dt_m$HEIGHT)
        dt_f<-dt_f[HEIGHT>=m1_f-(2*sd1_f)]
        dt_m<-dt_m[HEIGHT>=m1_m-(2*sd1_m)]
        PGS2_HRS_afr[[N]]<-rbind(dt_f, dt_m)
}
nrow(PGS2_HRS_afr[[1]])

lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~SEX, X))-> lm0_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~PGS, X))-> lm1_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~AGE, X))-> lm2_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~AGE2, X))-> lm3_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~EUR_ANC, X))-> lm4_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~PGS+AGE, X))-> lm5_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~PGS+AGE2,X))-> lm6_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~SEX+AGE+AGE2+EUR_ANC, X))-> lm7_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~SEX+AGE+AGE2+EUR_ANC+PGS, X))-> lm8_HRS_afr

partial_R2<-lapply(nam, function(X) partial.R2(lm7_HRS_afr[[X]], lm8_HRS_afr[[X]]))
names(partial_R2)<-nam
#######################

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
        R_sq=c(unlist(r2_HRS_afr[[I]]), partial_R2[[I]]),
        Med_Eur_Anc=c(unique(PGS3_HRS_afr[[I]][Quantile==a1[[I]][1]][,Med_Eur_Anc]), unique(PGS3_HRS_afr[[I]][Quantile==a1[[I]][2]][,Med_Eur_Anc]), median(PGS3_HRS_afr[[I]][, EUR_ANC])))
        B_HRS_afr[[I]][,N:=c(nrow(PGS3_HRS_afr[[I]][Quantile==a1[[I]][1]]), nrow(PGS3_HRS_afr[[I]][Quantile==a1[[I]][2]]),nrow(PGS3_HRS_afr[[I]]))]
        B_HRS_afr[[I]][,K:=1] #number of predictors. Need to check later if this is correct.
        B_HRS_afr[[I]][, LCL:=CI.Rsq(R_sq, k=K, n=N)[3]]
        B_HRS_afr[[I]][, UCL:=CI.Rsq(R_sq, k=K, n=N)[4]]
}

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
saveRDS(PGS3_HRS_afr, file='~/height_prediction/imputed/output/PGS3_HRS_afr.Rds')
saveRDS(results.HRS_afr, file='~/height_prediction/imputed/output/results.HRS_afr.Rds')
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
        labs(title = "HRS_AFR") + ylab(expression(Partial~R^2)) + xlab("European Ancestry Proportion")
        print(myp)
        ggsave(paste0('~/height_prediction/imputed/figs/HRS_afr_error_bars_', I, '.png'))
}

saveRDS(B_HRS_afr, file="~/height_prediction/imputed/output/B_HRS_afr.Rds")
########################################
##now HRS_eur
PGS_HRS_eur<-vector('list', 35)
names(PGS_HRS_eur)<-nam
for(N in nam){
	readRDS(paste0('~/height_prediction/imputed/output/PGS2', N, '.Rds'))[[2]]-> PGS_HRS_eur[[N]]
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

for(N in nam){
	data.table(ID=names(PGS_HRS_eur[[N]]), PGS=unlist(PGS_HRS_eur[[N]]))-> PGS2_HRS_eur[[N]]
	setkey(PGS2_HRS_eur[[N]], ID)
	PGS2_HRS_eur[[N]][Pheno_HRS_eur, nomatch=0]-> PGS2_HRS_eur[[N]]
	PGS2_HRS_eur[[N]][,EUR_ANC:=1]
	PGS2_HRS_eur[[N]][,AGE2:=AGE^2]
	PGS2_HRS_eur[[N]][,HEIGHT:=HEIGHT*100]
	PGS2_HRS_eur[[N]]$SEX<-as.factor(PGS2_HRS_eur[[N]]$SEX)
        dt_f<-PGS2_HRS_eur[[N]][SEX==2]
        dt_m<-PGS2_HRS_eur[[N]][SEX==1]
        sd1_f<-sd(dt_f$HEIGHT)
        m1_f<-mean(dt_f$HEIGHT)
        sd1_m<-sd(dt_m$HEIGHT)
        m1_m<-mean(dt_m$HEIGHT)
        dt_f<-dt_f[HEIGHT>=m1_f-(2*sd1_f)]
        dt_m<-dt_m[HEIGHT>=m1_m-(2*sd1_m)]
        PGS2_HRS_eur[[N]]<-rbind(dt_f, dt_m)
}
nrow(PGS2_HRS_eur[[N]])

lapply(PGS2_HRS_eur, function(X) lm(HEIGHT~SEX, X))-> lm1_HRS_eur
lapply(PGS2_HRS_eur, function(X) lm(HEIGHT~PGS, X))-> lm2_HRS_eur
lapply(PGS2_HRS_eur, function(X) lm(HEIGHT~AGE, X))-> lm3_HRS_eur
lapply(PGS2_HRS_eur, function(X) lm(HEIGHT~AGE2,X))-> lm4_HRS_eur
lapply(PGS2_HRS_eur, function(X) lm(HEIGHT~PGS+AGE, X))-> lm5_HRS_eur
lapply(PGS2_HRS_eur, function(X) lm(HEIGHT~PGS+AGE2,X))-> lm6_HRS_eur
lapply(PGS2_HRS_eur, function(X) lm(HEIGHT~SEX+AGE+AGE2, X))-> lm7_HRS_eur
lapply(PGS2_HRS_eur, function(X) lm(HEIGHT~SEX+AGE+AGE2+PGS, X))-> lm8_HRS_eur

partial_R2_eur<-lapply(nam, function(X) partial.R2(lm7_HRS_eur[[X]],lm8_HRS_eur[[X]])) #
names(partial_R2_eur)<-nam

lapply(PGS2_HRS_eur, function(X) X[, Quantile:='total'][,Med_Eur_Anc:=1])-> PGS3_HRS_eur


B_HRS_eur<-vector('list', length(PGS3_HRS_eur))
names(B_HRS_eur)<-names(PGS3_HRS_eur)

for (I in names(PGS3_HRS_eur)){
        B_HRS_eur[[I]]<-data.table(Quant="total",R_sq=partial_R2_eur[[I]],Med_Eur_Anc=1)
        B_HRS_eur[[I]][,N:=nrow(PGS3_HRS_eur[[I]])]
        B_HRS_eur[[I]][, K:=1] #number of predictors. Need to check later if this is correct.
        B_HRS_eur[[I]][, LCL:=CI.Rsq(R_sq, k=K, n=N)[3]]
        B_HRS_eur[[I]][, UCL:=CI.Rsq(R_sq, k=K, n=N)[4]]
}

### add confidence intervals calculated with bootstrap: https://www.statmethods.net/advstats/bootstrapping.html
results.HRS_eur<-vector('list', length(PGS3_HRS_eur))
names(results.HRS_eur)<-names(PGS3_HRS_eur)


for (I in names(PGS3_HRS_eur)){
        results.HRS_eur[[I]]<- boot(data=PGS2_HRS_eur[[I]], statistic=rsq.R2, R=999, formula1=HEIGHT~SEX+AGE+AGE2, formula2=HEIGHT~SEX+AGE+AGE2+PGS)
        results.HRS_eur[[I]]
        #names(results.HRS_eur[[I]])<-c(a1[[I]], "total")
        cat(I, ' done\n')
}


#95% confidence intervals.
saveRDS(results.HRS_eur, file='~/height_prediction/imputed/output/results.HRS_eur.Rds')
saveRDS(PGS3_HRS_eur, file='~/height_prediction/imputed/output/PGS3_HRS_eur.Rds')

boots.ci.HRS_eur<-lapply(results.HRS_eur, function(X)  boot.ci(X, type = c("norm", 'basic', "perc")))
names(boots.ci.HRS_eur)<-names(results.HRS_eur)

for (I in names(PGS3_HRS_eur)){
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

saveRDS(B_HRS_eur, file="~/height_prediction/imputed/output/B_HRS_eur.Rds")


###################################################
#combine all in a table

readRDS('~/height_prediction/gwas/HRS_afr/output/Nr_SNPs_HRS_afr.Rds')[Name %in% paste0('phys_', nam)]-> hrs_afr
readRDS('~/height_prediction/gwas/HRS_eur/output/Nr_SNPs_HRS.Rds')[Name %in% paste0('phys_', nam)]-> hrs_eur

setkey(hrs_afr, Name)
setkey(hrs_eur, Name)
dt<-data.table(Name=paste0('phys_',nam), HRS_afr_imp=unlist(partial_R2), HRS_eur_imp=unlist(partial_R2_eur), Nr_imp=unlist(lapply(nam, function(X) nrow(do.call(rbind,readRDS(paste0('~/height_prediction/imputed/output/vec_all_', X,'.Rds')))))))
setkey(dt, Name)
setkey(dt, Name)[hrs_afr][hrs_eur]-> dt2
dt2[, HRS_afr:=Part_R2]
dt2[, HRS_eur:=i.Part_R2]
dt2[,i.Nr:=NULL][,Part_R2:=NULL][,i.Part_R2:=NULL]
dt2[,eur_diff:=HRS_eur_imp-HRS_eur]
dt2[,afr_diff:=HRS_afr_imp-HRS_afr]
dt2[,Nr_diff:=Nr_imp-Nr]
saveRDS(dt2,'~/height_prediction/imputed/output/comparison.Rds')
txtStop()
