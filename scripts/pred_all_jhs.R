#!/usr/bin/env Rscript
##preable###
###########
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

a<-c('_p1.0000e-08','_p1.0000e-07','_p1.0000e-06','_p1.0000e-05','_p3.0000e-05','_p1.0000e-04','_p3.0000e-04','_p1.0000e-03','_p3.0000e-03','_p1.0000e-02', '_p3.0000e-02', '_p1.0000e-01', '_p3.0000e-01', '_p1.0000e+00')

a1<-paste0('~/height_prediction/ldpred/output/JHS.score_P+T_r0.20',a, '.txt')

b<-c('_p1.0000e-03', '_p3.0000e-03','_p1.0000e-02', '_p3.0000e-02', '_p1.0000e-01', '_p3.0000e-01','_p1.0000e+00', '-inf')

b1<-paste0('~/height_prediction/ldpred/output/JHS.score_LDpred', b, '.txt')

d<-c(unlist(a1), unlist(b1))
lapply(d , function(X) fread(X))-> PGS_JHS

#read in phenotype data
fread('~/height_prediction/input/JHS/JHS_phenotypes.txt')-> Pheno_JHS

#a partial R2 function
source('~/height_prediction/strat_prs/scripts/Rsq_R2.R')
##############
my_names<-c(paste0('pt', a), paste0('gibbs', b))
names(PGS_JHS)<-my_names

#add PGS to Pheno table in order to be able to make multiple analyses

#Pheno_JHS[, SUBJID:=paste0("0_", as.character(Pheno_JHS[, SUBJID]))]
setkey(Pheno_JHS, SUBJID)
#add ancestry
ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/JHS/rfmix_anc_chr', X, '.txt'))))

anc_JHS<-ancestry %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual
anc_JHS$SUBJID<-substr(anc_JHS[,SUBJID],3,9)
setkey(anc_JHS, SUBJID)

for (I in my_names){
        colnames(PGS_JHS[[I]])[1]<-'SUBJID'
	setkey(PGS_JHS[[I]], SUBJID)
        PGS_JHS[[I]][Pheno_JHS, nomatch=0]-> PGS_JHS[[I]]
        PGS_JHS[[I]][anc_JHS, nomatch=0]-> PGS_JHS[[I]]
        PGS_JHS[[I]][,age2:=age_baseline^2]
        PGS_JHS[[I]][,age:=age_baseline][, age_baseline:=NULL]
        PGS_JHS[[I]][AFR_ANC>=0.05]-> PGS_JHS[[I]]#filter out individuals that are not african...
        PGS_JHS[[I]][-which(is.na(PGS_JHS[[I]][,height_baseline])),]-> PGS_JHS[[I]]
        PGS_JHS[[I]][,HEIGHTX:=height_baseline]
        PGS_JHS[[I]][,height_baseline:=NULL]
        PGS_JHS[[I]][,sex:=sex_selfreport]
        PGS_JHS[[I]][,sex_selfreport:=NULL]
}
lapply(PGS_JHS, function(X) lm(HEIGHTX~PRS, X))-> lm1_JHS
lapply(PGS_JHS, function(X) lm(HEIGHTX~age, X))-> lm2_JHS
lapply(PGS_JHS, function(X) lm(HEIGHTX~age2, X))-> lm3_JHS
lapply(PGS_JHS, function(X) lm(HEIGHTX~EUR_ANC, X))-> lm4_JHS
lapply(PGS_JHS, function(X) lm(HEIGHTX~sex+age, X))-> lm5_JHS
lapply(PGS_JHS, function(X) lm(HEIGHTX~sex+age+age2, X))-> lm6_JHS
lapply(PGS_JHS, function(X) lm(HEIGHTX~sex+age+age2+EUR_ANC, X))-> lm7_JHS
lapply(PGS_JHS, function(X) lm(HEIGHTX~sex+age+age2+EUR_ANC+PRS,X))-> lm8_JHS

partial_r2_JHS<-lapply(1:length(PGS_JHS), function(X) partial.R2(lm7_JHS[[X]], lm8_JHS[[X]])) #min 1.8, max 4.3
names(partial_r2_JHS)<- names(PGS_JHS)


lapply(PGS_JHS, function(X) X[, Quantile:= cut(EUR_ANC,
                                breaks=quantile(EUR_ANC, probs=seq(0,1, by=1/2), na.rm=TRUE),
                                include.lowest=TRUE)])-> PGS2_JHS

lapply(1:length(PGS2_JHS), function(X) PGS2_JHS[[X]][,Med_Eur_Anc:=median(EUR_ANC),by=Quantile])

lapply(1:length(PGS2_JHS), function(X) as.character(unique((PGS2_JHS[[X]]$Quantile))))-> a

lapply(a, function(X) c(X[2], X[1]))-> a1

names(a1)<-names(PGS2_JHS)

r2_JHS<-vector('list', length(PGS2_JHS))
names(r2_JHS)<-names(PGS2_JHS)

for(I in names(r2_JHS)){
        r2_JHS[[I]]<-vector('list', length(a1[[I]]))
        names(r2_JHS[[I]])<-a1[[I]]
        for(i in a1[[I]]){
        r2_JHS[[I]][[i]]<-partial.R2(lm(HEIGHTX~sex+age+age2+EUR_ANC, PGS2_JHS[[I]][Quantile==i]),lm(HEIGHTX~sex+age+age2+EUR_ANC+PRS, PGS2_JHS[[I]][Quantile==i]))
        }
}


B_JHS<-vector('list', length(r2_JHS))
names(B_JHS)<-names(r2_JHS)

for (I in names(r2_JHS)){
        B_JHS[[I]]<-data.table(Quant=c(a1[[I]], "total"),
        R_sq=c(unlist(r2_JHS[[I]]), partial_r2_JHS[[I]]),
        Med_Eur_Anc=c(unique(PGS2_JHS[[I]][Quantile==a1[[I]][1]][,Med_Eur_Anc]), unique(PGS2_JHS[[I]][Quantile==a1[[I]][2]][,Med_Eur_Anc]), median(PGS2_JHS[[I]][, EUR_ANC])))
        B_JHS[[I]][,N:=c(nrow(PGS2_JHS[[I]][Quantile==a1[[I]][1]]), nrow(PGS2_JHS[[I]][Quantile==a1[[I]][2]]),nrow(PGS2_JHS[[I]]))]
        B_JHS[[I]][,K:=1] #number of predictors. Need to check later if this is correct.
        B_JHS[[I]][, LCL:=CI.Rsq(R_sq, k=K, n=N)[3]]
        B_JHS[[I]][, UCL:=CI.Rsq(R_sq, k=K, n=N)[4]]
}

### add confidence intervals calculated with bootstrap: https://www.statmethods.net/advstats/bootstrapping.html
results.JHS<-vector('list', length(PGS2_JHS))
names(results.JHS)<-names(PGS2_JHS)
for (I in names(PGS2_JHS)){
        results.JHS[[I]]<-vector('list', length(a1[[I]])+1)
        names(results.JHS[[I]])<-c(a1[[I]], "total")
        lapply(a1[[I]], function(i) boot(data=PGS2_JHS[[I]][Quantile==i], statistic=rsq.R2,R=1000, formula1=HEIGHTX~sex+age+age2+EUR_ANC, formula2=HEIGHTX~sex+age+age2+EUR_ANC+PRS))-> results.JHS[[I]]
        cat(I)
        cat(' done\n')
}
#saveRDS(PGS3_JHS, file='/project/mathilab/bbita/gwas_admix/new_height/JHS/PGS3_JHS.Rds')

#saveRDS(results.JHS, file='/project/mathilab/bbita/gwas_admix/new_height/JHS/results.JHS.Rds')

for (I in names(PGS2_JHS)){
        tp <- boot(data=PGS_JHS[[I]], statistic=rsq.R2, R=1000, formula1=HEIGHTX~sex+age+age2+EUR_ANC, formula2=HEIGHTX~sex+age+age2+EUR_ANC+PRS)
        list.append(results.JHS[[I]], tp)-> results.JHS[[I]]
        names(results.JHS[[I]])<-c(a1[[I]], "total")
        cat(I)
        cat(' done\n')
}
saveRDS(PGS2_JHS, file='~/height_prediction/ldpred/output/PGS3_JHS.Rds')
saveRDS(results.JHS, file='~/height_prediction/ldpred/output/results.JHS.Rds')

#confidence intervals


boots.ci.JHS<-lapply(results.JHS, function(Y) lapply(Y, function(X) boot.ci(X, type = c("norm", 'basic', "perc"))))
names(boots.ci.JHS)<-names(results.JHS)

for (I in names(PGS2_JHS)){
        B_JHS[[I]][1:2,]-> a
        B_JHS[[I]][3,]-> b
        a[,HVB_L:=sapply(a$Quant, function(X) as.numeric(gsub("\\]","",gsub("\\(","",gsub("\\[","",strsplit(X,",")[[1]])))))[1,]]
        a[,HVB_U:=sapply(a$Quant, function(X) as.numeric(gsub("\\]","",gsub("\\(","",gsub("\\[","",strsplit(X,",")[[1]])))))[2,]]
        b[,HVB_L:=1]
        b[,HVB_U:=1]
        rbind(a,b)->B_JHS[[I]]
        B_JHS[[I]][, Dataset:='JHS_AA']
        B_JHS[[I]][, boots_norm_L:=sapply(1:3, function(X) boots.ci.JHS[[I]][[X]]$normal[2])]
        B_JHS[[I]][, boots_norm_U:=sapply(1:3, function(X) boots.ci.JHS[[I]][[X]]$normal[3])]
        B_JHS[[I]][, boots_perc_L:=sapply(1:3, function(X) boots.ci.JHS[[I]][[X]]$perc[4])]
        B_JHS[[I]][, boots_perc_U:=sapply(1:3, function(X) boots.ci.JHS[[I]][[X]]$perc[5])]
        B_JHS[[I]][, boots_basic_L:=sapply(1:3, function(X) boots.ci.JHS[[I]][[X]]$basic[4])]
        B_JHS[[I]][, boots_basic_U:=sapply(1:3, function(X) boots.ci.JHS[[I]][[X]]$basic[5])]
}

saveRDS(B_JHS, file="~/height_prediction/ldpred/output/B_JHS.Rds")
