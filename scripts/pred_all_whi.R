#!/usr/bin/env Rscript
##preable###
###########
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
options(scipen=999)

a<-c('_p1.0000e-08','_p1.0000e-07','_p1.0000e-06','_p1.0000e-05','_p3.0000e-05','_p1.0000e-04','_p3.0000e-04','_p1.0000e-03','_p3.0000e-03','_p1.0000e-02', '_p3.0000e-02', '_p1.0000e-01', '_p3.0000e-01', '_p1.0000e+00')

a1<-paste0('~/height_prediction/ldpred/output/WHI.score_P+T_r0.20',a, '.txt')

b<-c('_p1.0000e-03', '_p3.0000e-03','_p1.0000e-02', '_p3.0000e-02', '_p1.0000e-01', '_p3.0000e-01','_p1.0000e+00', '-inf')

b1<-paste0('~/height_prediction/ldpred/output/WHI.score_LDpred', b, '.txt')

d<-c(unlist(a1), unlist(b1))
lapply(d , function(X) fread(X))-> PGS_WHI
lapply(PGS_WHI , function(X) X[,IID:=paste0('0_', IID)])


#read in phenotype data
fread('~/height_prediction/input/WHI/WHI_phenotypes.txt')-> Pheno_WHI

#a partial R2 function
source('~/height_prediction/strat_prs/scripts/Rsq_R2.R')
#a#############
my_names<-c(paste0('pt', a), paste0('gibbs', b))
names(PGS_WHI)<-my_names

#add PGS to Pheno table in order to be able to make multiple analyses

Pheno_WHI[, SUBJID:=paste0("0_", as.character(Pheno_WHI[, SUBJID]))]
setkey(Pheno_WHI, SUBJID)

#add ancestry
ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/WHI/rfmix_anc_chr', X, '.txt'))))

anc_WHI<-ancestry %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual

#admixture (obsolete)
#anc_WHI<-cbind(fread('~/height_prediction/input/WHI/WHI_b37_strand_prune_include.2.Q'), fread('~/height_prediction/input/WHI/WHI_b37_strand_prune_include.fam')[,V2])
#colnames(anc_WHI)<-c("AFR_ANC","EUR_ANC","SUBJID")
#anc_WHI<-data.table(AFR_ANC=unlist(indiv2), EUR_ANC=1-unlist(indiv2), SUBJID=as.integer(unique(gsub("_[A|B]", "",gsub("0:","",ind[,1])))))
#anc_WHI[,SUBJID:=paste0("0_", SUBJID)]
setkey(anc_WHI, SUBJID)

for (I in my_names){
        colnames(PGS_WHI[[I]])[1]<-'SUBJID'
	setkey(PGS_WHI[[I]], SUBJID)
        PGS_WHI[[I]][Pheno_WHI, nomatch=0]-> PGS_WHI[[I]]
        PGS_WHI[[I]][anc_WHI, nomatch=0]-> PGS_WHI[[I]]
        PGS_WHI[[I]][,age2:=AGE^2]
        PGS_WHI[[I]][AFR_ANC>=0.05]-> PGS_WHI[[I]] #filter out individuals that are not african...
        PGS_WHI[[I]][-which(is.na(PGS_WHI[[I]][,HEIGHTX])),]-> PGS_WHI[[I]]
	m1<-mean(PGS_WHI[[I]]$HEIGHTX)
	sd1<-sd(PGS_WHI[[I]]$HEIGHTX)
	PGS_WHI[[I]]<-PGS_WHI[[I]][HEIGHTX>=m1-(2*sd1) & HEIGHTX<=m1+(2*sd1)]
}

lapply(PGS_WHI, function(X) lm(HEIGHTX~PRS, X))-> lm1_WHI
lapply(PGS_WHI, function(X) lm(HEIGHTX~AGE, X))-> lm2_WHI
lapply(PGS_WHI, function(X) lm(HEIGHTX~age2, X))-> lm3_WHI
lapply(PGS_WHI, function(X) lm(HEIGHTX~EUR_ANC, X))-> lm4_WHI
lapply(PGS_WHI, function(X) lm(HEIGHTX~PRS+AGE, X))-> lm5_WHI
lapply(PGS_WHI, function(X) lm(HEIGHTX~PRS+age2, X))-> lm6_WHI
lapply(PGS_WHI, function(X) lm(HEIGHTX~AGE+age2+EUR_ANC, X))-> lm7_WHI
lapply(PGS_WHI, function(X) lm(HEIGHTX~AGE+age2+EUR_ANC+PRS, X))-> lm8_WHI

partial_r2_WHI<-lapply(1:length(PGS_WHI), function(X) partial.R2(lm7_WHI[[X]], lm8_WHI[[X]])) 
names(partial_r2_WHI)<- names(PGS_WHI)

###########################
lapply(PGS_WHI, function(X) X[, Quantile:= cut(EUR_ANC,
                                breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.25), na.rm=TRUE),
                                include.lowest=TRUE)])-> PGS2_WHI

lapply(1:length(PGS2_WHI), function(X) PGS2_WHI[[X]][,Med_Eur_Anc:=median(EUR_ANC),by=Quantile])

lapply(1:length(PGS2_WHI), function(X) as.character(unique((PGS2_WHI[[X]]$Quantile))))-> a

lapply(a, function(X) c(X[4],X[2], X[1], X[3]))-> a1

names(a1)<-names(PGS2_WHI)

r2_WHI<-vector('list', length(PGS2_WHI))
names(r2_WHI)<-names(PGS2_WHI)

for(I in names(r2_WHI)){
        r2_WHI[[I]]<-vector('list', length(a1[[I]]))
        names(r2_WHI[[I]])<-a1[[I]]
        for(i in a1[[I]]){
        r2_WHI[[I]][[i]]<-partial.R2(lm(HEIGHTX~AGE+age2+EUR_ANC, PGS2_WHI[[I]][Quantile==i]),lm(HEIGHTX~AGE+age2+EUR_ANC+PRS, PGS2_WHI[[I]][Quantile==i]))
        }
}

B_WHI<-vector('list', length(r2_WHI))
names(B_WHI)<-names(r2_WHI)


for (I in names(r2_WHI)){
	B_WHI[[I]]<-data.table(Quant=c(a1[[I]], "total"),
        R_sq=c(unlist(r2_WHI[[I]]), partial_r2_WHI[[I]]),
        Med_Eur_Anc=c(unique(PGS2_WHI[[I]][Quantile==a1[[I]][1]][,Med_Eur_Anc]), unique(PGS2_WHI[[I]][Quantile==a1[[I]][2]][,Med_Eur_Anc]),unique(PGS2_WHI[[I]][Quantile==a1[[I]][3]][,Med_Eur_Anc]),unique(PGS2_WHI[[I]][Quantile==a1[[I]][4]][,Med_Eur_Anc]), median(PGS2_WHI[[I]][, EUR_ANC])))
        B_WHI[[I]][,N:=c(nrow(PGS2_WHI[[I]][Quantile==a1[[I]][1]]), nrow(PGS2_WHI[[I]][Quantile==a1[[I]][2]]), nrow(PGS2_WHI[[I]][Quantile==a1[[I]][3]]), nrow(PGS2_WHI[[I]][Quantile==a1[[I]][4]]),nrow(PGS2_WHI[[I]]))]
        B_WHI[[I]][,K:=1] #number of predictors. Need to check later if this is correct.
        B_WHI[[I]][, LCL:=CI.Rsq(R_sq, k=K, n=N)[3]]
        B_WHI[[I]][, UCL:=CI.Rsq(R_sq, k=K, n=N)[4]]
}

### add confidence intervals calculated with bootstrap: https://www.statmethods.net/advstats/bootstrapping.html

results.WHI<-vector('list', length(PGS2_WHI))

for (I in names(PGS2_WHI)){
plot_grid(A_plot, B_plot, nrow=2, labels=c('A', 'B'))
	results.WHI[[I]]<-vector('list', length(a1[[I]])+1)
	names(results.WHI[[I]])<-c(a1[[I]], "total")
        lapply(a1[[I]], function(i) boot(data=PGS2_WHI[[I]][Quantile==i], statistic=rsq.R2,R=1000, formula1=HEIGHTX~AGE+age2+EUR_ANC, formula2=HEIGHTX~AGE+age2+EUR_ANC+PRS))-> results.WHI[[I]]
        cat(I)
        cat(' done\n')
}

for (I in names(PGS2_WHI)){
	tp <- boot(data=PGS2_WHI[[I]], statistic=rsq.R2, R=1000, formula1=HEIGHTX~AGE+age2+EUR_ANC, formula2=HEIGHTX~AGE+age2+EUR_ANC+PRS)
	list.append(results.WHI[[I]], tp)-> results.WHI[[I]]
	names(results.WHI[[I]])<-c(a1[[I]], "total")
        cat(I)
	cat(' done\n')
}
saveRDS(PGS2_WHI, file='~/height_prediction/ldpred/output/PGS3_WHI.Rds')
saveRDS(results.WHI, file='~/height_prediction/ldpred/output/results.WHI.Rds')

#confidence intervals
boots.ci.WHI<-lapply(results.WHI, function(Y) lapply(Y, function(X) boot.ci(X, type = c("norm", 'basic', "perc"))))
names(boots.ci.WHI)<-names(results.WHI)

for (I in names(PGS2_WHI)){
        B_WHI[[I]][1:4,]-> a
        B_WHI[[I]][5,]-> b
        a[,HVB_L:=sapply(a$Quant, function(X) as.numeric(gsub("\\]","",gsub("\\(","",gsub("\\[","",strsplit(X,",")[[1]])))))[1,]]
        a[,HVB_U:=sapply(a$Quant, function(X) as.numeric(gsub("\\]","",gsub("\\(","",gsub("\\[","",strsplit(X,",")[[1]])))))[2,]]
        b[,HVB_L:=1]
        b[,HVB_U:=1]
        rbind(a,b)->B_WHI[[I]]
        B_WHI[[I]][, Dataset:='WHI_AA']
        B_WHI[[I]][, boots_norm_L:=sapply(1:5, function(X) boots.ci.WHI[[I]][[X]]$normal[2])]
        B_WHI[[I]][, boots_norm_U:=sapply(1:5, function(X) boots.ci.WHI[[I]][[X]]$normal[3])]
        B_WHI[[I]][, boots_perc_L:=sapply(1:5, function(X) boots.ci.WHI[[I]][[X]]$perc[4])]
        B_WHI[[I]][, boots_perc_U:=sapply(1:5, function(X) boots.ci.WHI[[I]][[X]]$perc[5])]
        B_WHI[[I]][, boots_basic_L:=sapply(1:5, function(X) boots.ci.WHI[[I]][[X]]$basic[4])]
        B_WHI[[I]][, boots_basic_U:=sapply(1:5, function(X) boots.ci.WHI[[I]][[X]]$basic[5])]
}
saveRDS(B_WHI, file="~/height_prediction/ldpred/output/B_WHI.Rds")
