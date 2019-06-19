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
library(readr)
#############
#############

#read in PGS scores
readRDS('~/height_prediction/gwas/WHI/output/PGS_WHI.Rds')-> PGS_WHI

#read in phenotype data
fread('~/height_prediction/input/WHI/WHI_phenotypes.txt')-> Pheno_WHI

#a partial R2 function
source('~/height_prediction/strat_prs/scripts/Rsq_R2.R')
##############

#add PGS to Pheno table in order to be able to make multiple analyses

Pheno_WHI[, SUBJID:=paste0("0_", as.character(Pheno_WHI[, SUBJID]))]
setkey(Pheno_WHI, SUBJID)

ancestry<-vector('list', 22)
#add ancestry
for(chr in 1:22){
	what <- paste0("~/height_prediction/input/WHI/WHI_b37_strand_include_kgCY_chr", chr)  #pop1 is AFR
	ancestry[[chr]] <- read_fwf(paste0(what, "_rfmix_out.0.Viterbi.txt"), fwf_empty(paste0(what, "_rfmix_out.0.Viterbi.txt")))
	gc()
}

ind<-read.table(paste0(what, ".phind"), as.is=TRUE)
#Remove reference samples
samples.to.include <- ind[,3]=="Case"
ind <- ind[samples.to.include,]

indiv<-vector('list', ncol(ancestry[[22]]))
for(I in 1:ncol(ancestry[[22]])){	
	cat(I,"\r")
	indiv[[I]]<-(sum(ancestry[[1]][,I]==1)+ sum(ancestry[[2]][,I]==1)+ sum(ancestry[[3]][,I]==1)+sum(ancestry[[4]][,I]==1)+sum(ancestry[[5]][,I]==1)+sum(ancestry[[6]][,I]==1)+sum(ancestry[[7]][,I]==1)+sum(ancestry[[8]][,I]==1)+sum(ancestry[[9]][,I]==1)+sum(ancestry[[10]][,I]==1)+sum(ancestry[[11]][,I]==1)+sum(ancestry[[12]][,I]==1)+sum(ancestry[[13]][,I]==1)+sum(ancestry[[14]][,I]==1)+sum(ancestry[[15]][,I]==1)+sum(ancestry[[16]][,I]==1)+sum(ancestry[[17]][,I]==1)+sum(ancestry[[18]][,I]==1)+sum(ancestry[[19]][,I]==1)+sum(ancestry[[20]][,I]==1)+sum(ancestry[[21]][,I]==1)+ sum(ancestry[[22]][,I]==1))/sum(unlist(lapply(1:22, function(X) nrow(ancestry[[X]]))))
}

indiv2<-vector('list', ncol(ancestry[[22]])/2)

i=0
for(I in seq(from=1,to=length(indiv), by=2)){
	i=i+1	
	indiv2[[i]]<-sum(indiv[[I]],indiv[[I+1]])/2
	cat(i,"\r")
	}

#admixture (obsolete)
#anc_WHI<-cbind(fread('~/height_prediction/input/WHI/WHI_b37_strand_prune_include.2.Q'), fread('~/height_prediction/input/WHI/WHI_b37_strand_prune_include.fam')[,V2])
#colnames(anc_WHI)<-c("AFR_ANC","EUR_ANC","SUBJID")
anc_WHI<-data.table(AFR_ANC=unlist(indiv2), EUR_ANC=1-unlist(indiv2), SUBJID=as.integer(unique(gsub("_[A|B]", "",gsub("0:","",ind[,1])))))
anc_WHI[,SUBJID:=paste0("0_", SUBJID)]
setkey(anc_WHI, SUBJID)


PGS2_WHI<-vector('list', length(PGS_WHI))
names(PGS2_WHI)<-names(PGS_WHI)

for (I in names(PGS2_WHI)){
        data.table(SUBJID=names(PGS_WHI[[I]]), PGS=unlist(PGS_WHI[[I]]))-> PGS2_WHI[[I]]
        setkey(PGS2_WHI[[I]], SUBJID)
        PGS2_WHI[[I]][Pheno_WHI, nomatch=0]-> PGS2_WHI[[I]]
        PGS2_WHI[[I]][anc_WHI, nomatch=0]-> PGS2_WHI[[I]]
        PGS2_WHI[[I]][,age2:=AGE^2]
        PGS2_WHI[[I]][AFR_ANC>=0.05]-> PGS2_WHI[[I]] #filter out individuals that are not african...
        PGS2_WHI[[I]][-which(is.na(PGS2_WHI[[I]][,HEIGHTX])),]-> PGS2_WHI[[I]]
}

lapply(PGS2_WHI, function(X) lm(HEIGHTX~PGS, X))-> lm1_WHI
lapply(PGS2_WHI, function(X) lm(HEIGHTX~AGE, X))-> lm2_WHI
lapply(PGS2_WHI, function(X) lm(HEIGHTX~age2, X))-> lm3_WHI
lapply(PGS2_WHI, function(X) lm(HEIGHTX~EUR_ANC, X))-> lm4_WHI
lapply(PGS2_WHI, function(X) lm(HEIGHTX~PGS+AGE, X))-> lm5_WHI
lapply(PGS2_WHI, function(X) lm(HEIGHTX~PGS+age2, X))-> lm6_WHI
lapply(PGS2_WHI, function(X) lm(HEIGHTX~AGE+age2+EUR_ANC, X))-> lm7_WHI
lapply(PGS2_WHI, function(X) lm(HEIGHTX~AGE+age2+EUR_ANC+PGS, X))-> lm8_WHI


partial.R2(lm7_WHI[[67]],lm8_WHI[[67]]) #2.1%

partial_r2_WHI<-lapply(1:length(PGS2_WHI), function(X) partial.R2(lm7_WHI[[X]], lm8_WHI[[X]])) #min 2.27, max 4.83%
names(partial_r2_WHI)<- names(PGS2_WHI)

summary(unlist(partial_r2_WHI)) #min 1 ma 2.4%

which.max(unlist(partial_r2_WHI)) #

Nr_SNPs<-rep(NA, length(PGS2_WHI))
names(Nr_SNPs)<- names(PGS2_WHI)

for(I in names(Nr_SNPs)){
	readRDS(paste0('~/height_prediction/gwas/WHI/output/hei_', I, '.Rds'))-> smtg
	sum(sapply(smtg, function(X) nrow(X)))-> Nr_SNPs[I]
	cat(I, ' \n')
}
data.table(Nr=unlist(Nr_SNPs), Name=names(Nr_SNPs), Part_R2=unlist(partial_r2_WHI))-> A_table
saveRDS(A_table, file='~/height_prediction/gwas/WHI/utput/Nr_SNPs_WHI.Rds')

cor(A_table$Part_R2, A_table$Nr) #-0.589 #strong negative correlation
summary(lm(Part_R2~Nr,data=A_table))$r.squared #0.3375
#

###########################
######### Plots ###########
###########################

for(I in 1:length(PGS2_WHI)){
        A<-ggpairs(PGS2_WHI[[I]][,.(HEIGHTX, PGS, AGE, age2,EUR_ANC,WAISTX,WEIGHTX)])
        png(filename=paste0('~/height_prediction/gwas/WHI/figs/WHI_ggpairs_', names(PGS2_WHI)[I], '.png'))
        print(A)
        cat(names(PGS2_WHI)[I])
        cat(' done\n')
        dev.off()
}


lapply(PGS2_WHI, function(X) X[, Quantile:= cut(EUR_ANC,
                                breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.25), na.rm=TRUE),
                                include.lowest=TRUE)])-> PGS3_WHI

lapply(1:length(PGS3_WHI), function(X) PGS3_WHI[[X]][,Med_Eur_Anc:=median(EUR_ANC),by=Quantile])

lapply(1:length(PGS3_WHI), function(X) as.character(unique((PGS3_WHI[[X]]$Quantile))))-> a

lapply(a, function(X) c(X[3],X[4], X[1], X[2]))-> a1

names(a1)<-names(PGS3_WHI)

r2_WHI<-vector('list', length(PGS3_WHI))
names(r2_WHI)<-names(PGS3_WHI)

for(I in names(r2_WHI)){

        r2_WHI[[I]]<-vector('list', length(a1[[I]]))
        names(r2_WHI[[I]])<-a1[[I]]
        for(i in a1[[I]]){
        r2_WHI[[I]][[i]]<-partial.R2(lm(HEIGHTX~AGE+age2+EUR_ANC, PGS3_WHI[[I]][Quantile==i]),lm(HEIGHTX~AGE+age2+EUR_ANC+PGS, PGS3_WHI[[I]][Quantile==i]))
        }
}

B_WHI<-vector('list', length(r2_WHI))
names(B_WHI)<-names(r2_WHI)


for (I in names(r2_WHI)){
	B_WHI[[I]]<-data.table(Quant=c(a1[[I]], "total"),
        R_sq=c(unlist(r2_WHI[[I]]), partial_r2_WHI[[I]]),
        Med_Eur_Anc=c(unique(PGS3_WHI[[I]][Quantile==a1[[I]][1]][,Med_Eur_Anc]), unique(PGS3_WHI[[I]][Quantile==a1[[I]][2]][,Med_Eur_Anc]),unique(PGS3_WHI[[I]][Quantile==a1[[I]][3]][,Med_Eur_Anc]),unique(PGS3_WHI[[I]][Quantile==a1[[I]][4]][,Med_Eur_Anc]), median(PGS3_WHI[[I]][, EUR_ANC])))
        B_WHI[[I]][,N:=c(nrow(PGS3_WHI[[I]][Quantile==a1[[I]][1]]), nrow(PGS3_WHI[[I]][Quantile==a1[[I]][2]]), nrow(PGS3_WHI[[I]][Quantile==a1[[I]][3]]), nrow(PGS3_WHI[[I]][Quantile==a1[[I]][4]]),nrow(PGS3_WHI[[I]]))]
        B_WHI[[I]][,K:=1] #number of predictors. Need to check later if this is correct.
        B_WHI[[I]][, LCL:=CI.Rsq(R_sq, k=K, n=N)[3]]
        B_WHI[[I]][, UCL:=CI.Rsq(R_sq, k=K, n=N)[4]]
}

### add confidence intervals calculated with bootstrap: https://www.statmethods.net/advstats/bootstrapping.html
results.WHI<-vector('list', length(PGS3_WHI))

for (I in names(PGS3_WHI)){
	results.WHI[[I]]<-vector('list', length(a1[[I]])+1)
	names(results.WHI[[I]])<-c(a1[[I]], "total")
        lapply(a1[[I]], function(i) boot(data=PGS3_WHI[[I]][Quantile==i], statistic=rsq.R2,R=999, formula1=HEIGHTX~AGE+age2+EUR_ANC, formula2=HEIGHTX~AGE+age2+EUR_ANC+PGS))-> results.WHI[[I]]
        cat(I)
        cat(' done\n')
}

for (I in names(PGS3_WHI)){
	tp <- boot(data=PGS2_WHI[[I]], statistic=rsq.R2, R=999, formula1=HEIGHTX~AGE+age2+EUR_ANC, formula2=HEIGHTX~AGE+age2+EUR_ANC+PGS)
	list.append(results.WHI[[I]], tp)-> results.WHI[[I]]
	names(results.WHI[[I]])<-c(a1[[I]], "total")
        cat(I)
	cat(' done\n')
}
results.WHI<-results.WHI[81:160]
saveRDS(PGS3_WHI, file='~/height_prediction/gwas/WHI/output/PGS3_WHI.Rds')
saveRDS(results.WHI, file='~/height_prediction/gwas/WHI/output/results.WHI.Rds')

#confidence intervals

boots.ci.WHI<-lapply(results.WHI, function(Y) lapply(Y, function(X) boot.ci(X, type = c("norm", 'basic', "perc"))))
names(boots.ci.WHI)<-names(results.WHI)

for (I in names(PGS3_WHI)){
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

for (I in names(PGS3_WHI)){
        myp<-ggplot(B_WHI[[I]][1:4,], aes(x=Med_Eur_Anc, y=R_sq)) +
        geom_point(size=1.5, shape=21, fill="white") +
        geom_errorbar(aes(x=Med_Eur_Anc, ymin=boots_perc_L, ymax=boots_perc_U), width=0.05, size=0.8, color="cornflowerblue") +
        geom_line(color='lightgray')+
        geom_errorbarh(aes(x=Med_Eur_Anc, xmin=HVB_L, xmax=HVB_U), width=0.05, size=0.5, color="cornflowerblue") +
        labs(title = "WHI_AA") + ylab("R-squared") + xlab("European Ancestry Proportion")
        print(myp)
        ggsave(paste0('~/height_prediction/gwas/WHI/figs/error_bars_', I, '.png'))
}

saveRDS(B_WHI, file="~/height_prediction/gwas/WHI/output/B_WHI.Rds")
