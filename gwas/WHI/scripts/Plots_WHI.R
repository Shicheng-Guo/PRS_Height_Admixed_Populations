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

#add ancestry
ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/WHI/rfmix_anc_chr', X, '.txt'))))

anc_WHI<-ancestry %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual

#admixture (obsolete)
#anc_WHI<-cbind(fread('~/height_prediction/input/WHI/WHI_b37_strand_prune_include.2.Q'), fread('~/height_prediction/input/WHI/WHI_b37_strand_prune_include.fam')[,V2])
#colnames(anc_WHI)<-c("AFR_ANC","EUR_ANC","SUBJID")
#anc_WHI<-data.table(AFR_ANC=unlist(indiv2), EUR_ANC=1-unlist(indiv2), SUBJID=as.integer(unique(gsub("_[A|B]", "",gsub("0:","",ind[,1])))))
#anc_WHI[,SUBJID:=paste0("0_", SUBJID)]
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


partial.R2(lm7_WHI[[67]],lm8_WHI[[67]]) #4.7

partial_r2_WHI<-lapply(1:length(PGS2_WHI), function(X) partial.R2(lm7_WHI[[X]], lm8_WHI[[X]])) 
names(partial_r2_WHI)<- names(PGS2_WHI)

summary(unlist(partial_r2_WHI)) ##min 2.27, max 5.6%

which.max(unlist(partial_r2_WHI)) #77, or LD_100000_0.01_0.5. Excluding the LD ones, it's 68, phys_5000_0.0005

Nr_SNPs<-rep(NA, length(PGS2_WHI))
names(Nr_SNPs)<- names(PGS2_WHI)

for(I in names(Nr_SNPs)){
	readRDS(paste0('~/height_prediction/gwas/WHI/output/hei_', I, '.Rds'))-> smtg
	sum(sapply(smtg, function(X) nrow(X)))-> Nr_SNPs[I]
	cat(I, ' \n')
}
data.table(Nr=unlist(Nr_SNPs), Name=names(Nr_SNPs), Part_R2=unlist(partial_r2_WHI))-> A_table
saveRDS(A_table, file='~/height_prediction/gwas/WHI/output/Nr_SNPs_WHI.Rds')

cor(A_table$Part_R2, A_table$Nr) #0.89
summary(lm(Part_R2~Nr,data=A_table))$r.squared #0.79
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

lapply(a, function(X) c(X[4],X[2], X[1], X[3]))-> a1

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
saveRDS(PGS3_WHI, file='~/height_prediction/gwas/WHI/output/PGS3_WHI.Rds')
saveRDS(results.WHI, file='~/height_prediction/gwas/WHI/output/results.WHI.Rds')

#confidence intervals

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
