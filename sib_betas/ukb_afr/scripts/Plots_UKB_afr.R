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

txtStart(paste0("~/height_prediction/sib_betas/ukb_afr/plots_out.txt"))
#read in PGS scores
readRDS('~/height_prediction/sib_betas/ukb_afr/output/PGS_ukb_afr.Rds')-> PGS_UKB_afr
#read in phenotype data
fread('~/height_prediction/input/ukb_afr/UKB_AFR_pheno.txt')-> Pheno_UKB_afr
#a partial R2 function
source('~/height_prediction/strat_prs/scripts/Rsq_R2.R')

##############

#add PGS to Pheno table in order to be able to make multiple analyses
#Pheno_UKB_afr[,ID:=paste0(ID, "_", ID)]
for(j in 1:length(PGS_UKB_afr)){
	gsub("[0-9]+_","",names(PGS_UKB_afr[[j]]))-> names(PGS_UKB_afr[[j]])
}
as.character(Pheno_UKB_afr$ID)-> Pheno_UKB_afr$ID
setkey(Pheno_UKB_afr, ID)

#add ancestry
ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/ukb_afr/rfmix_anc_chr', X, '.txt'))))

anc_ukb_afr<-ancestry %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual

anc_ukb_afr[,ID:=SUBJID][,SUBJID:=NULL]
as.character(anc_ukb_afr$ID)-> anc_ukb_afr$ID
setkey(anc_ukb_afr, ID)

PGS2_UKB_afr<-vector('list', length(PGS_UKB_afr))
names(PGS2_UKB_afr)<-names(PGS_UKB_afr)


for (I in names(PGS_UKB_afr)){
        data.table(ID=names(PGS_UKB_afr[[I]]), PGS=unlist(PGS_UKB_afr[[I]]))-> PGS2_UKB_afr[[I]]
        setkey(PGS2_UKB_afr[[I]], ID)
        PGS2_UKB_afr[[I]][Pheno_UKB_afr, nomatch=0]-> PGS2_UKB_afr[[I]]
        PGS2_UKB_afr[[I]][anc_ukb_afr, nomatch=0]-> PGS2_UKB_afr[[I]]
        PGS2_UKB_afr[[I]][,age2:=Age^2]
	PGS2_UKB_afr[[I]][AFR_ANC>=0.05]-> PGS2_UKB_afr[[I]]
        PGS2_UKB_afr[[I]][-which(is.na(PGS2_UKB_afr[[I]][,Height])),]-> PGS2_UKB_afr[[I]]
	PGS2_UKB_afr[[I]][ID!="6007195"]-> PGS2_UKB_afr[[I]]
}

lapply(PGS2_UKB_afr, function(X) lm(Height~Sex, X))-> lm0_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~PGS, X))-> lm1_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~Age, X))-> lm2_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~age2, X))-> lm3_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~EUR_ANC, X))-> lm4_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~PGS+Age, X))-> lm5_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~PGS+age2, X))-> lm6_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~Sex+Age+age2+EUR_ANC, X))-> lm7_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~Sex+Age+age2+EUR_ANC+PGS, X))-> lm8_UKB_afr

partial.R2(lm7_UKB_afr[[63]],lm8_UKB_afr[[63]]) #0.0317
partial.R2(lm7_UKB_afr[[67]],lm8_UKB_afr[[67]]) #0.042

partial_r2_UKB_afr<-lapply(1:length(PGS2_UKB_afr), function(X) partial.R2(lm7_UKB_afr[[X]], lm8_UKB_afr[[X]])) 
names(partial_r2_UKB_afr)<-names(PGS2_UKB_afr)


Nr_SNPs<-rep(NA, length(PGS2_UKB_afr))
names(Nr_SNPs)<- names(PGS2_UKB_afr)

for(I in names(Nr_SNPs)){
        readRDS(paste0('~/height_prediction/sib_betas/ukb_afr/output/hei_', I, '.Rds'))-> smtg
        sum(sapply(smtg, function(X) nrow(X)))-> Nr_SNPs[I]
	remove(smtg)
	gc()
        cat(I, ' \n')
}

data.table(Nr=unlist(Nr_SNPs), Name=names(Nr_SNPs), Part_R2=unlist(partial_r2_UKB_afr))-> A_table
saveRDS(A_table, file='~/height_prediction/sib_betas/ukb_afr/output/Nr_SNPs_UKB_afr.Rds')

partial_r2_UKB_afr[which.max(partial_r2_UKB_afr)]
#
cor.test(unlist(Nr_SNPs), unlist(partial_r2_UKB_afr))# 0.77
summary(lm(Part_R2~Nr,data=A_table))$r.squared #0.59

for(I in 1:length(PGS2_UKB_afr)){
        A<-ggpairs(PGS2_UKB_afr[[I]][,.(Height, Sex, PGS, Age, age2,EUR_ANC)])
        png(filename=paste0('~/height_prediction/sib_betas/ukb_afr/figs/UKB_afr_ggpairs_', names(PGS2_UKB_afr)[I], '.png'))
        print(A)
        cat(names(PGS2_UKB_afr)[I])
        cat(' done\n')
        dev.off()
}


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
        r2_UKB_afr[[I]][[i]]<-partial.R2(lm(Height~Sex+Age+age2+EUR_ANC, PGS3_UKB_afr[[I]][Quantile==i]),lm(Height~Sex+Age+age2+EUR_ANC+PGS, PGS3_UKB_afr[[I]][Quantile==i]))
        }
}

B_UKB_afr<-vector('list', length(r2_UKB_afr))
names(B_UKB_afr)<-names(r2_UKB_afr)

for (I in names(r2_UKB_afr)){
        B_UKB_afr[[I]]<-data.table(Quant=c(a1[[I]], "total"),
        R_sq=c(unlist(r2_UKB_afr[[I]]), partial_r2_UKB_afr[[I]]),
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
                R=999, formula1=Height~Sex+Age+age2+EUR_ANC, formula2=Height~Sex+Age+age2+EUR_ANC+PGS))-> results.UKB_afr[[I]]
        cat(I)
        cat(' done\n')
}

for (I in names(PGS3_UKB_afr)){
        tp <- boot(data=PGS2_UKB_afr[[I]], statistic=rsq.R2, R=999, formula1=Height~Sex+Age+age2+EUR_ANC, formula2=Height~Sex+Age+age2+EUR_ANC+PGS)
        list.append(results.UKB_afr[[I]], tp)-> results.UKB_afr[[I]]
        names(results.UKB_afr[[I]])<-c(a1[[I]], "total")
        cat(' done\n')
}
saveRDS(PGS3_UKB_afr, file='~/height_prediction/sib_betas/ukb_afr/output/PGS3_UKB_afr.Rds')
saveRDS(results.UKB_afr, file='~/height_prediction/sib_betas/ukb_afr/output/results.UKB_afr.Rds')

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


for (I in names(PGS3_UKB_afr)){
       myp<-ggplot(B_UKB_afr[[I]][1:4,], aes(x=Med_Eur_Anc, y=R_sq)) +
        geom_point(size=1.5, shape=21, fill="white") +
        geom_errorbar(aes(x=Med_Eur_Anc, ymin=boots_perc_L, ymax=boots_perc_U), width=0.05, size=0.8, color="cornflowerblue") +
        geom_line(color='lightgray')+
        geom_errorbarh(aes(x=Med_Eur_Anc, xmin=HVB_L, xmax=HVB_U), width=0.05, size=0.5, color="cornflowerblue") +
        labs(title = "UKB_AFR") + ylab("R-squared") + xlab("European Ancestry Proportion")
        print(myp)
        ggsave(paste0('~/height_prediction/sib_betas/ukb_afr/figs/UKB_afr_error_bars_', I, '.png'))
}

saveRDS(B_UKB_afr, file="~/height_prediction/sib_betas/ukb_afr/output/B_UKB_afr.Rds")
txtStop()
