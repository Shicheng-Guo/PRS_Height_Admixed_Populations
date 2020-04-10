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
library(TeachingDemos)

txtStart(paste0("~/height_prediction/unweighted_prs/JHS/plots_out.txt"))

#read in PGS scores
readRDS('~/height_prediction/unweighted_prs/output/PGS_JHS.Rds')-> PGS_JHS

#read in phenotype data
fread('~/height_prediction/input/JHS/JHS_phenotypes.txt')-> Pheno_JHS

#a partial R2 function
source('~/height_prediction/strat_prs/scripts/Rsq_R2.R')
##############

#add PGS to Pheno table in order to be able to make multiple analyses

#Pheno_JHS[, SUBJID:=paste0("0_", as.character(Pheno_JHS[, SUBJID]))]
setkey(Pheno_JHS, SUBJID)

a<-unlist(lapply(names(PGS_JHS[[1]]), function(X) substr(X, 1,7)))
for(I in 1:length(PGS_JHS)){

names(PGS_JHS[[I]])<-a
}


#add ancestry
ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/JHS/rfmix_anc_chr', X, '.txt'))))

anc_JHS<-ancestry %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual
anc_JHS$SUBJID<-substr(anc_JHS[,SUBJID],3,9)
setkey(anc_JHS, SUBJID)


PGS2_JHS<-vector('list', length(PGS_JHS))
names(PGS2_JHS)<-names(PGS_JHS)

for (I in names(PGS2_JHS)){
        data.table(SUBJID=names(PGS_JHS[[I]]), PGS=unlist(PGS_JHS[[I]]))-> PGS2_JHS[[I]]
        setkey(PGS2_JHS[[I]], SUBJID)
        PGS2_JHS[[I]][Pheno_JHS, nomatch=0]-> PGS2_JHS[[I]]
        PGS2_JHS[[I]][anc_JHS, nomatch=0]-> PGS2_JHS[[I]]
        PGS2_JHS[[I]][,age2:=age_baseline^2]
        PGS2_JHS[[I]][,age:=age_baseline][, age_baseline:=NULL]
        PGS2_JHS[[I]][AFR_ANC>=0.05]-> PGS2_JHS[[I]]#filter out individuals that are not african...
        PGS2_JHS[[I]][-which(is.na(PGS2_JHS[[I]][,height_baseline])),]-> PGS2_JHS[[I]]
        PGS2_JHS[[I]][,HEIGHTX:=height_baseline]
        PGS2_JHS[[I]][,height_baseline:=NULL]
        PGS2_JHS[[I]][, sex:=sex_selfreport]
        PGS2_JHS[[I]][,sex_selfreport:=NULL]
}


lapply(PGS2_JHS, function(X) lm(HEIGHTX~sex,X))-> lm0_JHS
lapply(PGS2_JHS, function(X) lm(HEIGHTX~PGS, X))-> lm1_JHS
lapply(PGS2_JHS, function(X) lm(HEIGHTX~age, X))-> lm2_JHS
lapply(PGS2_JHS, function(X) lm(HEIGHTX~age2, X))-> lm3_JHS
lapply(PGS2_JHS, function(X) lm(HEIGHTX~EUR_ANC, X))-> lm4_JHS
lapply(PGS2_JHS, function(X) lm(HEIGHTX~sex+age, X))-> lm5_JHS
lapply(PGS2_JHS, function(X) lm(HEIGHTX~sex+age+age2, X))-> lm6_JHS
lapply(PGS2_JHS, function(X) lm(HEIGHTX~sex+age+age2+EUR_ANC, X))-> lm7_JHS
lapply(PGS2_JHS, function(X) lm(HEIGHTX~sex+age+age2+EUR_ANC+PGS,X))-> lm8_JHS

partial.R2(lm7_JHS[[35]],lm8_JHS[[35]])  

partial_r2_JHS<-lapply(1:length(PGS2_JHS), function(X) partial.R2(lm7_JHS[[X]], lm8_JHS[[X]])) 
names(partial_r2_JHS)<- names(PGS2_JHS)


Nr_SNPs<-rep(NA, length(PGS2_JHS))
names(Nr_SNPs)<- names(PGS2_JHS)

for(I in names(Nr_SNPs)){
        readRDS(paste0('~/height_prediction/gwas/JHS/output/hei_', I, '.Rds'))-> smtg
        sum(sapply(smtg, function(X) nrow(X)))-> Nr_SNPs[I]
        cat(I, ' \n')

}
data.table(Nr=unlist(Nr_SNPs), Name=names(Nr_SNPs), Part_R2=unlist(partial_r2_JHS))-> A_table
saveRDS(A_table, file='~/height_prediction/unweighted_prs/output/Nr_SNPs_JHS.Rds')

cor.test(unlist(Nr_SNPs), unlist(partial_r2_JHS))# 
summary(lm(Part_R2~Nr,data=A_table))$r.squared #
#

#interestingly, the LD pruned sets scrw things up. If we remvoe them


for(I in 1:length(PGS2_JHS)){
        A<-ggpairs(PGS2_JHS[[I]][,.(HEIGHTX, PGS, sex, age,EUR_ANC,weight_baseline)])
        png(filename=paste0('~/height_prediction/unweighted_prs/figs/JHS_ggpairs_', names(PGS2_JHS)[I], '.png'))
        print(A)
        cat(names(PGS2_JHS)[I])
        cat(' done\n')
        dev.off()
}



lapply(PGS2_JHS, function(X) X[, Quantile:= cut(EUR_ANC,
                                breaks=quantile(EUR_ANC, probs=seq(0,1, by=1/2), na.rm=TRUE),
                                include.lowest=TRUE)])-> PGS3_JHS

lapply(1:length(PGS3_JHS), function(X) PGS3_JHS[[X]][,Med_Eur_Anc:=median(EUR_ANC),by=Quantile])

lapply(1:length(PGS3_JHS), function(X) as.character(unique((PGS3_JHS[[X]]$Quantile))))-> a

lapply(a, function(X) c(X[2], X[1]))-> a1

names(a1)<-names(PGS3_JHS)

r2_JHS<-vector('list', length(PGS3_JHS))
names(r2_JHS)<-names(PGS3_JHS)

for(I in names(r2_JHS)){

        r2_JHS[[I]]<-vector('list', length(a1[[I]]))
        names(r2_JHS[[I]])<-a1[[I]]
        for(i in a1[[I]]){
        r2_JHS[[I]][[i]]<-partial.R2(lm(HEIGHTX~sex+age+age2+EUR_ANC, PGS3_JHS[[I]][Quantile==i]),lm(HEIGHTX~sex+age+age2+EUR_ANC+PGS, PGS3_JHS[[I]][Quantile==i]))
        }
}


B_JHS<-vector('list', length(r2_JHS))
names(B_JHS)<-names(r2_JHS)

for (I in names(r2_JHS)){
        B_JHS[[I]]<-data.table(Quant=c(a1[[I]], "total"),
        R_sq=c(unlist(r2_JHS[[I]]), partial_r2_JHS[[I]]),
        Med_Eur_Anc=c(unique(PGS3_JHS[[I]][Quantile==a1[[I]][1]][,Med_Eur_Anc]), unique(PGS3_JHS[[I]][Quantile==a1[[I]][2]][,Med_Eur_Anc]), median(PGS3_JHS[[I]][, EUR_ANC])))
        B_JHS[[I]][,N:=c(nrow(PGS3_JHS[[I]][Quantile==a1[[I]][1]]), nrow(PGS3_JHS[[I]][Quantile==a1[[I]][2]]),nrow(PGS3_JHS[[I]]))]
        B_JHS[[I]][,K:=1] #number of predictors. Need to check later if this is correct.
        B_JHS[[I]][, LCL:=CI.Rsq(R_sq, k=K, n=N)[3]]
        B_JHS[[I]][, UCL:=CI.Rsq(R_sq, k=K, n=N)[4]]
}

### add confidence intervals calculated with bootstrap: https://www.statmethods.net/advstats/bootstrapping.html
results.JHS<-vector('list', length(PGS3_JHS))
names(results.JHS)<-names(PGS3_JHS)
for (I in names(PGS3_JHS)){
        results.JHS[[I]]<-vector('list', length(a1[[I]])+1)
        names(results.JHS[[I]])<-c(a1[[I]], "total")
        lapply(a1[[I]], function(i) boot(data=PGS3_JHS[[I]][Quantile==i], statistic=rsq.R2,R=999, formula1=HEIGHTX~sex+age+age2+EUR_ANC, formula2=HEIGHTX~sex+age+age2+EUR_ANC+PGS))-> results.JHS[[I]]
        cat(I)
        cat(' done\n')
}
#saveRDS(PGS3_JHS, file='/project/mathilab/bbita/gwas_admix/new_height/JHS/PGS3_JHS.Rds')

#saveRDS(results.JHS, file='/project/mathilab/bbita/gwas_admix/new_height/JHS/results.JHS.Rds')

for (I in names(PGS3_JHS)){
        tp <- boot(data=PGS2_JHS[[I]], statistic=rsq.R2, R=999, formula1=HEIGHTX~sex+age+age2+EUR_ANC, formula2=HEIGHTX~sex+age+age2+EUR_ANC+PGS)
        list.append(results.JHS[[I]], tp)-> results.JHS[[I]]
        names(results.JHS[[I]])<-c(a1[[I]], "total")
        cat(I)
        cat(' done\n')
}
saveRDS(PGS3_JHS, file='~/height_prediction/unweighted_prs/output/PGS3_JHS.Rds')
saveRDS(results.JHS, file='~/height_prediction/unweighted_prs/output/results.JHS.Rds')

#confidence intervals


boots.ci.JHS<-lapply(results.JHS, function(Y) lapply(Y, function(X) boot.ci(X, type = c("norm", 'basic', "perc"))))
names(boots.ci.JHS)<-names(results.JHS)

for (I in names(PGS3_JHS)){
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

for (I in names(PGS3_JHS)){
        myp<-ggplot(B_JHS[[I]][1:2,], aes(x=Med_Eur_Anc, y=R_sq)) +
        geom_point(size=1.5, shape=21, fill="white") +
        geom_errorbar(aes(x=Med_Eur_Anc, ymin=boots_perc_L, ymax=boots_perc_U), width=0.05, size=0.8, color="cornflowerblue") +
        geom_line(color='lightgray')+
        geom_errorbarh(aes(x=Med_Eur_Anc, xmin=HVB_L, xmax=HVB_U), width=0.05, size=0.5, color="cornflowerblue") +
        labs(title = "JHS_AA") + ylab("R-squared") + xlab("European Ancestry Proportion")
        print(myp)
        ggsave(paste0('~/height_prediction/unweighted_prs/figs/JHS_error_bars_', I, '.png'))
}

saveRDS(B_JHS, file="~/height_prediction/unweighted_prs/output/B_JHS.Rds")
txtStop()
