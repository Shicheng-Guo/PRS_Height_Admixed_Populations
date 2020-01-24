##preable
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
readRDS('~/height_prediction/sib_betas/HRS_eur/output/PGS_HRS_eur.Rds')-> PGS_HRS_eur
#read in phenotype data
fread('~/height_prediction/input/HRS_eur/HRS_EUR_phenotypes.txt')-> Pheno_HRS_eur
###########
source('~/height_prediction/strat_prs/scripts/Rsq_R2.R')
##############

#add PGS to Pheno table in order to be able to make multiple analyses

Pheno_HRS_eur[,ID:=paste0(ID, "_", ID)]
setkey(Pheno_HRS_eur, ID)

PGS2_HRS_eur<-vector('list', length(PGS_HRS_eur))


names(PGS2_HRS_eur)<-names(PGS_HRS_eur)

for (I in names(PGS_HRS_eur)){
        data.table(ID=names(PGS_HRS_eur[[I]]), PGS=unlist(PGS_HRS_eur[[I]]))-> PGS2_HRS_eur[[I]]
        setkey(PGS2_HRS_eur[[I]], ID)
        PGS2_HRS_eur[[I]][Pheno_HRS_eur, nomatch=0]-> PGS2_HRS_eur[[I]]
	PGS2_HRS_eur[[I]][, EUR_ANC:=1]
        PGS2_HRS_eur[[I]][,AGE2:=AGE^2]
	PGS2_HRS_eur[[I]][,HEIGHT:=HEIGHT*100]
	#PGS2_UKB_eur[[I]][EUR_ANC>=0.05]->PGS2_UKB_eur[[I]]
	PGS2_HRS_eur[[I]]$SEX<-as.factor(PGS2_HRS_eur[[I]]$SEX)
        #PGS2_HRS_eur[[I]][-which(is.na(PGS2_HRS_eur[[I]][,HEIGHT])),]-> PGS2_HRS_eur[[I]]
}

lapply(PGS2_HRS_eur, function(X) lm(HEIGHT~SEX, X))-> lm1_HRS_eur
lapply(PGS2_HRS_eur, function(X) lm(HEIGHT~PGS, X))-> lm2_HRS_eur
lapply(PGS2_HRS_eur, function(X) lm(HEIGHT~AGE, X))-> lm3_HRS_eur
lapply(PGS2_HRS_eur, function(X) lm(HEIGHT~AGE2, X))-> lm4_HRS_eur
lapply(PGS2_HRS_eur, function(X) lm(HEIGHT~PGS+AGE, X))-> lm5_HRS_eur
lapply(PGS2_HRS_eur, function(X) lm(HEIGHT~PGS+AGE2, X))-> lm6_HRS_eur
lapply(PGS2_HRS_eur, function(X) lm(HEIGHT~SEX+AGE+AGE2, X))-> lm7_HRS_eur
lapply(PGS2_HRS_eur, function(X) lm(HEIGHT~SEX+AGE+AGE2+PGS, X))-> lm8_HRS_eur

partial.R2(lm7_HRS_eur[[63]],lm8_HRS_eur[[63]]) #12%

partial_r2_HRS_eur<-lapply(1:length(PGS2_HRS_eur), function(X) partial.R2(lm7_HRS_eur[[X]], lm8_HRS_eur[[X]])) #min 7.9% max 14.8%
names(partial_r2_HRS_eur)<-names(PGS2_HRS_eur)
which.max(unlist(partial_r2_HRS_eur)) #71, genet_0.3_0.0005
Nr_SNPs<-rep(NA, length(PGS2_HRS_eur))
names(Nr_SNPs)<- names(PGS2_HRS_eur)

for(I in names(Nr_SNPs)){
        readRDS(paste0('~/height_prediction/sib_betas/HRS_eur/output/hei_', I, '.Rds'))-> smtg
        sum(sapply(smtg, function(X) nrow(X)))-> Nr_SNPs[I]
	remove(smtg)
	gc()
        cat(I, ' \n')
}

data.table(Nr=unlist(Nr_SNPs), Name=names(Nr_SNPs), Part_R2=unlist(partial_r2_HRS_eur))-> A_table
saveRDS(A_table, file='~/height_prediction/sib_betas/HRS_eur/output/Nr_SNPs_HRS.Rds')


cor.test(unlist(Nr_SNPs), unlist(partial_r2_HRS_eur)) #43%
summary(lm(Part_R2~Nr,data=A_table))$r.squared #18%
#

#this weird negative corrrelation os due to the LD pruned sets: this is ebcause the LD blocks have few SNPs and explain A LOT and the LD_kb sets have a lot of SNPs and explain very little. 

#cor.test(A_table[-grep("LD",Name)]$Part_R2, A_table[-grep("LD",Name)]$Nr)


for(I in 1:length(PGS2_HRS_eur)){
        A<-ggpairs(PGS2_HRS_eur[[I]][,.(HEIGHT, SEX, PGS, AGE, AGE2)])
        png(filename=paste0('~/height_prediction/sib_betas/HRS_eur/figs/UKB_eur_ggpairs_', names(PGS2_HRS_eur)[I], '.png'))
        print(A)
        cat(names(PGS2_HRS_eur)[I])
        cat(' done\n')
        dev.off()
}

#lapply(PGS2_HRS_eur, function(X) X[, Quantile:= cut(EUR_ANC,
#                                breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.25), na.rm=TRUE),
#                                include.lowest=TRUE)])-> PGS3_HRS_eur

#names(PGS3_HRS_eur)<-names(PGS2_HRS_eur)

#lapply(1:length(PGS3_HRS_eur), function(X) PGS3_HRS_eur[[X]][,Med_Eur_Anc:=median(EUR_ANC),by=Quantile])

#lapply(1:length(PGS3_HRS_eur), function(X) as.character(unique((PGS3_HRS_eur[[X]]$Quantile))))-> a

#lapply(a, function(X) c(X[4],X[1], X[2], X[3]))-> a1 #check

#names(a1)<-names(PGS3_HRS_eur)

lapply(PGS2_HRS_eur, function(X) X[, Quantile:='total'][,Med_Eur_Anc:=1])-> PGS3_HRS_eur
#r2_HRS_eur<-vector('list', length(PGS3_HRS_eur))
#names(r2_HRS_eur)<-names(PGS3_HRS_eur)

#for(I in names(r2_HRS_eur)){
#        r2_HRS_eur[[I]]<-vector('list', length(a1[[I]]))
#        names(r2_HRS_eur[[I]])<-a1[[I]]
#        for(i in a1[[I]]){
#        r2_HRS_eur[[I]][[i]]<-partial.R2(lm(Height~Sex+Age+age2, PGS3_HRS_eur[[I]][Quantile==i]),lm(Height~Sex+Age+age2+PGS, PGS3_HRS_eur[[I]][Quantile==i]))
#        }
#}


B_HRS_eur<-vector('list', length(PGS3_HRS_eur))
names(B_HRS_eur)<-names(PGS3_HRS_eur)

for (I in names(PGS3_HRS_eur)){
        B_HRS_eur[[I]]<-data.table(Quant="total",R_sq=partial_r2_HRS_eur[[I]],Med_Eur_Anc=1)
        B_HRS_eur[[I]][,N:=nrow(PGS3_HRS_eur[[I]])]
        B_HRS_eur[[I]][, K:=1] #number of predictors. Need to check later if this is correct.
        B_HRS_eur[[I]][, LCL:=CI.Rsq(R_sq, k=K, n=N)[3]]
        B_HRS_eur[[I]][, UCL:=CI.Rsq(R_sq, k=K, n=N)[4]]
}

### add confidence intervals calculated with bootstrap: https://www.statmethods.net/advstats/bootstrapping.html
results.HRS_eur<-vector('list', length(PGS3_HRS_eur))
names(results.HRS_eur)<-names(PGS3_HRS_eur)

#for (I in names(PGS3_HRS_eur)){
#        results.HRS_eur[[I]]<-vector('list', length(a1[[I]])+1)
#        names(results.HRS_eur[[I]])<-c(a1[[I]], "total")
#       lapply(a1[[I]], function(i) boot(data=PGS3_HRS_eur[[I]][Quantile==i], statistic=rsq.R2,R=999, formula1=Height~Sex+Age+age2, formula2=Height~Sex+Age+age2+PGS))-> results.HRS_eur[[I]]
#        cat(I)
#        cat(' done\n')
#}

for (I in names(PGS3_HRS_eur)){
        results.HRS_eur[[I]]<- boot(data=PGS2_HRS_eur[[I]], statistic=rsq.R2, R=999, formula1=HEIGHT~SEX+AGE+AGE2, formula2=HEIGHT~SEX+AGE+AGE2+PGS)
        results.HRS_eur[[I]]
        #names(results.HRS_eur[[I]])<-c(a1[[I]], "total")
        cat(I, ' done\n')
}

#95% confidence intervals.
saveRDS(results.HRS_eur, file='~/height_prediction/sib_betas/HRS_eur/output/results.HRS_eur.Rds')
saveRDS(PGS3_HRS_eur, file='~/height_prediction/sib_betas/HRS_eur/output/PGS3_HRS_eur.Rds')


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

saveRDS(B_HRS_eur, file="~/height_prediction/sib_betas/HRS_eur/output/B_HRS_eur.Rds")


A_table[grep("phys", Name),]->dt_phys
A_table[grep("genet", Name),]->dt_genet
A_table[grep("LD",    Name),]->dt_LD

#factor(dt$Set)-> dt$Setp
factor(dt_phys$Name)-> dt_phys$Name
factor(dt_genet$Name)-> dt_genet$Name
factor(dt_LD$Name)-> dt_LD$Name
factor(dt_LD$Name, levels(dt_LD$Name)[c(4,5,3,1,2)])-> dt_LD$Name
factor(dt_genet$Name, levels(dt_genet$Name)[c(5,4,3,2,1,10,9,8,7,6,15,14,13,12,11,20,19,18,17,16,25,24,23,22,21, 30, 29, 28, 27, 26, 35,34,33,32,31)])-> dt_genet$Name
factor(dt_phys$Name, levels(dt_phys$Name)[c(24,27,30,33,35,4,7,10,13,15,16,17,18,19,20,22,25,28,31,34,36,37,38,39,40,2,5,8,11,14,21,23,26,29,32,1,3,6,9,12)])-> dt_phys$Name

melt(dt_LD)-> dt_LD
rbind(dt_LD[grep("_250000_", dt_LD$Set)][, window:=250000], dt_LD[grep("_100000_", dt_LD$Set)][, window:=100000], dt_LD[grep("_50000_", dt_LD$Set)][, window:=50000], dt_LD[grep("block", dt_LD$Set)][, window:="-"])-> dt_LD
dt_LD[, method:=gsub("LD_50000_0.01_0.5", "LD_0.01_0.5", gsub("LD_100000_0.01_0.5", "LD_0.01_0.5", gsub("LD_250000_0.01_0.5", "LD_0.01_0.5", gsub("LD_block_0_0_AFR", "LD_block_AFR", gsub("LD_block_0_0_EUR", "LD_block_EUR", dt_LD[,Set])))))]

as.factor(dt_LD$window)-> dt_LD$window
factor(dt_LD$window, levels(dt_LD$window)[c(1,4,2,3)])-> dt_LD$window

ggplot(dt_LD,aes(x=method, y=value, colour=window, shape=variable)) + geom_point(size=2.5, alpha=1)
ggsave('~/height_prediction/sib_betas/figs/reg_rsq_eur_anc_LD.png')
#

A_table[, Method:=NA]
A_table[, Window:=NA]
A_table$Method[grep("phys", A_table$Name)]<-'phys'
A_table$Method[grep("genet", A_table$Name)]<-'genet'
A_table$Method[grep("LD_block", A_table$Name)]<-'LD_block'
A_table$Window[grep("phys_5000_", A_table$Name)]<-"5000"
A_table$Window[grep("phys_10000_", A_table$Name)]<-"10000"
A_table$Window[grep("phys_25000_", A_table$Name)]<-"25000"
A_table$Window[grep("phys_50000_", A_table$Name)]<-"50000"
A_table$Window[grep("phys_75000_", A_table$Name)]<-"75000"
A_table$Window[grep("phys_100000_", A_table$Name)]<-"100000"
A_table$Window[grep("phys_500000_", A_table$Name)]<-"500000"
A_table$Window[grep("phys_1000000_", A_table$Name)]<-"1000000"
A_table$Window[grep("genet_1_", A_table$Name)]<-"1cM"
A_table$Window[grep("genet_0.5_", A_table$Name)]<-"0.5cM"
A_table$Window[grep("genet_0.3_", A_table$Name)]<-"0.3cM"
A_table$Window[grep("genet_0.25_", A_table$Name)]<-"0.25cM"
A_table$Window[grep("genet_0.2_", A_table$Name)]<-"0.2cM"
A_table$Window[grep("genet_0.15_", A_table$Name)]<-"0.15cM"
A_table$Window[grep("genet_0.1_", A_table$Name)]<-"0.1cM"
A_table$Window[grep("LD_250000_", A_table$Name)]<-"LD"
A_table$Window[grep("LD_100000_", A_table$Name)]<-"LD"
A_table$Window[grep("LD_50000_", A_table$Name)]<-"LD"

A_table$Method[which(is.na(A_table$Method))]<-"LD"
ggplot(A_table, aes(x=Nr, y=Part_R2, colour=Method, group=Window)) + geom_point(size=2) + geom_line()+ geom_text(A_table[Name=='phys_100000_0.0005'], aes(label=Name)) +
theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),  axis.text.x=element_text(size=9), axis.text.y=element_text(size=9))
ggsave('~/height_prediction/sib_betas/HRS_eur/figs/ashg_like.png')

