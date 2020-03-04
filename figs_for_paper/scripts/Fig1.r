library("optparse")
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(wesanderson)
library(rlist)
library(asbio)
library(GGally)
library(tidyr)
library(hexbin)
library(cowplot)
library(psychometric)
library(boot)
library(RColorBrewer)
options(scipen=999)

i<-'phys_100000_0.0005'
args<-'gwas'
readRDS(paste0('~/height_prediction/', args[1], '/WHI/output/results.WHI.Rds'))[[i]]-> results.WHI
readRDS(paste0('~/height_prediction/',args[1], '/JHS/output/results.JHS.Rds'))[[i]]-> results.JHS
readRDS(paste0('~/height_prediction/',args[1], '/ukb_afr/output/results.UKB_afr.Rds'))[[i]]-> results.UKB_afr
readRDS(paste0('~/height_prediction/',args[1], '/HRS_eur/output/results.HRS_eur.Rds'))[[i]]-> results.HRS_eur
readRDS(paste0('~/height_prediction/',args[1], '/HRS_afr/output/results.HRS_afr.Rds'))[[i]]-> results.HRS_afr
readRDS(paste0('~/height_prediction/',args[1], '/WHI/output/B_WHI.Rds'))[[i]]->B_WHI
readRDS(paste0('~/height_prediction/',args[1], '/JHS/output/B_JHS.Rds'))[[i]]->B_JHS
readRDS(paste0('~/height_prediction/',args[1], '/ukb_afr/output/B_UKB_afr.Rds'))[[i]]-> B_UKB_afr
readRDS(paste0('~/height_prediction/',args[1], '/HRS_eur/output/B_HRS_eur.Rds'))[[i]]-> B_HRS_eur
readRDS(paste0('~/height_prediction/',args[1],'/HRS_afr/output/B_HRS_afr.Rds'))[[i]]-> B_HRS_afr
readRDS(paste0('~/height_prediction/',args[1],'/ukb_afr/output/PGS3_UKB_afr.Rds'))[[i]]-> PGS3_UKB_afr
readRDS(paste0('~/height_prediction/',args[1],'/WHI/output/PGS3_WHI.Rds'))[[i]]-> PGS3_WHI
readRDS(paste0('~/height_prediction/',args[1],'/JHS/output/PGS3_JHS.Rds'))[[i]]-> PGS3_JHS
readRDS(paste0('~/height_prediction/',args[1],'/HRS_eur/output/PGS3_HRS_eur.Rds'))[[i]]-> PGS3_HRS_eur
readRDS(paste0('~/height_prediction/',args[1],'/HRS_afr/output/PGS3_HRS_afr.Rds'))[[i]]-> PGS3_HRS_afr

ALL<-rbind(B_JHS[1:2,][, Dataset:='JHS_afr'], B_WHI[1:4,][, Dataset:='WHI_afr'], B_UKB_afr[1:4,][,Dataset:='UKB_afr'],B_HRS_afr[1:2,][, Dataset:='HRS_afr'], _HRS_eur)

tmp<-1/c(var(results.JHS[[1]]$t), var(results.JHS[[2]]$t), var(results.WHI[[1]]$t), var(results.WHI[[2]]$t), var(results.WHI[[3]]$t), var(results.WHI[[4]]$t), var(results.UKB_afr[[1]]$t),var(results.UKB_afr[[2]]$t), var(results.UKB_afr[[3]]$t), var(results.UKB_afr[[4]]$t), var(results.HRS_afr[[1]]$t), var(results.HRS_afr[[2]]$t), var(results.HRS_eur$t))  #weighing lm by boostrap replicates.

ALL[,W:=tmp]
ALL$Dataset<-factor(ALL$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "HRS_eur"))
my_colrs<-c(brewer.pal(4, 'Set1'),"#101010")


A_plot<-ggplot(ALL2, aes(x=Med_Eur_Anc, y=R_sq,colour=Dataset)) +
	geom_point(aes(shape=Dataset), size=1.5, fill="white", alpha=0.8) + stat_smooth(data=ALL2,method = "lm", mapping = aes(weight = W), col='black') +
        geom_errorbar(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset,ymin=boots_perc_L, ymax=boots_perc_U), width=0.05, size=0.8) +
        #geom_line(color='lightgray')+
        geom_errorbarh(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset, xmin=HVB_L, xmax=HVB_U), width=0.05, size=0.5) +
        scale_color_manual(values=my_colrs) +
        ylab(expression(paste("Partial R"^"2"))) + xlab("European Ancestry Proportion") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.key=element_blank(), legend.background=element
_blank(), legend.title=element_blank(), legend.text=element_text(size=15),legend.position = c(0.15,0.85))
        print(A_plot)


#####################3
#Second plot
##Load other PRSs
all_datasets<-vector('list', 5)
names(all_datasets)<-c('WHI_afr', 'JHS_afr', 'UKB_afr', 'HRS_afr', 'HRS_eur')
prs<-vector('list', 22)
for(chr in 1:22){
	prs[[chr]]<-readRDS(paste0('~/height_prediction/loc_anc_analysis/output/chr', chr,'_0_prs_WHI_EuR.Rds'))
}
samp_names<-names(prs[[1]])
dt<-data.table(SUBJID=samp_names, PRS_0=NA)

prs2<-vector('list', length(samp_names))
names(prs2)<-samp_names


for (S in samp_names){
	prs2[[S]]<-sum(prs[[1]][[S]],prs[[2]][[S]],prs[[3]][[S]],prs[[4]][[S]],prs[[5]][[S]],prs[[6]][[S]], prs[[7]][[S]], prs[[8]][[S]], prs[[9]][[S]], prs[[10]][[S]],prs[[11]][[S]], prs[[12]][[S]], prs[[13]][[S]],prs[[14]][[S]], prs[[15]][[S]], prs[[16]][[S]], prs[[17]][[S]], prs[[18]][[S]], prs[[19]][[S]], prs[[20]][[S]], prs[[21]][[S]],prs[[22]][[S]], na.rm=T)
}

dt[, PRS_0:=scale(unlist(prs2))]

#phenotype
fread('~/height_prediction/input/WHI/WHI_phenotypes.txt')-> Pheno_WHI

Pheno_WHI[, SUBJID:=paste0("0_", as.character(Pheno_WHI[, SUBJID]))]
setkey(Pheno_WHI, SUBJID)
#add ancestry
ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/WHI/rfmix_anc_chr', X, '.txt'))))
anc_WHI<-ancestry %>% dplyr::group_by(SUBJID) %>% dplyr::summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual
setkey(anc_WHI, SUBJID)
##

#PRS_EUR<-data.table(PRS_EUR=unlist(readRDS('~/height_prediction/gwas/WHI/output/PGS_WHI_phys_100000_0.0005.Rds')),SUBJID=names(unlist(readRDS('~/height_prediction/gwas/WHI/output/PGS_WHI_phys_100000_0.0005.Rds'))))
#a<-merge(dt, PRS_EUR, by="SUBJID")
setkey(dt, SUBJID)
setkey(Pheno_WHI, SUBJID)
dt[Pheno_WHI][anc_WHI]-> final
final2<-final[,c("SUBJID","PRS_0","AGE", "EUR_ANC", "HEIGHTX")]

final2[,AGE2:=AGE^2]

#
r2_whi_afr<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0, data=final2)) # 1%, compared to 3.87% for when using all SNPs and alfa=0
#
saveRDS(final2, '~/height_prediction/loc_anc_analysis/output/all_PRS_WHI_Eur.Rds')
all_datasets[['WHI_afr']]<-final2
remove(final2,dt)

all_datasets[['WHI_afr']][, Quantile:= cut(EUR_ANC,
                        breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.25), na.rm=TRUE),
                        include.lowest=TRUE)]

all_datasets[['WHI_afr']][,Med_Eur_Anc:=median(EUR_ANC),by=Quantile]

as.character(unique((all_datasets[['WHI_afr']]$Quantile)))-> a

c(a[4],a[2], a[1], a[3])-> a1

r2_WHI<-vector('list', length(a1))
names(r2_WHI)<-a1
for(i in a1){
        r2_WHI[[i]]<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC,all_datasets[['WHI_afr']][Quantile==i]),lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0, all_datasets[['WHI_afr']][Quantile==i]))
}
#

B_WHI_afr<-data.table(Quant=c(a1, "total"),R_sq=c(unlist(r2_WHI),r2_whi_afr),
Med_Eur_Anc=c(unique(all_datasets[['WHI_afr']][Quantile==a1[1]][,Med_Eur_Anc]), unique(all_datasets[['WHI_afr']][Quantile==a1[2]][,Med_Eur_Anc]), unique(all_datasets[['WHI_afr']][Quantile==a1[3]][,Med_Eur_Anc]), unique(all_datasets[['WHI_afr']][Quantile==a1[4]][,Med_Eur_Anc]), median(all_datasets[['WHI_afr']][, EUR_ANC])))
B_WHI_afr[,N:=c(nrow(all_datasets[['WHI_afr']][Quantile==a1[1]]), nrow(all_datasets[['WHI_afr']][Quantile==a1[2]]), nrow(all_datasets[['WHI_afr']][Quantile==a1[3]]), nrow(all_datasets[['WHI_afr']][Quantile==a1[4]]),nrow(all_datasets[['WHI_afr']]))]
B_WHI_afr[,K:=1] #number of predictors. Need to check later if this is correct.
B_WHI_afr[, LCL:=CI.Rsq(R_sq, k=K, n=N)[3]]
B_WHI_afr[, UCL:=CI.Rsq(R_sq, k=K, n=N)[4]]

### add confidence intervals calculated with bootstrap: https://www.statmethods.net/advstats/bootstrapping.html
results.WHI_afr<-vector('list', length(a1)+1)
names(results.WHI_afr)<-c(a1, "total")
for(i in a1){
        boot(data=all_datasets[['WHI_afr']][Quantile==i], statistic=rsq.R2,R=999, formula1=HEIGHTX~AGE+AGE2+EUR_ANC, formula2=HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0)-> results.WHI_afr[[i]]
}

results.WHI_afr[['total']] <- boot(data=all_datasets[['WHI_afr']], statistic=rsq.R2, R=999, formula1=HEIGHTX~AGE+AGE2+EUR_ANC, formula2=HEIGHTX~AGE+AGE2+EUR_ANC+PRS_0)
#confidence intervals
boots.ci.WHI_afr<-lapply(results.WHI_afr, function(X) boot.ci(X, type = c("norm", 'basic', "perc")))
names(boots.ci.WHI_afr)<-names(results.WHI_afr)

B_WHI_afr[1:4,]-> a
B_WHI_afr[5,]-> b
a[,HVB_L:=sapply(a$Quant, function(X) as.numeric(gsub("\\]","",gsub("\\(","",gsub("\\[","",strsplit(X,",")[[1]])))))[1,]]
a[,HVB_U:=sapply(a$Quant, function(X) as.numeric(gsub("\\]","",gsub("\\(","",gsub("\\[","",strsplit(X,",")[[1]])))))[2,]]
b[,HVB_L:=1]
b[,HVB_U:=1]
rbind(a,b)->B_WHI_afr
B_WHI_afr[, Dataset:='WHI_afr']
B_WHI_afr[, boots_norm_L:=sapply(1:5, function(X) boots.ci.WHI_afr[[X]]$normal[2])]
B_WHI_afr[, boots_norm_U:=sapply(1:5, function(X) boots.ci.WHI_afr[[X]]$normal[3])]
B_WHI_afr[, boots_perc_L:=sapply(1:5, function(X) boots.ci.WHI_afr[[X]]$perc[4])]
B_WHI_afr[, boots_perc_U:=sapply(1:5, function(X) boots.ci.WHI_afr[[X]]$perc[5])]
B_WHI_afr[, boots_basic_L:=sapply(1:5, function(X) boots.ci.WHI_afr[[X]]$basic[4])]
B_WHI_afr[, boots_basic_U:=sapply(1:5, function(X) boots.ci.WHI_afr[[X]]$basic[5])]
#####################################################################################
#####################################################################################
##JHS

####
prs<-vector('list', 22)
for(chr in 1:22){
        prs[[chr]]<-readRDS(paste0('~/height_prediction/loc_anc_analysis/output/chr', chr,'_0_prs_JHS_EuR.Rds'))
}
samp_names<-names(prs[[1]])
dt<-data.table(SUBJID=samp_names, PRS_0=NA)

prs2<-vector('list', length(samp_names))
names(prs2)<-samp_names


for (S in samp_names){
        prs2[[S]]<-sum(prs[[1]][[S]],prs[[2]][[S]],prs[[3]][[S]],prs[[4]][[S]],prs[[5]][[S]],prs[[6]][[S]], prs[[7]][[S]], prs[[8]][[S]], prs[[9]][[S]], prs[[10]][[S]],prs[[11]][[S]], prs[[12]][[S]], prs[[13]][[S]],prs[[14]][[S]], prs[[15]][[S]], prs[[16]][[S]], prs[[17]][[S]], prs[[18]][[S]], prs[[19]][[S]], prs[[20]][[S]], prs[[21]][[S]],prs[[22]][[S]], na.rm=T)
}

dt[, PRS_0:=scale(unlist(prs2))]
dt$SUBJID<-substr(dt$SUBJID, 1,7)
#phenotype
fread('~/height_prediction/input/JHS/JHS_phenotypes.txt')-> Pheno_JHS
Pheno_JHS[,AGE2:=age_baseline^2]
setkey(Pheno_JHS, SUBJID)
#add ancestry
anc_JHS<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/JHS/rfmix_anc_chr', X, '.txt'))))
anc_JHS$SUBJID<-substr(anc_JHS[,SUBJID],3,9)
anc_JHS<-anc_JHS %>% dplyr::group_by(SUBJID) %>% dplyr::summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual
setkey(anc_JHS, SUBJID)
##
setkey(dt, SUBJID)
setkey(Pheno_JHS, SUBJID)
dt[Pheno_JHS][anc_JHS]-> final
final[, AGE:=age_baseline]
final[, HEIGHTX:=height_baseline]
final2<-final[,c("SUBJID","PRS_0","EUR_ANC", "SEX", "HEIGHTX", "AGE2", "AGE")]

r2_jhs_afr<-partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_0, data=final2)) # 0.75%
#
saveRDS(final2, '~/height_prediction/loc_anc_analysis/output/all_PRS_JHS_Eur.Rds')
all_datasets[['JHS_afr']]<-final2
remove(final2,dt)
all_datasets[['JHS_afr']][, Quantile:= cut(EUR_ANC,
                        breaks=quantile(EUR_ANC, probs=seq(0,1, by=1/2), na.rm=TRUE),
                        include.lowest=TRUE)]

all_datasets[['JHS_afr']][,Med_Eur_Anc:=median(EUR_ANC),by=Quantile]

as.character(unique((all_datasets[['JHS_afr']]$Quantile)))-> a

c(a[2], a[1])-> a1

r2_JHS_afr<-vector('list', length(a1))
names(r2_JHS_afr)<-a1
for(i in a1){
        r2_JHS_afr[[i]]<-partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC,all_datasets[['JHS_afr']][Quantile==i]),lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_0, all_datasets[['JHS_afr']][Quantile==i]))
}

B_JHS_afr<-data.table(Quant=c(a1, "total"),R_sq=c(unlist(r2_JHS_afr),r2_jhs_afr),
Med_Eur_Anc=c(unique(all_datasets[['JHS_afr']][Quantile==a1[1]][,Med_Eur_Anc]), unique(all_datasets[['JHS_afr']][Quantile==a1[2]][,Med_Eur_Anc]), median(all_datasets[['JHS_afr']][, EUR_ANC])))
B_JHS_afr[,N:=c(nrow(all_datasets[['JHS_afr']][Quantile==a1[1]]), nrow(all_datasets[['JHS_afr']][Quantile==a1[2]]),nrow(all_datasets[['JHS_afr']]))]
B_JHS_afr[,K:=1] #number of predictors. Need to check later if this is correct.
B_JHS_afr[, LCL:=CI.Rsq(R_sq, k=K, n=N)[3]]
B_JHS_afr[, UCL:=CI.Rsq(R_sq, k=K, n=N)[4]]

### add confidence intervals calculated with bootstrap: https://www.statmethods.net/advstats/bootstrapping.html
results.JHS_afr<-vector('list', length(a1)+1)
names(results.JHS_afr)<-c(a1, "total")
for(i in a1){
        boot(data=all_datasets[['JHS_afr']][Quantile==i], statistic=rsq.R2,R=999, formula1=HEIGHTX~SEX+AGE+AGE2+EUR_ANC, formula2=HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_0)-> results.JHS_afr[[i]]
}

results.JHS_afr[['total']] <- boot(data=all_datasets[['JHS_afr']], statistic=rsq.R2, R=999, formula1=HEIGHTX~SEX+AGE+AGE2+EUR_ANC, formula2=HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_0)
#confidence intervals
boots.ci.JHS_afr<-lapply(results.JHS_afr, function(X) boot.ci(X, type = c("norm", 'basic', "perc")))
names(boots.ci.JHS_afr)<-names(results.JHS_afr)


B_JHS_afr[1:2,]-> a
B_JHS_afr[3,]-> b
a[,HVB_L:=sapply(a$Quant, function(X) as.numeric(gsub("\\]","",gsub("\\(","",gsub("\\[","",strsplit(X,",")[[1]])))))[1,]]
a[,HVB_U:=sapply(a$Quant, function(X) as.numeric(gsub("\\]","",gsub("\\(","",gsub("\\[","",strsplit(X,",")[[1]])))))[2,]]
b[,HVB_L:=1]
b[,HVB_U:=1]
rbind(a,b)->B_JHS_afr
B_JHS_afr[, Dataset:='JHS_afr']
B_JHS_afr[, boots_norm_L:=sapply(1:3, function(X) boots.ci.JHS_afr[[X]]$normal[2])]
B_JHS_afr[, boots_norm_U:=sapply(1:3, function(X) boots.ci.JHS_afr[[X]]$normal[3])]
B_JHS_afr[, boots_perc_L:=sapply(1:3, function(X) boots.ci.JHS_afr[[X]]$perc[4])]
B_JHS_afr[, boots_perc_U:=sapply(1:3, function(X) boots.ci.JHS_afr[[X]]$perc[5])]
B_JHS_afr[, boots_basic_L:=sapply(1:3, function(X) boots.ci.JHS_afr[[X]]$basic[4])]
B_JHS_afr[, boots_basic_U:=sapply(1:3, function(X) boots.ci.JHS_afr[[X]]$basic[5])]

#####################################################################################
#####################################################################################
#UKB_afr
prs<-vector('list', 22)
for(chr in 1:22){
        prs[[chr]]<-readRDS(paste0('~/height_prediction/loc_anc_analysis/output/chr', chr,'_0_prs_UKB_afr_EuR.Rds'))
}
samp_names<-names(prs[[1]])
dt<-data.table(SUBJID=samp_names, PRS_0=NA)

prs2<-vector('list', length(samp_names))
names(prs2)<-samp_names


for (S in samp_names){
        prs2[[S]]<-sum(prs[[1]][[S]],prs[[2]][[S]],prs[[3]][[S]],prs[[4]][[S]],prs[[5]][[S]],prs[[6]][[S]], prs[[7]][[S]], prs[[8]][[S]], prs[[9]][[S]], prs[[10]][[S]],prs[[11]][[S]], prs[[12]][[S]], prs[[13]][[S]],prs[[14]][[S]], prs[[15]][[S]], prs[[16]][[S]], prs[[17]][[S]], prs[[18]][[S]], prs[[19]][[S]], prs[[20]][[S]], prs[[21]][[S]],prs[[22]][[S]], na.rm=T)
}

dt[, PRS_0:=scale(unlist(prs2))]

#phenotype
fread('~/height_prediction/input/ukb_afr/UKB_AFR_pheno.txt')-> Pheno_UKB_afr

Pheno_UKB_afr[, SUBJID:=paste0(as.character(Pheno_UKB_afr[, ID]),"_", as.character(Pheno_UKB_afr[,ID]))]
setkey(Pheno_UKB_afr, SUBJID)
#add ancestry
ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/ukb_afr/rfmix_anc_chr', X, '.txt'))))
anc_UKB_afr<-ancestry %>% dplyr::group_by(SUBJID) %>% dplyr::summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual
anc_UKB_afr[, SUBJID:=paste0(as.character(anc_UKB_afr[, SUBJID]),"_", as.character(anc_UKB_afr[,SUBJID]))]
setkey(anc_UKB_afr, SUBJID)
##

setkey(dt, SUBJID)
setkey(Pheno_UKB_afr, SUBJID)
dt[Pheno_UKB_afr][anc_UKB_afr]-> final
final[, HEIGHTX:=Height][,AGE:=Age]
final[,AGE2:=AGE^2][, SEX:=Sex]

final2<-final[,c("SUBJID","PRS_0","AGE", "EUR_ANC", "HEIGHTX", 'AGE2', "SEX")]
#c
r2_ukb_afr<-partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_0, data=final2)) #0.75%
saveRDS(final2, '~/height_prediction/loc_anc_analysis/output/all_PRS_UKB_afr_Eur.Rds')
all_datasets[['UKB_afr']]<-final2
remove(final2,dt)
all_datasets[['UKB_afr']][, Quantile:= cut(EUR_ANC,
                        breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.25), na.rm=TRUE),
                        include.lowest=TRUE)]

all_datasets[['UKB_afr']][,Med_Eur_Anc:=median(EUR_ANC),by=Quantile]

as.character(unique((all_datasets[['UKB_afr']]$Quantile)))-> a

c(a[2],a[4], a[3], a[1])-> a1

r2_UKB_afr<-vector('list', length(a1))
names(r2_UKB_afr)<-a1
for(i in a1){
        r2_UKB_afr[[i]]<-partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC,all_datasets[['UKB_afr']][Quantile==i]),lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_0, all_datasets[['UKB_afr']][Quantile==i]))
}

B_UKB_afr<-data.table(Quant=c(a1, "total"),R_sq=c(unlist(r2_UKB_afr),r2_ukb_afr),
Med_Eur_Anc=c(unique(all_datasets[['UKB_afr']][Quantile==a1[1]][,Med_Eur_Anc]), unique(all_datasets[['UKB_afr']][Quantile==a1[2]][,Med_Eur_Anc]), unique(all_datasets[['UKB_afr']][Quantile==a1[3]][,Med_Eur_Anc]), unique(all_datasets[['UKB_afr']][Quantile==a1[4]][,Med_Eur_Anc]), median(all_datasets[['UKB_afr']][, EUR_ANC])))
B_UKB_afr[,N:=c(nrow(all_datasets[['UKB_afr']][Quantile==a1[1]]), nrow(all_datasets[['UKB_afr']][Quantile==a1[2]]),nrow(all_datasets[['UKB_afr']][Quantile==a1[3]]), nrow(all_datasets[['UKB_afr']][Quantile==a1[4]]), nrow(all_datasets[['UKB_afr']]))]
B_UKB_afr[,K:=1] #number of predictors. Need to check later if this is correct.
B_UKB_afr[, LCL:=CI.Rsq(R_sq, k=K, n=N)[3]]
B_UKB_afr[, UCL:=CI.Rsq(R_sq, k=K, n=N)[4]]

### add confidence intervals calculated with bootstrap: https://www.statmethods.net/advstats/bootstrapping.html
results.UKB_afr<-vector('list', length(a1)+1)
names(results.UKB_afr)<-c(a1, "total")
for(i in a1){
        boot(data=all_datasets[['UKB_afr']][Quantile==i], statistic=rsq.R2,R=999, formula1=HEIGHTX~SEX+AGE+AGE2+EUR_ANC, formula2=HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_0)-> results.UKB_afr[[i]]
}

results.UKB_afr[['total']] <- boot(data=all_datasets[['UKB_afr']], statistic=rsq.R2, R=999, formula1=HEIGHTX~SEX+AGE+AGE2+EUR_ANC, formula2=HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_0)
#confidence intervals
boots.ci.UKB_afr<-lapply(results.UKB_afr, function(X) boot.ci(X, type = c("norm", 'basic', "perc")))
names(boots.ci.UKB_afr)<-names(results.UKB_afr)
B_UKB_afr[1:4,]-> a
B_UKB_afr[5,]-> b
a[,HVB_L:=sapply(a$Quant, function(X) as.numeric(gsub("\\]","",gsub("\\(","",gsub("\\[","",strsplit(X,",")[[1]])))))[1,]]
a[,HVB_U:=sapply(a$Quant, function(X) as.numeric(gsub("\\]","",gsub("\\(","",gsub("\\[","",strsplit(X,",")[[1]])))))[2,]]
b[,HVB_L:=1]
b[,HVB_U:=1]
rbind(a,b)->B_UKB_afr
B_UKB_afr[, Dataset:='UKB_afr']
B_UKB_afr[, boots_norm_L:=sapply(1:5, function(X) boots.ci.UKB_afr[[X]]$normal[2])]
B_UKB_afr[, boots_norm_U:=sapply(1:5, function(X) boots.ci.UKB_afr[[X]]$normal[3])]
B_UKB_afr[, boots_perc_L:=sapply(1:5, function(X) boots.ci.UKB_afr[[X]]$perc[4])]
B_UKB_afr[, boots_perc_U:=sapply(1:5, function(X) boots.ci.UKB_afr[[X]]$perc[5])]
B_UKB_afr[, boots_basic_L:=sapply(1:5, function(X) boots.ci.UKB_afr[[X]]$basic[4])]
B_UKB_afr[, boots_basic_U:=sapply(1:5, function(X) boots.ci.UKB_afr[[X]]$basic[5])]
#####################################################################################
#####################################################################################
##HRS_afr
prs<-vector('list', 22)
for(chr in 1:22){
        prs[[chr]]<-readRDS(paste0('~/height_prediction/loc_anc_analysis/output/chr', chr,'_0_prs_HRS_afr_EuR.Rds'))
}
samp_names<-names(prs[[1]])
dt<-data.table(SUBJID=samp_names, PRS_0=NA)

prs2<-vector('list', length(samp_names))
names(prs2)<-samp_names


for (S in samp_names){
        prs2[[S]]<-sum(prs[[1]][[S]],prs[[2]][[S]],prs[[3]][[S]],prs[[4]][[S]],prs[[5]][[S]],prs[[6]][[S]], prs[[7]][[S]], prs[[8]][[S]], prs[[9]][[S]], prs[[10]][[S]],prs[[11]][[S]], prs[[12]][[S]], prs[[13]][[S]],prs[[14]][[S]], prs[[15]][[S]], prs[[16]][[S]], prs[[17]][[S]], prs[[18]][[S]], prs[[19]][[S]], prs[[20]][[S]], prs[[21]][[S]],prs[[22]][[S]], na.rm=T)
}

dt[, PRS_0:=scale(unlist(prs2))]
#phenotype
fread('~/height_prediction/input/HRS_afr/HRS_AFR_phenotypes.txt', fill=T)-> Pheno_HRS_afr

Pheno_HRS_afr[, SUBJID:=paste0(as.character(Pheno_HRS_afr[, ID]),"_", as.character(Pheno_HRS_afr[,ID]))]
setkey(Pheno_HRS_afr, SUBJID)
#add ancestry
ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/HRS_afr/rfmix_anc_chr', X, '.txt'))))
anc_HRS_afr<-ancestry %>% dplyr::group_by(SUBJID) %>% dplyr::summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual
anc_HRS_afr[, SUBJID:=paste0(as.character(anc_HRS_afr[, SUBJID]),"_", as.character(anc_HRS_afr[,SUBJID]))]
setkey(anc_HRS_afr, SUBJID)
##

setkey(dt, SUBJID)
setkey(Pheno_HRS_afr, SUBJID)
dt[Pheno_HRS_afr][anc_HRS_afr]-> final
final[, HEIGHTX:=HEIGHT]
final[,AGE2:=AGE^2]

final2<-final[,c("SUBJID","PRS_0","AGE", "EUR_ANC", "HEIGHTX", 'AGE2', "SEX")]

#c
r2_hrs_afr<-partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_0, data=final2)) #0.87 %
saveRDS(final2, '~/height_prediction/loc_anc_analysis/output/all_PRS_HRS_afr_Eur.Rds')
all_datasets[['HRS_afr']]<-final2
remove(final2,dt)

all_datasets[['HRS_afr']][, Quantile:= cut(EUR_ANC,
                        breaks=quantile(EUR_ANC, probs=seq(0,1, by=1/2), na.rm=TRUE),
                        include.lowest=TRUE)]

all_datasets[['HRS_afr']][,Med_Eur_Anc:=median(EUR_ANC),by=Quantile]

as.character(unique((all_datasets[['HRS_afr']]$Quantile)))-> a

c(a[2], a[1])-> a1

r2_HRS_afr<-vector('list', length(a1))
names(r2_HRS_afr)<-a1
for(i in a1){
        r2_HRS_afr[[i]]<-partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC,all_datasets[['HRS_afr']][Quantile==i]),lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_0, all_datasets[['HRS_afr']][Quantile==i]))
}

B_HRS_afr<-data.table(Quant=c(a1, "total"),R_sq=c(unlist(r2_HRS_afr),r2_hrs_afr),
Med_Eur_Anc=c(unique(all_datasets[['HRS_afr']][Quantile==a1[1]][,Med_Eur_Anc]), unique(all_datasets[['HRS_afr']][Quantile==a1[2]][,Med_Eur_Anc]), median(all_datasets[['HRS_afr']][, EUR_ANC])))
B_HRS_afr[,N:=c(nrow(all_datasets[['HRS_afr']][Quantile==a1[1]]), nrow(all_datasets[['HRS_afr']][Quantile==a1[2]]),nrow(all_datasets[['HRS_afr']]))]
B_HRS_afr[,K:=1] #number of predictors. Need to check later if this is correct.
B_HRS_afr[, LCL:=CI.Rsq(R_sq, k=K, n=N)[3]]
B_HRS_afr[, UCL:=CI.Rsq(R_sq, k=K, n=N)[4]]

### add confidence intervals calculated with bootstrap: https://www.statmethods.net/advstats/bootstrapping.html
results.HRS_afr<-vector('list', length(a1)+1)
names(results.HRS_afr)<-c(a1, "total")
for(i in a1){
	boot(data=all_datasets[['HRS_afr']][Quantile==i], statistic=rsq.R2,R=999, formula1=HEIGHTX~SEX+AGE+AGE2+EUR_ANC, formula2=HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_0)-> results.HRS_afr[[i]]
}

results.HRS_afr[['total']] <- boot(data=all_datasets[['HRS_afr']], statistic=rsq.R2, R=999, formula1=HEIGHTX~SEX+AGE+AGE2+EUR_ANC, formula2=HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_0)
#confidence intervals
boots.ci.HRS_afr<-lapply(results.HRS_afr, function(X) boot.ci(X, type = c("norm", 'basic', "perc")))
names(boots.ci.HRS_afr)<-names(results.HRS_afr)

B_HRS_afr[1:2,]-> a
B_HRS_afr[3,]-> b
a[,HVB_L:=sapply(a$Quant, function(X) as.numeric(gsub("\\]","",gsub("\\(","",gsub("\\[","",strsplit(X,",")[[1]])))))[1,]]
a[,HVB_U:=sapply(a$Quant, function(X) as.numeric(gsub("\\]","",gsub("\\(","",gsub("\\[","",strsplit(X,",")[[1]])))))[2,]]
b[,HVB_L:=1]
b[,HVB_U:=1]
rbind(a,b)->B_HRS_afr
B_HRS_afr[, Dataset:='HRS_afr']
B_HRS_afr[, boots_norm_L:=sapply(1:3, function(X) boots.ci.HRS_afr[[X]]$normal[2])]
B_HRS_afr[, boots_norm_U:=sapply(1:3, function(X) boots.ci.HRS_afr[[X]]$normal[3])]
B_HRS_afr[, boots_perc_L:=sapply(1:3, function(X) boots.ci.HRS_afr[[X]]$perc[4])]
B_HRS_afr[, boots_perc_U:=sapply(1:3, function(X) boots.ci.HRS_afr[[X]]$perc[5])]
B_HRS_afr[, boots_basic_L:=sapply(1:3, function(X) boots.ci.HRS_afr[[X]]$basic[4])]
B_HRS_afr[, boots_basic_U:=sapply(1:3, function(X) boots.ci.HRS_afr[[X]]$basic[5])]
#####################################################################################
#####################################################################################
readRDS('~/height_prediction/gwas/HRS_eur/output/B_HRS_eur.Rds')[[i]]-> B_HRS_eur
readRDS('~/height_prediction/gwas/HRS_eur/output/results.HRS_eur.Rds')[[i]]->results.HRS_eur
ALL2<-rbind(B_JHS_afr[1:2,], B_WHI_afr[1:4,], B_UKB_afr[1:4,], B_HRS_afr[1:2,], B_HRS_eur)
ALL2$Dataset<-factor(ALL2$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "HRS_eur"))
tmp<-1/c(var(results.JHS_afr[[1]]$t), var(results.JHS_afr[[2]]$t), var(results.WHI_afr[[1]]$t), var(results.WHI_afr[[2]]$t), var(results.WHI_afr[[3]]$t), var(results.WHI_afr[[4]]$t), var(results.UKB_afr[[1]]$t),var(results.UKB_afr[[2]]$t), var(results.UKB_afr[[3]]$t), var(results.UKB_afr[[4]]$t), var(results.HRS_afr[[1]]$t), var(results.HRS_afr[[2]]$t),var(results.HRS_eur$t))  #weighing lm by boostrap replicates.

ALL2[,W:=tmp]

        B_plot<-ggplot(ALL2, aes(x=Med_Eur_Anc, y=R_sq,colour=Dataset)) +
        geom_point(aes(shape=Dataset), size=1.5, fill="white", alpha=0.8) + stat_smooth(data=ALL2,method = "lm", mapping = aes(weight = W), col='black') +
        geom_errorbar(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset,ymin=boots_perc_L, ymax=boots_perc_U), width=0.05, size=0.8) +
        geom_errorbarh(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset, xmin=HVB_L, xmax=HVB_U), width=0.05, size=0.5) +
        scale_color_manual(values=c(brewer.pal(4, 'Set1'),"#101010")) +
        ylab(expression(paste("Partial R"^"2"))) + xlab("European Ancestry Proportion") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y
= element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.key=element_blank(), legend.background=
element_blank(), legend.title=element_blank(), legend.text=element_text(size=15),legend.position = c(0.15,0.85))
        print(my_plot)
ggsave('~/height_prediction/loc_anc_analysis/figs/r2_plot_eur_only.pdf')
summary(lm(ALL$R_sq~ALL$Med_Eur_Anc))


plot_grid(A_plot, B_plot, labels=c("A","B"), nrow=1)

png('~/height_prediction/figs_for_paper/figs/Fig1.png')
