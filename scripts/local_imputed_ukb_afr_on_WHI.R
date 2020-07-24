#!/usr/bin/env Rscript
library(readr)
library(data.table)
library(reshape2)
library(asbio)
library(ggplot2)
library(dplyr)
library(parallel)
#########################
source('~/height_prediction/strat_prs/scripts/Rsq_R2.R')
source('~/height_prediction/gwas/WHI/scripts/PolygenicScore_v3.R')
source('~/height_prediction/gwas/WHI/scripts/short_fun_imputed.R')
library(TeachingDemos)

txtStart(paste0("~/height_prediction/gwas/WHI/local_imputed_out.txt"))
########################
cat('checkpoint number 1\n')
final_plink<-readRDS('~/height_prediction/loc_anc_analysis/output/final_plink.Rds')
setkey(final_plink, CHR, POS)
gc()
cat('checkpoint number 2\n')
######## skip skip
#res_all<-vector('list', 22)
#for(I in 1:22){
#	res_all[[I]]<-short_fun_imputed(args=I)
#	saveRDS(res_all[[I]], file=paste0("~/height_prediction/imputed/output/chr", I, "_prs.Rds"))
#	cat('Chr ', I, ' done\n')
#}

#saveRDS(res_all, file="~/height_prediction/imputed/output/all_prs.Rds")
#cat('checkpoint number X\n')

#res_all_HRS<-vector('list', 22)
#for(I in 1:22){
#        res_all_HRS[[I]]<-short_fun_imputed_v2(args=I)
#        saveRDS(res_all_HRS[[I]], file=paste0("~/height_prediction/imputed/output/chr", I, "_prs_HRS.Rds"))
#        cat('Chr ', I, ' done\n')
# }
#saveRDS(res_all_HRS, file="~/height_prediction/imputed/output/all_prs_HRS.Rds")

#res_all_JHS<-vector('list', 22)
#for(I in 1:22){
#        res_all_JHS[[I]]<-short_fun_imputed_v3(args=I)
#        saveRDS(res_all_JHS[[I]], file=paste0("~/height_prediction/imputed/output/chr", I, "_prs_JHS.Rds"))
#        cat('Chr ', I, ' done\n')
# }
#saveRDS(res_all_JHS, file="~/height_prediction/imputed/output/all_prs_JHS.Rds")
######## skip skip
readRDS("~/height_prediction/imputed/output/all_prs.Rds")-> res_all
names(res_all)<-c("PLINK", "PLINK_Tstat1")
readRDS("~/height_prediction/imputed/output/all_prs_HRS.Rds")-> res_all_HRS
names(res_all_HRS)<-c("PLINK", "PLINK_Tstat1")
readRDS("~/height_prediction/imputed/output/all_prs_JHS.Rds")-> res_all_JHS
names(res_all_JHS)<-c("PLINK", "PLINK_Tstat1")

cat('checkpoint number 3\n')
data.table(SUBJID=names(res_all[[1]][[1]]), PRS_plink=(unlist(res_all[[1]][[1]])+unlist(res_all[[2]][[1]])+unlist(res_all[[3]][[1]])+unlist(res_all[[4]][[1]])+unlist(res_all[[5]][[1]])+unlist(res_all[[6]][[1]])+unlist(res_all[[7]][[1]])+unlist(res_all[[8]][[1]])+unlist(res_all[[9]][[1]])+unlist(res_all[[10]][[1]])+unlist(res_all[[11]][[1]])+unlist(res_all[[12]][[1]])+unlist(res_all[[13]][[1]])+unlist(res_all[[14]][[1]])+unlist(res_all[[15]][[1]])+unlist(res_all[[16]][[1]])+unlist(res_all[[17]][[1]])+ unlist(res_all[[18]][[1]])+unlist(res_all[[19]][[1]])+unlist(res_all[[20]][[1]])+unlist(res_all[[21]][[1]])+unlist(res_all[[22]][[1]])),PRS_plink_tstat_1=(unlist(res_all[[1]][[2]])+unlist(res_all[[2]][[2]])+unlist(res_all[[3]][[2]])+unlist(res_all[[4]][[2]])+unlist(res_all[[5]][[2]])+unlist(res_all[[6]][[2]])+unlist(res_all[[7]][[2]])+unlist(res_all[[8]][[2]])+unlist(res_all[[9]][[2]])+unlist(res_all[[10]][[2]])+unlist(res_all[[11]][[2]])+unlist(res_all[[12]][[2]])+unlist(res_all[[13]][[2]])+unlist(res_all[[14]][[2]])+unlist(res_all[[15]][[2]])+unlist(res_all[[16]][[2]])+unlist(res_all[[17]][[2]])+ unlist(res_all[[18]][[2]])+unlist(res_all[[19]][[2]])+unlist(res_all[[20]][[2]])+unlist(res_all[[21]][[2]])+unlist(res_all[[22]][[2]])),PRS_EUR=unlist(readRDS('~/height_prediction/gwas/WHI/output/PGS_WHI_phys_100000_0.0005.Rds')))-> a


#a[, PRS_plink:=scale(PRS_plink)]
#a[, PRS_plink_tstat_1:=scale(PRS_plink_tstat_1)]

remove(res_all)
gc()
data.table(SUBJID=names(res_all_HRS[[1]][[1]]),PRS_plink=(unlist(res_all_HRS[[1]][[1]])+unlist(res_all_HRS[[2]][[1]])+unlist(res_all_HRS[[3]][[1]])+unlist(res_all_HRS[[4]][[1]])+unlist(res_all_HRS[[5]][[1]])+unlist(res_all_HRS[[6]][[1]])+unlist(res_all_HRS[[7]][[1]])+unlist(res_all_HRS[[8]][[1]])+unlist(res_all_HRS[[9]][[1]])+unlist(res_all_HRS[[10]][[1]])+unlist(res_all_HRS[[11]][[1]])+unlist(res_all_HRS[[12]][[1]])+unlist(res_all_HRS[[13]][[1]])+unlist(res_all_HRS[[14]][[1]])+unlist(res_all_HRS[[15]][[1]])+unlist(res_all_HRS[[16]][[1]])+unlist(res_all_HRS[[17]][[1]])+ unlist(res_all_HRS[[18]][[1]])+unlist(res_all_HRS[[19]][[1]])+unlist(res_all_HRS[[20]][[1]])+unlist(res_all_HRS[[21]][[1]])+unlist(res_all_HRS[[22]][[1]])),PRS_plink_tstat_1=(unlist(res_all_HRS[[1]][[2]])+unlist(res_all_HRS[[2]][[2]])+unlist(res_all_HRS[[3]][[2]])+unlist(res_all_HRS[[4]][[2]])+unlist(res_all_HRS[[5]][[2]])+unlist(res_all_HRS[[6]][[2]])+unlist(res_all_HRS[[7]][[2]])+unlist(res_all_HRS[[8]][[2]])+unlist(res_all_HRS[[9]][[2]])+unlist(res_all_HRS[[10]][[2]])+unlist(res_all_HRS[[11]][[2]])+unlist(res_all_HRS[[12]][[2]])+unlist(res_all_HRS[[13]][[2]])+unlist(res_all_HRS[[14]][[2]])+unlist(res_all_HRS[[15]][[2]])+unlist(res_all_HRS[[16]][[2]])+unlist(res_all_HRS[[17]][[2]])+ unlist(res_all_HRS[[18]][[2]])+unlist(res_all_HRS[[19]][[2]])+unlist(res_all_HRS[[20]][[2]])+unlist(res_all_HRS[[21]][[2]])+unlist(res_all_HRS[[22]][[2]])),PRS_EUR=unlist(readRDS('~/height_prediction/gwas/HRS_afr/output/PGS_HRS_afr_phys_100000_0.0005.Rds')))-> a_HRS

#a_HRS[, PRS_EUR:=scale(PRS_EUR)]
#a_HRS[, PRS_plink:=scale(PRS_plink)]
#a_HRS[, PRS_plink_tstat_1:=scale(PRS_plink_tstat_1)]

remove(res_all_HRS)
gc()
data.table(SUBJID=names(res_all_JHS[[1]][[1]]),PRS_plink=(unlist(res_all_JHS[[1]][[1]])+unlist(res_all_JHS[[2]][[1]])+unlist(res_all_JHS[[3]][[1]])+unlist(res_all_JHS[[4]][[1]])+unlist(res_all_JHS[[5]][[1]])+unlist(res_all_JHS[[6]][[1]])+unlist(res_all_JHS[[7]][[1]])+unlist(res_all_JHS[[8]][[1]])+unlist(res_all_JHS[[9]][[1]])+unlist(res_all_JHS[[10]][[1]])+unlist(res_all_JHS[[11]][[1]])+unlist(res_all_JHS[[12]][[1]])+unlist(res_all_JHS[[13]][[1]])+unlist(res_all_JHS[[14]][[1]])+unlist(res_all_JHS[[15]][[1]])+unlist(res_all_JHS[[16]][[1]])+unlist(res_all_JHS[[17]][[1]])+ unlist(res_all_JHS[[18]][[1]])+unlist(res_all_JHS[[19]][[1]])+unlist(res_all_JHS[[20]][[1]])+unlist(res_all_JHS[[21]][[1]])+unlist(res_all_JHS[[22]][[1]])),PRS_plink_tstat_1=(unlist(res_all_JHS[[1]][[2]])+unlist(res_all_JHS[[2]][[2]])+unlist(res_all_JHS[[3]][[2]])+unlist(res_all_JHS[[4]][[2]])+unlist(res_all_JHS[[5]][[2]])+unlist(res_all_JHS[[6]][[2]])+unlist(res_all_JHS[[7]][[2]])+unlist(res_all_JHS[[8]][[2]])+unlist(res_all_JHS[[9]][[2]])+unlist(res_all_JHS[[10]][[2]])+unlist(res_all_JHS[[11]][[2]])+unlist(res_all_JHS[[12]][[2]])+unlist(res_all_JHS[[13]][[2]])+unlist(res_all_JHS[[14]][[2]])+unlist(res_all_JHS[[15]][[2]])+unlist(res_all_JHS[[16]][[2]])+unlist(res_all_JHS[[17]][[2]])+ unlist(res_all_JHS[[18]][[2]])+unlist(res_all_JHS[[19]][[2]])+unlist(res_all_JHS[[20]][[2]])+unlist(res_all_JHS[[21]][[2]])+unlist(res_all_JHS[[22]][[2]])),PRS_EUR=unlist(readRDS('~/height_prediction/gwas/JHS/output/PGS_JHS_phys_100000_0.0005.Rds')))-> a_JHS

#a_JHS[, PRS_plink:=scale(PRS_plink)]
#a_JHS[, PRS_plink_tstat_1:=scale(PRS_plink_tstat_1)]
remove(res_all_JHS)

#phenotype
fread('~/height_prediction/input/WHI/WHI_phenotypes.txt')-> Pheno_WHI

Pheno_WHI[, SUBJID:=paste0("0_", as.character(Pheno_WHI[, SUBJID]))]
setkey(Pheno_WHI, SUBJID)
#add ancestry
ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/WHI/rfmix_anc_chr', X, '.txt'))))
anc_WHI<-ancestry %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual
setkey(anc_WHI, SUBJID)
##
setkey(a, SUBJID)
setkey(Pheno_WHI, SUBJID)
a[Pheno_WHI][anc_WHI]-> final
m1<-mean(final$HEIGHTX,na.rm=T)
sd1<-sd(final$HEIGHTX, na.rm=T)
final<-final[HEIGHTX>=m1-(2*sd1) & HEIGHTX<=m1+(2*sd1)][AFR_ANC>=0.05]
#
#HRS
fread('~/height_prediction/input/HRS_afr/HRS_AFR_phenotypes.txt', fill=T)[,1:4]-> Pheno_HRS
Pheno_HRS[, SUBJID:=paste0(ID, "_", ID)]
setkey(Pheno_HRS, SUBJID)
#add ancestry
ancestry_HRS<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/HRS_afr/rfmix_anc_chr', X, '.txt'))))
anc_HRS<-ancestry_HRS %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual
anc_HRS[, SUBJID:=paste0(SUBJID, "_", SUBJID)]
setkey(anc_HRS, SUBJID)
##
setkey(a_HRS, SUBJID)
setkey(Pheno_HRS, SUBJID)
a_HRS[Pheno_HRS][anc_HRS]-> final_HRS
final_HRS[,SEX:=ifelse(SEX==1, "Male", "Female")]
final_HRS[, HEIGHTX:=100*HEIGHT]
dt_f<-final_HRS[SEX=='Female']
dt_m<-final_HRS[SEX=='Male']
sd1_f<-sd(dt_f$HEIGHTX)
m1_f<-mean(dt_f$HEIGHTX)
sd1_m<-sd(dt_m$HEIGHTX)
m1_m<-mean(dt_m$HEIGHTX)
dt_f<-dt_f[HEIGHTX>=m1_f-(2*sd1_f)]
dt_m<-dt_m[HEIGHTX>=m1_m-(2*sd1_m)]
final_HRS<-rbind(dt_f, dt_m)
##JHS

fread('~/height_prediction/input/JHS/JHS_phenotypes.txt', fill=T)-> Pheno_JHS
#Pheno_JHS[, SUBJID:=paste0("0_", SUBJID, "_", FAM)]
setkey(Pheno_JHS, SUBJID)
#add ancestry
ancestry_JHS<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/JHS/rfmix_anc_chr', X, '.txt'))))
anc_JHS<-ancestry_JHS %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual
#anc_JHS[,SUBJID:=gsub(":", "_", SUBJID)]
anc_JHS$SUBJID<-substr(anc_JHS[,SUBJID],3,9)
setkey(anc_JHS, SUBJID)
##
a_JHS[,SUBJID:=substr(a_JHS[,SUBJID],1,7)]
setkey(a_JHS, SUBJID)
setkey(Pheno_JHS, SUBJID)
a_JHS[Pheno_JHS, nomatch=0][anc_JHS, nomatch=0]-> final_JHS
final_JHS[,HEIGHTX:=height_baseline]
final_JHS[,height_baseline:=NULL]
final_JHS[,AGE:=age_baseline]
final_JHS[AFR_ANC>=0.05]
final_JHS<-final_JHS[-which(!complete.cases(final_JHS$HEIGHTX))]

final[,c("SUBJID", "PRS_plink","PRS_EUR", "PRS_plink_tstat_1", "SEX", "HEIGHTX", "EUR_ANC","AFR_ANC","AGE")]-> final2
final_JHS[,c("SUBJID","PRS_plink","PRS_EUR", "PRS_plink_tstat_1",  "SEX", "HEIGHTX", "EUR_ANC","AFR_ANC","AGE")]-> final2_JHS
final_HRS[,c("SUBJID", "PRS_plink","PRS_EUR", "PRS_plink_tstat_1",  "SEX", "HEIGHTX", "EUR_ANC","AFR_ANC","AGE")]-> final2_HRS
cat('another checkpoint\n')

final2[,AGE2:=AGE^2][, PRS_eur_afr_plink:=(EUR_ANC*PRS_EUR)+(AFR_ANC*PRS_plink)] #changed this 03.03.2020
final2_JHS[,AGE2:=AGE^2][, PRS_eur_afr_plink:=(EUR_ANC*PRS_EUR)+(AFR_ANC*PRS_plink)]
final2_HRS[,AGE2:=AGE^2][, PRS_eur_afr_plink:=(EUR_ANC*PRS_EUR)+(AFR_ANC*PRS_plink)]# changed this 03.03.2020

fwrite(final2, file="~/height_prediction/imputed/output/all_prs_whi.txt", sep="\t")
fwrite(final2_JHS, file="~/height_prediction/imputed/output/all_prs_jhs.txt", sep="\t")
fwrite(final2_HRS, file="~/height_prediction/imputed/output/all_prs_hrs.txt", sep="\t")
#WHI
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2))*100 #3.578924 o#same as PRS1_WHI(alpha=0) and PRS2_WHI(alpha=0)
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_plink, data=final2))*100 #1.290213 #same as PRS1_WHI(alpha=1)
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_plink_tstat_1, data=final2))*100 #1.138397
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_eur_afr_plink, data=final2))*100 # #2.194672 same as PRS2_WHI(alpha=1)?
#here we see that there's potential in including AFR betas::
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2),lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR+PRS_eur_afr_plink, data=final2))*100 #0.3077749 #gain from PRS2
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2),lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR+PRS_plink, data=final2))*100 #0.2631669 #gain from PRS1
#
#JHS
partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC, data=final2_JHS), lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_EUR, data=final2_JHS))*100 #3.825324 same as PRS1(alpha=0) and PRS2(alpha=0)
partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC, data=final2_JHS), lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_plink, data=final2_JHS))*100 #1.484073 #same as PRS1(alpha=1)
partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC, data=final2_JHS), lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_plink_tstat_1, data=final2_JHS))*100 #1.666295
partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC, data=final2_JHS), lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_eur_afr_plink, data=final2_JHS))*100 # #1.966994 same as PRS2(alpha=1)
#here we see that there's potential in including AFR betas::
partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_EUR, data=final2_JHS),lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_EUR+PRS_eur_afr_plink, data=final2_JHS))*100 #0.1824423 #gain from PRS2
partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_EUR, data=final2_JHS),lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_EUR+PRS_plink, data=final2_JHS))*100 #0.2590213 #gain from PRS1
#HRS
partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC, data=final2_HRS), lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_EUR, data=final2_HRS))*100 #3.103462 same as PRS1(alpha=0) and PRS2(alpha-0)
partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC, data=final2_HRS), lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_plink, data=final2_HRS))*100 #0.7241404 same as PRS1(alpha=1)
partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC, data=final2_HRS), lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_plink_tstat_1, data=final2_HRS))*100 #0.818256
partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC, data=final2_HRS), lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_eur_afr_plink, data=final2_HRS))*100 # #1.24982 same as PRS2(alpha=1)
#here we see that there's potential in including AFR betas::
partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_EUR, data=final2_HRS),lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_EUR+PRS_eur_afr_plink, data=final2_HRS))*100 # 0.1121853 #gain from PRS2
partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_EUR, data=final2_HRS),lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_EUR+PRS_plink, data=final2_HRS))*100 #0.09691087 #gain from PRS1
###################
##################a
a_vec<-seq(0,1,0.01)
##################
PRS1_WHI<-function(alpha=0.9, dt='final2'){
        my_dt<-get(dt)
        my_dt[,PRS_comb:=((1-alpha)*PRS_EUR)+(alpha*PRS_plink)]
        res2<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=my_dt), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_comb,data=my_dt ))
        return(res2)
}
PRS1_JHS<-function(alpha=0.9, dt='final2_JHS'){
        my_dt<-get(dt)
        my_dt[,PRS_comb:=((1-alpha)*PRS_EUR)+(alpha*PRS_plink)]
        res2<-partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC, data=my_dt), lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_comb, data=my_dt))
        return(res2)
}
PRS1_HRS<-function(alpha=0.9, dt='final2_HRS'){
        my_dt<-get(dt)
        my_dt[,PRS_comb:=((1-alpha)*PRS_EUR)+(alpha*PRS_plink)]
        res2<-partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC, data=my_dt), lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_comb, data=my_dt))
        return(res2)
}
lapply(a_vec, function(X) PRS1_WHI(alpha=X))-> PRS1_WHI_res
lapply(a_vec, function(X) PRS1_JHS(alpha=X))-> PRS1_JHS_res
lapply(a_vec, function(X) PRS1_HRS(alpha=X))-> PRS1_HRS_res
rbind(data.table(Part_R2=unlist(PRS1_WHI_res), Alpha=a_vec, Dataset='WHI_afr'), data.table(Part_R2=unlist(PRS1_JHS_res), Alpha=a_vec, Dataset='JHS_afr'), data.table(Part_R2=unlist(PRS1_HRS_res), Alpha=a_vec, Dataset='HRS_afr'))-> PRS1_all_res
PRS1_all_res %>% group_by(Dataset) %>% summarise(Alpha=Alpha[which.max(Part_R2)], max(Part_R2)) %>% as.data.table
#
optimize(PRS1_WHI, interval=c(0,1), maximum=T, tol = 0.0001) #one liner for max 0.2099902//0.03832673
optimize(PRS1_JHS, interval=c(0,1), maximum=T, tol = 0.0001) #one liner for max 0.2145054//0.04074437
optimize(PRS1_HRS, interval=c(0,1), maximum=T, tol = 0.0001) #one liner for max 0.1315083//0.03197365
#
PRS1_all_res[, Test:='PRS1']
saveRDS(PRS1_all_res, file='~/height_prediction/loc_anc_analysis/output/PRS1_all_res.Rds')
ggplot(PRS1_all_res, aes(x=Alpha, y=Part_R2, colour=Dataset)) + geom_point(size=0.8)+ geom_line() +
#geom_vline(xintercept=optimize(PRS1_WHI, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, col='red', lty=2) +
#geom_vline(xintercept=optimize(PRS1_HRS_JHS, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, col='red', lty=2) +
geom_hline(yintercept=optimize(PRS1_WHI, interval=c(0,1), maximum=T, tol = 0.0001)$objective, col='red', lty=2) +
geom_hline(yintercept=optimize(PRS1_JHS, interval=c(0,1), maximum=T, tol = 0.0001)$objective, col='red', lty=2) +
geom_hline(yintercept=optimize(PRS1_HRS, interval=c(0,1), maximum=T, tol = 0.0001)$objective, col='red', lty=2)+
coord_cartesian(ylim = c(0, 0.048), xlim=c(0,0.6))
ggsave('~/height_prediction/loc_anc_analysis/figs/PRS1.pdf')
###################
PRS2_WHI<-function(alpha, dt='final2'){
        my_dt<-get(dt)
	my_dt[, PRS_comb:=(alpha*(1-EUR_ANC)*PRS_plink)+((1-alpha)+(alpha*EUR_ANC))*PRS_EUR]
        res<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=my_dt), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_comb, data=my_dt))
        return(res)
}

PRS2_JHS<-function(alpha, dt='final2_JHS'){
        my_dt<-get(dt)
	my_dt[, PRS_comb:=(alpha*(1-EUR_ANC)*PRS_plink)+((1-alpha)+(alpha*EUR_ANC))*PRS_EUR]
        res<-partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC, data=my_dt), lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_comb, data=my_dt))
        return(res)
}
PRS2_HRS<-function(alpha, dt='final2_HRS'){
        my_dt<-get(dt)
	my_dt[, PRS_comb:=(alpha*(1-EUR_ANC)*PRS_plink)+((1-alpha)+(alpha*EUR_ANC))*PRS_EUR]
	res<-partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC, data=my_dt), lm(HEIGHTX~SEX+AGE+AGE2+EUR_ANC+PRS_comb, data=my_dt))
        return(res)
}

PRS2_WHI_res<-lapply(a_vec, function(X) PRS2_WHI(alpha=X))
PRS2_JHS_res<-lapply(a_vec, function(X) PRS2_JHS(alpha=X))
PRS2_HRS_res<-lapply(a_vec, function(X) PRS2_HRS(alpha=X))

PRS2_all_res<-rbind(data.table(Part_R2=unlist(PRS2_WHI_res), Alpha=a_vec, Dataset='WHI_afr'), data.table(Part_R2=unlist(PRS2_JHS_res), Alpha=a_vec, Dataset='JHS_afr'),data.table(Part_R2=unlist(PRS2_HRS_res), Alpha=a_vec, Dataset='HRS_afr'))
PRS2_all_res %>% group_by(Dataset) %>% summarise(Alpha=Alpha[which.max(Part_R2)], max(Part_R2)) %>% as.data.table

optimize(PRS2_WHI, interval=c(0,1), maximum=T, tol = 0.0001) # 0.2919605//0.03875684
optimize(PRS2_JHS, interval=c(0,1), maximum=T, tol = 0.0001) # 0.2248384//0.04000787
optimize(PRS2_HRS, interval=c(0,1), maximum=T, tol = 0.0001) # 0.174482//0.03212165

PRS2_all_res[, Test:='PRS2']
saveRDS(PRS2_all_res, file='~/height_prediction/loc_anc_analysis/output/PRS2_all_res.Rds')
ggplot(PRS2_all_res, aes(x=Alpha, y=Part_R2, colour=Dataset)) + geom_point(size=0.8)+ geom_line() +
#geom_vline(xintercept=optimize(PRS2_WHI, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, col='red', lty=2) +
#geom_vline(xintercept=optimize(PRS2_HRS_JHS, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, col='red', lty=2) +
geom_hline(yintercept=optimize(PRS2_WHI, interval=c(0,1), maximum=T, tol = 0.0001)$objective, col='red', lty=2) +
geom_hline(yintercept=optimize(PRS2_JHS, interval=c(0,1), maximum=T, tol = 0.0001)$objective, col='red', lty=2) +
geom_hline(yintercept=optimize(PRS2_HRS, interval=c(0,1), maximum=T, tol = 0.0001)$objective, col='red', lty=2) +
coord_cartesian(ylim = c(0, 0.048), xlim=c(0,0.6))
ggsave('~/height_prediction/loc_anc_analysis/figs/PRS2.pdf')

rbind(PRS1_all_res, PRS2_all_res)->PRS1_2_all_res
saveRDS(PRS1_2_all_res, file='~/height_prediction/loc_anc_analysis/output/PRS1_2_all_res.Rds')
#The End
####################################################################################
txtStop()
