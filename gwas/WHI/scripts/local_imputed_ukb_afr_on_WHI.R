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
#
res_all<-vector('list', 22)
for(I in 1:22){
	res_all[[I]]<-short_fun_imputed(args=I)
	saveRDS(res_all[[I]], file=paste0("~/height_prediction/imputed/output/chr", I, "_prs.Rds"))
	cat('Chr ', I, ' done\n')
}


saveRDS(res_all, file="~/height_prediction/imputed/output/all_prs.Rds")
cat('checkpoint number X\n')

res_all_HRS<-vector('list', 22)
for(I in 1:22){
        res_all_HRS[[I]]<-short_fun_imputed_v2(args=I)
        saveRDS(res_all_HRS[[I]], file=paste0("~/height_prediction/imputed/output/chr", I, "_prs_HRS.Rds"))
        cat('Chr ', I, ' done\n')

 }
saveRDS(res_all_HRS, file="~/height_prediction/imputed/output/all_prs_HRS.Rds")

res_all_JHS<-vector('list', 22)
for(I in 1:22){
        res_all_JHS[[I]]<-short_fun_imputed_v3(args=I)
        saveRDS(res_all_JHS[[I]], file=paste0("~/height_prediction/imputed/output/chr", I, "_prs_JHS.Rds"))
        cat('Chr ', I, ' done\n')
 }
saveRDS(res_all_JHS, file="~/height_prediction/imputed/output/all_prs_JHS.Rds")

readRDS("~/height_prediction/imputed/output/all_prs.Rds")-> res_all
names(res_all)<-c("PLINK", "PLINK_Tstat1")
readRDS("~/height_prediction/imputed/output/all_prs_HRS.Rds")-> res_all_HRS
names(res_all_HRS)<-c("PLINK", "PLINK_Tstat1")
readRDS("~/height_prediction/imputed/output/all_prs_JHS.Rds")-> res_all_JHS
names(res_all_JHS)<-c("PLINK", "PLINK_Tstat1")

cat('checkpoint number 3\n')
data.table(SUBJID=names(res_all[[1]][[1]]), PRS_plink=(unlist(res_all[[1]][[1]])+unlist(res_all[[2]][[1]])+unlist(res_all[[3]][[1]])+unlist(res_all[[4]][[1]])+unlist(res_all[[5]][[1]])+unlist(res_all[[6]][[1]])+unlist(res_all[[7]][[1]])+unlist(res_all[[8]][[1]])+unlist(res_all[[9]][[1]])+unlist(res_all[[10]][[1]])+unlist(res_all[[11]][[1]])+unlist(res_all[[12]][[1]])+unlist(res_all[[13]][[1]])+unlist(res_all[[14]][[1]])+unlist(res_all[[15]][[1]])+unlist(res_all[[16]][[1]])+unlist(res_all[[17]][[1]])+ unlist(res_all[[18]][[1]])+unlist(res_all[[19]][[1]])+unlist(res_all[[20]][[1]])+unlist(res_all[[21]][[1]])+unlist(res_all[[22]][[1]])),PRS_plink_tstat_1=(unlist(res_all[[1]][[2]])+unlist(res_all[[2]][[2]])+unlist(res_all[[3]][[2]])+unlist(res_all[[4]][[2]])+unlist(res_all[[5]][[2]])+unlist(res_all[[6]][[2]])+unlist(res_all[[7]][[2]])+unlist(res_all[[8]][[2]])+unlist(res_all[[9]][[2]])+unlist(res_all[[10]][[2]])+unlist(res_all[[11]][[2]])+unlist(res_all[[12]][[2]])+unlist(res_all[[13]][[2]])+unlist(res_all[[14]][[2]])+unlist(res_all[[15]][[2]])+unlist(res_all[[16]][[2]])+unlist(res_all[[17]][[2]])+ unlist(res_all[[18]][[2]])+unlist(res_all[[19]][[2]])+unlist(res_all[[20]][[2]])+unlist(res_all[[21]][[2]])+unlist(res_all[[22]][[2]])),PRS_EUR=unlist(readRDS('~/height_prediction/gwas/WHI/output/PGS_WHI_phys_100000_0.0005.Rds')))-> a


a[, PRS_plink:=scale(PRS_plink)]
a[, PRS_plink_tstat_1:=scale(PRS_plink_tstat_1)]

remove(res_all)
gc()
data.table(SUBJID=names(res_all_HRS[[1]][[1]]),PRS_plink=(unlist(res_all_HRS[[1]][[1]])+unlist(res_all_HRS[[2]][[1]])+unlist(res_all_HRS[[3]][[1]])+unlist(res_all_HRS[[4]][[1]])+unlist(res_all_HRS[[5]][[1]])+unlist(res_all_HRS[[6]][[1]])+unlist(res_all_HRS[[7]][[1]])+unlist(res_all_HRS[[8]][[1]])+unlist(res_all_HRS[[9]][[1]])+unlist(res_all_HRS[[10]][[1]])+unlist(res_all_HRS[[11]][[1]])+unlist(res_all_HRS[[12]][[1]])+unlist(res_all_HRS[[13]][[1]])+unlist(res_all_HRS[[14]][[1]])+unlist(res_all_HRS[[15]][[1]])+unlist(res_all_HRS[[16]][[1]])+unlist(res_all_HRS[[17]][[1]])+ unlist(res_all_HRS[[18]][[1]])+unlist(res_all_HRS[[19]][[1]])+unlist(res_all_HRS[[20]][[1]])+unlist(res_all_HRS[[21]][[1]])+unlist(res_all_HRS[[22]][[1]])),PRS_plink_tstat_1=(unlist(res_all_HRS[[1]][[2]])+unlist(res_all_HRS[[2]][[2]])+unlist(res_all_HRS[[3]][[2]])+unlist(res_all_HRS[[4]][[2]])+unlist(res_all_HRS[[5]][[2]])+unlist(res_all_HRS[[6]][[2]])+unlist(res_all_HRS[[7]][[2]])+unlist(res_all_HRS[[8]][[2]])+unlist(res_all_HRS[[9]][[2]])+unlist(res_all_HRS[[10]][[2]])+unlist(res_all_HRS[[11]][[2]])+unlist(res_all_HRS[[12]][[2]])+unlist(res_all_HRS[[13]][[2]])+unlist(res_all_HRS[[14]][[2]])+unlist(res_all_HRS[[15]][[2]])+unlist(res_all_HRS[[16]][[2]])+unlist(res_all_HRS[[17]][[2]])+ unlist(res_all_HRS[[18]][[2]])+unlist(res_all_HRS[[19]][[2]])+unlist(res_all_HRS[[20]][[2]])+unlist(res_all_HRS[[21]][[2]])+unlist(res_all_HRS[[22]][[2]])),PRS_EUR=unlist(readRDS('~/height_prediction/gwas/HRS_afr/output/PGS_HRS_afr_phys_100000_0.0005.Rds')))-> a_HRS

#a_HRS[, PRS_EUR:=scale(PRS_EUR)]
a_HRS[, PRS_plink:=scale(PRS_plink)]
a_HRS[, PRS_plink_tstat_1:=scale(PRS_plink_tstat_1)]

remove(res_all_HRS)
gc()
data.table(SUBJID=names(res_all_JHS[[1]][[1]]),PRS_plink=(unlist(res_all_JHS[[1]][[1]])+unlist(res_all_JHS[[2]][[1]])+unlist(res_all_JHS[[3]][[1]])+unlist(res_all_JHS[[4]][[1]])+unlist(res_all_JHS[[5]][[1]])+unlist(res_all_JHS[[6]][[1]])+unlist(res_all_JHS[[7]][[1]])+unlist(res_all_JHS[[8]][[1]])+unlist(res_all_JHS[[9]][[1]])+unlist(res_all_JHS[[10]][[1]])+unlist(res_all_JHS[[11]][[1]])+unlist(res_all_JHS[[12]][[1]])+unlist(res_all_JHS[[13]][[1]])+unlist(res_all_JHS[[14]][[1]])+unlist(res_all_JHS[[15]][[1]])+unlist(res_all_JHS[[16]][[1]])+unlist(res_all_JHS[[17]][[1]])+ unlist(res_all_JHS[[18]][[1]])+unlist(res_all_JHS[[19]][[1]])+unlist(res_all_JHS[[20]][[1]])+unlist(res_all_JHS[[21]][[1]])+unlist(res_all_JHS[[22]][[1]])),PRS_plink_tstat_1=(unlist(res_all_JHS[[1]][[2]])+unlist(res_all_JHS[[2]][[2]])+unlist(res_all_JHS[[3]][[2]])+unlist(res_all_JHS[[4]][[2]])+unlist(res_all_JHS[[5]][[2]])+unlist(res_all_JHS[[6]][[2]])+unlist(res_all_JHS[[7]][[2]])+unlist(res_all_JHS[[8]][[2]])+unlist(res_all_JHS[[9]][[2]])+unlist(res_all_JHS[[10]][[2]])+unlist(res_all_JHS[[11]][[2]])+unlist(res_all_JHS[[12]][[2]])+unlist(res_all_JHS[[13]][[2]])+unlist(res_all_JHS[[14]][[2]])+unlist(res_all_JHS[[15]][[2]])+unlist(res_all_JHS[[16]][[2]])+unlist(res_all_JHS[[17]][[2]])+ unlist(res_all_JHS[[18]][[2]])+unlist(res_all_JHS[[19]][[2]])+unlist(res_all_JHS[[20]][[2]])+unlist(res_all_JHS[[21]][[2]])+unlist(res_all_JHS[[22]][[2]])),PRS_EUR=unlist(readRDS('~/height_prediction/gwas/JHS/output/PGS_JHS_phys_100000_0.0005.Rds')))-> a_JHS

a_JHS[, PRS_plink:=scale(PRS_plink)]
a_JHS[, PRS_plink_tstat_1:=scale(PRS_plink_tstat_1)]
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
final_JHS[,SEX:=ifelse(SEX==1, "Male", "Female")]
final[,c("SUBJID", "PRS_plink","PRS_EUR", "PRS_plink_tstat_1", "SEX", "HEIGHTX", "EUR_ANC","AFR_ANC","AGE")]-> final2
final_JHS[,c("SUBJID","PRS_plink","PRS_EUR", "PRS_plink_tstat_1",  "SEX", "HEIGHTX", "EUR_ANC","AFR_ANC","AGE")]-> final2_JHS
final_HRS[,c("SUBJID", "PRS_plink","PRS_EUR", "PRS_plink_tstat_1",  "SEX", "HEIGHTX", "EUR_ANC","AFR_ANC","AGE")]-> final2_HRS
final2_HRS_JHS<-rbind(final2_HRS[, Dt:='HRS'], final2_JHS[,Dt:='JHS'])
#ok, so for ukb_afr and WHI (both), POP1 is AFR
cat('another checkpoint\n')

final2[,AGE2:=AGE^2][, PRS_eur_afr_plink:=(mean(EUR_ANC)*(PRS_EUR))+(mean(AFR_ANC)*PRS_plink)]
final2_HRS_JHS[,AGE2:=AGE^2][, PRS_eur_afr_plink:=(mean(EUR_ANC)*(PRS_EUR))+(mean(AFR_ANC)*PRS_plink)]

fwrite(final2, file="~/height_prediction/imputed/all_prs_whi.txt", sep="\t")
fwrite(final2_HRS_JHS, file="~/height_prediction/imputed/output/all_prs_hrs_jsh.txt", sep="\t")

partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2))*100 #4.12%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_plink, data=final2))*100 #1.68
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_plink_tstat_1, data=final2))*100 #1.55%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_eur_afr_plink, data=final2))*100 # #4.47%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2),lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR+PRS_eur_afr_plink, data=final2))*100 #0.40
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2),lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR+PRS_plink, data=final2))*100 #0.40
#
partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+Dt+EUR_ANC, data=final2_HRS_JHS), lm(HEIGHTX~SEX+AGE+AGE2+Dt+EUR_ANC+PRS_EUR, data=final2_HRS_JHS))*100 #2.97%
partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+Dt+EUR_ANC, data=final2_HRS_JHS), lm(HEIGHTX~SEX+AGE+AGE2+Dt+EUR_ANC+PRS_plink, data=final2_HRS_JHS))*100 #1.05%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2_HRS_JHS), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_eur_afr_plink, data=final2_HRS_JHS))*100 #2.41
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2_HRS_JHS),lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR+PRS_eur_afr_plink, data=final2_HRS_JHS))*100 #0.07%

###################
##################a
a_vec<-seq(0,1,0.001)
##################
my_alpha_v2<-function(alpha=0.9, dt='final2'){
        my_dt<-get(dt)
        my_dt[,PRS_comb:=((1-alpha)*PRS_EUR)+(alpha*PRS_plink)]
        res2<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=my_dt), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_comb,data=my_dt ))
        return(res2)
}
my_alpha_conc_v2<-function(alpha=0.9, dt='final2_HRS_JHS'){
        my_dt<-get(dt)
        my_dt[,PRS_comb:=((1-alpha)*PRS_EUR)+(alpha*PRS_plink)]
        res2<-partial.R2(lm(HEIGHTX~SEX+AGE+AGE2+Dt+EUR_ANC, data=my_dt), lm(HEIGHTX~SEX+AGE+AGE2+Dt+EUR_ANC+PRS_comb, data=my_dt))
        return(res2)
}
lapply(a_vec, function(X) my_alpha_v2(alpha=X))-> wei_PRS_v2
lapply(a_vec, function(X) my_alpha_conc_v2(alpha=X))-> wei_PRS_conc_v2
rbind(data.table(part_R2=unlist(wei_PRS_v2), alfa=a_vec, Dataset='WHI'), data.table(part_R2=unlist(wei_PRS_conc_v2), alfa=a_vec, Dataset='JHS+HRS'))-> wanna_plot_v2

#
saveRDS(wanna_plot_v2, file='~/height_prediction/imputed/output/temp_dt.Rds')
ggplot(wanna_plot_v2, aes(x=alfa, y=part_R2, colour=Dataset)) + geom_point(size=0.8)+ geom_line() +
geom_vline(xintercept=optimize(my_alpha_v2, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, col='red', lty=2) +
geom_vline(xintercept=optimize(my_alpha_conc_v2, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, col='red', lty=2) +
geom_hline(yintercept=optimize(my_alpha_v2, interval=c(0,1), maximum=T, tol = 0.0001)$objective, col='red', lty=2) +
geom_hline(yintercept=optimize(my_alpha_conc_v2, interval=c(0,1), maximum=T, tol = 0.0001)$objective, col='red', lty=2) +
coord_cartesian(ylim = c(0, 0.048), xlim=c(0,0.6))
ggsave('~/height_prediction/imputed/figs/alfa_plink.pdf')
optimize(my_alpha_v2, interval=c(0,1), maximum=T, tol = 0.0001) #one liner for max 0.7094478 //0.04503568
optimize(my_alpha_conc_v2, interval=c(0,1), maximum=T, tol = 0.0001) #one liner for max  0.6716053//  0.03176889

##################
a_vec2<-seq(from=0, to=1, by=0.001)

###################
max_alpha_v2<-function(alpha, dt='final2'){
        my_dt<-get(dt)
        my_dt[,PRS_comb:=(AFR_ANC*alpha*PRS_plink)+((1-(alpha*AFR_ANC))*PRS_EUR)]
        res<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=my_dt), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_comb, data=my_dt))
        return(res)
}

max_alpha_conc_v2<-function(alpha, dt='final2_HRS_JHS'){
        my_dt<-get(dt)
        my_dt[,PRS_comb:=(AFR_ANC*alpha*PRS_plink)+((1-(alpha*AFR_ANC))*PRS_EUR)]
        res<-partial.R2(lm(HEIGHTX~Dt+SEX+AGE+AGE2+EUR_ANC, data=my_dt), lm(HEIGHTX~Dt+SEX+AGE+AGE2+EUR_ANC+PRS_comb, data=my_dt))
        return(res)
}

all_prs_v2<-lapply(a_vec, function(X) max_alpha_v2(alpha=X))
all_prs_conc_v2<-lapply(a_vec, function(X) max_alpha_conc_v2(alpha=X))

optimize(max_alpha_v2, interval=c(0,1), maximum=T, tol = 0.0001) # 0.7749016//0.04514398
optimize(max_alpha_conc_v2, interval=c(0,1), maximum=T, tol = 0.0001) #  0.6501442 // 0.03099347

temp_dt_v2<-rbind(data.table(part_R2=unlist(all_prs_v2), alfa=a_vec, Dataset='WHI'), data.table(part_R2=unlist(all_prs_conc_v2), alfa=a_vec, Dataset='JHS+HRS'))

saveRDS(temp_dt_v2, file='~/height_prediction/imputed/output/temp_dt_v2.Rds')
ggplot(temp_dt_v2, aes(x=alfa, y=part_R2, colour=Dataset)) + geom_point(size=0.8)+ geom_line() +
geom_vline(xintercept=optimize(max_alpha_v2, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, col='red', lty=2) +
geom_vline(xintercept=optimize(max_alpha_conc_v2, interval=c(0,1), maximum=T, tol = 0.0001)$maximum, col='red', lty=2) +
geom_hline(yintercept=optimize(max_alpha_v2, interval=c(0,1), maximum=T, tol = 0.0001)$objective, col='red', lty=2) +
geom_hline(yintercept=optimize(max_alpha_conc_v2, interval=c(0,1), maximum=T, tol = 0.0001)$objective, col='red', lty=2) +
coord_cartesian(ylim = c(0, 0.048), xlim=c(0,0.6))
ggsave('~/height_prediction/imputed/figs/alfa_plink_withAnc.pdf')





#The End

#playground
LA<-readRDS('~/height_prediction/loc_anc_analysis/output/part_R2_WHI.Rds')[, Dataset:='WHI']
rbind(LA, data.table(Alfa=LA$Alfa, part_R2=NA,Dataset='JHS+HRS'))-> LA #for now
LA[, alpha:=Alfa][,Alfa:=NULL]
LA[, PRS3:=part_R2][,part_R2:=NULL]
test4<-cbind(wanna_plot_v2, temp_dt_v2)
test4[,alfa:=NULL]
test4[, alpha:=alfa]
test4[,alfa:=NULL]
colnames(test4)[1]<-"PRS1"
colnames(test4)[3]<-"PRS2"
test4[,Dataset:=NULL]

merge(test4, LA,by=c('alpha', 'Dataset'))-> test4
melt(test4,id=c('Dataset', 'alpha'))-> test5

 
test5[, Dataset2:=ifelse(Dataset=='WHI', "Women's Health Initiative", "Jackson Heart Study + \nHealth and Retirement Study")]
test5[alpha<=0.5]-> test5
ggplot(test5, aes(x=alpha, y=value,colour=variable)) + facet_wrap(~Dataset2) +
geom_line(size=1.2) + coord_cartesian(xlim=c(0,0.5), ylim=c(0.02, 0.048)) + scale_color_manual(values=c("darkseagreen4", "darkslateblue", "darkorange3")) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.key=element_blank(), legend.background=element_blank(),legend.title=element_blank(), legend.text=element_text(size=18)) + labs(x=expression(alpha), y=expression(Partial~R^2)) 

ggsave('~/height_prediction/imputed/figs/test_this_plink.pdf')
###

#inlcuyding local ancestry PRS


test4<-cbind(wanna_plot_v2, temp_dt_v2)
test4[,alfa:=NULL]
test4[, alpha:=alfa]
test4[,alfa:=NULL]
colnames(test4)[1]<-"PRS1"
colnames(test4)[3]<-"PRS2"
test4[,Dataset:=NULL]
##test4[Dataset=='WHI']-> test4

#test4[,PRS3:=readRDS('~/height_prediction/loc_anc_analysis/output/temp_dt.Rds')$part_R2]
#test4[,PRS4:=readRDS('~/height_prediction/loc_anc_analysis/output/temp_dt_v2.Rds')$part_R2]

melt(test4,id=c('Dataset', 'alpha'))-> test5


#test5[, Dataset2:=ifelse(Dataset=='WHI', "Women's Health Initiative", "Jackson Heart Study + \nHealth and Retirement Study")]
ggplot(test5, aes(x=alpha, y=value,colour=variable))  + facet_wrap(~Dataset2) + 
geom_line(size=1.2) + coord_cartesian(xlim=c(0,0.5), ylim=c(0.02, 0.048)) + scale_color_manual(values=c("darkseagreen4", "darkslateblue", "deeppink4", "gray7")) +
theme(strip.text.x = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.key=element_blank(), legend.background=element_blank(),legend.title=element_blank(), legend.text=element_text(size=18)) + labs(x=expression(alpha), y=expression(Partial~R^2))
ggsave('~/height_prediction/loc_anc_analysis/figs/multi_prs.pdf')

txtStop()
