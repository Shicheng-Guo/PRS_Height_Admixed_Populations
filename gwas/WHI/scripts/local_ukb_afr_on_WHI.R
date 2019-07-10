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
source('~/height_prediction/gwas/WHI/scripts/short_fun.R')
########################
cat('checkpoint number 1\n')

plink<-fread('~/height_prediction/runSmartpCA-master/UKB_AFR/association_v3.Res.Height.glm.linear.adjusted', header=T, fill=T)
colnames(plink)[2]<-'MarkerName'
colnames(plink)[1]<-'CHR'
N<-8816*2
#
plink2<-fread('~/height_prediction/runSmartpCA-master/UKB_AFR/test3.txt', fill=T)
colnames(plink2)<-c("CHR","POS", "MarkerName","REF","ALT","A1","TEST","OBS_CT","PLINK", "SE","T_STAT", "UNADJ") #plink is BETA
setkey(plink, MarkerName, CHR, UNADJ)
setkey(plink2, MarkerName, CHR, UNADJ)
plink[plink2, nomatch=0]-> final_plink
final_plink$CHR<-as.numeric(final_plink$CHR)
arrange(final_plink, CHR,POS) %>% as.data.table -> final_plink
setkey(final_plink, MarkerName, CHR, POS)

cat('checkpoint number 2\n')

res_all<-vector('list', 22)
for(I in 1:22){
	short_fun(args=I)-> res
	saveRDS(res, file=paste0("~/height_prediction/gwas/WHI/output/chr", I, "_prs.Rds"))
	res_all[[I]]<-res
	cat('Chr ', I, ' done\n')
}

saveRDS(res_all, file="~/height_prediction/gwas/WHI/output/all_prs.Rds")

cat('checkpoint number 3\n')
data.table(SUBJID=names(res_all[[1]][[1]]), 
PRS_POP1=(unlist(res_all[[1]][[1]])+unlist(res_all[[2]][[1]])+unlist(res_all[[3]][[1]])+unlist(res_all[[4]][[1]])+unlist(res_all[[5]][[1]])+unlist(res_all[[6]][[1]])+unlist(res_all[[7]][[1]])+unlist(res_all[[8]][[1]])+unlist(res_all[[9]][[1]])+unlist(res_all[[10]][[1]])+unlist(res_all[[11]][[1]])+unlist(res_all[[12]][[1]])+unlist(res_all[[13]][[1]])+unlist(res_all[[14]][[1]])+unlist(res_all[[15]][[1]])+unlist(res_all[[16]][[1]])+unlist(res_all[[17]][[1]])+ unlist(res_all[[18]][[1]])+unlist(res_all[[19]][[1]])+unlist(res_all[[20]][[1]])+unlist(res_all[[21]][[1]])+unlist(res_all[[22]][[1]])), 
PRS_POP2=(unlist(res_all[[1]][[2]])+unlist(res_all[[2]][[2]])+unlist(res_all[[3]][[2]])+unlist(res_all[[4]][[2]])+unlist(res_all[[5]][[2]])+unlist(res_all[[6]][[2]])+unlist(res_all[[7]][[2]])+unlist(res_all[[8]][[2]])+unlist(res_all[[9]][[2]])+unlist(res_all[[10]][[2]])+unlist(res_all[[11]][[2]])+unlist(res_all[[12]][[2]])+unlist(res_all[[13]][[2]])+unlist(res_all[[14]][[2]])+unlist(res_all[[15]][[2]])+unlist(res_all[[16]][[2]])+unlist(res_all[[17]][[2]])+ unlist(res_all[[18]][[2]])+unlist(res_all[[19]][[2]])+unlist(res_all[[20]][[2]])+unlist(res_all[[21]][[2]])+unlist(res_all[[22]][[2]])), 
PRS_all=(unlist(res_all[[1]][[3]])+unlist(res_all[[2]][[3]])+unlist(res_all[[3]][[3]])+unlist(res_all[[4]][[3]])+unlist(res_all[[5]][[3]])+unlist(res_all[[6]][[3]])+unlist(res_all[[7]][[3]])+unlist(res_all[[8]][[3]])+unlist(res_all[[9]][[3]])+unlist(res_all[[10]][[3]])+unlist(res_all[[11]][[3]])+unlist(res_all[[12]][[3]])+unlist(res_all[[13]][[3]])+unlist(res_all[[14]][[3]])+unlist(res_all[[15]][[3]])+unlist(res_all[[16]][[3]])+unlist(res_all[[17]][[3]])+ unlist(res_all[[18]][[3]])+unlist(res_all[[19]][[3]])+unlist(res_all[[20]][[3]])+unlist(res_all[[21]][[3]])+unlist(res_all[[22]][[3]])), 
PRS_all_tstat_1=(unlist(res_all[[1]][[4]])+unlist(res_all[[2]][[4]])+unlist(res_all[[3]][[4]])+unlist(res_all[[4]][[4]])+unlist(res_all[[5]][[4]])+unlist(res_all[[6]][[4]])+unlist(res_all[[7]][[4]])+unlist(res_all[[8]][[4]])+unlist(res_all[[9]][[4]])+unlist(res_all[[10]][[4]])+unlist(res_all[[11]][[4]])+unlist(res_all[[12]][[4]])+unlist(res_all[[13]][[4]])+unlist(res_all[[14]][[4]])+unlist(res_all[[15]][[4]])+unlist(res_all[[16]][[4]])+unlist(res_all[[17]][[4]])+ unlist(res_all[[18]][[4]])+unlist(res_all[[19]][[4]])+unlist(res_all[[20]][[4]])+unlist(res_all[[21]][[4]])+unlist(res_all[[22]][[4]])),
PRS_plink=(unlist(res_all[[1]][[5]])+unlist(res_all[[2]][[5]])+unlist(res_all[[3]][[5]])+unlist(res_all[[4]][[5]])+unlist(res_all[[5]][[5]])+unlist(res_all[[6]][[5]])+unlist(res_all[[7]][[5]])+unlist(res_all[[8]][[5]])+unlist(res_all[[9]][[5]])+unlist(res_all[[10]][[5]])+unlist(res_all[[11]][[5]])+unlist(res_all[[12]][[5]])+unlist(res_all[[13]][[5]])+unlist(res_all[[14]][[5]])+unlist(res_all[[15]][[5]])+unlist(res_all[[16]][[5]])+unlist(res_all[[17]][[5]])+ unlist(res_all[[18]][[5]])+unlist(res_all[[19]][[5]])+unlist(res_all[[20]][[5]])+unlist(res_all[[21]][[5]])+unlist(res_all[[22]][[5]])),
PRS_plink_tstat_1=(unlist(res_all[[1]][[6]])+unlist(res_all[[2]][[6]])+unlist(res_all[[3]][[6]])+unlist(res_all[[4]][[6]])+unlist(res_all[[5]][[6]])+unlist(res_all[[6]][[6]])+unlist(res_all[[7]][[6]])+unlist(res_all[[8]][[6]])+unlist(res_all[[9]][[6]])+unlist(res_all[[10]][[6]])+unlist(res_all[[11]][[6]])+unlist(res_all[[12]][[6]])+unlist(res_all[[13]][[6]])+unlist(res_all[[14]][[6]])+unlist(res_all[[15]][[6]])+unlist(res_all[[16]][[6]])+unlist(res_all[[17]][[6]])+ unlist(res_all[[18]][[6]])+unlist(res_all[[19]][[6]])+unlist(res_all[[20]][[6]])+unlist(res_all[[21]][[6]])+unlist(res_all[[22]][[6]])),
PRS_EUR=unlist(readRDS('~/height_prediction/gwas/WHI/output/PGS_WHI_phys_100000_0.0005.Rds')))-> a

a[, PRS_POP1:=scale(PRS_POP1)]
a[, PRS_POP2:=scale(PRS_POP2)]
a[, PRS_all:=scale(PRS_all)]
a[, PRS_EUR:=scale(PRS_EUR)]
a[, PRS_plink:=scale(PRS_plink)]
a[, PRS_plink_tstat_1:=scale(PRS_plink_tstat_1)]
a[, PRS_all_tstat_1:=scale(PRS_all_tstat_1)]

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

final[,c("SUBJID","PRS_POP1","PRS_POP2","PRS_all", "PRS_plink","PRS_EUR","PRS_all_tstat_1", "PRS_plink_tstat_1",  "SEX", "HEIGHTX", "EUR_ANC","AFR_ANC","AGE")]-> final2

#ok, so for ukb_afr and WHI (both), POP1 is AFR
cat('another checkpoint\n')

final2[,AGE2:=AGE^2][, PRS_eur_afr_local:=(mean(EUR_ANC)*(PRS_EUR))+(mean(AFR_ANC)*PRS_all)]
final2[,AGE2:=AGE^2][, PRS_eur_afr_plink:=(mean(EUR_ANC)*(PRS_EUR))+(mean(AFR_ANC)*PRS_plink)]
fwrite(final2, file="~/height_prediction/gwas/WHI/output/all_prs_whi.txt", sep="\t")

partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2))*100 #4.1%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_all, data=final2))*100 #0.35%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_plink, data=final2))*100 #0.39%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_POP1, data=final2))*100 #0.09 %
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_POP2, data=final2))*100 #0.1896%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_plink_tstat_1, data=final2))*100 #0.34%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_all_tstat_1, data=final2))*100 #0.31%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_eur_afr_local, data=final2))*100 #1.197%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_eur_afr_plink, data=final2))*100 # #1.289
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2),lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR+PRS_eur_afr_plink, data=final2))*100 #0.07
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2),lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR+PRS_eur_afr_local, data=final2))*100 #0.05

###################
##################
a_vec<-seq(0,1,0.001)
wei_PRS<-vector('list', length(a_vec))
names(wei_PRS)<-a_vec
part_r2<-c()
for(i in 1:length(a_vec)){
	cat(i, '\n')
	wei_PRS[[i]]<-((1-a_vec[i])*final2$PRS_EUR)+(a_vec[i]*final2$PRS_all)
	part_r2[i]<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(final2$HEIGHTX~final2$AGE+final2$AGE2+final2$EUR_ANC+wei_PRS[[i]]))
}

data.table(part_R2=part_r2, alfa=a_vec)-> wanna_plot

ggplot(wanna_plot, aes(x=alfa, y=part_R2)) + geom_point(size=0.8)+ geom_line() + geom_vline(xintercept=wanna_plot[which.max(wanna_plot$part_R2),alfa], col='red') + geom_hline(yintercept=partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2)), col='red') + coord_cartesian(ylim = c(0, 0.048))
ggsave('~/height_prediction/gwas/WHI/figs/alfa_part_R2_WHI_all.pdf')
###################
##################
a_vec<-seq(0,1,0.001)
wei_PRS<-vector('list', length(a_vec))
names(wei_PRS)<-a_vec
part_r2<-c()
for(i in 1:length(a_vec)){
	cat(i, '\n')
	wei_PRS[[i]]<-((1-a_vec[i])*final2$PRS_EUR)+(a_vec[i]*final2$PRS_plink)
	part_r2[i]<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(final2$HEIGHTX~final2$AGE+final2$AGE2+final2$EUR_ANC+wei_PRS[[i]]))
}

data.table(part_R2=part_r2, alfa=a_vec)-> wanna_plot2

ggplot(wanna_plot2, aes(x=alfa, y=part_R2)) + geom_point(size=0.8)+ geom_line() + geom_vline(xintercept=wanna_plot[which.max(wanna_plot$part_R2),alfa], col='red') + geom_hline(yintercept=partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2)), col='red')

ggsave('~/height_prediction/gwas/WHI/figs/alfa_part_R2_WHI_plink.pdf')
###################
##################
a_vec2<-seq(from=0, to=1, by=0.001)

all_indivs<-vector('list', nrow(final2))
names(all_indivs)<-final2[,SUBJID]
counter<-0
system.time(for(S in final2[,SUBJID]){
	wei_PRS_2<-vector('list', length(a_vec2))
	names(wei_PRS_2)<-a_vec2
	wei_PRS_2<-lapply(1:length(a_vec2), function(i)( final2[SUBJID==S]$AFR_ANC*a_vec2[i]*final2[SUBJID==S]$PRS_all[[1]])+(final2[SUBJID==S]$AFR_ANC*(1-a_vec2[i])*final2[SUBJID==S]$PRS_EUR[[1]]))
	counter+1-> counter
        cat(counter, "\r")
	all_indivs[[S]]<-wei_PRS_2
	gc()
}
) #5240 seconds
temp_dt<-matrix(ncol=2, nrow=length(a_vec2))
colnames(temp_dt)<-c('alpha', 'part_r2')
temp_dt[,1]<-a_vec2
for (j in 1:length(a_vec2)){
	temp_dt[j,2]<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(final2$HEIGHTX~final2$AGE+final2$AGE2+final2$EUR_ANC+unlist(lapply(all_indivs, function(X) X[[j]]))))
	}
as.data.frame(temp_dt)-> temp_dt
temp_dt[which.max(temp_dt$part_r2),]
saveRDS(temp_dt, file='~/height_prediction/gwas/WHI/output/temp_dt_LA_all.Rds')
ggplot(temp_dt, aes(x=alpha, y=part_r2)) + geom_point(size=0.8)+ geom_line()
ggsave('~/height_prediction/gwas/WHI/figs/indiv_LA_alfa_part_R2_WHI.pdf')
temp_dt<-readRDS('~/height_prediction/gwas/WHI/output/temp_dt_LA_all.Rds')
###################
##################
a_vec2<-seq(from=0, to=1, by=0.001)

all_indivs<-vector('list', nrow(final2))
names(all_indivs)<-final2[,SUBJID]
counter<-0
system.time(for(S in final2[,SUBJID]){
        wei_PRS_2<-vector('list', length(a_vec2))
        names(wei_PRS_2)<-a_vec2
        wei_PRS_2<-lapply(1:length(a_vec2), function(i)( final2[SUBJID==S]$AFR_ANC*a_vec2[i]*final2[SUBJID==S]$PRS_plink[[1]])+(final2[SUBJID==S]$AFR_ANC*(1-a_vec2[i])*final2[SUBJID==S]$PRS_EUR[[1]]))
        counter+1-> counter
        cat(counter, "\r")
        all_indivs[[S]]<-wei_PRS_2
        gc()
}
) #5240 seconds
temp_dt2<-matrix(ncol=2, nrow=length(a_vec2))
colnames(temp_dt2)<-c('alpha', 'part_r2')
temp_dt2[,1]<-a_vec2
for (j in 1:length(a_vec2)){
        temp_dt2[j,2]<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(final2$HEIGHTX~final2$AGE+final2$AGE2+final2$EUR_ANC+unlist(lapply(all_indivs, function(X) X[[j]]))))
        }
as.data.frame(temp_dt2)-> temp_dt2
temp_dt2[which.max(temp_dt2$part_r2),]
saveRDS(temp_dt2, file='~/height_prediction/gwas/WHI/output/temp_dt_plink.Rds')
ggplot(temp_dt2, aes(x=alpha, y=part_r2)) + geom_point(size=0.8)+ geom_line()
ggsave('~/height_prediction/gwas/WHI/figs/indiv_plink_alfa_part_R2_WHI.pdf')
temp_dt2<-readRDS('~/height_prediction/gwas/WHI/output/temp_dt_plink.Rds')

#The End

#playground

test<-cbind(wanna_plot,temp_dt)
test4<-cbind(wanna_plot2, temp_dt2)
test[,alfa:=NULL]
test4[,alfa:-NULL]

colnames(test)[1]<-"Without Ancestry"
colnames(test)[3]<-"With Ancestry"

colnames(test4)[1]<-"Without Ancestry"
colnames(test4)[3]<-"With Ancestry"

melt(test, id='alpha')-> test2
melt(test4, id='alpha')-> test5
colnames(test2)[3]<-'Part_R2'
colnames(test5)[3]<-'Part_R2'

ggplot(test2, aes(x=alpha, y=Part_R2, colour=variable)) + geom_point(size=0.8)+ geom_line()
ggsave('~/height_prediction/gwas/WHI/figs/test_this_all.pdf')

ggplot(test5, aes(x=alpha, y=Part_R2, colour=variable)) + geom_point(size=0.8)+ geom_line()
ggsave('~/height_prediction/gwas/WHI/figs/test_this_plink.pdf')


