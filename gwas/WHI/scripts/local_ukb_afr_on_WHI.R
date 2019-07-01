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
########################
cat('checkpoint number 1\n')

short_fun<-function(args=22, data=prun2){
        what <- paste0("~/height_prediction/input/ukb_afr/UKB_kgCY_chr", args)
        betas<-fread(paste0('~/height_prediction/gwas/ukb_afr/output/AS_Beta_chr', args,'example.txt'))
        prun<-setDT(readRDS('~/height_prediction/gwas/WHI/output/hei_phys_100000_0.0005_v2.Rds')[[args]])
        fread(paste0(what, '.phsnp'))-> snps
        colnames(snps)<-c("ID", "CHR", "V3", "POS","REF","ALT")
        tag='phys_100000_0.0005'
        cbind(betas, snps)-> betas2
        data<-betas2[which(betas2[, POS] %in% prun[, POS]),]
        prun[POS %in% betas2[, POS]]-> prun2
	prun[POS %in% betas2[abs(Tstat_all)>1][, POS]]-> prun_tstat_all_1
	#prun[POS %in% betas2[abs(Tstat_all)>2][, POS]]-> prun_tstat_all_2
        #data$POS==prun2$POS
        prun2[, POP1:=data$POP1]
        prun2[, POP2:=data$POP2]
        prun2[, ALL:=data$ALL]
	prun_tstat_all_1[, POP1:=data$POP1]
        prun_tstat_all_1[, POP2:=data$POP2]
        prun2_tstat_all_1[, ALL:=data$ALL]
	
	final_plink[CHR==args]-> plink
	setkey(plink, CHR, POS)
	setkey(prun2, CHR,POS)
	setkey(prun_tstat_all_1, CHR,POS)
	prun2[plink, nomatch=0]-> prun2
	prun_tstat_all_1[plink, nomatch=0]-> prun_tstat_all_1
	prun2[plink, nomatch=0][abs(T_STAT)>1]-> prun2_TSTAT_plink_1
	#prun2[plink, nomatch=0][abs(T_STAT)>2]-> prun2_TSTAT_plink_2
	prun2[,i.MarkerName.1:=NULL][,UNADJ:=NULL][,GC:=NULL][,BONF:=NULL][,HOLM:=NULL][,SIDAK_SS:=NULL][,SIDAK_SD:=NULL][,FDR_BH:=NULL][,FDR_BY:=NULL][,i.REF:=NULL][,i.ALT:=NULL][,TEST:=NULL][,OBS_CT:=NULL][,i.SE:=NULL][,T_STAT:=NULL][,i.MarkerName:=NULL]
	prun2_TSTAT_plink_1[,i.MarkerName.1:=NULL][,UNADJ:=NULL][,GC:=NULL][,BONF:=NULL][,HOLM:=NULL][,SIDAK_SS:=NULL][,SIDAK_SD:=NULL][,FDR_BH:=NULL][,FDR_BY:=NULL][,i.REF:=NULL][,i.ALT:=NULL][,TEST:=NULL][,OBS_CT:=NULL][,i.SE:=NULL][,T_STAT:=NULL][,i.MarkerName:=NULL]
	prun_tstat_all_1[,i.MarkerName.1:=NULL][,UNADJ:=NULL][,GC:=NULL][,BONF:=NULL][,HOLM:=NULL][,SIDAK_SS:=NULL][,SIDAK_SD:=NULL][,FDR_BH:=NULL][,FDR_BY:=NULL][,i.REF:=NULL][,i.ALT:=NULL][,TEST:=NULL][,OBS_CT:=NULL][,i.SE:=NULL][,T_STAT:=NULL][,i.MarkerName:=NULL]
        test<-PolScore3(beta="POP1",data=prun2)	
        test2<-PolScore3(beta="POP2",data=prun2)
        test3<-PolScore3(beta='ALL',data=prun2)
	test3b<-PolScore3(beta='ALL',data=prun_tstat_all_1)
	test4<-PolScore3(beta='PLINK',data=prun2) #plink
	test4b<-PolScore3(beta='PLINK',data=prun2_TSTAT_plink_1)
	return(list(test,test2,test3,test3b,test4, test4b))
}

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


#saveRDS(res, file="~/height_prediction/gwas/WHI/output/all_prs.Rds")

#saveRDS(res, file="~/height_prediction/gwas/WHI/output/all_prs_v2.Rds")
cat('checkpoint number 3\n')
data.table(SUBJID=names(res_all[[1]][[1]]), 
PRS_POP1=(unlist(res_all[[1]][[1]])+unlist(res_all[[2]][[1]])+unlist(res_all[[3]][[1]])+unlist(res_all[[4]][[1]])+unlist(res_all[[5]][[1]])+unlist(res_all[[6]][[1]])+unlist(res_all[[7]][[1]])+unlist(res_all[[8]][[1]])+unlist(res_all[[9]][[1]])+unlist(res_all[[10]][[1]])+unlist(res_all[[11]][[1]])+unlist(res_all[[12]][[1]])+unlist(res_all[[13]][[1]])+unlist(res_all[[14]][[1]])+unlist(res_all[[15]][[1]])+unlist(res_all[[16]][[1]])+unlist(res_all[[17]][[1]])+ unlist(res_all[[18]][[1]])+unlist(res_all[[19]][[1]])+unlist(res_all[[20]][[1]])+unlist(res_all[[21]][[1]])+unlist(res_all[[22]][[1]])), 
PRS_POP2=(unlist(res_all[[1]][[2]])+unlist(res_all[[2]][[2]])+unlist(res_all[[3]][[2]])+unlist(res_all[[4]][[2]])+unlist(res_all[[5]][[2]])+unlist(res_all[[6]][[2]])+unlist(res_all[[7]][[2]])+unlist(res_all[[8]][[2]])+unlist(res_all[[9]][[2]])+unlist(res_all[[10]][[2]])+unlist(res_all[[11]][[2]])+unlist(res_all[[12]][[2]])+unlist(res_all[[13]][[2]])+unlist(res_all[[14]][[2]])+unlist(res_all[[15]][[2]])+unlist(res_all[[16]][[2]])+unlist(res_all[[17]][[2]])+ unlist(res_all[[18]][[2]])+unlist(res_all[[19]][[2]])+unlist(res_all[[20]][[2]])+unlist(res_all[[21]][[2]])+unlist(res_all[[22]][[2]])), 
PRS_all=(unlist(res_all[[1]][[3]])+unlist(res_all[[2]][[3]])+unlist(res_all[[3]][[3]])+unlist(res_all[[4]][[3]])+unlist(res_all[[5]][[3]])+unlist(res_all[[6]][[3]])+unlist(res_all[[7]][[3]])+unlist(res_all[[8]][[3]])+unlist(res_all[[9]][[3]])+unlist(res_all[[10]][[3]])+unlist(res_all[[11]][[3]])+unlist(res_all[[12]][[3]])+unlist(res_all[[13]][[3]])+unlist(res_all[[14]][[3]])+unlist(res_all[[15]][[3]])+unlist(res_all[[16]][[3]])+unlist(res_all[[17]][[3]])+ unlist(res_all[[18]][[3]])+unlist(res_all[[19]][[3]])+unlist(res_all[[20]][[3]])+unlist(res_all[[21]][[3]])+unlist(res_all[[22]][[3]])), 
PRS_plink=(unlist(res_all[[1]][[4]])+unlist(res_all[[2]][[4]])+unlist(res_all[[3]][[4]])+unlist(res_all[[4]][[4]])+unlist(res_all[[5]][[4]])+unlist(res_all[[6]][[4]])+unlist(res_all[[7]][[4]])+unlist(res_all[[8]][[4]])+unlist(res_all[[9]][[4]])+unlist(res_all[[10]][[4]])+unlist(res_all[[11]][[4]])+unlist(res_all[[12]][[4]])+unlist(res_all[[13]][[4]])+unlist(res_all[[14]][[4]])+unlist(res_all[[15]][[4]])+unlist(res_all[[16]][[4]])+unlist(res_all[[17]][[4]])+ unlist(res_all[[18]][[4]])+unlist(res_all[[19]][[4]])+unlist(res_all[[20]][[4]])+unlist(res_all[[21]][[4]])+unlist(res_all[[22]][[4]])),
PRS_EUR=unlist(readRDS('~/height_prediction/gwas/WHI/output/PGS_WHI_phys_100000_0.0005.Rds')))-> a

a[, PRS_POP1:=scale(PRS_POP1)]
a[, PRS_POP2:=scale(PRS_POP2)]
a[, PRS_all:=scale(PRS_all)]
a[, PRS_EUR:=scale(PRS_EUR)]
a[, PRS_plink:=scale(PRS_plink)]

fread('~/height_prediction/input/WHI/WHI_phenotypes.txt')-> Pheno_WHI

Pheno_WHI[, SUBJID:=paste0("0_", as.character(Pheno_WHI[, SUBJID]))]
setkey(Pheno_WHI, SUBJID)
#add ancestry
#admixture, obsolete

ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/WHI/rfmix_anc_chr', X, '.txt'))))

anc_WHI<-ancestry %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual

setkey(anc_WHI, SUBJID)

##
setkey(a, SUBJID)
setkey(Pheno_WHI, SUBJID)

a[Pheno_WHI][anc_WHI]-> final

final[,c("SUBJID","PRS_POP1","PRS_POP2","PRS_all", "PRS_plink","PRS_EUR", "SEX", "HEIGHTX", "EUR_ANC","AFR_ANC","AGE")]-> final2

#ok, so for ukb_afr and WHI (both), POP1 is AFR
cat('another checkpoint\n')

final2[,AGE2:=AGE^2][, PRS_eur_afr_local:=(mean(EUR_ANC)*(PRS_EUR))+(mean(AFR_ANC)*PRS_all)]
final2[,AGE2:=AGE^2][, PRS_eur_afr_plink:=(mean(EUR_ANC)*(PRS_EUR))+(mean(AFR_ANC)*PRS_plink)]
fwrite(final2, file="~/height_prediction/gwas/WHI/output/all_prs_whi.txt", sep="\t")

partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2))*100 #4.1%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_all, data=final2))*100 #0.35%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_plink, data=final2))*100 #0.39%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_POP1, data=final2))*100 #0.0019 %
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_POP2, data=final2))*100 #0.04%

partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_eur_afr_local, data=final2))*100 #0.5%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_eur_afr_plink, data=final2))*100 # #0.54


#read in plink betas

combo_prs<-readRDS('~/height_prediction/gwas/ukb_afr/output/combo_prs.Rds')


a_vec<-seq(0,1,0.001)
wei_PRS<-vector('list', length(a_vec))
names(wei_PRS)<-a_vec
part_r2<-c()
for(i in 1:length(a_vec)){

cat(i, '\n')
wei_PRS[[i]]<-(a_vec[i]*final2$PRS_EUR)+((1-a_vec[i])*final2$PRS_all)
part_r2[i]<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(final2$HEIGHTX~final2$AGE+final2$AGE2+final2$EUR_ANC+wei_PRS[[i]]))
}

data.table(part_R2=part_r2, alfa=a_vec)-> wanna_plot

ggplot(wanna_plot, aes(x=alfa, y=part_R2)) + geom_point(size=0.8)+ geom_line() + geom_vline(xintercept=wanna_plot[which.max(wanna_plot$part_R2),alfa], col='red') + geom_hline(yintercept=partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2)), col='red') + coord_cartesian(ylim = c(0, 0.048))
ggsave('~/height_prediction/gwas/WHI/figs/alfa_part_R2_WHI_all.pdf')
##

a_vec<-seq(0,1,0.001)
wei_PRS<-vector('list', length(a_vec))
names(wei_PRS)<-a_vec
part_r2<-c()
for(i in 1:length(a_vec)){

cat(i, '\n')
wei_PRS[[i]]<-(a_vec[i]*final2$PRS_EUR)+((1-a_vec[i])*final2$PRS_plink)
part_r2[i]<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(final2$HEIGHTX~final2$AGE+final2$AGE2+final2$EUR_ANC+wei_PRS[[i]]))
}

data.table(part_R2=part_r2, alfa=a_vec)-> wanna_plot

ggplot(wanna_plot, aes(x=alfa, y=part_R2)) + geom_point(size=0.8)+ geom_line() + geom_vline(xintercept=wanna_plot[which.max(wanna_plot$part_R2),alfa], col='red') + geom_hline(yintercept=partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2)), col='red')

ggsave('~/height_prediction/gwas/WHI/figs/alfa_part_R2_WHI_plink.pdf')

a_vec2<-seq(from=0, to=1, by=0.001)

all_indivs<-vector('list', nrow(final2))
names(all_indivs)<-final2[,SUBJID]
counter<-0
for(S in final2[,SUBJID]){
	wei_PRS_2<-vector('list', length(a_vec2))
	names(wei_PRS)<-a_vec
	part_r2_2<-c()
	for(i in 1:length(a_vec2)){
		cat(i, '\r')
		wei_PRS_2[[i]]<-(final2[SUBJID==S]$AFR_ANC*a_vec2[i]*final2[SUBJID==S]$PRS_all[[1]])+(final2[SUBJID==S]$EUR_ANC*final2[SUBJID==S]$PRS_EUR[[1]])
	}
	all_indivs[[S]]<-wei_PRS_2
	counter+1-> counter
	cat(counter, "\n")
}


	temp_dt<-matrixc(ncol=2, nrow=length(a_vec2))
	for (j in 1:length(a_vec2)){
		part_r2_2[j,2]<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(final2$HEIGHTX~final2$AGE+final2$AGE2+final2$EUR_ANC+unlist(lapply(all_indivs, function(X) X[[j]]))))
	}
	setDT(temp_dt)
	temp_dt[which.max(temp_dt$part_R2),]

ggplot(temp_dt, aes(x=alpha, y=part_R2)) + geom_point(size=0.8)+ geom_line()

ggsave('~/height_prediction/gwas/WHI/figs/indiv_LA_alfa_part_R2_WHI.pdf')


