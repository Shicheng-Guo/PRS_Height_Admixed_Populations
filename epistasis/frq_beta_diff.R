#!/usr/bin/env Rscript
############################
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
}
library(data.table)
library(dplyr)
library(ggplot2)
library(EnvStats)
options(scipen=999)
W<-as.numeric(args[1])#window size
#read in frequencies
readRDS('~/height_prediction/epistasis/output/table_HRS_eur.Rds')->hrs_eur
readRDS('~/height_prediction/epistasis/output/table_HRS_afr.Rds')->hrs_afr
readRDS('~/height_prediction/epistasis/output/table_ukb_eur.Rds')->ukb_eur
readRDS('~/height_prediction/epistasis/output/table_ukb_afr.Rds')->ukb_afr
#calculate frequency difference per SNP

setkey(hrs_eur, CHR, POS)
setkey(hrs_afr, CHR, POS)
setkey(ukb_afr, CHR, POS)
setkey(ukb_eur, CHR, POS)
comb<-hrs_eur[hrs_afr, nomatch=0]
comb_ukb<-ukb_eur[ukb_afr, nomatch=0]
comb$Freq_effect_allele<-as.numeric(comb$Freq_effect_allele)
comb$i.Freq_effect_allele<-as.numeric(comb$i.Freq_effect_allele)
comb[,Freq_Diff:=Freq_effect_allele-i.Freq_effect_allele]
comb_ukb$Freq_effect_allele<-as.numeric(comb_ukb$Freq_effect_allele)
comb_ukb$i.Freq_effect_allele<-as.numeric(comb_ukb$i.Freq_effect_allele)
comb_ukb[,Freq_Diff:=(Freq_effect_allele-i.Freq_effect_allele)]
comb[, POS1:=POS]
comb[, POS2:=POS]
comb_ukb[, POS1:=POS]
comb_ukb[, POS2:=POS]
#betas_ukb<-readRDS('~/height_prediction/gwas/ukb_afr/output/betas_phys_100000_0.0005_20000.Rds')
final_plink<-readRDS('~/height_prediction/loc_anc_analysis/output/final_plink.Rds')
select(do.call(rbind, readRDS('~/height_prediction/gwas/ukb_afr/output/hei_phys_100000_0.0005_v2.Rds')), CHR, POS, b, SE, p, N)-> prs_snps
select(do.call(rbind, readRDS('~/height_prediction/gwas/HRS_afr/output/hei_phys_100000_0.0005_v2.Rds')), CHR, POS, b, SE, p, N)-> prs_snps_hrs
#coords<-betas_ukb[,.(CHR,POS)][, POS1:=POS-(W/2)][,POS2:=POS+(W/2)]
coords<-prs_snps[,.(CHR,POS)][, POS1:=POS-(W/2)][,POS2:=POS+(W/2)]
coords<-split(coords, by='CHR')
lapply(coords, function(X) setkey(X, CHR, POS1, POS2))

coords_hrs<-prs_snps_hrs[,.(CHR,POS)][, POS1:=POS-(W/2)][,POS2:=POS+(W/2)]
coords_hrs<-split(coords_hrs, by='CHR')
lapply(coords_hrs, function(X) setkey(X, CHR, POS1, POS2))

RES<-vector('list', 22)
RES_UKB<-vector('list', 22)

for(chr in 22:1){
	tp<-comb[CHR==chr]
	setkey(tp, CHR, POS1, POS2)
	coords_hrs[[chr]]$CHR<-as.integer(coords_hrs[[chr]]$CHR)
	setkey(coords_hrs[[chr]], CHR,POS1, POS2)
	foverlaps(tp, coords_hrs[[chr]], type='within', nomatch=0)[, Win:=paste0(POS1, "|",POS2)]-> res
	res[,MeanFreqDiff:=mean(Freq_Diff^2, na.rm=T), by="Win"]
	RES[[chr]]<-res
}
for(chr in 22:1){
	tp<-comb_ukb[CHR==chr]
	setkey(tp, CHR, POS1, POS2)
        coords[[chr]]$CHR<-as.integer(coords[[chr]]$CHR)
        setkey(coords[[chr]], CHR,POS1, POS2)
        foverlaps(tp, coords[[chr]], type='within', nomatch=0)[, Win:=paste0(POS1, "|",POS2)]-> res
        res[,MeanFreqDiff:=mean(Freq_Diff^2, na.rm=T), by="Win"]
	RES_UKB[[chr]]<-res
	cat(chr,'\r')
}
do.call(rbind, RES)-> RES2
do.call(rbind, RES_UKB)-> RES2_UKB
RES2[, Data:='HRS']
RES2_UKB[, Data:='UKB']

RES2<-unique(RES2, by=c('CHR','Win'))

#calculate beta diff
merge(final_plink, prs_snps_hrs, by=c('CHR', 'POS'))-> betas_plink_hrs
betas_plink_hrs[, Beta_Diff:=b-PLINK]
betas_plink_hrs[, Beta_Diff_Chisq:=(Beta_Diff/sqrt(((SE^2)+(SE_plink^2))))^2]
betas_plink_hrs$CHR<-as.integer(betas_plink_hrs$CHR)
setkey(betas_plink_hrs,CHR, POS)
setkey(RES2, CHR, POS)
dt2<-betas_plink_hrs[RES2, nomatch=0]
#check correlation between freq_diff and beta_diff
RES2_UKB<-unique(RES2_UKB, by=c('CHR','Win'))
#add betas (POP 1 and UKB_eur)
#calculate beta diff
merge(final_plink, prs_snps, by=c('CHR', 'POS'))-> betas_plink
betas_plink[, Beta_Diff:=b-PLINK]
betas_plink[, Beta_Diff_Chisq:=(Beta_Diff/sqrt(((SE^2)+(SE_plink^2))))^2]
betas_plink$CHR<-as.integer(betas_plink$CHR)
setkey(betas_plink,CHR, POS)
setkey(RES2_UKB, CHR, POS)
dt2_UKB<-betas_plink[RES2_UKB, nomatch=0]

rbind(dt2, dt2_UKB)-> my_dt

ggplot(my_dt, aes(x=MeanFreqDiff, y=sqrt(Beta_Diff_Chisq))) + geom_point(cex=0.5, col="lightgray") + geom_smooth(method='lm', se=F, col='black') + labs(x="Mean Squared Frequency Difference", y=expression(chi^2), cex=15) + facet_wrap(~Data, nrow=2) + theme_bw()
ggsave(paste0('~/height_prediction/epistasis/figs/epistasis_', W, '.pdf'))

#p-values
summary(lm(Beta_Diff_Chisq~MeanFreqDiff,data=dt2_UKB))
summary(lm(Beta_Diff_Chisq~MeanFreqDiff,data=dt2))

attr(dt2_UKB$Beta_Diff_Chisq,"scaled:scale")<-NULL
attr(dt2_UKB$Beta_Diff_Chisq,"scaled:center")<-NULL
model_ukb<-lm(Beta_Diff_Chisq~MeanFreqDiff,data=dt2_UKB)
attr(dt2$Beta_Diff_Chisq,"scaled:scale")<-NULL
attr(dt2$Beta_Diff_Chisq,"scaled:center")<-NULL
model_hrs<-lm(Beta_Diff_Chisq~MeanFreqDiff,data=dt2)


pdf(paste0('figs/boxcoxhrs_', W, '.pdf'))
plot(boxcox(model_hrs))
dev.off()

pdf(paste0('figs/boxcoxukb_', W, '.pdf'))
plot(boxcox(model_ukb))
dev.off()
