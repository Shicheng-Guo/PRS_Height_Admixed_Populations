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
library(ggpubr)
options(scipen=999)
W<-as.numeric(args[1])#window size
#read in frequencies
library(TeachingDemos)

txtStart(paste0('epistasis/outout_', W, '.txt'))
readRDS('output/table_HRS_eur.Rds')->hrs_eur
readRDS('output/table_HRS_afr.Rds')->hrs_afr
#calculate frequency difference per SNP

setkey(hrs_eur, CHR, POS)
setkey(hrs_afr, CHR, POS)
comb<-hrs_eur[hrs_afr, nomatch=0]
comb$Freq_effect_allele<-as.numeric(comb$Freq_effect_allele)
comb$i.Freq_effect_allele<-as.numeric(comb$i.Freq_effect_allele)
comb[,Freq_Diff:=Freq_effect_allele-i.Freq_effect_allele]
comb[, POS1:=POS]
comb[, POS2:=POS]
#betas_ukb<-readRDS('~/height_prediction/gwas/ukb_afr/output/betas_phys_100000_0.0005_20000.Rds')
final_plink<-readRDS('~/height_prediction/loc_anc_analysis/output/final_plink.Rds')
select(do.call(rbind, readRDS('~/height_prediction/gwas/HRS_afr/output/hei_phys_100000_0.0005_v2.Rds')), CHR, POS, b, SE, p, N)-> prs_snps_hrs
coords_hrs<-prs_snps_hrs[,.(CHR,POS)][, POS1:=POS-(W/2)][,POS2:=POS+(W/2)]
coords_hrs<-split(coords_hrs, by='CHR')
lapply(coords_hrs, function(X) setkey(X, CHR, POS1, POS2))

RES<-vector('list', 22)

for(chr in 22:1){
	tp<-comb[CHR==chr]
	setkey(tp, CHR, POS1, POS2)
	coords_hrs[[chr]]$CHR<-as.integer(coords_hrs[[chr]]$CHR)
	setkey(coords_hrs[[chr]], CHR,POS1, POS2)
	foverlaps(tp, coords_hrs[[chr]], type='within', nomatch=0)[, Win:=paste0(POS1, "|",POS2)]-> res
	res[,MeanFreqDiff:=mean(Freq_Diff^2, na.rm=T), by="Win"]
	RES[[chr]]<-res
}
do.call(rbind, RES)-> RES2
RES2[, Data:='HRS']

RES2<-unique(RES2, by=c('CHR','Win'))

#calculate beta diff
merge(final_plink, prs_snps_hrs, by=c('CHR', 'POS'))-> betas_plink_hrs
betas_plink_hrs[, Beta_Diff:=b-PLINK]
betas_plink_hrs[, Beta_Diff_Chisq:=(Beta_Diff/sqrt(((SE^2)+(SE_plink^2))))^2]
betas_plink_hrs$CHR<-as.integer(betas_plink_hrs$CHR)
setkey(betas_plink_hrs,CHR, POS)
setkey(RES2, CHR, POS)
dt2<-betas_plink_hrs[RES2, nomatch=0]
dt2[,Quantile:=cut(MeanFreqDiff, breaks=quantile(MeanFreqDiff, probs=seq(0,1, by=0.2), na.rm=T), include.lowest=T)]
dt2[, MeanBetaDiffChisq:=mean(Beta_Diff_Chisq, na.rm=T), by=Quantile]
dt2[, MedianMeanFreqDiff:=median(MeanFreqDiff, na.rm=T), by=Quantile]
saveRDS(dt2, paste0('~/height_prediction/epistasis/output/dt2_', W, '.Rds'))
#check correlation between freq_diff and beta_diff
#calculate beta diff
nrow(dt2)
my_dt<-dt2[Beta_Diff_Chisq<=15]
nrow(my_dt)
model_hrs<-lm(Beta_Diff_Chisq~MeanFreqDiff,data=dt2)
require(broom)
glance(model_hrs)
pval<-glance(model_hrs)$p.value
ggplot(my_dt, aes(x=MeanFreqDiff, y=Beta_Diff_Chisq)) + geom_point(cex=0.5, col="lightgray") + 
geom_smooth(method='lm', se=T, col='black') + 
geom_point(aes(x=MedianMeanFreqDiff, y=MeanBetaDiffChisq), col='red', cex=0.5) +
labs(x="Mean Squared Frequency Difference", y=expression(chi[diff]^2), cex=18) + 
annotate("text", x=0.3, y=0.15, label=paste("p=", round(pval,4))) +
#facet_wrap(~Data, nrow=2) + 
theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18))
ggsave(paste0('~/height_prediction/epistasis/figs/epistasis_', W, '.pdf'))

#boxcox
attr(dt2$Beta_Diff_Chisq,"scaled:scale")<-NULL
attr(dt2$Beta_Diff_Chisq,"scaled:center")<-NULL
model_hrs<-lm(Beta_Diff_Chisq~MeanFreqDiff,data=dt2)
summary(model_hrs)


#pdf(paste0('figs/boxcoxhrs_', W, '.pdf'))
#plot(boxcox(model_hrs))
#dev.off()

txtStop()
