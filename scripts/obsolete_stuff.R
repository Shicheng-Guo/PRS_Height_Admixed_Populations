#!/usr/bin/env Rscript
library(data.table)
library(snpStats)
library(ggplot2)
library(dplyr)
library(parallel)
library(mgcViz)
library(SOAR)
library(grid)
#source('~/height_prediction/gwas/ukb_afr/scripts/compare_betas.R') #Run the analysis
##
combo<-fread('~/height_prediction/gwas/ukb_afr/output/plink_local_anc_effect_sizes.txt')
combo_prs<-fread('~/height_prediction/gwas/ukb_afr/output/combo_prs.txt')
pdf('~/height_prediction/figs_for_paper/figs/qqplot_test.pdf')
qq.chisq(combo$Beta_Diff_Chisq, df=2,  main="QQ plot ALL SNPs",  xlab="Expected", ylab="Observed")
qq.chisq(combo_prs$Beta_Diff_Chisq, df=2,  main="QQ plot PRS",  xlab="Expected", ylab="Observed")
dev.off()



##QQ plot for plink UKB

observed <- sort(combo$UNADJ)
lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

observed2 <- sort(combo$GC)
lobs2 <- -(log10(observed2))

expected2 <- c(1:length(observed2))
lexp2 <- -(log10(expected2 / (length(expected2)+1)))

df<-data.table(Obs=c(lobs, lobs2), Exp=c(lexp, lexp2), Class=c(rep("UNADJ", length(observed)),rep("GC", length(observed2))))

ggplot(df, aes(x=Exp, y=Obs, colour=Class)) + geom_point() +   geom_abline(intercept = 0, slope = 1, col="black") + ggtitle("Plink GWAS for UKB admixed inviduals")
ggsave('~/height_prediction/figs_for_paper/figs/qqplot_plink_UKB.pdf')

##manhattan for plink
library(qqman)
combo[, P_all:=2*(pt(abs(Tstat_all), df=N-13, lower.tail=F))]
combo[, P_pop1:=2*(pt(abs(Tstat_pop1), df=N-13, lower.tail=F))]
combo[, P_pop2:=2*(pt(abs(Tstat_pop2), df=N-13, lower.tail=F))]
lamb<- median(qchisq(combo$P_all, df=1), na.rm=T)/qchisq(0.5, df=1)
combo[, P_all_GC:= pchisq(qchisq(P_all, df=1)/lamb, df=1)]
lamb1<- median(qchisq(combo$P_pop1, df=1), na.rm=T)/qchisq(0.5, df=1)
combo[, P_pop1_GC:= pchisq(qchisq(P_pop1, df=1)/lamb1, df=1)]
lamb2<- median(qchisq(combo$P_pop2, df=1), na.rm=T)/qchisq(0.5, df=1)
combo[, P_pop2_GC:= pchisq(qchisq(P_pop2, df=1)/lamb2, df=1)]
saveRDS(combo_prs, file="~/height_prediction/gwas/ukb_afr/output/combo_prs.Rds")
saveRDS(combo, file="~/height_prediction/gwas/ukb_afr/output/combo_local.Rds")

df2<-combo[,.(CHR,POS,MarkerName,UNADJ)]
df2[,BP:=POS][,POS:=NULL]
df2[,P:=UNADJ][,UNADJ:=NULL]
df2[,SNP:=MarkerName][,MarkerName:=NULL]
df2[order(CHR, BP)]-> df2

png('~/height_prediction/figs_for_paper/figs/manhattan_plink_UKB.png')
manhattan(df2, genomewideline=F, suggestiveline=F)
dev.off()


df2<-combo[,.(CHR,POS,MarkerName,P_all)]
df2[,P:=P_all][,P_all:=NULL]
df2[,SNP:=MarkerName][,MarkerName:=NULL]
df2[,BP:=POS][,POS:=NULL]
df2[order(CHR, BP)]-> df2

png('~/height_prediction/figs_for_paper/figs/manhattan_LA_UKB.png')
manhattan(df2, genomewideline=F, suggestiveline=F)
dev.off()


#

#QQ plot for LA UKB

observed <- sort(combo$P_all)
lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

observed_GC <- sort(combo$P_all_GC)
lobs_GC <- -(log10(observed_GC))

expected_GC <- c(1:length(observed_GC))
lexp_GC <- -(log10(expected_GC / (length(expected_GC)+1)))

observed2 <- sort(combo$P_pop1)
lobs2 <- -(log10(observed2))

expected2 <- c(1:length(observed2))
lexp2 <- -(log10(expected2 / (length(expected2)+1)))

observed2_GC <- sort(combo$P_pop1_GC)
lobs2_GC <- -(log10(observed2_GC))

expected2_GC <- c(1:length(observed2_GC))
lexp2_GC <- -(log10(expected2_GC / (length(expected2_GC)+1)))

observed3 <- sort(combo$P_pop2)
lobs3 <- -(log10(observed3))

expected3 <- c(1:length(observed3))
lexp3 <- -(log10(expected3 / (length(expected3)+1)))

observed3_GC <- sort(combo$P_pop2_GC)
lobs3_GC<- -log10(observed3_GC)


expected3_GC <- c(1:length(observed3_GC))
lexp3_GC <- -(log10(expected3_GC / (length(expected3_GC)+1)))

df<-data.table(Obs=c(lobs, lobs_GC, lobs2,lobs2_GC, lobs3, lobs3_GC), Exp=c(lexp, lexp_GC,lexp2,lexp2_GC, lexp3, lexp3_GC), Class=c(rep("ALL", length(observed)), rep("ALL_GC", length(observed_GC)),rep("POP1", length(observed2)), rep("POP1_GC", length(observed2_GC)), rep('POP2',length(observed3)), rep("POP2_GC", length(observed3_GC))))

pdf('~/height_prediction/figs_for_paper/figs/qqplot_LA_UKB.pdf')
ggplot(df, aes(x=Exp, y=Obs, colour=Class)) + geom_point(cex=0.5) +   geom_abline(intercept = 0, slope = 1, col="black") + ggtitle("Local ancestry GWAS for UKB admixed inviduals")
dev.off()



