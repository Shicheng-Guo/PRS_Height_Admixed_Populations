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
library(cowplot)
options(scipen=999)

i<-'gibbs-inf'
readRDS(paste0('~/height_prediction/ldpred/output/results.WHI.Rds'))[[i]]-> results.WHI
readRDS(paste0('~/height_prediction/ldpred/output/results.JHS.Rds'))[[i]]-> results.JHS
readRDS(paste0('~/height_prediction/ldpred/output/results.UKB_afr.Rds'))[[22]]-> results.UKB_afr
readRDS(paste0('~/height_prediction/ldpred/output/results.HRS_eur.Rds'))[[i]]-> results.HRS_eur
readRDS(paste0('~/height_prediction/ldpred/output/results.HRS_afr.Rds'))[[i]]-> results.HRS_afr
readRDS(paste0('~/height_prediction/ldpred/output/B_WHI.Rds'))[[i]]->B_WHI
readRDS(paste0('~/height_prediction/ldpred/output/B_JHS.Rds'))[[i]]->B_JHS
readRDS(paste0('~/height_prediction/ldpred/output/B_UKB_afr_all.Rds'))[[i]]-> B_UKB_afr
readRDS(paste0('~/height_prediction/ldpred/output/B_HRS_eur.Rds'))[[i]]-> B_HRS_eur
readRDS(paste0('~/height_prediction/ldpred/output/B_HRS_afr.Rds'))[[i]]-> B_HRS_afr
readRDS(paste0('~/height_prediction/ldpred/output/PGS2_UKB_afr.Rds'))[[i]]-> PGS3_UKB_afr
readRDS(paste0('~/height_prediction/ldpred/output/PGS3_WHI.Rds'))[[i]]-> PGS3_WHI
readRDS(paste0('~/height_prediction/ldpred/output/PGS3_JHS.Rds'))[[i]]-> PGS3_JHS
readRDS(paste0('~/height_prediction/ldpred/output/PGS3_HRS_eur.Rds'))[[i]]-> PGS3_HRS_eur
readRDS(paste0('~/height_prediction/ldpred/output/PGS3_HRS_afr.Rds'))[[i]]-> PGS3_HRS_afr

ALL<-rbind(B_JHS[1:2,][, Dataset:='JHS_afr'], B_WHI[1:4,][, Dataset:='WHI_afr'], B_UKB_afr[1:4,][,Dataset:='UKB_afr'],B_HRS_afr[1:2,][, Dataset:='HRS_afr'], B_HRS_eur)

tmp<-1/c(var(results.JHS[[1]]$t), var(results.JHS[[2]]$t), var(results.WHI[[1]]$t), var(results.WHI[[2]]$t), var(results.WHI[[3]]$t), var(results.WHI[[4]]$t), var(results.UKB_afr[[1]]$t),var(results.UKB_afr[[2]]$t), var(results.UKB_afr[[3]]$t), var(results.UKB_afr[[4]]$t), var(results.HRS_afr[[1]]$t), var(results.HRS_afr[[2]]$t), var(results.HRS_eur$t))  #weighing lm by boostrap replicates.

ALL[,W:=tmp]
ALL$Dataset<-factor(ALL$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "HRS_eur"))
my_colrs<-c(brewer.pal(4, 'Set1'),"#101010")
shapes<-c(16,17,18,6,15)

A_plot<-ggplot(ALL, aes(x=Med_Eur_Anc, y=R_sq,colour=Dataset)) +
	geom_point(aes(shape=Dataset), size=5, fill="white", alpha=0.8) + stat_smooth(data=ALL,method = "lm", mapping = aes(weight = sqrt(W)), col='darkgray', lty=2, lwd=1) +
        geom_errorbar(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset,ymin=boots_perc_L, ymax=boots_perc_U), width=0.05, size=0.8) +
        #geom_line(color='lightgray')+
#        geom_errorbarh(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset, xmin=HVB_L, xmax=HVB_U), width=0.05, size=0.5) +
        scale_color_manual(values=my_colrs) + coord_cartesian(ylim = c(0, 0.16)) +
	scale_shape_manual(values=shapes)+
        ylab(expression(paste("Partial R"^"2"))) + xlab("European Ancestry Proportion") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=15), 
	axis.text.y=element_text(size=15), legend.key=element_blank(), legend.background=element_blank(), legend.title=element_blank(), legend.text=element_text(size=15), legend.position=c(0.10,0.85))


png('~/height_prediction/figs_for_paper/figs/ldpred-inf.png', res=300, unit="in", height=12, width=8)
print(A_plot)
dev.off()

