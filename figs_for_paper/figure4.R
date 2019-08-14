library(data.table)
library(reshape)
library(ggplot2)
library(dplyr)
library(cowplot)

source('~/height_prediction/strat_prs/scripts/Rsq_R2.R')
##panel figure

#R2 vs ancestry for HRS_eur, HRS_afr, UKB_afr imputed data





##
#barplot eur_afr, imp, non-imp, prunings
comp<-readRDS('~/height_prediction/imputed/output/comparison.Rds')

dt<-comp[Name %in% c("phys_500000_0.000005", "phys_100000_0.0005")][, .(Name, HRS_afr_imp, HRS_afr, HRS_eur_imp, HRS_eur)]
melt(dt)-> dt2
dt2$Name<-rep(c("100 Kb p<0.005","500 Kb p<0.000005"),4)
dt2[, Dataset:=c(rep('HRS_afr',4), rep('HRS_eur',4))]
dt2[,Method:=rep(c(rep("Imputed",2), rep("Non-imputed",2)),2)]

ggplot(dt2, aes(x=Name, y=value, fill=Method)) + facet_wrap(~Dataset) + geom_bar(stat="identity",position=position_dodge(), alpha=0.7) + scale_fill_brewer(name="Dataset", type="div", palette='Dark2') + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position='bottom', legend.title=element_blank()) + labs(y=expression(Partial~R^2)) + scale_alpha_manual(values=c(0.1, 0.5))
ggsave('~/height_prediction/figs_for_paper/test.pdf')


comp_ukb<-readRDS('~/height_prediction/imputed/output/comparison_ukb.Rds')
dt_ukb<-comp_ukb[Name %in% c("phys_500000_0.000005", "phys_100000_0.0005")][, .(Name, UKB_afr_imp, UKB_afr)]
melt(dt_ukb)-> dt2_ukb

dt2_ukb$Name<-rep(c("100 Kb p<0.005","500 Kb p<0.000005"),2)
dt2_ukb[, Dataset:=rep('UKB_afr',4)]
dt2_ukb[,Method:=c(rep("Imputed",2), rep("Non-imputed",2))]


rbind(dt2, dt2_ukb)

rbind(dt2_ukb, dt2)-> dt_final
dt_final$Dataset<-factor(dt_final$Dataset, levels=c("UKB_afr", "HRS_afr", "HRS_eur"))

ggplot(dt_final, aes(x=Name, y=value, fill=Method)) + facet_wrap(~Dataset) + geom_bar(stat="identity",position=position_dodge(), alpha=0.7) + scale_fill_brewer(name="Dataset", type="div", palette='Dark2') + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position='bottom', legend.title=element_blank(), axis.text.x = element_text(angle = 45)) + labs(y=expression(Partial~R^2)) + scale_alpha_manual(values=c(0.1, 0.5))
ggsave('~/height_prediction/figs_for_paper/Fig4_B.pdf')


#now do R2 vs anc
##############################################################
#combine all

##
readRDS(paste0('~/height_prediction/imputed/output/results.UKB_afr.Rds'))-> results.UKB_afr
readRDS(paste0('~/height_prediction/imputed//output/results.HRS_eur.Rds'))-> results.HRS_eur
readRDS(paste0('~/height_prediction/imputed/output/results.HRS_afr.Rds'))-> results.HRS_afr
readRDS(paste0('~/height_prediction/imputed/output/B_UKB_afr.Rds'))-> B_UKB_afr
readRDS(paste0('~/height_prediction/imputed/output/B_HRS_eur.Rds'))-> B_HRS_eur
readRDS(paste0('~/height_prediction/imputed/output/B_HRS_afr.Rds'))-> B_HRS_afr
readRDS(paste0('~/height_prediction/imputed/output/PGS3_UKB_afr.Rds'))-> PGS3_UKB_afr
readRDS(paste0('~/height_prediction/imputed/output/PGS3_HRS_eur.Rds'))-> PGS3_HRS_eur
readRDS(paste0('~/height_prediction/imputed/output/PGS3_HRS_afr.Rds'))-> PGS3_HRS_afr

I<-which(names(B_UKB_afr)=="100000_0.0005")
ALL2<-rbind(B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_HRS_eur[[I]])
tmp<-1/c(var(results.UKB_afr[[I]][[1]]$t),var(results.UKB_afr[[I]][[2]]$t), var(results.UKB_afr[[I]][[3]]$t), var(results.UKB_afr[[I]][[4]]$t), var(results.HRS_afr[[I]][[1]]$t), var(results.HRS_afr[[I]][[2]]$t),var(results.HRS_eur[[I]]$t))
cbind(ALL2, W=tmp)-> ALL2
my_plotA<-ggplot(ALL2, aes(x=Med_Eur_Anc, y=R_sq,colour=Dataset)) +
        geom_point(aes(shape=Dataset), size=1.5, fill="white") + stat_smooth(data=ALL2,method = "lm", mapping = aes(weight = W), col='black') +
        geom_errorbar(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset,ymin=boots_perc_L, ymax=boots_perc_U), width=0.05, size=0.8) +
        geom_errorbarh(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset, xmin=HVB_L, xmax=HVB_U), width=0.05, size=0.5) +  scale_color_brewer(palette="Dark2")+
        ylab(expression(paste("Partial R"^"2"))) + xlab("European Ancestry Proportion") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),axis.text.x=element_text(size=9), axis.text.y=element_text(size=9), legend.key=element_blank(), legend.background=element_blank(), legend.title=element_blank(), legend.position="bottom")


my_plotB<-ggplot(dt_final, aes(x=Name, y=value, fill=Method)) + facet_wrap(~Dataset, ncol = 1) + geom_bar(stat="identity",position=position_dodge(), alpha=0.7) + scale_fill_brewer(name="Dataset", type="div", palette='Dark2') + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position='bottom', legend.title=element_blank(), axis.title.x=element_blank()) + labs(y=expression(Partial~R^2), xlab="") + scale_alpha_manual(values=c(0.1, 0.5))


png('~/height_prediction/figs_for_paper/Figure4.png',width = 900, height = 500, units = "px")
plot_grid(my_plotA, my_plotB, labels=c("A","B"), nrow=1)
dev.off()

ggsave('~/height_prediction/figs_for_paper/Figure4.pdf')



