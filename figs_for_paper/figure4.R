library(data.table)
library(reshape)
library(ggplot2)
library(dplyr)

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

ggplot(dt2, aes(x=Name, y=value, fill=Method)) + facet_wrap(~Dataset) + geom_bar(stat="identity",position=position_dodge(), alpha=0.7) + scale_fill_brewer(name="Dataset", type="div", palette='Dark2') + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),gend.position='bottom', legend.title=element_blank()) + labs(y=expression(Partial~R^2)) + scale_alpha_manual(values=c(0.1, 0.5))
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
ggsave('~/height_prediction/figs_for_paper/Fig4_A.pdf')


#now do R2 vs anc




