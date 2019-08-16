#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
}
library(ggplot2)
library(reshape2)
library(data.table)
library(reshape)
library(RColorBrewer)
library(MASS)
library(cowplot)

dtset<-args[1]
source('~/height_prediction/strat_prs/scripts/fancy_scientific.R')
a_list<-vector('list', 2)
names(a_list)<-c("AA","CEU")

a_list[['AA']]<-vector('list', 5)
a_list[['CEU']]<-vector('list', 5)
names(a_list[['AA']])<-c("HRS_eur","JHS","WHI", "ukb_afr", "HRS_afr")
names(a_list[['CEU']])<-c("HRS_eur","JHS","WHI","ukb_afr", "HRS_afr")

for(I in c("WHI","JHS","ukb_afr","HRS_eur", "HRS_afr")){
	a_list[['AA']][[I]]<-vector('list', 3)
	a_list[['CEU']][[I]]<-vector('list', 3)
	names(a_list[['AA']][[I]])<-c("5000","10000","20000")
 	names(a_list[['CEU']][[I]])<-c("5000","10000","20000")
	for(J in names(a_list[['AA']][[I]])){	
		a_list[['AA']][[I]][[J]]<-readRDS(paste0('~/height_prediction/strat_prs/output/part_R2_', I,"_",dtset, "_", J, "_AA.Rds"))
		a_list[['CEU']][[I]][[J]]<-readRDS(paste0('~/height_prediction/strat_prs/output/part_R2_', I, "_", dtset,"_", J, "_CEU.Rds"))
	}
}

if(dtset=='gwas'){
	r2_vec<-c(0.1210222,0.03587304,0.04100905,0.03776804,0.02376103)
} else if (dtset=='sib_betas'){
	r2_vec<-c(0.07511196,0.00967814,0.01993696,0.02717107,0.01060191)
}
df1<-data.table(Quantile=rep(c("q1","q2","q3","q4"),3),
Win=c(rep(5000,4), rep(10000,4), rep(20000,4)),
HRS_eur=unlist(a_list[['AA']])[1:12]/r2_vec[1], #updated
JHS_afr=unlist(a_list[['AA']])[13:24]/r2_vec[2], #updated
WHI_afr=unlist(a_list[['AA']])[25:36]/r2_vec[3], #updated
UKB_afr=unlist(a_list[['AA']])[37:48]/r2_vec[4], #updated
HRS_afr=unlist(a_list[['AA']])[49:60]/r2_vec[5] #updaed
)

df3<-data.table(Quantile=rep(c("q1","q2","q3","q4"),3),
Win=c(rep(5000,4), rep(10000,4), rep(20000,4)),
HRS_eur=unlist(a_list[['CEU']])[1:12]/r2_vec[1], #updated
JHS_afr=unlist(a_list[['CEU']])[13:24]/r2_vec[2], #updated
WHI_afr=unlist(a_list[['CEU']])[25:36]/r2_vec[3], #updated
UKB_afr=unlist(a_list[['CEU']])[37:48]/r2_vec[4], #updated
HRS_afr=unlist(a_list[['CEU']])[49:60]/r2_vec[5] #updaed
)


df1[, Map:='AA_Map']
df3[, Map:='CEU_Map']
df4<-rbind(df1,df3)


melt(df4, id=c("Quantile","Win", "Map"))-> df5

plot1<-ggplot(df5[Win==10000], aes(x=Quantile, y=value, fill=variable)) + geom_bar(stat='identity', position='dodge', alpha=0.7) + facet_grid(. ~Map) + labs(y=expression(paste("Relative partial R"^"2")),x="cM") +  scale_fill_brewer(name="Dataset", type="div", palette='Dark2') + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom", legend.title=element_blank()) + scale_x_discrete(labels=c("Low", expression(symbol('\256')), expression(symbol('\256')), "High"))

print(plot1)

ggsave(paste0('~/height_prediction/strat_prs/figs/barplot_AA_CEU_', dtset,'.pdf'))


#

df1<-data.table(Quantile=rep(c("q1","q2","q3","q4"),3),
Win=c(rep(5000,4), rep(10000,4), rep(20000,4)),
HRS_eur=unlist(a_list[['AA']])[1:12]/c(readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_HRS_eur_',dtset, '_5000_AA.Rds')),readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_HRS_eur_', dtset, '_10000_AA.Rds')),readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_HRS_eur_', dtset, '_20000_AA.Rds'))),
JHS_afr=unlist(a_list[['AA']])[13:24]/c(readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_JHS_', dtset, '_5000_AA.Rds')),readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_JHS_', dtset, '_10000_AA.Rds')),readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_JHS_', dtset,'_20000_AA.Rds'))),
WHI_afr=unlist(a_list[['AA']])[25:36]/c(readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_WHI_', dtset, '_5000_AA.Rds')),readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_WHI_', dtset, '_10000_AA.Rds')), readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_WHI_',dtset,'_20000_AA.Rds'))),
UKB_afr=unlist(a_list[['AA']])[37:48]/c(readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_ukb_afr_', dtset,'_5000_AA.Rds')),readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_ukb_afr_', dtset, '_10000_AA.Rds')),readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_ukb_afr_',dtset,'_20000_AA.Rds'))),
HRS_afr=unlist(a_list[['AA']])[49:60]/c(readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_HRS_afr_',dtset,'_5000_AA.Rds')),readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_HRS_afr_', dtset, '_10000_AA.Rds')),readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_HRS_afr_', dtset,'_20000_AA.Rds')))
)

df3<-data.table(Quantile=rep(c("q1","q2","q3","q4"),3),
Win=c(rep(5000,4), rep(10000,4), rep(20000,4)),
HRS_eur=unlist(a_list[['CEU']])[1:12]/c(readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_HRS_eur_',dtset, '_5000_AA.Rds')),readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_HRS_eur_', dtset, '_10000_AA.Rds')),readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_HRS_eur_', dtset, '_20000_AA.Rds'))),
JHS_afr=unlist(a_list[['CEU']])[13:24]/c(readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_JHS_', dtset, '_5000_AA.Rds')),readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_JHS_', dtset, '_10000_AA.Rds')),readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_JHS_', dtset,'_20000_AA.Rds'))),
WHI_afr=unlist(a_list[['CEU']])[25:36]/c(readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_WHI_', dtset, '_5000_AA.Rds')),readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_WHI_', dtset, '_10000_AA.Rds')), readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_WHI_',dtset,'_20000_AA.Rds'))),
UKB_afr=unlist(a_list[['CEU']])[37:48]/c(readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_ukb_afr_', dtset,'_5000_AA.Rds')),readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_ukb_afr_', dtset, '_10000_AA.Rds')),readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_ukb_afr_',dtset,'_20000_AA.Rds'))),
HRS_afr=unlist(a_list[['CEU']])[49:60]/c(readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_HRS_afr_',dtset,'_5000_AA.Rds')),readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_HRS_afr_', dtset, '_10000_AA.Rds')),readRDS(paste0('~/height_prediction/strat_prs/output/Nr_SNPs_HRS_afr_', dtset,'_20000_AA.Rds')))
)

df1[, Map:='AA_Map']
df3[, Map:='CEU_Map']
df4<-rbind(df1,df3)

melt(df4, id=c("Quantile","Win","Map"))-> df5


plot2<-ggplot(df5[Win==10000], aes(x=Quantile, y=value, fill=variable)) +  scale_y_continuous(labels=fancy_scientific) +
geom_bar(stat='identity', position='dodge', alpha=0.7) + facet_grid(. ~Map) + labs(y=expression('Partial R'^2*'/Number of SNPs'), x="cM") + scale_fill_brewer(name="Dataset", type="div", palette='Dark2') + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position='bottom', legend.title=element_blank()) + scale_x_discrete(labels=c("Low", expression(symbol('\256')), expression(symbol('\256')), "High"))
print(plot2)
ggsave(paste0('~/height_prediction/strat_prs/figs/barplot_AA_CEU_v2_', dtset,'.pdf'))




###



plot_grid(plot1, plot2,labels = c("A", "B"), nrow=2, align="v")

ggsave(paste0('~/height_prediction/strat_prs/figs/barplot_ALL_',dtset,'.pdf'))



##

if(dtset=='gwas'){
	beta<-readRDS(paste0('~/height_prediction/', dtset, '/ukb_afr/output/betas_phys_100000_0.0005_10000.Rds'))
} else {
	beta<-do.call(rbind, readRDS(paste0('~/height_prediction/', dtset, '/ukb_afr/output/betas_phys_100000_0.0005_10000.Rds')))
}
plot3<-ggplot(beta[AA.rate<=0.05], aes(x=AA.rate, y=Beta_Diff_Chisq)) + geom_point(cex=0.5, col='light gray') + geom_smooth(method='lm') + labs(y=expression(beta[eur]~-beta[afr]), x="cM (AA_Map)") + labels=c("0.1", "10","1000"))

#plot2b<-plot2 + theme(legend.text=element_text(size=8)) + theme(legend.title=element_blank())
plot4<-ggplot(beta[CEU.rate<=0.05], aes(x=CEU.rate, y=Beta_Diff_Chisq)) + geom_point(cex=0.5, col='light gray') + geom_smooth(method='lm') + labs(y=expression(beta[eur]~-beta[afr]), x="cM (CEU_Map)") +  scale_y_log10(breaks=c(0.1,10,1000), labels=c("0.1", "10","1000"))



#plot3<-ggplot(beta[AA.rate<=0.05], aes(x=AA.rate, y=y2)) + geom_point(cex=0.5, col='light gray') + geom_smooth(method='lm') + labs(y='Beta Difference', x="Recomb. Rate  (AA Map)") + scale_y_log10()

#plot2b<-plot2 + theme(legend.text=element_text(size=8)) + theme(legend.title=element_blank())
#plot4<-ggplot(beta[CEU.rate<=0.05], aes(x=CEU.rate, y=y2)) + geom_point(cex=0.5, col='light gray') + geom_smooth(method='lm') + labs(y='Beta Difference', x="Recomb. Rate (CEU Map)") + scale_y_log10()

plot1a<-plot1 + guides(fill=FALSE)
plot_grid(plot1a, plot3,plot2,plot4,labels = c("A", "C", "B","D"), nrow=2, align="v")
ggsave(paste0('~/height_prediction/strat_prs/figs/panel_', args[1], '.pdf'))


