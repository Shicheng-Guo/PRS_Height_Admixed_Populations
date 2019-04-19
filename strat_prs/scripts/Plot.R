library(ggplot2)
library(reshape2)
library(data.table)
library(reshape)
library(RColorBrewer)
library(MASS)
library(cowplot)

source('fancy_scientific.R')
a_list<-vector('list', 2)
names(a_list)<-c("AA","CEU")

a_list[['AA']]<-vector('list', 6)
a_list[['CEU']]<-vector('list', 6)
names(a_list[['AA']])<-c("WHI","JHS","pennBB_afr","ukb_afr","pennBB_eur","UKB_eur")
names(a_list[['CEU']])<-c("WHI","JHS","pennBB_afr","ukb_afr","pennBB_eur","UKB_eur")

for(I in c("WHI","JHS","pennBB_afr","ukb_afr","pennBB_eur","UKB_eur")){
	a_list[['AA']][[I]]<-vector('list', 3)
	a_list[['CEU']][[I]]<-vector('list', 3)
	names(a_list[['AA']][[I]])<-c("5000","10000","20000")
 	names(a_list[['CEU']][[I]])<-c("5000","10000","20000")
	for(J in names(a_list[['AA']][[I]])){	
		a_list[['AA']][[I]][[J]]<-readRDS(paste0('/project/mathilab/bbita/gwas_admix/height_prediction/strat_prs/output/part_R2_', I, "_", J, "_AA.Rds"))
		a_list[['CEU']][[I]][[J]]<-readRDS(paste0('/project/mathilab/bbita/gwas_admix/height_prediction/strat_prs/output/part_R2_', I, "_", J, "_CEU.Rds"))
	}
}

df1<-data.table(Quantile=rep(c("q1","q2","q3","q4"),3),
Win=c(rep(5000,4), rep(10000,4), rep(20000,4)),
JHS=c(unlist(a_list[['AA']][[2]][[1]]), unlist(a_list[['AA']][[2]][[2]]), unlist(a_list[['AA']][[2]][[3]]))/0.03561384,
WHI=c(unlist(a_list[['AA']][[1]][[1]]),  unlist(a_list[['AA']][[1]][[2]]),  unlist(a_list[['AA']][[1]][[3]]))/0.04081778,
PennBB_afr=c(unlist(a_list[['AA']][[3]][[1]]),unlist(a_list[['AA']][[3]][[2]]),unlist(a_list[['AA']][[3]][[3]]))/0.03432268,
Ukb_afr=c(unlist(a_list[['AA']][[4]][[1]]), unlist(a_list[['AA']][[4]][[2]]), unlist(a_list[['AA']][[4]][[3]]))/0.0362177,
PennBB_eur=c(unlist(a_list[['AA']][[5]][[1]]), unlist(a_list[['AA']][[5]][[2]]), unlist(a_list[['AA']][[5]][[3]]))/0.1500369,
Ukb_eur=c(unlist(a_list[['AA']][[6]][[1]]), unlist(a_list[['AA']][[6]][[2]]), unlist(a_list[['AA']][[6]][[3]]))/0.225137
)

df3<-data.table(Quantile=rep(c("q1","q2","q3","q4"),3),
Win=c(rep(5000,4), rep(10000,4), rep(20000,4)),
JHS=c(unlist(a_list[['CEU']][[2]][[1]]), unlist(a_list[['CEU']][[2]][[2]]), unlist(a_list[['CEU']][[2]][[3]]))/0.03561384,
WHI=c(unlist(a_list[['CEU']][[1]][[1]]), unlist(a_list[['CEU']][[1]][[2]]), unlist(a_list[['CEU']][[1]][[3]]))/0.04081778,
PennBB_afr=c(unlist(a_list[['CEU']][[3]][[1]]), unlist(a_list[['CEU']][[3]][[2]]), unlist(a_list[['CEU']][[3]][[3]]))/0.03432268,
Ukb_afr=c(unlist(a_list[['CEU']][[4]][[1]]), unlist(a_list[['CEU']][[4]][[2]]), unlist(a_list[['CEU']][[4]][[3]]))/0.0362177,
PennBB_eur=c(unlist(a_list[['CEU']][[5]][[1]]),unlist(a_list[['CEU']][[5]][[2]]), unlist(a_list[['CEU']][[5]][[3]]))/0.1500369,
Ukb_eur=c(unlist(a_list[['CEU']][[6]][[1]]),unlist(a_list[['CEU']][[6]][[2]]), unlist(a_list[['CEU']][[6]][[3]]))/0.225137
)

df1[, Map:='AA']
df3[, Map:='CEU']
df4<-rbind(df1,df3)



melt(df4, id=c("Quantile","Win", "Map"))-> df5

plot1<-ggplot(df5[Win==10000], aes(x=Quantile, y=value, fill=variable)) + geom_bar(stat='identity', position='dodge', alpha=0.7) + facet_grid(. ~Map) + labs(y=expression(paste("Relative R"^"2")),x="Recombination Rate") +  scale_fill_brewer(name="Dataset", type="div", palette='Dark2') + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none")


ggsave('/project/mathilab/bbita/gwas_admix/height_prediction/strat_prs/figs/barplot_AA_CEU.pdf')


#


df1<-data.table(Quantile=rep(c("q1","q2","q3","q4"),3),
Win=c(rep(5000,4), rep(10000,4), rep(20000,4)),
JHS=c(unlist(a_list[['AA']][[2]][[1]])/readRDS('../output/Nr_SNPs_JHS_5000_AA.Rds'), unlist(a_list[['AA']][[2]][[2]])/readRDS('../output/Nr_SNPs_JHS_10000_AA.Rds'), unlist(a_list[['AA']][[2]][[3]])/readRDS('../output/Nr_SNPs_JHS_20000_AA.Rds')),
WHI=c(unlist(a_list[['AA']][[1]][[1]])/readRDS('../output/Nr_SNPs_WHI_5000_AA.Rds'), unlist(a_list[['AA']][[1]][[2]])/readRDS('../output/Nr_SNPs_WHI_10000_AA.Rds'), unlist(a_list[['AA']][[1]][[3]])/readRDS('../output/Nr_SNPs_WHI_20000_AA.Rds')),
PennBB_afr=c(unlist(a_list[['AA']][[3]][[1]])/readRDS('../output/Nr_SNPs_pennBB_afr_5000_AA.Rds'),unlist(a_list[['AA']][[3]][[2]])/readRDS('../output/Nr_SNPs_pennBB_afr_10000_AA.Rds'),unlist(a_list[['AA']][[3]][[3]])/readRDS('../output/Nr_SNPs_pennBB_afr_20000_AA.Rds')),
Ukb_afr=c(unlist(a_list[['AA']][[4]][[1]])/readRDS('../output/Nr_SNPs_ukb_afr_5000_AA.Rds'),unlist(a_list[['AA']][[4]][[2]])/readRDS('../output/Nr_SNPs_ukb_afr_10000_AA.Rds'),unlist(a_list[['AA']][[4]][[3]])/readRDS('../output/Nr_SNPs_ukb_afr_20000_AA.Rds')),
PennBB_eur=c(unlist(a_list[['AA']][[5]][[1]])/readRDS('../output/Nr_SNPs_pennBB_eur_5000_AA.Rds'), unlist(a_list[['AA']][[5]][[2]])/readRDS('../output/Nr_SNPs_pennBB_eur_10000_AA.Rds'),unlist(a_list[['AA']][[5]][[3]])/readRDS('../output/Nr_SNPs_pennBB_eur_20000_AA.Rds')),
Ukb_eur=c(unlist(a_list[['AA']][[6]][[1]])/readRDS('../output/Nr_SNPs_ukb_eur_5000_AA.Rds'), unlist(a_list[['AA']][[6]][[2]])/readRDS('../output/Nr_SNPs_ukb_eur_10000_AA.Rds'), unlist(a_list[['AA']][[6]][[3]])/readRDS('../output/Nr_SNPs_ukb_eur_20000_AA.Rds'))
)

df3<-data.table(Quantile=rep(c("q1","q2","q3","q4"),3),
Win=c(rep(5000,4), rep(10000,4), rep(20000,4)),
JHS=c(unlist(a_list[['CEU']][[2]][[1]])/readRDS('../output/Nr_SNPs_JHS_5000_CEU.Rds'), unlist(a_list[['CEU']][[2]][[2]])/readRDS('../output/Nr_SNPs_JHS_10000_CEU.Rds'), unlist(a_list[['CEU']][[2]][[3]])/readRDS('../output/Nr_SNPs_JHS_20000_CEU.Rds')),
WHI=c(unlist(a_list[['CEU']][[1]][[1]])/readRDS('../output/Nr_SNPs_WHI_5000_CEU.Rds'), unlist(a_list[['CEU']][[1]][[2]])/readRDS('../output/Nr_SNPs_WHI_10000_CEU.Rds'), unlist(a_list[['CEU']][[1]][[3]])/readRDS('../output/Nr_SNPs_WHI_20000_CEU.Rds')),
PennBB_afr=c(unlist(a_list[['CEU']][[3]][[1]])/readRDS('../output/Nr_SNPs_pennBB_afr_5000_CEU.Rds'),unlist(a_list[['CEU']][[3]][[2]])/readRDS('../output/Nr_SNPs_pennBB_afr_10000_CEU.Rds'),unlist(a_list[['CEU']][[3]][[3]])/readRDS('../output/Nr_SNPs_pennBB_afr_20000_CEU.Rds')),
Ukb_afr=c(unlist(a_list[['CEU']][[4]][[1]])/readRDS('../output/Nr_SNPs_ukb_afr_5000_CEU.Rds'),unlist(a_list[['CEU']][[4]][[2]])/readRDS('../output/Nr_SNPs_ukb_afr_10000_CEU.Rds'),unlist(a_list[['CEU']][[4]][[3]])/readRDS('../output/Nr_SNPs_ukb_afr_20000_CEU.Rds')),
PennBB_eur=c(unlist(a_list[['CEU']][[5]][[1]])/readRDS('../output/Nr_SNPs_pennBB_eur_5000_CEU.Rds'), unlist(a_list[['CEU']][[5]][[2]])/readRDS('../output/Nr_SNPs_pennBB_eur_10000_CEU.Rds'),unlist(a_list[['CEU']][[5]][[3]])/readRDS('../output/Nr_SNPs_pennBB_eur_20000_CEU.Rds')),
Ukb_eur=c(unlist(a_list[['CEU']][[6]][[1]])/readRDS('../output/Nr_SNPs_ukb_eur_5000_CEU.Rds'), unlist(a_list[['CEU']][[6]][[2]])/readRDS('../output/Nr_SNPs_ukb_eur_10000_CEU.Rds'), unlist(a_list[['CEU']][[6]][[3]])/readRDS('../output/Nr_SNPs_ukb_eur_20000_CEU.Rds'))
)

df1[, Map:='AA']
df3[, Map:='CEU']
df4<-rbind(df1,df3)

melt(df4, id=c("Quantile","Win","Map"))-> df5


plot2<-ggplot(df5[Win==10000], aes(x=Quantile, y=value, fill=variable)) +  scale_y_continuous(labels=fancy_scientific) +
geom_bar(stat='identity', position='dodge', alpha=0.7) + facet_grid(. ~Map) + labs(y=expression('R'^2*'/Number of SNPs'), x="Recombination Rate") + scale_fill_brewer(name="Dataset", type="div", palette='Dark2') + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position='bottom') 

ggsave('/project/mathilab/bbita/gwas_admix/height_prediction/strat_prs/figs/barplot_AA_CEU_v2.pdf')


plot_grid(plot1, plot2,labels = c("A", "B"), nrow=2, align="v")

ggsave('/project/mathilab/bbita/gwas_admix/height_prediction/strat_prs/figs/barplot_ALL.pdf')


