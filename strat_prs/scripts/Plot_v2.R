
#if (length(args)==0) {
#  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
#}
library(ggplot2)
library(reshape2)
library(data.table)
library(reshape)
library(RColorBrewer)
library(MASS)
library(cowplot)
library(dplyr)
options(scipen=999)

args<-'gwas'
dtset<-args[1]
source('~/height_prediction/strat_prs/scripts/fancy_scientific.R')
a_list<-vector('list', 2)
ci<-vector('list', 2)
names(ci)<-c('AA','CEU')
names(a_list)<-c("AA","CEU")

a_list[['AA']]<-vector('list', 5)
a_list[['CEU']]<-vector('list', 5)
ci[['AA']]<-vector('list', 5)
ci[['CEU']]<-vector('list', 5)
names(a_list[['AA']])<-c("ukb_afr", "WHI", "JHS", "HRS_afr",  "HRS_eur")
names(ci[['AA']])<-c("ukb_afr", "WHI", "JHS", "HRS_afr",  "HRS_eur")
for(I in c("WHI","JHS","ukb_afr","HRS_eur", "HRS_afr")){
	a_list[['AA']][[I]]<-vector('list', 6)
	ci[['AA']][[I]]<-vector('list', 6)
	names(a_list[['AA']][[I]])<-c("3000","6000","10000","20000", "40000", "100000")
	names(ci[['AA']][[I]])<-c("3000","6000","10000","20000", "40000", "100000")
	for(J in names(a_list[['AA']][[I]])){	
		a_list[['AA']][[I]][[J]]<-readRDS(paste0('~/height_prediction/strat_prs/output/part_R2_', I,"_",dtset, "_", J, "_AA_v2.Rds"))
		ci[['AA']][[I]][[J]]<-readRDS(paste0('~/height_prediction/strat_prs/output/results_', I,"_",dtset, "_", J, "_AA_v2.Rds"))
	}
}

if(dtset=='gwas'){ ##UKB_afr, WHI_afr, JHS_afr, HRS_afr,HRS_eur
	r2_vec<-c(0.041,0.041,0.038,0.025, 0.125)
	r2_l_vec<-c(0.032,0.033, 0.024,0.015,0.112)
	r2_u_vec<-c(0.049,0.051,0.058, 0.039,0.138)
} else if (dtset=='sib_betas'){
#	r2_vec<-c(0.07511196,0.00967814,0.01993696,0.02717107,0.01060191)
}

df1<-rbind(data.table(Quantile=rep(c("q1","q2","q3","q4"),6),
Dataset='HRS_eur',Win=c(rep(3000, 4), rep(6000,4), rep(10000,4), rep(20000,4), rep(40000,4), rep(100000,4)),
R2=unlist(a_list[['AA']][['HRS_eur']])/r2_vec[5],
Perc_L=c(unlist(lapply(ci[['AA']][['HRS_eur']][['3000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['HRS_eur']][['6000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['HRS_eur']][['10000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['HRS_eur']][['20000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['HRS_eur']][['40000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['HRS_eur']][['100000']], function(X) X$percent[4]))),
Perc_U=c(unlist(lapply(ci[['AA']][['HRS_eur']][['3000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['HRS_eur']][['6000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['HRS_eur']][['10000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['HRS_eur']][['20000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['HRS_eur']][['40000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['HRS_eur']][['100000']], function(X) X$percent[5])))),
data.table(Quantile=rep(c("q1","q2","q3","q4"),6),
Dataset='JHS_afr',Win=c(rep(3000, 4), rep(6000,4), rep(10000,4), rep(20000,4), rep(40000,4), rep(100000,4)),
R2=unlist(a_list[['AA']][['JHS']])/r2_vec[3],
Perc_L=c(unlist(lapply(ci[['AA']][['JHS']][['3000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['JHS']][['6000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['JHS']][['10000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['JHS']][['20000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['JHS']][['40000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['JHS']][['100000']], function(X) X$percent[4]))),
Perc_U=c(unlist(lapply(ci[['AA']][['JHS']][['3000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['JHS']][['6000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['JHS']][['10000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['JHS']][['20000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['JHS']][['40000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['JHS']][['100000']], function(X) X$percent[5])))),
data.table(Quantile=rep(c("q1","q2","q3","q4"),6),
Dataset='WHI_afr',Win=c(rep(3000, 4), rep(6000,4), rep(10000,4), rep(20000,4), rep(40000,4), rep(100000,4)),
R2=unlist(a_list[['AA']][['WHI']])/r2_vec[2],
Perc_L=c(unlist(lapply(ci[['AA']][['WHI']][['3000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['WHI']][['6000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['WHI']][['10000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['WHI']][['20000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['WHI']][['40000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['WHI']][['100000']], function(X) X$percent[4]))),
Perc_U=c(unlist(lapply(ci[['AA']][['WHI']][['3000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['WHI']][['6000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['WHI']][['10000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['WHI']][['20000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['WHI']][['40000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['WHI']][['100000']], function(X) X$percent[5])))),
data.table(Quantile=rep(c("q1","q2","q3","q4"),6),
Dataset='UKB_afr',Win=c(rep(3000, 4), rep(6000,4), rep(10000,4), rep(20000,4), rep(40000,4), rep(100000,4)),
R2=unlist(a_list[['AA']][['ukb_afr']])/r2_vec[1],
Perc_L=c(unlist(lapply(ci[['AA']][['ukb_afr']][['3000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['ukb_afr']][['6000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['ukb_afr']][['10000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['ukb_afr']][['20000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['ukb_afr']][['40000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['ukb_afr']][['100000']], function(X) X$percent[4]))),
Perc_U=c(unlist(lapply(ci[['AA']][['ukb_afr']][['3000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['ukb_afr']][['6000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['ukb_afr']][['10000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['ukb_afr']][['20000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['ukb_afr']][['40000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['ukb_afr']][['100000']], function(X) X$percent[5])))),
data.table(Quantile=rep(c("q1","q2","q3","q4"),6),
Dataset='HRS_afr',Win=c(rep(3000, 4), rep(6000,4), rep(10000,4), rep(20000,4), rep(40000,4), rep(100000,4)),
R2=unlist(a_list[['AA']][['HRS_afr']])/r2_vec[4],
Perc_L=c(unlist(lapply(ci[['AA']][['HRS_afr']][['3000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['HRS_afr']][['6000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['HRS_afr']][['10000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['HRS_afr']][['20000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['HRS_afr']][['40000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['HRS_afr']][['100000']], function(X) X$percent[4]))),
Perc_U=c(unlist(lapply(ci[['AA']][['HRS_afr']][['3000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['HRS_afr']][['6000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['HRS_afr']][['10000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['HRS_afr']][['20000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['HRS_afr']][['40000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['HRS_afr']][['100000']], function(X) X$percent[5])))))

####
df2<-rbind(data.table(Quantile=rep(c("q1","q2","q3","q4"),6),
Dataset='HRS_eur',Win=c(rep(3000, 4), rep(6000,4), rep(10000,4), rep(20000,4), rep(40000,4), rep(100000,4)),
R2=unlist(a_list[['AA']][['HRS_eur']]),
Perc_L=c(unlist(lapply(ci[['AA']][['HRS_eur']][['3000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['HRS_eur']][['6000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['HRS_eur']][['10000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['HRS_eur']][['20000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['HRS_eur']][['40000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['HRS_eur']][['100000']], function(X) X$percent[4]))),
Perc_U=c(unlist(lapply(ci[['AA']][['HRS_eur']][['3000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['HRS_eur']][['6000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['HRS_eur']][['10000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['HRS_eur']][['20000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['HRS_eur']][['40000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['HRS_eur']][['100000']], function(X) X$percent[5])))),
data.table(Quantile=rep(c("q1","q2","q3","q4"),6),
Dataset='JHS_afr',Win=c(rep(3000, 4), rep(6000,4), rep(10000,4), rep(20000,4), rep(40000,4), rep(100000,4)),
R2=unlist(a_list[['AA']][['JHS']]),
Perc_L=c(unlist(lapply(ci[['AA']][['JHS']][['3000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['JHS']][['6000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['JHS']][['10000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['JHS']][['20000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['JHS']][['40000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['JHS']][['100000']], function(X) X$percent[4]))),
Perc_U=c(unlist(lapply(ci[['AA']][['JHS']][['3000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['JHS']][['6000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['JHS']][['10000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['JHS']][['20000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['JHS']][['40000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['JHS']][['100000']], function(X) X$percent[5])))),
data.table(Quantile=rep(c("q1","q2","q3","q4"),6),
Dataset='WHI_afr',Win=c(rep(3000, 4), rep(6000,4), rep(10000,4), rep(20000,4), rep(40000,4), rep(100000,4)),
R2=unlist(a_list[['AA']][['WHI']]),
Perc_L=c(unlist(lapply(ci[['AA']][['WHI']][['3000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['WHI']][['6000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['WHI']][['10000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['WHI']][['20000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['WHI']][['40000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['WHI'
]][['100000']], function(X) X$percent[4]))),
Perc_U=c(unlist(lapply(ci[['AA']][['WHI']][['3000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['WHI']][['6000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['WHI']][['10000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['WHI']][['20000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['WHI']][['40000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['WHI']][['100000']], function(X) X$percent[5])))),
data.table(Quantile=rep(c("q1","q2","q3","q4"),6),
Dataset='UKB_afr',Win=c(rep(3000, 4), rep(6000,4), rep(10000,4), rep(20000,4), rep(40000,4), rep(100000,4)),
R2=unlist(a_list[['AA']][['ukb_afr']]),
Perc_L=c(unlist(lapply(ci[['AA']][['ukb_afr']][['3000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['ukb_afr']][['6000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['ukb_afr']][['10000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['ukb_afr']][['20000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['ukb_afr']][['40000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['ukb_afr']][['100000']], function(X) X$percent[4]))),
Perc_U=c(unlist(lapply(ci[['AA']][['ukb_afr']][['3000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['ukb_afr']][['6000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['ukb_afr']][['10000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['ukb_afr']][['20000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['ukb_afr']][['40000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['ukb_afr']][['100000']], function(X) X$percent[5])))),
data.table(Quantile=rep(c("q1","q2","q3","q4"),6),
Dataset='HRS_afr',Win=c(rep(3000, 4), rep(6000,4), rep(10000,4), rep(20000,4), rep(40000,4), rep(100000,4)),
R2=unlist(a_list[['AA']][['HRS_afr']]),
Perc_L=c(unlist(lapply(ci[['AA']][['HRS_afr']][['3000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['HRS_afr']][['6000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['HRS_afr']][['10000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['HRS_afr']][['20000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['HRS_afr']][['40000']], function(X) X$percent[4])), unlist(lapply(ci[['AA']][['HRS_afr']][['100000']], function(X) X$percent[4]))),
Perc_U=c(unlist(lapply(ci[['AA']][['HRS_afr']][['3000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['HRS_afr']][['6000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['HRS_afr']][['10000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['HRS_afr']][['20000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['HRS_afr']][['40000']], function(X) X$percent[5])), unlist(lapply(ci[['AA']][['HRS_afr']][['100000']], function(X) X$percent[5])))))

#
df1[, Set:='Relative'] #order ##UKB_afr, WHI_afr, JHS_afr, HRS_afr,HRS_eur
df1[, Perc_U:=Perc_U/c(rep(r2_u_vec[5],24), rep(r2_u_vec[3], 24), rep(r2_u_vec[2], 24), rep(r2_u_vec[1], 24), rep(r2_u_vec[4], 24))]
df1[, Perc_L:=Perc_L/c(rep(r2_l_vec[5],24), rep(r2_l_vec[3], 24), rep(r2_l_vec[2], 24), rep(r2_l_vec[1], 24), rep(r2_l_vec[4], 24))]
df2[, Set:='Absolute']

df4<-rbind(df1, df2)
#melt(df4, id=c("Quantile","Win", "Map", "Perc_L", "Perc_U", "Dataset"))-> df5
df4$Dataset<-factor(df4$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "HRS_eur"))

df<-df4[Win==20000]
df$Set<-factor(df$Set, levels=c("Absolute", "Relative"))
pd <- position_dodge(0.5)
plot1<-ggplot(df, aes(x=Quantile, y=R2, colour=Dataset)) + 
#geom_bar(stat='identity', position='dodge', alpha=0.8) + 
geom_point(position=pd) +
geom_errorbar(aes(ymin=Perc_L, ymax=Perc_U), position = pd) +
facet_wrap(. ~Set, scales='free_y') + 
labs(y=expression(paste("Partial R"^"2")), x="Recombination Rate")+  
scale_colour_manual(values=c(brewer.pal(4, 'Set1'),"#101010")) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.title=element_blank(), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),legend.text=element_text(size=12),legend.position = "bottom") + 
scale_x_discrete(labels=c("Low", expression(symbol('\256')), expression(symbol('\256')), "High"))

#axis.line = element_line(colour = "black"), legend.position = "none", legend.title=element_blank(), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=15),axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.text=element_text(size=10)) + scale_x_discrete(labels=c("Low", expression(symbol('\256')), expression(symbol('\256')), "High"))

print(plot1)
ggsave(paste0('~/height_prediction/strat_prs/figs/v2_barplot_AA_', dtset,'.pdf'))

#
df2$Dataset<-factor(df2$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "HRS_eur"))

df2<-df2[Win==20000]

pd <- position_dodge(0.5)
plot2<-ggplot(df2, aes(x=Quantile, y=R2, colour=Dataset)) +
#geom_bar(stat='identity', position='dodge', alpha=0.8) +
geom_point(position=pd) +
geom_errorbar(aes(ymin=Perc_L, ymax=Perc_U), position = pd) +
#faet_grid(. ~Map) +
labs(y=expression(paste("Partial R"^"2")),x="Recombination Rate") +
scale_colour_manual(values=c(brewer.pal(4, 'Set1'),"#101010")) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "right", legend.direction='vertical', legend.title=element_blank(), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.text=element_text(size=10)) + scale_x_discrete(labels=c("Low", expression(symbol('\256')), expression(symbol('\256')), "High"))


#

#########################
########################
########################

if(dtset=='gwas'){
	#beta<-readRDS(paste0('~/height_prediction/', dtset, '/ukb_afr/output/betas_phys_100000_0.0005_20000.Rds'))  #THIS NEEDS TO BE PLINK EFFECT SIZES NEED TO FIX
	beta<-select(do.call(rbind,readRDS('~/height_prediction/gwas/WHI/output/hei_phys_100000_0.0005_v2.Rds')), CHR, POS, b, MarkerName, Allele1, SE) %>% as.data.table
	beta1<-readRDS('~/height_prediction/loc_anc_analysis/output/final_plink.Rds')
} else {
	beta<-do.call(rbind, readRDS(paste0('~/height_prediction/', dtset, '/ukb_afr/output/betas_phys_100000_0.0005_20000.Rds')))
}

merge(beta, beta1, by=c('CHR', 'POS'))-> beta2
#pause


BETA<-vector('list', 22)
prun<-"phys_100000_0.0005"
hei<-readRDS(paste0('~/height_prediction/gwas/WHI/output/hei_', prun, '_v2.Rds'))
for(chr in 22:1){
cat('start chr ')
cat(chr)
cat('\n')
rate.dist<-20000
betas<-beta2[CHR==chr]
maps<-fread(paste0('zcat /project/mathilab/data/maps_b37/maps_chr.', chr, '.gz'))
snps <- read.table(paste0("~/height_prediction/input/WHI/WHI_b37_strand_include_kgCY_chr", chr, ".phsnp"), as.is=TRUE)
colnames(snps) <- c("ID", "CHR", "Map", "POS", "REF", "ALT")
betas<-merge(snps, betas, by=c('CHR', 'POS'))
setDT(betas)
cat('checkpoint \n')
AA.rate <- approxfun(maps$Physical_Pos, maps$AA_Map)
YRI.rate <- approxfun(maps$Physical_Pos, maps$YRI_LD)
CEU.rate <- approxfun(maps$Physical_Pos, maps$CEU_LD)
COMBINED.rate<-approxfun(maps$Physical_Pos,maps$COMBINED_LD)
betas$AA.rate <- 0
betas$CEU_YRI_diff.rate <- 0
betas$CEU.rate <- 0
betas$YRI.rate <- 0
betas$COMBINED.rate<-0
cat('checkpoint \n')
for(i in 1:NROW(betas)){
    cat(i, '\r')
    AA.x <- AA.rate(betas$POS[i]+(rate.dist/2))-AA.rate(betas$POS[i]-(rate.dist/2))
    betas$AA.rate[i] <- AA.x
    CEU.x <- CEU.rate(betas$POS[i]+(rate.dist/2))-CEU.rate(betas$POS[i]-(rate.dist/2))
    betas$CEU.rate[i] <- CEU.x
    YRI.x <- YRI.rate(betas$POS[i]+(rate.dist/2))-YRI.rate(betas$POS[i]-(rate.dist/2))
    betas$YRI.rate[i]<-YRI.x
    COMBINED.x<-COMBINED.rate(betas$POS[i]+(rate.dist/2))-COMBINED.rate(betas$POS[i]-(rate.dist/2))
    betas$COMBINED.rate[i]<-COMBINED.x
    betas$CEU_YRI_diff.rate[i] <- abs(CEU.x-YRI.x)
}
betas[,Beta_Diff:=b-PLINK]
betas[,Beta_Diff_Chisq:=(Beta_Diff/sqrt(((SE^2)+(SE_plink^2))))^2]
#At this point you will restrict the betas to the SNPS that you are using in the PRS.
BETA[[chr]]<-betas
cat('chr ')
cat(chr)
cat(' done\n')
}
cat('sleep\n')
Sys.sleep(30)
#saveRDS(beta, 'betas.Rds')
#pdf(paste0("chr", chr, "_rate_against_AA_map.pdf"))
do.call(rbind, BETA)-> BETA
setDT(BETA)
BETA[order(CHR, POS)]-> BETA
as.factor(BETA$CHR)-> BETA$CHR

saveRDS(BETA, file=paste0('~/height_prediction/gwas/WHI/output/plink_whi.Rds'))

#########################
########################
########################
#########################
########################
########################
BETA[,Quantile:=cut(AA.rate, breaks=quantile(AA.rate, probs=seq(0,1, by=0.05), na.rm=T), include.lowest=T)]
BETA[, MeanBetaDiffChisq:=mean(Beta_Diff_Chisq, na.rm=T), by=Quantile]
BETA[, MedianRecRate:=median(AA.rate, na.rm=T), by=Quantile]
lm_test0<-lm(Beta_Diff_Chisq~AA.rate, data=BETA)
require(broom)
glance(lm_test0)
pval<-glance(lm_test0)$p.value
plot3<-ggplot(BETA, aes(x=AA.rate, y=Beta_Diff_Chisq)) + geom_point(cex=0.5, col='light gray') + geom_smooth(method='lm', se=T, lwd=1, col="black") + labs(y=expression(chi[diff]^2), x="Recombination Rate") + geom_point(aes(x=MedianRecRate, y=MeanBetaDiffChisq, col="red"), cex=0.5) + 
annotate("text", x=0.4, y=5, label=paste("p=", round(pval,4))) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none", legend.title=element_blank(), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.text=element_text(size=15))
plot1a<-plot1 + guides(fill=FALSE)
cat('CHECKPOINT\n')
ld<-do.call(rbind, lapply(1:22, function(X) fread(paste0('zcat ~/height_prediction/figs_for_paper/eur_w_ld_chr/', X,'.l2.ldscore.gz'))))
#beta<-readRDS('~/height_prediction/gwas/ukb_afr/output/betas_phys_100000_0.0005_20000.Rds')
colnames(ld)[3]<-'POS'
colnames(ld)[2]<-'MarkerName.y'
setkey(ld, CHR, POS)
#beta$CHR<-as.integer(beta$CHR)
setkey(BETA, CHR, POS)
#beta[ld,nomatch=0]-> test
test<-merge(BETA,ld, by=c('CHR', 'POS','MarkerName.y'))

test[,Quantile:=cut(L2, breaks=quantile(L2, probs=seq(0,1, by=0.05), na.rm=T), include.lowest=T)]
test[, MeanBetaDiffChisq:=mean(Beta_Diff_Chisq, na.rm=T), by=Quantile]
test[, MedianL2:=median(L2, na.rm=T), by=Quantile]

#plot4<-ggplot(test, aes(x=MedianL2, y=MeanBetaDiffChisq)) + geom_point(cex=0.5, col='light gray') + geom_smooth(method='lm', se=F) + labs(y=expression(Mean~chi^2), x="LD Score" )
lm_test<-lm(Beta_Diff_Chisq~L2, data=test)
require(broom)
glance(lm_test)
pval<-glance(lm_test)$p.value
plot4<-ggplot(test, aes(x=L2, y=Beta_Diff_Chisq)) + geom_point(cex=0.5, col='light gray') + geom_smooth(method='lm', se=T, col='black') + 
labs(y=expression(chi[diff]^2), x="LD Score" ) + 
geom_point(aes(x=MedianL2, y=MeanBetaDiffChisq, col="red"), cex=0.5) + 
annotate("text", x=250, y=5, label=paste("p=", round(pval,4))) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none", legend.title=element_blank(), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.text=element_text(size=15))

plot_grid(plot1,plot_grid(plot3,plot4,labels = c("C", "D"), nrow=1), nrow=2, labels="A", align="v")
ggsave(paste0('~/height_prediction/strat_prs/figs/panel_', args[1], '_v2.pdf'))
#The End
