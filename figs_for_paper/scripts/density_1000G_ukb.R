#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)
#if (length(args)==0) {
#  stop("At least one argument must be supplied (input file).n", call.=FALSE)
#}
#Load packages
library("optparse")
library(data.table)
library(dplyr)
library(ggplot2);library(reshape2); library(wesanderson)
library(rlist)
library(cowplot)
source('~/height_prediction/scripts/PolygenicScore_v2.R')
args<-"phys_100000_0.0005"
#********************************
#* Read in 1000G samples data ***
#********************************
samples<-fread('/project/mathilab/data/1kg/20130502_phase3_final/integrated_call_samples_v3.20130502.ALL.panel', fill=T)[,.(sample, pop, super_pop)]
eur<-unique(samples[super_pop=='EUR']$pop)
afr<-unique(samples[super_pop=='AFR']$pop)
amr<-unique(samples[super_pop=='AMR']$pop)
eas<-unique(samples[super_pop=='EAS']$pop)
sas<-unique(samples[super_pop=='SAS']$pop)
unique(samples[, pop])-> pops
unique(samples[,super_pop])-> superpops
allpops<-c(afr, eur, sas, eas, amr, superpops)
#*********************************
#* Calculate Polygenic Scores ******
#*********************************
#Run PRS for 1000 Genomes using UKBiobank	
ukb_height<-fread('zcat ~/height_prediction/gwas/input/50.assoc.tsv.gz') #read in GWAS summary statistics for height from the Uk Biobank
ukb_height[,c("CHR", "POS","Allele2","Allele1") := tstrsplit(variant, ":", fixed=TRUE)][,variant:=NULL]  #fix columns. In the UKB, variants are listed as CHR:POS: REF:ALT, where ALT is the effect allele. So, in order to be compatilble with my scripts (where allele 1 is effect all
ukb_height[, MarkerName:=rsid][, N:=nCompleteSamples][, AC:=NULL][, b:=beta][,p:=pval]
ukb_height[,rsid:=NULL][,nCompleteSamples:=NULL][, beta:=NULL][, pval:=NULL][, SE:=se]
ukb_height[,.(MarkerName,Allele1,Allele2, b, SE, p, N, CHR, POS)]-> ukb_height

ukb_height[,CHR:=as.integer(CHR)]
ukb_height[,POS:=as.integer(POS)]

ukb_height[order(CHR,POS)]-> ukb_height

cat('Got SNP IDs, CHR and POS.\n')

hei<-vector('list', 22)
for(Z in 1:22){
	fread(paste0('zcat /project/mathilab/bbita/gwas_admix/height/LD_prun/hei_SNPs_phys_100000_0.0005_chr', Z, '.vcf.gz'))-> hei[[Z]]
        colnames(hei[[Z]])[1:3]<-c('CHR','POS', 'MarkerName')
        cat(paste0('chr ', Z, ' done\n'))
        setkey(hei[[Z]], CHR, POS)
        unique(hei[[Z]])-> hei[[Z]]
        hei[[Z]][REF %in% c("A","C","G","T")]-> hei[[Z]]
        setkey(ukb_height, CHR, POS)
        setkey(hei[[Z]], CHR, POS)
        hei[[Z]][ukb_height, nomatch=0]-> hei[[Z]]
}

names(hei)<-seq(1:22)
saveRDS(hei, file=paste0('~/height_prediction/gwas/ukb_eur/output/hei_phys_100000_0.0005_1000g.Rds'))
allpops<-c(unique(samples$pop), unique(samples$super_pop))
PGS<-vector('list', length(allpops))
names(PGS)<-allpops
pops<-c(unique(samples$pop))
superpops<-unique(samples$super_pop)


for (P in pops){
	cat(paste0("Population is ", P, '\n'))
	PGS[[P]]<-vector('list', 22)
        for (CR in 1:22){
		print(paste0("Chromosome is ", CR))
                PolScore_1000G(POP=P, CHR=CR, superpop=F)-> PGS[[P]][[CR]]
                cat(paste0(CR, '  done\n'))
        }
        cat (paste0(P, ' done\n'))
        gc()
}

for (P in superpops){
	PGS[[P]]<-vector('list', 22)
	cat(paste0("Super-population is ", P, '\n'))
        for(CR in 1:22){
		print(paste0("Chromosome is ", CR))
                PolScore_1000G(CHR=CR, POP=P, superpop=T)-> PGS[[P]][[CR]]
                cat(paste0(CR, '  done\n'))
	}
        cat(paste0(P, ' done\n'))
        gc()
}
vector('list', length(allpops))-> PGS2
names(PGS2)<-allpops
for (P in 1:length(allpops)){ #create list with names of individuals for each population. Each individual has a polygenic score.
	PGS2[[P]]<-vector('list', length(names(PGS[[P]][[21]])))
        names(PGS2[[P]])<-names(PGS[[P]][[21]])
}
for (P in allpops){
	samps<-names(PGS[[P]][[21]]) #sum PGS across chromosomes.
	for (S in samps){
        	PGS[[P]][[1]][[S]] + PGS[[P]][[2]][[S]] + PGS[[P]][[3]][[S]] + PGS[[P]][[4]][[S]] + PGS[[P]][[5]][[S]] + PGS[[P]][[6]][[S]] + PGS[[P]][[7]][[S]] + PGS[[P]][[8]][[S]] + PGS[[P]][[9]][[S]] + PGS[[P]][[10]][[S]] + PGS[[P]][[11]][[S]] +
                PGS[[P]][[12]][[S]]+PGS[[P]][[13]][[S]]+PGS[[P]][[14]][[S]]+PGS[[P]][[15]][[S]]+PGS[[P]][[16]][[S]]+PGS[[P]][[17]][[S]]+PGS[[P]][[18]][[S]]+PGS[[P]][[19]][[S]]+PGS[[P]][[20]][[S]]+PGS[[P]][[21]][[S]]+PGS[[P]][[22]][[S]]->PGS2[[P]][[S]]
        }
}

saveRDS(PGS2, file=paste0('~/height_prediction/gwas/ukb_eur/output/PGS_1000G_phys_100000_0.0005.Rds'))


#*******************************
	if(dt=='final2'){
#* Prepare data for plotting ***
#*******************************
path<-paste0('~/height_prediction/figs_for_paper/figs/')
PGS<-readRDS('~/height_prediction/gwas/ukb_eur/output/PGS_1000G_phys_100000_0.0005.Rds')

PGS2<- lapply(PGS, function(X) unlist(X))
dat <- lapply(PGS2, function(x) cbind(x = seq_along(x), y = x))
list.names <- names(dat)
lns <- sapply(dat, nrow)
as.data.table(do.call(rbind, dat))-> dat
dat$group <- rep(list.names, lns)
as.factor(dat$group)-> dat$group
rbind(dat[group %in% afr], dat[group %in% eur], dat[group %in% sas], dat[group %in% eas], dat[group %in% amr], dat[group %in% superpops])-> dat
dat$group<-factor(dat$group, levels=allpops)
 mean(dat[group=='CEU'][,y])->m1 #normalize PGS by CEU PGS values
sd(dat[group=='CEU'][,y])->sd1
dat[, y1:=scale(as.matrix(y), center=m1, scale=sd1), by=group]

##Plotting

ggplot(dat[group %in% pops],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "All populations (1000G)") + scale_color_brewer(palette="PRGn")
ggsave(paste0(path,'PGS_1000G_pops_UKB.png'))
ggplot(dat[group %in% superpops],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "Superpopulations (1000G)") + scale_color_brewer(palette="PRGn")
ggsave(paste0(path, 'PGS_1000G_superpops_UKB.png'))
ggplot(dat[group %in% c(afr,eur)],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "African and European Populations (1000G)") + scale_color_brewer(palette="PRGn")
ggsave(paste0(path,'PGS_1000G_EUR_AFR_UKB.png'))
ggplot(dat[group %in% eur],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "European Populations (1000G)") + scale_color_brewer(palette="PRGn")
ggsave(paste0(path,'PGS_1000G_EUR_UKB.png'))
ggplot(dat[group %in% afr],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "African Populations (1000G)") + scale_color_brewer(palette="PRGn")
ggsave(paste0(path, 'PGS_1000G_AFR_UKB.png'))
ggplot(dat[group %in% amr],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "American Populations (1000G)") + scale_color_brewer(palette="PRGn")
ggsave(paste0(path, 'PGS_1000G_AMR_UKB.png'))
ggplot(dat[group %in% sas],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "South Asian Populations (1000G)") + scale_color_brewer(palette="PRGn")
ggsave(paste0(path, 'PGS_1000G_SAS_UKB.png'))
ggplot(dat[group %in% eas],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "East Asian Populations (1000G)") + scale_color_brewer(palette="PRGn")
ggsave(paste0(path,'PGS_1000G_EAS_UKB.png'))
print('Done making figures')

#
paste0('/project/mathilab/bbita/gwas_admix/height/PGS_1000G_',args[1],'.Rds')-> sm
readRDS(sm)->PGS_giant
PGS2_giant<- lapply(PGS_giant, function(X) unlist(X))
dat_giant <- lapply(PGS2_giant, function(x) cbind(x = seq_along(x), y = x))
list.names_giant <- names(dat_giant)
lns_giant <- sapply(dat_giant, nrow)
as.data.table(do.call(rbind, dat_giant))-> dat_giant
dat_giant$group <- rep(list.names_giant, lns_giant)
as.factor(dat_giant$group)-> dat_giant$group
rbind(dat_giant[group %in% afr], dat_giant[group %in% eur], dat_giant[group %in% sas], dat_giant[group %in% eas], dat_giant[group %in% amr], dat_giant[group %in% superpops])-> dat_giant
dat_giant$group<-factor(dat_giant$group, levels=allpops)
mean(dat_giant[group=='CEU'][,y])->m1_giant #normalize PGS by CEU PGS values
sd(dat_giant[group=='CEU'][,y])->sd1_giant
dat_giant[, y1:=scale(as.matrix(y), center=m1_giant, scale=sd1_giant), by=group]



ggplot(dat_giant[group %in% pops],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "All populations (1000G)") + scale_color_brewer(palette="PRGn")
ggsave(paste0(path,'PGS_1000G_pops_GIANT.png'))
ggplot(dat_giant[group %in% superpops],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "Superpopulations (1000G)") + scale_color_brewer(palette="PRGn")
ggsave(paste0(path, 'PGS_1000G_superpops_GIANT.png'))
ggplot(dat_giant[group %in% c(afr,eur)],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "African and European Populations (1000G)") + scale_color_brewer(palette="PRGn")
ggsave(paste0(path,'PGS_1000G_EUR_AFR_GIANT.png'))
ggplot(dat_giant[group %in% eur],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "European Populations (1000G)") + scale_color_brewer(palette="PRGn")
ggsave(paste0(path,'PGS_1000G_EUR_GIANT.png'))
ggplot(dat_giant[group %in% afr],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "African Populations (1000G)") + scale_color_brewer(palette="PRGn")
ggsave(paste0(path, 'PGS_1000G_AFR_GIANT.png'))
ggplot(dat_giant[group %in% amr],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "American Populations (1000G)") + scale_color_brewer(palette="PRGn")
ggsave(paste0(path, 'PGS_1000G_AMR_GIANT.png'))
ggplot(dat_giant[group %in% sas],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "South Asian Populations (1000G)") + scale_color_brewer(palette="PRGn")
ggsave(paste0(path, 'PGS_1000G_SAS_GIANT.png'))
ggplot(dat_giant[group %in% eas],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "East Asian Populations (1000G)") + scale_color_brewer(palette="PRGn")
ggsave(paste0(path,'PGS_1000G_EAS_GIANT.png'))
print('Done making figures')


##
my_plotA<-ggplot(dat[group %in% superpops],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "UK Biobank") + scale_color_brewer(palette="PRGn") + coord_cartesian(xlim=c(-9,3))
my_plotB<-ggplot(dat_giant[group %in% superpops],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title ="GIANT") + scale_color_brewer(palette="PRGn") + coord_cartesian(xlim=c(-9,3))
pdf('~/height_prediction/figs_for_paper/figs/SM_S1.pdf', width=9, height=7)
plot_grid(my_plotB, my_plotA, labels=c("A","B"), nrow=2)
dev.off()

