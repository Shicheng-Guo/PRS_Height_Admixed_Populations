#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
#Load packages
library("optparse")
library(data.table)
library(dplyr)
library(ggplot2);library(reshape2); library(wesanderson)
library(rlist)
args<-"LD_prun"
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
#* Read in Polygenic Scores ******
#*********************************
if(args[1]=='small' | args[1]=='LD_prun'){
	paste0('/project/mathilab/bbita/gwas_admix/height/ukbiobank/PGS_1000G_',args[1],'.Rds')-> sm
	readRDS(sm)-> PGS
} else if (args[1]=='ALL'){
	readRDS('/project/mathilab/bbita/gwas_admix/height/ukbiobank/PGS_1000G_big_1_2.Rds')-> chr1_2
        readRDS('/project/mathilab/bbita/gwas_admix/height/ukbiobank/PGS_1000G_big_3_4.Rds')-> chr3_4
        readRDS('/project/mathilab/bbita/gwas_admix/height/ukbiobank/PGS_1000G_big_5_7.Rds')-> chr5_7
        readRDS('/project/mathilab/bbita/gwas_admix/height/ukbiobank/PGS_1000G_big_8_10.Rds')->chr8_10
        readRDS('/project/mathilab/bbita/gwas_admix/height/ukbiobank/PGS_1000G_big_11_13.Rds')-> chr11_13
        readRDS('/project/mathilab/bbita/gwas_admix/height/ukbiobank/PGS_1000G_big_14_16.Rds')-> chr14_16
        readRDS('/project/mathilab/bbita/gwas_admix/height/ukbiobank/PGS_1000G_big_17_19.Rds')-> chr17_19
        readRDS('/project/mathilab/bbita/gwas_admix/height/ukbiobank/PGS_1000G_big_20_22.Rds')-> chr20_22
	PGS<-vector('list', length(allpops))
	names(PGS)<-allpops
	for(P in allpops){
		samps<-names(chr20_22[[P]][[21]])	
		vector('list', length(samps))-> PGS[[P]]
		names(PGS[[P]])<-samps
		for (S in samps){
       		chr1_2[[P]][[1]][[S]]+chr1_2[[P]][[2]][[S]]+ chr3_4[[P]][[3]][[S]] + chr3_4[[P]][[4]][[S]]+chr5_7[[P]][[5]][[S]]+chr5_7[[P]][[6]][[S]]+chr5_7[[P]][[7]][[S]]+chr8_10[[P]][[8]][[S]]+chr8_10[[P]][[9]][[S]]+chr8_10[[P]][[10]][[S]]+chr11_13[[P]][[11]][[S]] + 
		chr11_13[[P]][[12]][[S]]+chr11_13[[P]][[13]][[S]]+chr14_16[[P]][[14]][[S]]+chr14_16[[P]][[15]][[S]]+chr14_16[[P]][[16]][[S]]+chr17_19[[P]][[17]][[S]]+chr17_19[[P]][[18]][[S]]+chr17_19[[P]][[19]][[S]]+chr20_22[[P]][[20]][[S]]+chr20_22[[P]][[21]][[S]]+chr20_22[[P]][[22]][[S]]->PGS[[P]][[S]]	
		}
	}
	saveRDS(PGS, file='/project/mathilab/bbita/gwas_admix/height/ukbiobank/PGS_1000G_big_all.Rds')
}

#*******************************
#* Prepare data for plotting ***
#*******************************
if(args[1]=='small' | args[1]=='ALL'){
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
} else if(args[1]=='LD_prun'){
	lapply(PGS, function(X) lapply(X, function(Y) unlist(Y)))-> PGS2
	names(PGS2)<-names(PGS)
	dat <- lapply(PGS2, function(Q) lapply(Q, function(Z) cbind(Z = seq_along(Z), y = Z)))
	list.names <- lapply(dat, function(X) names(X))
	lns <- lapply(dat, function(X) sapply(X, nrow))
	lapply(dat, function(X) as.data.table(do.call(rbind, X)))-> dat
	for(I in 1:length(dat)){
		dat[[I]]$group <- rep(list.names[[I]], lns[[I]])
		as.factor(dat[[I]]$group)-> dat[[I]]$group
		rbind(dat[[I]][group %in% afr], dat[[I]][group %in% eur], dat[[I]][group %in% sas], dat[[I]][group %in% eas], dat[[I]][group %in% amr], dat[[I]][group %in% superpops])-> dat[[I]]
		dat[[I]]$group<-factor(dat[[I]]$group, levels=allpops)
		mean(dat[[I]][group=='CEU'][,y])->m1
		sd(dat[[I]][group=='CEU'][,y])->sd1
		dat[[I]][, y1:=scale(as.matrix(y), center=m1, scale=sd1), by=group]
	}
}
#***************
#* Plotting ***
#***************
if(args[1]=='small' | args[1]=='ALL'){
	try(system(paste0('mkdir -m 777 /project/mathilab/bbita/gwas_admix/height/ukbiobank/', args[1], '/figs/')))
	path<-paste0('/project/mathilab/bbita/gwas_admix/height/ukbiobank/', args[1], '/figs/')
	ggplot(dat[group %in% pops],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "All populations (1000G)") + scale_color_brewer(palette="PRGn")
	ggsave(paste0(path,'PGS_1000G_pops.png'))
	ggplot(dat[group %in% superpops],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "Superpopulations (1000G)") + scale_color_brewer(palette="PRGn")
	ggsave(paste0(path, 'PGS_1000G_superpops.png'))
	ggplot(dat[group %in% c(afr,eur)],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "African and European Populations (1000G)") + scale_color_brewer(palette="PRGn")
	ggsave(paste0(path,'PGS_1000G_EUR_AFR.png'))
	ggplot(dat[group %in% eur],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "European Populations (1000G)") + scale_color_brewer(palette="PRGn")
	ggsave(paste0(path,'PGS_1000G_EUR.png'))
	ggplot(dat[group %in% afr],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "African Populations (1000G)") + scale_color_brewer(palette="PRGn")
	ggsave(paste0(path, 'PGS_1000G_AFR.png'))
	ggplot(dat[group %in% amr],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "American Populations (1000G)") + scale_color_brewer(palette="PRGn")
	ggsave(paste0(path, 'PGS_1000G_AMR.png'))
	ggplot(dat[group %in% sas],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "South Asian Populations (1000G)") + scale_color_brewer(palette="PRGn")
	ggsave(paste0(path, 'PGS_1000G_SAS.png'))
	ggplot(dat[group %in% eas],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "East Asian Populations (1000G)") + scale_color_brewer(palette="PRGn")
	ggsave(paste0(path,'PGS_1000G_EAS.png'))
	print('Done making figures')
} else if(args[1]=='LD_prun'){
	for(I in names(PGS2)){
		path<-paste0('/project/mathilab/bbita/gwas_admix/height/ukbiobank/figs/', I, '/')
		try(system(paste0('mkdir -m 777 /project/mathilab/bbita/gwas_admix/height/ukbiobank/figs/')))
		try(system(paste0('mkdir -m 777 ', path)))	
		ggplot(dat[[I]][group %in% pops],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "All populations (1000G)") + scale_color_brewer(palette="PRGn")
		ggsave(paste0(path,'PGS_1000G_pops.png'))
		ggplot(dat[[I]][group %in% superpops],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "Superpopulations (1000G)") + scale_color_brewer(palette="PRGn")
		ggsave(paste0(path, 'PGS_1000G_superpops.png'))
		ggplot(dat[[I]][group %in% c(afr,eur)],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "African and European Populations (1000G)") + scale_color_brewer(palette="PRGn")
		ggsave(paste0(path,'PGS_1000G_EUR_AFR.png'))
		ggplot(dat[[I]][group %in% eur],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "European Populations (1000G)") + scale_color_brewer(palette="PRGn")
		ggsave(paste0(path,'PGS_1000G_EUR.png'))
		ggplot(dat[[I]][group %in% afr],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "African Populations (1000G)") + scale_color_brewer(palette="PRGn")
		ggsave(paste0(path, 'PGS_1000G_AFR.png'))
		ggplot(dat[[I]][group %in% amr],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "American Populations (1000G)") + scale_color_brewer(palette="PRGn")
		ggsave(paste0(path, 'PGS_1000G_AMR.png'))
		ggplot(dat[[I]][group %in% sas],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "South Asian Populations (1000G)") + scale_color_brewer(palette="PRGn")
		ggsave(paste0(path, 'PGS_1000G_SAS.png'))
		ggplot(dat[[I]][group %in% eas],aes(x=y1, fill=group)) + geom_density(alpha=0.30) + labs(x = "Polygenic Score") + labs(title = "East Asian Populations (1000G)") + scale_color_brewer(palette="PRGn")
		ggsave(paste0(path,'PGS_1000G_EAS.png'))
		print('Done making figures')
	}	
}

#***************
# Other Plots **
#***************

if(args[1]=='LD_prun'){
	data.table(diff_EUR_AFR=unlist(lapply(dat, function(X) mean(X[group=='EUR'][,y1])-mean(X[group=='AFR'][,y1]))), Set=names(dat))-> dt
	dt[grep("phys",dt$Set),]->dt_phys
	dt[grep("genet",dt$Set),]->dt_genet
	dt[grep("LD",dt$Set),]->dt_LD
#use code from fig5 from NCD paper: https://github.com/bbitarello/NCV_dir_package/blob/master/scripts/enrich_plot.R  
# https://raw.githubusercontent.com/bbitarello/NCD-Statistics/master/Figures_main/Fig5.tiff
		factor(dt$Set)-> dt$Set
		factor(dt_phys$Set)-> dt_phys$Set
		factor(dt_genet$Set)-> dt_genet$Set
		factor(dt_LD$Set)-> dt_LD$Set
		factor(dt_LD$Set, levels(dt_LD$Set)[c(4,5,1,2,3)])-> dt_LD$Set
		factor(dt_genet$Set, levels(dt_genet$Set)[c(5,4,3,2,1,10,9,8,7,6,15,14,13,12,11,20,19,18,17,16,25,24,23,22,21)])-> dt_genet$Set
		factor(dt_phys$Set, levels(dt_phys$Set)[c(25,23,21,19,17,10,8,6,4,2,15,14,13,12,11,24,22,20,18,16,30,29,28,27,26, 9,7,5,3,1)])-> dt_phys$Set
		h=0
		a="E(EUR)-E(AFR)"
		ggplot(dt,aes(x=Set, y=diff_EUR_AFR, colour=Set)) + geom_point(size=2.5, alpha = 5/10) + geom_hline(yintercept=h, color="darkgray", linetype="dashed") + theme(legend.position="none",axis.text.x=element_text(angle=45, hjust=1))
		ggsave("figs/Set_comps.png")
		ggplot(dt_phys,aes(x=Set, y=diff_EUR_AFR, colour=Set)) + geom_point(size=2.5, alpha = 5/10) + geom_hline(yintercept=h, color="darkgray", linetype="dashed") + theme(legend.position="none",axis.text.x=element_text(angle=45, hjust=1))
		ggsave("figs/Set_comps_phys.png")
		ggplot(dt_genet,aes(x=Set, y=diff_EUR_AFR, colour=Set)) + geom_point(size=2.5, alpha = 5/10) + geom_hline(yintercept=h, color="darkgray", linetype="dashed") + theme(legend.position="none",axis.text.x=element_text(angle=45, hjust=1))
                ggsave("figs/Set_comps_genet.png")
		ggplot(dt_LD,aes(x=Set, y=diff_EUR_AFR, colour=Set)) + geom_point(size=2.5, alpha = 5/10) + geom_hline(yintercept=h, color="darkgray", linetype="dashed") + theme(legend.position="none",axis.text.x=element_text(angle=45, hjust=1))
                ggsave("figs/Set_comps_LD.png")
}

#add another figure, with diff_EUR_AFR vs Nr.SNPs:




dt[,Set:=as.character(Set)]
fwrite(data.table(Set=c("|---", dt$Set), diff_EUR_AFR=c("---|", dt$diff_EUR_AFR)), file='test_table_2.txt', col.names=T, quote=F, sep="|")

system("paste test_table.txt test_table_2.txt > test_table_3.txt")
system('cat test_table_3.txt >> README.md')
#******
#END *
#*****
