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
names(PGS)<-names(allpops)
pops<-c(unique(samples$pop))
superpops<-unique(samples$super_pop)

for (P in allpops){
	PGS[[P]]<-vector('list',22)
        names(PGS[[P]])<-c(1:22)
}

for (P in pops){
	cat(paste0("Population is ", P, '\n'))
        for (CR in 1:22){
		print(paste0("Chromosome is ", CR))
                PolScore_1000G(POP='CEU', CHR=CR, superpop=F)-> PGS[[P]][[CR]]
                cat(paste0(CR, '  done\n'))
        }
        cat (paste0(P, ' done\n'))
        gc()
}

for (P in superpops){
	cat(paste0("Super-population is ", P, '\n'))
        for(CR in 1:22){
		print(paste0("Chromosome is ", CR))
                PolScore_1000G(CHR=CR, POP='CEU', superpop=T)-> PGS[[P]][[CR]]
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

saveRDS(PGS2, file=paste0('~/heigth_prediction/gwas/ukb_eur/output/PGS_1000G_phys_100000_0.0005.Rds'))


#*******************************
#* Prepare data for plotting ***
#*******************************
path<-paste0('~/height_prediction/figs_for_paper/figs/')
PGS<-readRDS('~/heigth_prediction/gwas/ukb_eur/output/PGS_1000G_phys_100000_0.0005.Rds')

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
