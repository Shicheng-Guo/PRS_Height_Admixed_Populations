#!/usr/bin/env Rscript
############################
#args = commandArgs(trailingOnly=TRUE)
#if (length(args)==0) {
#  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
##preable
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
options(scipen=999)
library(TeachingDemos)
args<-'gwas'

txtStart(paste0("~/height_prediction/", args, "/outout.txt"))
# Your code
##############################################################
#combine all

##
readRDS(paste0('~/height_prediction/', args[1], '/WHI/output/results.WHI.Rds'))-> results.WHI
readRDS(paste0('~/height_prediction/',args[1], '/JHS/output/results.JHS.Rds'))-> results.JHS
readRDS(paste0('~/height_prediction/',args[1], '/ukb_afr/output/results.UKB_afr.Rds'))-> results.UKB_afr
readRDS(paste0('~/height_prediction/',args[1], '/HRS_eur/output/results.HRS_eur.Rds'))-> results.HRS_eur
readRDS(paste0('~/height_prediction/',args[1], '/HRS_afr/output/results.HRS_afr.Rds'))-> results.HRS_afr
readRDS(paste0('~/height_prediction/',args[1], '/WHI/output/B_WHI.Rds'))->B_WHI
readRDS(paste0('~/height_prediction/',args[1], '/JHS/output/B_JHS.Rds'))->B_JHS
readRDS(paste0('~/height_prediction/',args[1], '/ukb_afr/output/B_UKB_afr.Rds'))-> B_UKB_afr
readRDS(paste0('~/height_prediction/',args[1], '/HRS_eur/output/B_HRS_eur.Rds'))-> B_HRS_eur
readRDS(paste0('~/height_prediction/',args[1],'/HRS_afr/output/B_HRS_afr.Rds'))-> B_HRS_afr
readRDS(paste0('~/height_prediction/',args[1],'/ukb_afr/output/PGS3_UKB_afr.Rds'))-> PGS3_UKB_afr
readRDS(paste0('~/height_prediction/',args[1],'/WHI/output/PGS3_WHI.Rds'))-> PGS3_WHI
readRDS(paste0('~/height_prediction/',args[1],'/JHS/output/PGS3_JHS.Rds'))-> PGS3_JHS
readRDS(paste0('~/height_prediction/',args[1],'/HRS_eur/output/PGS3_HRS_eur.Rds'))-> PGS3_HRS_eur
readRDS(paste0('~/height_prediction/',args[1],'/HRS_afr/output/PGS3_HRS_afr.Rds'))-> PGS3_HRS_afr

if(args[1]=='sib_betas'){
	readRDS(paste0('~/height_prediction/', args[1], '/ukb_eur/output/results.UKB_eur.Rds'))-> results.UKB_eur
	readRDS(paste0('~/height_prediction/', args[1], '/ukb_eur/output/B_UKB_eur.Rds'))-> B_UKB_eur
	readRDS(paste0('~/height_prediction/', args[1], '/ukb_eur/output/PGS3_UKB_eur.Rds'))-> PGS3_UKB_eur
	lapply(B_UKB_eur, function(X) X[,Dataset:="UKB_eur"])
}
for(I in names(B_JHS)){ #JHS lacks the LD prunning methods
	if(args[1]=='sib_betas'){
		ALL<-rbind(B_JHS[[I]][1:2,][, Dataset:='JHS_afr'], B_WHI[[I]][1:4,][, Dataset:='WHI_afr'], B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'], B_UKB_eur[[I]], B_HRS_eur[[I]])
		ALL$Dataset<-factor(ALL$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "UKB_eur", "HRS_eur"))
		my_colrs<-c(brewer.pal(5, 'Set1'),"#101010")
	} else{
	ALL<-rbind(B_JHS[[I]][1:2,][, Dataset:='JHS_afr'], B_WHI[[I]][1:4,][, Dataset:='WHI_afr'], B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_HRS_eur[[I]])
	ALL$Dataset<-factor(ALL$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "HRS_eur"))
	my_colrs<-c(brewer.pal(4, 'Set1'),"#101010")
	}
	my_plot<-ggplot(ALL, aes(x=Med_Eur_Anc, y=R_sq, colour=Dataset, shape=Dataset)) +
	geom_point(size=1.5, fill="white", alpha=0.8) +
	geom_errorbar(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset,ymin=boots_perc_L, ymax=boots_perc_U), width=0.05, size=0.8) +
	geom_errorbarh(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset, xmin=HVB_L, xmax=HVB_U), width=0.05, size=0.5) +  scale_color_manual(values=my_colrs) +
	ylab(expression(paste("Partial R"^"2"))) + xlab("European Ancestry Proportion") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),legend.key=element_blank(),legend.background=element_blank(),legend.title=element_blank(), legend.text=element_text(size=12))
	print(my_plot)
	if(args=='sib_betas'){
	ggsave(paste0('~/height_prediction/figs/sib_error_bars_all_v2_', I, '.png'))
	} else{
	ggsave(paste0('~/height_prediction/figs/error_bars_all_v2_', I, '.png'))
	}
}
#
ALL2<-vector('list', length(names(B_JHS)))
names(ALL2)<-names(B_JHS)
for(I in names(B_JHS)){
	if(args[1]=='sib_betas'){
		ALL2[[I]]<-rbind(B_JHS[[I]][1:2,][, Dataset:='JHS_afr'], B_WHI[[I]][1:4,][, Dataset:='WHI_afr'], B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_UKB_eur[[I]], B_HRS_eur[[I]])
		tmp<-1/c(var(results.JHS[[I]][[1]]$t), var(results.JHS[[I]][[2]]$t), var(results.WHI[[I]][[1]]$t), var(results.WHI[[I]][[2]]$t), var(results.WHI[[I]][[3]]$t), var(results.WHI[[I]][[4]]$t), var(results.UKB_afr[[I]][[1]]$t),var(results.UKB_afr[[I]][[2]]$t), var(results.UKB_afr[[I]][[3]]$t), var(results.UKB_afr[[I]][[4]]$t), var(results.HRS_afr[[I]][[1]]$t), var(results.HRS_afr[[I]][[2]]$t), var(results.UKB_eur[[I]]$t),var(results.HRS_eur[[I]]$t))
		my_colrs<-c(brewer.pal(5, 'Set1'),"#101010")
		ALL2[[I]]$Dataset<-factor(ALL2[[I]]$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "UKB_eur","HRS_eur"))
	} else{
		ALL2[[I]]<-rbind(B_JHS[[I]][1:2,][, Dataset:='JHS_afr'], B_WHI[[I]][1:4,][, Dataset:='WHI_afr'], B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_HRS_eur[[I]])
		tmp<-1/c(var(results.JHS[[I]][[1]]$t), var(results.JHS[[I]][[2]]$t), var(results.WHI[[I]][[1]]$t), var(results.WHI[[I]][[2]]$t), var(results.WHI[[I]][[3]]$t), var(results.WHI[[I]][[4]]$t), var(results.UKB_afr[[I]][[1]]$t),var(results.UKB_afr[[I]][[2]]$t), var(results.UKB_afr[[I]][[3]]$t), var(results.UKB_afr[[I]][[4]]$t), var(results.HRS_afr[[I]][[1]]$t), var(results.HRS_afr[[I]][[2]]$t), var(results.HRS_eur[[I]]$t))  #weighing lm by boostrap replicates.
		my_colrs<-c(brewer.pal(4, 'Set1'),"#101010")
		ALL2[[I]]$Dataset<-factor(ALL2[[I]]$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "HRS_eur"))
	}
	cbind(ALL2[[I]], W=tmp)-> ALL2[[I]]
	my_plot2<-ggplot(ALL2[[I]], aes(x=Med_Eur_Anc, y=R_sq)) +
	geom_point(size=1.5, shape=21, fill="white", alpha=0.8) + stat_smooth(method = "lm", mapping = aes(weight = W), col='black') +
	scale_color_manual(values=my_colrs) +
	ylab(expression(paste("Partial R"^"2"))) + xlab("European Ancestry Proportion") +

	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.key=element_blank(), legend.background=element_blank(),legend.title=element_blank(), legend.text=element_text(size=12))
	print(my_plot2)
	if(args[1]=='sib_betas'){
	ggsave(paste0('~/height_prediction/figs/sib_error_bars_all_v3_', I, '.png'))
	} else{
	ggsave(paste0('~/height_prediction/figs/error_bars_all_v3_', I, '.png'))
	}
	cat(I, 'done\n')
}

a<-data.table(Name=names(B_JHS), Intercept=unlist(lapply(1:80, function(I) coef(lm(R_sq~Med_Eur_Anc, weights= W, data=ALL2[[I]]))[[1]])), Slope=unlist(lapply(1:80, function(I) coef(lm(R_sq~Med_Eur_Anc, weights= W, data=ALL2[[I]]))[[2]])), R_sq=unlist(lapply(1:80, function(I) summary(lm(R_sq~Med_Eur_Anc, weights= W, data=ALL2[[I]]))[9])), P=unlist(lapply(1:80, function(I) summary(lm(R_sq~Med_Eur_Anc, weights= W, data=ALL2[[I]]))$coefficients[8])))
a[Name=='phys_100000_0.0005']

fwrite(a, file=paste0("~/height_prediction/figs_for_paper/figs/SM_Table1_", args[1], ".txt", sep=","))
ALL2b<-vector('list', length(names(B_JHS)))
#add ancestry
names(ALL2b)<-names(B_JHS)

for(I in names(B_JHS)){
	if(args[1]=='sib_betas'){
       	ALL2b[[I]]<-rbind(B_JHS[[I]][1:2,][, Dataset:='JHS_afr'], B_WHI[[I]][1:4,][, Dataset:='WHI_afr'], B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_UKB_eur[[I]],B_HRS_eur[[I]])
              tmp<-1/c(var(results.JHS[[I]][[1]]$t), var(results.JHS[[I]][[2]]$t), var(results.WHI[[I]][[1]]$t), var(results.WHI[[I]][[2]]$t), var(results.WHI[[I]][[3]]$t), var(results.WHI[[I]][[4]]$t), var(results.UKB_afr[[I]][[1]]$t),var(results.UKB_afr[[I]][[2]]$t), var(results.UKB_afr[[I]][[3]]$t), var(results.UKB_afr[[I]][[4]]$t), var(results.HRS_afr[[I]][[1]]$t), var(results.HRS_afr[[I]][[2]]$t),  var(results.UKB_eur[[I]]$t),var(results.HRS_eur[[I]]$t))
		ALL2b[[I]]$Dataset<-factor(ALL2b[[I]]$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr", "UKB_eur", "HRS_eur"))
		my_colrs<-c(brewer.pal(5, 'Set1'),"#101010")
	} else{
        	ALL2b[[I]]<-rbind(B_JHS[[I]][1:2,][, Dataset:='JHS_afr'], B_WHI[[I]][1:4,][, Dataset:='WHI_afr'], B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_HRS_eur[[I]])
                tmp<-1/c(var(results.JHS[[I]][[1]]$t), var(results.JHS[[I]][[2]]$t), var(results.WHI[[I]][[1]]$t), var(results.WHI[[I]][[2]]$t), var(results.WHI[[I]][[3]]$t), var(results.WHI[[I]][[4]]$t), var(results.UKB_afr[[I]][[1]]$t),var(results.UKB_afr[[I]][[2]]$t), var(results.UKB_afr[[I]][[3]]$t), var(results.UKB_afr[[I]][[4]]$t), var(results.HRS_afr[[I]][[1]]$t), var(results.HRS_afr[[I]][[2]]$t), var(results.HRS_eur[[I]]$t))  #weighing lm by boostrap replicates.
		ALL2b[[I]]$Dataset<-factor(ALL2b[[I]]$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "HRS_eur"))
		my_colrs<-c(brewer.pal(4, 'Set1'),"#101010")
        }
        cbind(ALL2b[[I]], W=tmp)-> ALL2b[[I]]
	my_plot2<-ggplot(ALL2b[[I]], aes(x=Med_Eur_Anc, y=R_sq)) +
                geom_point(size=1.5, shape=21, fill="white", alpha=0.8) + stat_smooth(method = "lm", mapping = aes(weight = W), col='black') +
                scale_color_manual(values=my_colrs) +
		ylab(expression(paste("Partial R"^"2"))) + xlab("European Ancestry Proportion") +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.key=element_blank(), legend.background=element_blank(),legend.title=element_blank(), legend.text=element_text(size=12))
        print(my_plot2)
	if(args[1]=='sib_betas'){
	ggsave(paste0('~/height_prediction/figs/sib_error_bars_all_v3b_', I, '.png'))	
	} else{
        ggsave(paste0('~/height_prediction/figs/error_bars_all_v3b_', I, '.png'))
	}
	cat(I, 'done\n')
}

for(I in names(B_JHS)){
	if(args[1]=='sib_betas'){
	ALL2[[I]]$Dataset<-factor(ALL2[[I]]$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "UKB_eur","HRS_eur"))
	my_colrs<-c(brewer.pal(5, 'Set1'),"#101010")
	} else{
	ALL2[[I]]$Dataset<-factor(ALL2[[I]]$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "HRS_eur"))
	my_colrs<-c(brewer.pal(4, 'Set1'),"#101010")
	}
        my_plot<-ggplot(ALL2[[I]], aes(x=Med_Eur_Anc, y=R_sq,colour=Dataset)) +
        geom_point(aes(shape=Dataset), size=1.5, fill="white", alpha=0.8) + stat_smooth(data=ALL2[[I]],method = "lm", mapping = aes(weight = W), col='black') +
        geom_errorbar(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset,ymin=boots_perc_L, ymax=boots_perc_U), width=0.05, size=0.8) +
        #geom_line(color='lightgray')+
        geom_errorbarh(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset, xmin=HVB_L, xmax=HVB_U), width=0.05, size=0.5) + 
	scale_color_manual(values=my_colrs) + 
        ylab(expression(paste("Partial R"^"2"))) + xlab("European Ancestry Proportion") +
  	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.key=element_blank(), legend.background=element_blank(), legend.title=element_blank(), legend.text=element_text(size=15),legend.position = c(0.15,0.85)) 
	print(my_plot)
	if(args[1]=='sib_betas'){
	ggsave(paste0('~/height_prediction/figs/sib_error_bars_all_v4_', I, '.png'))
	} else{
        ggsave(paste0('~/height_prediction/figs/error_bars_all_v4_', I, '.png'))
	}
	cat(I, 'done\n')
}
#
lm(ALL2[[63]]$R_sq ~ ALL2[[63]]$Med_Eur_Anc, weights=ALL2[[63]]$W)
summary(lm( ALL2[[63]]$R_sq ~ ALL2[[63]]$Med_Eur_Anc, weights=ALL2[[63]]$W)) #p-value: 0.006, adj-r2=0.466 for sibs;0.0000001534, adj-r2=0.92 for gwas
#stop here 04/09/2019
ALL3<-vector('list', length(names(B_JHS)))
names(ALL3)<- names(B_JHS)

for (I in names(B_JHS)){
	if(args[1]=='sib_betas'){
		ALL3[[I]]<-rbind(B_JHS[[I]][1:2,][, Dataset:='JHS_afr'], B_WHI[[I]][1:4,][, Dataset:='WHI_afr'], B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_HRS_eur[[I]], B_UKB_eur[[I]])
	 tmp<-lm(R_sq~Med_Eur_Anc,weights=1/
        c(var(results.JHS[[I]][[1]]$t), var(results.JHS[[I]][[2]]$t), var(results.WHI[[I]][[1]]$t), var(results.WHI[[I]][[2]]$t), var(results.WHI[[I]][[3]]$t), var(results.WHI[[I]][[4]]$t), var(results.UKB_afr[[I]][[1]]$t),var(results.UKB_afr[[I]][[2]]$t), var(results.UKB_afr[[I]][[3]]$t), var(results.UKB_afr[[I]][[4]]$t), var(results.HRS_afr[[I]][[1]]$t), var(results.HRS_afr[[I]][[2]]$t),var(results.UKB_eur[[I]]$t),var(results.HRS_eur[[I]]$t)), data=ALL3[[I]])
	} else{
                ALL3[[I]]<-rbind(B_JHS[[I]][1:2,][, Dataset:='JHS_afr'], B_WHI[[I]][1:4,][, Dataset:='WHI_afr'], B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_HRS_eur[[I]])
	        tmp<-lm(R_sq~Med_Eur_Anc,weights=1/
        c(var(results.JHS[[I]][[1]]$t), var(results.JHS[[I]][[2]]$t), var(results.WHI[[I]][[1]]$t), var(results.WHI[[I]][[2]]$t), var(results.WHI[[I]][[3]]$t), var(results.WHI[[I]][[4]]$t), var(results.UKB_afr[[I]][[1]]$t),var(results.UKB_afr[[I]][[2]]$t), var(results.UKB_afr[[I]][[3]]$t), var(results.UKB_afr[[I]][[4]]$t), var(results.HRS_afr[[I]][[1]]$t), var(results.HRS_afr[[I]][[2]]$t),var(results.HRS_eur[[I]]$t)), data=ALL3[[I]])
	}
	ALL3[[I]][,Set:=I]
	readRDS(paste0('~/height_prediction/', args[1],'/WHI/output/Nr_SNPs_WHI.Rds'))[Name==I][, Nr]->a
	readRDS(paste0('~/height_prediction/', args[1],'/ukb_afr/output/Nr_SNPs_UKB_afr.Rds'))[Name==I][, Nr]->b
	readRDS(paste0('~/height_prediction/', args[1],'/JHS/output/Nr_SNPs_JHS.Rds'))[Name==I][, Nr]->d
	readRDS(paste0('~/height_prediction/', args[1],'/HRS_eur/output/Nr_SNPs_HRS.Rds'))[Name==I][, Nr]->f
	#readRDS('Nr_SNPs_pennBB_afr.Rds')[Name==I][, Nr]->f
	#readRDS('Nr_SNPs_pennBB_eur.Rds')[Name==I][, Nr]->g
	ALL3[[I]][,Intercept:=coef(tmp)[[1]]][,Slope:=coef(tmp)[[2]]]
	ALL3[[I]][,Slope_Intercept:=sum(coef(tmp))]
	ALL3[[I]][, Nr_SNPs_WHI:=a]
	ALL3[[I]][, Nr_SNPs_UKB:=b]
	ALL3[[I]][, Nr_SNPs_JHS:=d]
	ALL3[[I]][, Nr_SNPs_HRS_eur:=f]
	ALL3[[I]][, Nr_SNPs_HRS_afr:=f]
	cat(I, ' \n')
}
#if(args[1]=='gwas'){
do.call(rbind,ALL3)[,.(Quant,Set,Intercept,Slope_Intercept, Slope, Nr_SNPs_WHI, Nr_SNPs_HRS_eur, R_sq, Med_Eur_Anc)]->ALL4
#add nr of snps

ALL4[grep("phys",  ALL4$Set),][,.(Set,Intercept,Slope_Intercept)]->dt_phys
ALL4[grep("genet", ALL4$Set),][,.(Set,Intercept,Slope_Intercept)]->dt_genet
ALL4[grep("LD",    ALL4$Set),][,.(Set,Intercept,Slope_Intercept)]->dt_LD

#factor(dt$Set)-> dt$Setp
factor(dt_phys$Set)-> dt_phys$Set
factor(dt_genet$Set)-> dt_genet$Set
factor(dt_LD$Set)-> dt_LD$Set
factor(dt_LD$Set, levels(dt_LD$Set)[c(4,5,3,1,2)])-> dt_LD$Set
factor(dt_genet$Set, levels(dt_genet$Set)[c(5,4,3,2,1,10,9,8,7,6,15,14,13,12,11,20,19,18,17,16,25,24,23,22,21, 30, 29, 28, 27, 26, 35,34,33,32,31)])-> dt_genet$Set
factor(dt_phys$Set, levels(dt_phys$Set)[c(24,27,30,33,35,4,7,10,13,15,16,17,18,19,20,22,25,28,31,34,36,37,38,39,40,2,5,8,11,14,21,23,26,29,32,1,3,6,9,12)])-> dt_phys$Set

melt(dt_LD)-> dt_LD
rbind(dt_LD[grep("_250000_", dt_LD$Set)][, window:=250000], dt_LD[grep("_100000_", dt_LD$Set)][, window:=100000], dt_LD[grep("_50000_", dt_LD$Set)][, window:=50000], dt_LD[grep("block", dt_LD$Set)][, window:="-"])-> dt_LD
dt_LD[, method:=gsub("LD_50000_0.01_0.5", "LD_0.01_0.5", gsub("LD_100000_0.01_0.5", "LD_0.01_0.5", gsub("LD_250000_0.01_0.5", "LD_0.01_0.5", gsub("LD_block_0_0_AFR", "LD_block_AFR", gsub("LD_block_0_0_EUR", "LD_block_EUR", dt_LD[,Set])))))]

as.factor(dt_LD$window)-> dt_LD$window
factor(dt_LD$window, levels(dt_LD$window)[c(1,4,2,3)])-> dt_LD$window

ggplot(dt_LD,aes(x=method, y=value, colour=window, shape=variable)) + geom_point(size=2.5, alpha=1) + theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
ggsave('~/height_prediction/figs/reg_rsq_eur_anc_LD.png')
#
#
melt(dt_genet)->dt_genet
dt_genet[, c("method", "window","p") := tstrsplit(Set, "_")]
as.factor(dt_genet$window)-> dt_genet$window
ggplot() + geom_point(data=dt_genet,aes(x=p, y=value, colour=window, shape=variable), size=2.5, alpha = 1) + theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12)) 
ggsave('~/height_prediction/figs/reg_rsq_eur_anc_genet.png')
#
#
melt(dt_phys)-> dt_phys
dt_phys[, c("method", "window","p") := tstrsplit(Set, "_")]
as.factor(dt_phys$window)-> dt_phys$window
factor(dt_phys$window, levels(dt_phys$window)[c(3,7,2,8,6,4,1,5)])-> dt_phys$window
ggplot() + geom_point(data=dt_phys,aes(x=p, y=value, colour=window, shape=variable), size=2.5, alpha = 1) + theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
ggsave('~/height_prediction/figs/reg_rsq_eur_anc_phys.png')

#plot as function of Nr_SNPs

data.table(Set=unique(ALL4$Set), Nr_SNPs=unique(ALL4$Nr_SNPs_HRS_eur))->A1

ALL4[grep("phys",  ALL4$Set),][,.(Set,Intercept,Slope, Slope_Intercept, R_sq, Med_Eur_Anc)][Med_Eur_Anc==1][, Med_Eur_Anc:=NULL]->dt_phys
ALL4[grep("genet", ALL4$Set),][,.(Set,Intercept,Slope,Slope_Intercept, R_sq, Med_Eur_Anc)][Med_Eur_Anc==1][, Med_Eur_Anc:=NULL]->dt_genet
ALL4[grep("LD",    ALL4$Set),][,.(Set,Intercept,Slope, Slope_Intercept, R_sq, Med_Eur_Anc)][Med_Eur_Anc==1][, Med_Eur_Anc:=NULL]->dt_LD
dt_LD %>% dplyr::group_by(Set) %>% dplyr::mutate(R_sq=min(R_sq)) %>% as.data.table-> dt_LD #remove ukb_eur
dt_phys %>% dplyr::group_by(Set) %>% dplyr::mutate(R_sq=min(R_sq)) %>% as.data.table-> dt_phys
dt_genet %>% dplyr::group_by(Set) %>% dplyr::mutate(R_sq=min(R_sq)) %>% as.data.table-> dt_genet

#factor(dt$Set)-> dt$Setp
factor(dt_phys$Set)-> dt_phys$Set
factor(dt_genet$Set)-> dt_genet$Set
factor(dt_LD$Set)-> dt_LD$Set
factor(dt_LD$Set, levels(dt_LD$Set)[c(4,5,3,1,2)])-> dt_LD$Set
factor(dt_genet$Set, levels(dt_genet$Set)[c(5,4,3,2,1,10,9,8,7,6,15,14,13,12,11,20,19,18,17,16,25,24,23,22,21, 30, 29, 28, 27, 26, 35,34,33,32,31)])-> dt_genet$Set
factor(dt_phys$Set, levels(dt_phys$Set)[c(24,27,30,33,35,4,7,10,13,15,16,17,18,19,20,22,25,28,31,34,36,37,38,39,40,2,5,8,11,14,21,23,26,29,32,1,3,6,9,12)])-> dt_phys$Set


melt(dt_LD)-> dt_LD
rbind(dt_LD[grep("_250000_", dt_LD$Set)][, window:=250000], dt_LD[grep("_100000_", dt_LD$Set)][, window:=100000], dt_LD[grep("_50000_", dt_LD$Set)][, window:=50000], dt_LD[grep("block", dt_LD$Set)][, window:="-"])-> dt_LD
dt_LD[, method:=gsub("LD_50000_0.01_0.5", "LD_0.01_0.5", gsub("LD_100000_0.01_0.5", "LD_0.01_0.5", gsub("LD_250000_0.01_0.5", "LD_0.01_0.5", gsub("LD_block_0_0_AFR", "LD_block_AFR", gsub("LD_block_0_0_EUR", "LD_block_EUR", dt_LD[,Set])))))]

as.factor(dt_LD$window)-> dt_LD$window
factor(dt_LD$window, levels(dt_LD$window)[c(1,4,2,3)])-> dt_LD$window

setkey(dt_LD, Set)
dt_LD[A1,nomatch=0]-> dt_LD

melt(dt_phys)-> dt_phys
dt_phys[, c("method", "window","p") := tstrsplit(Set, "_")]
as.factor(dt_phys$window)-> dt_phys$window
factor(dt_phys$window, levels(dt_phys$window)[c(3,7,2,8,6,4,1,5)])-> dt_phys$window

setkey(A1, Set)
setkey(dt_phys, Set)
dt_phys[A1,nomatch=0]-> dt_phys

melt(dt_genet)->dt_genet
dt_genet[, c("method", "window","p") := tstrsplit(Set, "_")]
as.factor(dt_genet$window)-> dt_genet$window

setkey(dt_genet, Set)
dt_genet[A1,nomatch=0]-> dt_genet

#
gsub("Slope_Intercept", "Slope+Intercept",dt_phys$variable)-> dt_phys$variable
gsub("R_sq", "Partial R-squared",dt_phys$variable)-> dt_phys$variable
dt_phys[variable!="Slope"]-> dt_phys
plotA<-ggplot() + geom_point(data=dt_phys,aes(x=Nr_SNPs, y=value, colour=window, shape=variable), size=1.5, alpha = 0.7) + 
theme(axis.title.y = element_text(size = 12), axis.title.x=element_text(size=12), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=9), legend.title = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
guides(shape=FALSE,color=guide_legend(override.aes=list(shape=15))) + ylab("Value")
#ggsave('~/height_prediction/figs/reg_rsq_eur_anc_phys_v2.png')
gsub("Slope_Intercept", "Slope+Intercept",dt_genet$variable)-> dt_genet$variable
gsub("R_sq", "Partial R-squared",dt_genet$variable)-> dt_genet$variable
dt_genet[variable!="Slope"]-> dt_genet
plotB<-ggplot() + geom_point(data=dt_genet,aes(x=Nr_SNPs, y=value, colour=window, shape=variable), size=1.5, alpha = 0.7) + 
theme(axis.title.y = element_text(size = 12), axis.title.x=element_text(size=12), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=9), legend.title = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
guides(color=guide_legend(override.aes=list(shape=15))) + ylab("Value")
#ggsave('~/height_prediction/figs/reg_rsq_eur_anc_genet_v2.png')
gsub("Slope_Intercept", "Slope+Intercept",dt_LD$variable)-> dt_LD$variable
gsub("R_sq", "Partial R-squared",dt_LD$variable)-> dt_LD$variable
dt_LD[variable!="Slope"]-> dt_LD
plotC<-ggplot() + geom_point(data=dt_LD,aes(x=Nr_SNPs, y=value, colour=Set, shape=variable), size=1.5, alpha = 0.7) + 
theme(axis.title.y = element_text(size = 12), axis.title.x=element_text(size=12), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=9), legend.title = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape=FALSE,color=guide_legend(override.aes=list(shape=15)))  + ylab("Value")
#ggsave('~/height_prediction/figs/reg_rsq_eur_anc_LD_v2.png')

png('~/height_prediction/figs/panel_Nr_snps.png', units="in", height=11, width=7, res=600)
plot_grid(plotC, plotB,plotA,  labels = c("A", "B", "C"), nrow=3, align="v")
dev.off()

#try to use just one plot
dt_LD[, p:=0.01]
dt_LD[, .(Set, variable, value, method, window,p, Nr_SNPs)]-> dt_LD
rbind(dt_phys,dt_genet, dt_LD)-> dt


cat('STOP HERE STOP HERE\n')
ggplot() + geom_point(data=dt,aes(x=Nr_SNPs, y=value, colour=method, shape=p), size=1.2, alpha = 0.7) + theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text=element_text(size=12)) + facet_wrap(~variable, nrow=2, scales='free_y')
ggsave('~/height_prediction/figs/test.png')

ggplot(dt, aes(x=Nr_SNPs, y=value, colour=method, shape=p)) + geom_point(size=1.2, alpha = 0.7) + scale_shape_manual(values=c(16,3,15,0,17,12)) + geom_line(alpha=0.4) + theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text=element_text(size=12)) +  facet_wrap(~variable, nrow=2, scales='free_y')
ggsave('~/height_prediction/figs/test2.png')

ggplot(dt[variable=='R_sq'], aes(x=Nr_SNPs, y=value, colour=method, shape=p)) + geom_point(size=1.2, alpha = 0.7) + scale_shape_manual(values=c(16,3,15,0,17,12)) + geom_line(alpha=0.4) + theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text=element_text(size=12))
ggsave('~/height_prediction/figs/test2b.png')

ggplot(dt, aes(x=Nr_SNPs, y=value, colour=window, shape=method)) + geom_point(size=1.2, alpha = 0.7) + scale_shape_manual(values=c(16,3,13,0,17,12)) + geom_line(alpha=0.4) + theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text=element_text(size=12)) +  facet_wrap(~variable, nrow=2, scales='free_y')
ggsave('~/height_prediction/figs/test3.png')


#ggplot(dt, aes(x=Nr_SNPs, y=value, shape=p, colour=int)) + geom_point(size=1.2, alpha = 0.7) + geom_line(alpha=0.4) +  facet_wrap(~variable, nrow=2, scales='free_y')
#ggsave('figs/test4.png')

###
combo<-vector('list', length(PGS3_JHS))
names(combo)<-names(PGS3_JHS)

for (I in names(PGS3_JHS)){
	rbind(PGS3_WHI[[I]][,.(SUBJID,AGE, age2, HEIGHTX,PGS, SEX,EUR_ANC)][,SUBJ_ID:=SUBJID][, age:=AGE][, sex:='FEMALE'][, SUBJID:=NULL][, AGE:=NULL][, SEX:=NULL][,Dataset:='WHI_afr'][, Res.Height:=resid(lm(HEIGHTX~age+age2))],
	PGS3_HRS_afr[[I]][,SUBJ_ID:=ID][, age:=AGE][, age2:=AGE2][, HEIGHTX:=HEIGHT][, sex:=SEX][,.(SUBJ_ID,age, age2, HEIGHTX,PGS, sex, EUR_ANC)][,Dataset:='HRS_afr'][, Res.Height:=resid(lm(HEIGHTX~sex+age+age2))],
	PGS3_HRS_eur[[I]][,SUBJ_ID:=ID][, age:=AGE][, age2:=AGE2][, HEIGHTX:=HEIGHT][, sex:=SEX][,.(SUBJ_ID,age, age2, HEIGHTX,PGS, sex, EUR_ANC)][,Dataset:='HRS_eur'][, Res.Height:=resid(lm(HEIGHTX~sex+age+age2))],
	PGS3_UKB_afr[[I]][, sex:=Sex][, HEIGHTX:=Height][,SUBJ_ID:=ID][,.(SUBJ_ID, Age, age2, HEIGHTX, PGS, sex, EUR_ANC)][, age:=Age][, Age:=NULL][,Dataset:='UKB_afr'][, Res.Height:=resid(lm(HEIGHTX~sex+age+age2))],
	PGS3_JHS[[I]][,.(SUBJID, age, age2, HEIGHTX, PGS, sex, EUR_ANC)][,SUBJ_ID:=SUBJID][, SUBJID:=NULL][,Dataset:='JHS_afr'][, Res.Height:=resid(lm(HEIGHTX~sex+age+age2))])[,Prun_Set:=I][,Nr_SNPs:=0]-> combo[[I]]
	readRDS(paste0('~/height_prediction/', args[1],'/WHI/output/Nr_SNPs_WHI.Rds'))[Name==I][, Nr]->a
        readRDS(paste0('~/height_prediction/', args[1],'/ukb_afr/output/Nr_SNPs_UKB_afr.Rds'))[Name==I][, Nr]->b
        readRDS(paste0('~/height_prediction/', args[1],'/JHS/output/Nr_SNPs_JHS.Rds'))[Name==I][, Nr]->d
        readRDS(paste0('~/height_prediction/', args[1],'/HRS_eur/output/Nr_SNPs_HRS.Rds'))[Name==I][, Nr]->f
	combo[[I]][Dataset=='WHI_afr']$Nr_SNPs<-a
	combo[[I]][Dataset=='JHS_afr']$Nr_SNPs<-d
	combo[[I]][Dataset=='UKB_afr']$Nr_SNPs<-b
	combo[[I]][Dataset=='HRS_afr']$Nr_SNPs<-f
	combo[[I]][Dataset=='HRS_eur']$Nr_SNPs<-f
	combo[[I]][, Std.PRS:=scale(PGS), by=Dataset]-> combo[[I]]
	cat(I, ' done\n')
}
do.call(rbind,combo)-> combo2
#combo2[Dataset=='pennBB_EA', EUR_ANC:=1]
fwrite(combo2, file="~/height_prediction/output/WHI_JHS_UKB_HRS_summ.txt", sep="\t")
#gzip('WHI_pennBB_UKB_summ.txt',destname='WHI_pennBB_UKB_summ.txt.gz')
###
#Define F statistc	

F_x<-function(x,ind_index, data){
#p(x|y)
y<-1-x
z<-quantile(data[[ind_index]][,Res.Height], probs=x)[[1]]
z2<-quantile(data[[ind_index]][,PGS], probs=x)[[1]]
a<-nrow(data[[ind_index]][Std.PRS>=z2])
b<-nrow(data[[ind_index]][Std.PRS>=z2][Res.Height>=z])
res<-(b/a)/y
res<-res
return(res)
}

F2_x<-function(x, data){
#p(x|y)
y<-1-x
z<-quantile(data[,Res.Height], probs=x)[[1]]
z2<-quantile(data[,Std.PRS], probs=x)[[1]]
a<-nrow(data[Std.PRS>=z2])
b<-nrow(data[Std.PRS>=z2][Res.Height>=z])
res<-((b+0.5)/(a+0.5))/y #add 0.5 to b to avoid 0 in numerator. Haldane-Ascombe correction
#res<-log(res)
return(res)
}

F3_x<-function(x, data, data2){
#p(x|y)
y<-1-x
z<-quantile(data[,Res.Height], probs=x)[[1]]
z2<-quantile(data2[Dataset=='HRS_eur'][,Std.PRS], probs=x)[[1]]
a<-nrow(data[Std.PRS>=z2])
tp<-nrow(data[Std.PRS>=z2][Res.Height>=z])
fp<-nrow(data[Std.PRS<z2][Res.Height>=z])
tn<-nrow(data[Std.PRS<z2][Res.Height<z])
fn<-nrow(data[Std.PRS>=z2][Res.Height<z])
res<-list(F=((tp+0.5)/(a+0.5))/y, G=((fp+0.5)/(a+0.5))/y, f=tp/(tp+fn), g=fp/(fp+tn)) #add 0.5 to b to avoid 0 in numerator. Haldane-Ascombe correction
#res<-log(res)
return(res)
}

F4_x<-function(x, data){
#p(x|y)
y<-1-x
z<-quantile(data[,Res.Height], probs=x)[[1]]
z2<-quantile(data[,Std.PRS], probs=x)[[1]]
a<-nrow(data[Std.PRS>=z2])
tp<-nrow(data[Std.PRS>=z2][Res.Height>=z])
fp<-nrow(data[Std.PRS<z2][Res.Height>=z])
tn<-nrow(data[Std.PRS<z2][Res.Height<z])
fn<-nrow(data[Std.PRS>=z2][Res.Height<z])
res<-list(F=((tp+0.5)/(a+0.5))/y, G=((fp+0.5)/(a+0.5))/y, f=tp/(tp+fn), g=fp/(fp+tn)) #add 0.5 to b to avoid 0 in numerator. Haldane-Ascombe correction
#res<-log(res)
return(res)
}

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F3_X=F3_x(x=X,data=combo[[Y]][Dataset=='WHI_afr'],  data2=combo[[Y]]), Quantile=X))))-> AA
names(AA)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F2_X=F2_x(x=X,data=combo[[Y]][Dataset=='WHI_afr']), Quantile=X))))-> AAA
names(AAA)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F3_X=F3_x(x=X,data=combo[[Y]][Dataset=='JHS_afr'],  data2=combo[[Y]]), Quantile=X))))-> AJ
names(AJ)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F2_X=F2_x(x=X,data=combo[[Y]][Dataset=='JHS_afr']), Quantile=X))))-> AAJ
names(AAJ)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F3_X=F3_x(x=X,data=combo[[Y]][Dataset=='UKB_afr'],  data2=combo[[Y]]), Quantile=X))))-> AU
names(AU)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F2_X=F2_x(x=X,data=combo[[Y]][Dataset=='UKB_afr']), Quantile=X))))-> AAU
names(AAU)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F3_X=F3_x(x=X,data=combo[[Y]][Dataset=='HRS_eur'],  data2=combo[[Y]]), Quantile=X))))-> AP
names(AP)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F2_X=F2_x(x=X,data=combo[[Y]][Dataset=='HRS_eur']), Quantile=X))))-> AAP
names(AAP)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F3_X=F3_x(x=X,data=combo[[Y]][Dataset=='HRS_afr'],  data2=combo[[Y]]), Quantile=X))))-> APE
names(APE)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F2_X=F2_x(x=X,data=combo[[Y]][Dataset=='HRS_afr']), Quantile=X))))-> AAPE
names(AAPE)<-names(B_JHS)


I<-names(AA)[63]

png(paste0('~/height_prediction/figs/OR_WHI_', I,  ".png"))
ggplot(AA[[I]], aes(x=Quantile, y=F3_X.F)) +
geom_point(size=2) + labs(title="Odds ratio of P(>=Xth HEIGHT quantile|>=Xth PRS quantile)", y="OR") + geom_hline(yintercept=1, linetype="dashed", color = "orange") +
theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text=element_text(size=12))
dev.off()

tmp1<-AA[[I]]
tmp2<-AAA[[I]]
tmp3<-AJ[[I]]
tmp4<-AAJ[[I]]
tmp5<-AU[[I]]
tmp6<-AAU[[I]]
tmp7<-AP[[I]]
tmp8<-AAP[[I]]
tmp9<-APE[[I]]
tmp10<-AAPE[[I]]
#
setDT(tmp1)
tmp1[,OR:=F3_X.F]
setDT(tmp3)
tmp3[,OR:=F3_X.F]
setDT(tmp5)
tmp5[,OR:=F3_X.F]
setDT(tmp7)
tmp7[,OR:=F3_X.F]
setDT(tmp9)
tmp9[,OR:=F3_X.F]
tmp1[, F3_X.F:=NULL][, F3_X.G:=NULL][,F3_X.f:=NULL][,F3_X.g:=NULL][, Group:="WHI_afr_Eur"]
tmp3[, F3_X.F:=NULL][, F3_X.G:=NULL][,F3_X.f:=NULL][,F3_X.g:=NULL][, Group:="JHS_afr_Eur"]
tmp5[, F3_X.F:=NULL][, F3_X.G:=NULL][,F3_X.f:=NULL][,F3_X.g:=NULL][, Group:="UKB_AFR_Eur"]
tmp7[,F3_X.F:=NULL][, F3_X.G:=NULL][,F3_X.f:=NULL][,F3_X.g:=NULL][, Group:="HRS_afr_Eur"]
tmp9[,F3_X.F:=NULL][, F3_X.G:=NULL][,F3_X.f:=NULL][,F3_X.g:=NULL][, Group:="HRS_eur_Eur"]
#
setDT(tmp2)
tmp2[,OR:=F2_X][, F2_X:=NULL][, Group:='WHI_afr_matched']
setDT(tmp4)
tmp4[,OR:=F2_X][, F2_X:=NULL][, Group:="JHS_afr_matched"]
setDT(tmp6)
tmp6[,OR:=F2_X][, F2_X:=NULL][, Group:="UKB_afr_matched"]
setDT(tmp8)
tmp8[,OR:=F2_X][, F2_X:=NULL][, Group:="HRS_afr_matched"]
setDT(tmp10)
tmp10[,OR:=F2_X][, F2_X:=NULL][,Group:="HRS_eur_matched"]


all<-rbind(tmp1,tmp2,tmp3,tmp4,tmp5, tmp6,tmp7, tmp8, tmp9, tmp10)
all[, Group2:=c(rep("WHI_afr", 100), rep("JHS_afr", 100), rep("UKB_afr", 100), rep("HRS_afr",100), rep("HRS_eur",100))]
all[, Group3:=c(rep("Eur", 50), rep("Matched",50), rep("Eur",50), rep("Matched",50), rep("Eur",50), rep("Matched", 50), rep("Eur",50), rep("Matched", 50), rep("Eur", 50),  rep("Matched", 50))]
all[, Prs.Quantile:=Quantile]

pdf(paste0('~/height_prediction/figs/OR_WHI_JHS_', I,  ".pdf"))
ggplot(all[Prs.Quantile<=0.975], aes(x=Prs.Quantile, y=OR, group=Group, color=Group2, linetype=Group3)) + geom_line() +
geom_point(size=0.7) + labs(y="OR", x="Quantile") +
theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text=element_text(size=12))
dev.off()


pdf(paste0('~/height_prediction/figs/OR_WHI_JHS_', I,  "_v2.pdf"))
ggplot(all[Prs.Quantile<=0.975], aes(x=Prs.Quantile, y=OR, group=Group, color=Group2, linetype=Group3)) + geom_smooth() +
geom_point(size=0.7) + labs(y="OR", x="Quantile") + 
theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text=element_text(size=12))
dev.off()
#lapply(1:length(names(A)),function(X) combo[[X]][,Quantile:= cut(EUR_ANC,breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.2), na.rm=TRUE),include.lowest=TRUE)])


WOW<-vector('list', length(combo));test<-vector('list', length(combo))
names(WOW)<-names(B_WHI);names(test)<-names(B_WHI)
for(I in names(AA)){
	a<-combo[[I]][Dataset=='WHI_afr'][,Quantile:= cut(EUR_ANC,breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.25), na.rm=TRUE),include.lowest=TRUE)]
	b<-combo[[I]][Dataset=='HRS_afr'][,Quantile:= cut(EUR_ANC,breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.5), na.rm=TRUE),include.lowest=TRUE)]
	d<-combo[[I]][Dataset=='UKB_afr'][,Quantile:= cut(EUR_ANC,breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.25), na.rm=TRUE),include.lowest=TRUE)]
	f<-combo[[I]][Dataset=='HRS_eur'][,EUR_ANC:=1][,Quantile:=as.factor(as.character(1))]
	g<-combo[[I]][Dataset=='JHS_afr'][,Quantile:= cut(EUR_ANC,breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.5), na.rm=TRUE),include.lowest=TRUE)]
	test[[I]]<-rbind(a,b,d,f,g);test[[I]]$Quantile<-as.character(test[[I]]$Quantile)
	WOW[[I]]<-vector('list', 5);names(WOW[[I]])<- c('WHI_afr', 'HRS_afr','UKB_afr', 'HRS_eur', 'JHS_afr')
	lapply(levels(a$Quantile), function(X) F2_x(data=a[Quantile==X], x=0.975))-> WOW[[I]][['WHI_afr']][[1]];list.append(WOW[[I]][['WHI_afr']][[1]], F2_x(data=a, x=0.975))-> WOW[[I]][['WHI_afr']][[1]]	
	lapply(levels(a$Quantile), function(X) F2_x(data=a[Quantile==X], x=0.95))-> WOW[[I]][['WHI_afr']][[2]];list.append(WOW[[I]][['WHI_afr']][[2]], F2_x(data=a, x=0.95))-> WOW[[I]][['WHI_afr']][[2]]
	lapply(levels(a$Quantile), function(X) F2_x(data=a[Quantile==X], x=0.9))-> WOW[[I]][['WHI_afr']][[3]];list.append(WOW[[I]][['WHI_afr']][[3]], F2_x(data=a, x=0.9))-> WOW[[I]][['WHI_afr']][[3]]
	lapply(levels(a$Quantile), function(X) F2_x(data=a[Quantile==X], x=0.85))-> WOW[[I]][['WHI_afr']][[4]];list.append(WOW[[I]][['WHI_afr']][[4]], F2_x(data=a, x=0.85))-> WOW[[I]][['WHI_afr']][[4]]
	lapply(levels(a$Quantile), function(X) F2_x(data=a[Quantile==X], x=0.8))-> WOW[[I]][['WHI_afr']][[5]];list.append(WOW[[I]][['WHI_afr']][[5]], F2_x(data=a, x=0.8))-> WOW[[I]][['WHI_afr']][[5]]
	names(WOW[[I]][['WHI_afr']])<-c('0.975perc', '0.95perc', '0.9perc', '0.85perc', '0.8perc');names(WOW[[I]][['WHI_afr']][['0.9perc']])<-c(levels(a$Quantile), 'all')
	names(WOW[[I]][['WHI_afr']][['0.975perc']])<-c(levels(a$Quantile), 'all');names(WOW[[I]][['WHI_afr']][['0.95perc']])<-c(levels(a$Quantile), 'all')
	names(WOW[[I]][['WHI_afr']][['0.85perc']])<-c(levels(a$Quantile), 'all');names(WOW[[I]][['WHI_afr']][['0.8perc']])<-c(levels(a$Quantile),'all') 
	lapply(levels(b$Quantile), function(X) F2_x(data=b[Quantile==X], x=0.975))-> WOW[[I]][['HRS_afr']][[1]];list.append(WOW[[I]][['HRS_afr']][[1]], F2_x(data=b, x=0.975))-> WOW[[I]][['HRS_afr']][[1]]
	lapply(levels(b$Quantile), function(X) F2_x(data=b[Quantile==X], x=0.95))-> WOW[[I]][['HRS_afr']][[2]];list.append(WOW[[I]][['HRS_afr']][[2]], F2_x(data=b, x=0.95))-> WOW[[I]][['HRS_afr']][[2]]
	lapply(levels(b$Quantile), function(X) F2_x(data=b[Quantile==X], x=0.9))-> WOW[[I]][['HRS_afr']][[3]];list.append(WOW[[I]][['HRS_afr']][[3]], F2_x(data=b, x=0.9))-> WOW[[I]][['HRS_afr']][[3]]
	lapply(levels(b$Quantile), function(X) F2_x(data=b[Quantile==X], x=0.85))-> WOW[[I]][['HRS_afr']][[4]];list.append(WOW[[I]][['HRS_afr']][[4]], F2_x(data=b, x=0.85))-> WOW[[I]][['HRS_afr']][[4]]
        lapply(levels(b$Quantile), function(X) F2_x(data=b[Quantile==X], x=0.8))-> WOW[[I]][['HRS_afr']][[5]];list.append(WOW[[I]][['HRS_afr']][[5]], F2_x(data=b, x=0.8))-> WOW[[I]][['HRS_afr']][[5]]
	names(WOW[[I]][['HRS_afr']])<-c('0.975perc', '0.95perc', '0.9perc', '0.85perc', '0.8perc');names(WOW[[I]][['HRS_afr']][['0.9perc']])<-c(levels(b$Quantile), 'all')
	names(WOW[[I]][['HRS_afr']][['0.95perc']])<-c(levels(b$Quantile), 'all');names(WOW[[I]][['HRS_afr']][['0.975perc']])<-c(levels(b$Quantile), 'all')
	names(WOW[[I]][['HRS_afr']][['0.85perc']])<-c(levels(b$Quantile), 'all');names(WOW[[I]][['HRS_afr']][['0.8perc']])<-c(levels(b$Quantile), 'all')
	lapply(levels(d$Quantile), function(X) F2_x(data=d[Quantile==X], x=0.975))-> WOW[[I]][['UKB_afr']][[1]];list.append(WOW[[I]][['UKB_afr']][[1]], F2_x(data=d, x=0.975))-> WOW[[I]][['UKB_afr']][[1]]
	lapply(levels(d$Quantile), function(X) F2_x(data=d[Quantile==X], x=0.95))-> WOW[[I]][['UKB_afr']][[2]];list.append(WOW[[I]][['UKB_afr']][[2]], F2_x(data=d, x=0.95))-> WOW[[I]][['UKB_afr']][[2]]
	lapply(levels(d$Quantile), function(X) F2_x(data=d[Quantile==X], x=0.9))-> WOW[[I]][['UKB_afr']][[3]];list.append(WOW[[I]][['UKB_afr']][[3]], F2_x(data=d, x=0.9))-> WOW[[I]][['UKB_afr']][[3]]
	lapply(levels(d$Quantile), function(X) F2_x(data=d[Quantile==X], x=0.85))-> WOW[[I]][['UKB_afr']][[4]];list.append(WOW[[I]][['UKB_afr']][[4]], F2_x(data=d, x=0.85))-> WOW[[I]][['UKB_afr']][[4]]	
	lapply(levels(d$Quantile), function(X) F2_x(data=d[Quantile==X], x=0.8))-> WOW[[I]][['UKB_afr']][[5]];list.append(WOW[[I]][['UKB_afr']][[5]], F2_x(data=d, x=0.8))-> WOW[[I]][['UKB_afr']][[5]]
	names(WOW[[I]][['UKB_afr']])<-c('0.975perc', '0.95perc', '0.9perc','0.85perc', '0.8perc');names(WOW[[I]][['UKB_afr']][['0.9perc']])<-c(levels(d$Quantile), 'all')
	names(WOW[[I]][['UKB_afr']][['0.95perc']])<-c(levels(d$Quantile), 'all');names(WOW[[I]][['UKB_afr']][['0.975perc']])<-c(levels(d$Quantile), 'all')
	names(WOW[[I]][['UKB_afr']][['0.85perc']])<-c(levels(d$Quantile), 'all'); names(WOW[[I]][['UKB_afr']][['0.8perc']])<-c(levels(d$Quantile), 'all')
	lapply(unique(f$Quantile), function(X) F2_x(data=f[Quantile==X], x=0.975))-> WOW[[I]][['HRS_eur']][[1]]
	lapply(unique(f$Quantile), function(X) F2_x(data=f[Quantile==X], x=0.95))-> WOW[[I]][['HRS_eur']][[2]];lapply(unique(f$Quantile), function(X) F2_x(data=f[Quantile==X], x=0.9))-> WOW[[I]][['HRS_eur']][[3]]
	lapply(unique(f$Quantile), function(X) F2_x(data=f[Quantile==X], x=0.85))-> WOW[[I]][['HRS_eur']][[4]];lapply(unique(f$Quantile), function(X) F2_x(data=f[Quantile==X], x=0.8))-> WOW[[I]][['HRS_eur']][[5]]
	names(WOW[[I]][['HRS_eur']])<-c('0.975perc', '0.95perc', '0.9perc','0.85perc', '0.8perc');names(WOW[[I]][['HRS_eur']][['0.9perc']])<-unique(f$Quantile)
	names(WOW[[I]][['HRS_eur']][['0.95perc']])<-unique(f$Quantile);names(WOW[[I]][['HRS_eur']][['0.975perc']])<-unique(f$Quantile)
	names(WOW[[I]][['HRS_eur']][['0.85perc']])<-unique(f$Quantile);names(WOW[[I]][['HRS_eur']][['0.8perc']])<-unique(f$Quantile)
	lapply(levels(g$Quantile), function(X) F2_x(data=g[Quantile==X], x=0.975))-> WOW[[I]][['JHS_afr']][[1]];list.append(WOW[[I]][['JHS_afr']][[1]], F2_x(data=d, x=0.975))-> WOW[[I]][['JHS_afr']][[1]]
        lapply(levels(g$Quantile), function(X) F2_x(data=g[Quantile==X], x=0.95))-> WOW[[I]][['JHS_afr']][[2]];list.append(WOW[[I]][['JHS_afr']][[2]], F2_x(data=d, x=0.95))-> WOW[[I]][['JHS_afr']][[2]]
        lapply(levels(g$Quantile), function(X) F2_x(data=g[Quantile==X], x=0.9))-> WOW[[I]][['JHS_afr']][[3]];list.append(WOW[[I]][['JHS_afr']][[3]], F2_x(data=d, x=0.9))-> WOW[[I]][['JHS_afr']][[3]]
        lapply(levels(g$Quantile), function(X) F2_x(data=g[Quantile==X], x=0.85))-> WOW[[I]][['JHS_afr']][[4]];list.append(WOW[[I]][['JHS_afr']][[4]], F2_x(data=d, x=0.85))-> WOW[[I]][['JHS_afr']][[4]]
        lapply(levels(g$Quantile), function(X) F2_x(data=g[Quantile==X], x=0.8))-> WOW[[I]][['JHS_afr']][[5]];list.append(WOW[[I]][['JHS_afr']][[5]], F2_x(data=d, x=0.8))-> WOW[[I]][['JHS_afr']][[5]]
        names(WOW[[I]][['JHS_afr']])<-c('0.975perc', '0.95perc', '0.9perc','0.85perc', '0.8perc');names(WOW[[I]][['JHS_afr']][['0.9perc']])<-c(levels(g$Quantile), 'all')
        names(WOW[[I]][['JHS_afr']][['0.95perc']])<-c(levels(g$Quantile), 'all');names(WOW[[I]][['JHS_afr']][['0.975perc']])<-c(levels(g$Quantile), 'all');names(WOW[[I]][['JHS_afr']][['0.85perc']])<-c(levels(g$Quantile), 'all'); names(WOW[[I]][['JHS_afr']][['0.8perc']])<-c(levels(g$Quantile), 'all')
	cat(I);cat('\n')
}

WOW2<-vector('list', length(B_WHI));test2<-vector('list', length(B_WHI))
names(WOW2)<-names(B_WHI);names(test2)<-names(B_WHI)
for(I in names(AA)){
        a<-combo[[I]][Dataset=='WHI_afr'][,Quantile:= cut(EUR_ANC,breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.25), na.rm=TRUE),include.lowest=TRUE)]
        b<-combo[[I]][Dataset=='HRS_afr'][,Quantile:= cut(EUR_ANC,breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.5), na.rm=TRUE),include.lowest=TRUE)]
        d<-combo[[I]][Dataset=='UKB_afr'][,Quantile:= cut(EUR_ANC,breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.25), na.rm=TRUE),include.lowest=TRUE)]
	f<-combo[[I]][Dataset=='HRS_eur'][,EUR_ANC:=1][,Quantile:=as.factor(as.character(1))]
	g<-combo[[I]][Dataset=='JHS_afr'][,Quantile:= cut(EUR_ANC,breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.5), na.rm=TRUE),include.lowest=TRUE)]
        test2[[I]]<-rbind(a,b,d,f,g);test2[[I]]$Quantile<-as.character(test2[[I]]$Quantile)
        WOW2[[I]]<-vector('list', 5);names(WOW2[[I]])<- c('WHI_afr', 'HRS_afr','UKB_afr', 'HRS_eur', 'JHS_afr')
        lapply(levels(a$Quantile), function(X) F3_x(data=a[Quantile==X], x=0.975,  data2=combo[[I]]))-> WOW2[[I]][['WHI_afr']][[1]];list.append(WOW2[[I]][['WHI_afr']][[1]], F3_x(data=a, x=0.975,data2=combo[[I]]))-> WOW2[[I]][['WHI_afr']][[1]]
        lapply(levels(a$Quantile), function(X) F3_x(data=a[Quantile==X], x=0.95, data2=combo[[I]]))-> WOW2[[I]][['WHI_afr']][[2]];list.append(WOW2[[I]][['WHI_afr']][[2]], F3_x(data=a, x=0.95,data2=combo[[I]]))-> WOW2[[I]][['WHI_afr']][[2]]
        lapply(levels(a$Quantile), function(X) F3_x(data=a[Quantile==X], x=0.9, data2=combo[[I]]))-> WOW2[[I]][['WHI_afr']][[3]];list.append(WOW2[[I]][['WHI_afr']][[3]], F3_x(data=a, x=0.9,data2=combo[[I]]))-> WOW2[[I]][['WHI_afr']][[3]]
        lapply(levels(a$Quantile), function(X) F3_x(data=a[Quantile==X], x=0.85, data2=combo[[I]]))-> WOW2[[I]][['WHI_afr']][[4]];list.append(WOW2[[I]][['WHI_afr']][[4]], F3_x(data=a, x=0.85,data2=combo[[I]]))-> WOW2[[I]][['WHI_afr']][[4]]
	lapply(levels(a$Quantile), function(X) F3_x(data=a[Quantile==X], x=0.8, data2=combo[[I]]))-> WOW2[[I]][['WHI_afr']][[5]];list.append(WOW2[[I]][['WHI_afr']][[5]], F3_x(data=a, x=0.8,data2=combo[[I]]))-> WOW2[[I]][['WHI_afr']][[5]]
        names(WOW2[[I]][['WHI_afr']])<-c('0.975perc', '0.95perc', '0.9perc', '0.85perc', '0.8perc');names(WOW2[[I]][['WHI_afr']][['0.9perc']])<-c(levels(a$Quantile), 'all')
        names(WOW2[[I]][['WHI_afr']][['0.95perc']])<-c(levels(a$Quantile), 'all');names(WOW2[[I]][['WHI_afr']][['0.975perc']])<-c(levels(a$Quantile), 'all')
 	names(WOW2[[I]][['WHI_afr']][['0.85perc']])<-c(levels(a$Quantile), 'all');names(WOW2[[I]][['WHI_afr']][['0.8perc']])<-c(levels(a$Quantile), 'all')
        lapply(levels(b$Quantile), function(X) F3_x(data=b[Quantile==X], x=0.975,  data2=combo[[I]]))-> WOW2[[I]][['HRS_afr']][[1]];list.append(WOW2[[I]][['HRS_afr']][[1]], F3_x(data=b, x=0.975,data2=combo[[I]]))-> WOW2[[I]][['HRS_afr']][[1]]
        lapply(levels(b$Quantile), function(X) F3_x(data=b[Quantile==X], x=0.95, data2=combo[[I]]))-> WOW2[[I]][['HRS_afr']][[2]];list.append(WOW2[[I]][['HRS_afr']][[2]], F3_x(data=b, x=0.95,data2=combo[[I]]))-> WOW2[[I]][['HRS_afr']][[2]]
        lapply(levels(b$Quantile), function(X) F3_x(data=b[Quantile==X], x=0.9, data2=combo[[I]]))-> WOW2[[I]][['HRS_afr']][[3]];list.append(WOW2[[I]][['HRS_afr']][[3]], F3_x(data=b, x=0.9,data2=combo[[I]]))-> WOW2[[I]][['HRS_afr']][[3]]
        lapply(levels(b$Quantile), function(X) F3_x(data=b[Quantile==X], x=0.85, data2=combo[[I]]))-> WOW2[[I]][['HRS_afr']][[4]];list.append(WOW2[[I]][['HRS_afr']][[4]], F3_x(data=b, x=0.85,data2=combo[[I]]))-> WOW2[[I]][['HRS_afr']][[4]]
	lapply(levels(b$Quantile), function(X) F3_x(data=b[Quantile==X], x=0.8, data2=combo[[I]]))-> WOW2[[I]][['HRS_afr']][[5]];list.append(WOW2[[I]][['HRS_afr']][[5]], F3_x(data=b, x=0.8,data2=combo[[I]]))-> WOW2[[I]][['HRS_afr']][[5]]
        names(WOW2[[I]][['HRS_afr']])<-c('0.975perc', '0.95perc', '0.9perc', '0.85perc', '0.8perc')
        names(WOW2[[I]][['HRS_afr']][['0.9perc']])<-c(levels(b$Quantile), 'all')
        names(WOW2[[I]][['HRS_afr']][['0.95perc']])<-c(levels(b$Quantile), 'all')
        names(WOW2[[I]][['HRS_afr']][['0.975perc']])<-c(levels(b$Quantile), 'all')
        names(WOW2[[I]][['HRS_afr']][['0.85perc']])<-c(levels(b$Quantile), 'all')
	names(WOW2[[I]][['HRS_afr']][['0.8perc']])<-c(levels(b$Quantile), 'all')
        lapply(levels(d$Quantile), function(X) F3_x(data=d[Quantile==X], x=0.975,  data2=combo[[I]]))-> WOW2[[I]][['UKB_afr']][[1]];list.append(WOW2[[I]][['UKB_afr']][[1]], F3_x(data=d, x=0.975,data2=combo[[I]]))-> WOW2[[I]][['UKB_afr']][[1]]
        lapply(levels(d$Quantile), function(X) F3_x(data=d[Quantile==X], x=0.95, data2=combo[[I]]))-> WOW2[[I]][['UKB_afr']][[2]];list.append(WOW2[[I]][['UKB_afr']][[2]], F3_x(data=d, x=0.95,data2=combo[[I]]))-> WOW2[[I]][['UKB_afr']][[2]]
        lapply(levels(d$Quantile), function(X) F3_x(data=d[Quantile==X], x=0.9, data2=combo[[I]]))-> WOW2[[I]][['UKB_afr']][[3]];list.append(WOW2[[I]][['UKB_afr']][[3]], F3_x(data=d, x=0.9,data2=combo[[I]]))-> WOW2[[I]][['UKB_afr']][[3]]
        lapply(levels(d$Quantile), function(X) F3_x(data=d[Quantile==X], x=0.85, data2=combo[[I]]))-> WOW2[[I]][['UKB_afr']][[4]];list.append(WOW2[[I]][['UKB_afr']][[4]], F3_x(data=d, x=0.85,data2=combo[[I]]))-> WOW2[[I]][['UKB_afr']][[4]]
	lapply(levels(d$Quantile), function(X) F3_x(data=d[Quantile==X], x=0.8, data2=combo[[I]]))-> WOW2[[I]][['UKB_afr']][[5]];list.append(WOW2[[I]][['UKB_afr']][[5]], F3_x(data=d, x=0.8,data2=combo[[I]]))-> WOW2[[I]][['UKB_afr']][[5]]
        names(WOW2[[I]][['UKB_afr']])<-c('0.975perc', '0.95perc', '0.9perc','0.85perc',  '0.8perc')
        names(WOW2[[I]][['UKB_afr']][['0.9perc']])<-c(levels(d$Quantile), 'all')
        names(WOW2[[I]][['UKB_afr']][['0.95perc']])<-c(levels(d$Quantile), 'all')
        names(WOW2[[I]][['UKB_afr']][['0.975perc']])<-c(levels(d$Quantile), 'all')
        names(WOW2[[I]][['UKB_afr']][['0.85perc']])<-c(levels(d$Quantile), 'all')
	names(WOW2[[I]][['UKB_afr']][['0.8perc']])<-c(levels(d$Quantile), 'all')
        lapply(unique(f$Quantile), function(X) F3_x(data=f[Quantile==X], x=0.975,  data2=combo[[I]]))-> WOW2[[I]][['HRS_eur']][[1]];
        lapply(unique(f$Quantile), function(X) F3_x(data=f[Quantile==X], x=0.95, data2=combo[[I]]))-> WOW2[[I]][['HRS_eur']][[2]];
        lapply(unique(f$Quantile), function(X) F3_x(data=f[Quantile==X], x=0.9, data2=combo[[I]]))-> WOW2[[I]][['HRS_eur']][[3]];
        lapply(unique(f$Quantile), function(X) F3_x(data=f[Quantile==X], x=0.85, data2=combo[[I]]))-> WOW2[[I]][['HRS_eur']][[4]];
	lapply(unique(f$Quantile), function(X) F3_x(data=f[Quantile==X], x=0.8, data2=combo[[I]]))-> WOW2[[I]][['HRS_eur']][[5]];
        names(WOW2[[I]][['HRS_eur']])<-c('0.975perc', '0.95perc', '0.9perc','0.85perc', '0.8perc');names(WOW2[[I]][['HRS_eur']][['0.9perc']])<-unique(f$Quantile)
        names(WOW2[[I]][['HRS_eur']][['0.95perc']])<-unique(f$Quantile);names(WOW2[[I]][['HRS_eur']][['0.975perc']])<-unique(f$Quantile)
        names(WOW2[[I]][['HRS_eur']][['0.85perc']])<-unique(f$Quantile);names(WOW2[[I]][['HRS_eur']][['0.8perc']])<-unique(f$Quantile)
	lapply(levels(g$Quantile), function(X) F3_x(data=g[Quantile==X], x=0.975,data2=combo[[I]]))-> WOW2[[I]][['JHS_afr']][[1]];list.append(WOW2[[I]][['JHS_afr']][[1]], F3_x(data=d, x=0.975,data2=combo[[I]]))-> WOW2[[I]][['JHS_afr']][[1]]
        lapply(levels(g$Quantile), function(X) F3_x(data=g[Quantile==X], x=0.95,data2=combo[[I]]))-> WOW2[[I]][['JHS_afr']][[2]];list.append(WOW2[[I]][['JHS_afr']][[2]], F3_x(data=d, x=0.95,data2=combo[[I]]))-> WOW2[[I]][['JHS_afr']][[2]]
        lapply(levels(g$Quantile), function(X) F3_x(data=g[Quantile==X], x=0.9,data2=combo[[I]]))-> WOW2[[I]][['JHS_afr']][[3]];list.append(WOW2[[I]][['JHS_afr']][[3]], F3_x(data=d, x=0.9,data2=combo[[I]]))-> WOW2[[I]][['JHS_afr']][[3]]
        lapply(levels(g$Quantile), function(X) F3_x(data=g[Quantile==X], x=0.85,data2=combo[[I]]))-> WOW2[[I]][['JHS_afr']][[4]];list.append(WOW2[[I]][['JHS_afr']][[4]], F3_x(data=d, x=0.85,data2=combo[[I]]))-> WOW2[[I]][['JHS_afr']][[4]]
        lapply(levels(g$Quantile), function(X) F3_x(data=g[Quantile==X], x=0.8,data2=combo[[I]]))-> WOW2[[I]][['JHS_afr']][[5]];list.append(WOW2[[I]][['JHS_afr']][[5]], F3_x(data=d, x=0.8,data2=combo[[I]]))-> WOW2[[I]][['JHS_afr']][[5]]
        names(WOW2[[I]][['JHS_afr']])<-c('0.975perc', '0.95perc', '0.9perc','0.85perc', '0.8perc');names(WOW2[[I]][['JHS_afr']][['0.9perc']])<-c(levels(g$Quantile), 'all')
        names(WOW2[[I]][['JHS_afr']][['0.95perc']])<-c(levels(g$Quantile), 'all');names(WOW2[[I]][['JHS_afr']][['0.975perc']])<-c(levels(g$Quantile), 'all');names(WOW2[[I]][['JHS_afr']][['0.85perc']])<-c(levels(g$Quantile), 'all'); names(WOW2[[I]][['JHS_afr']][['0.8perc']])<-c(levels(g$Quantile), 'all')
        cat(I)
        cat('\n')
}

OR_table<-vector('list', length(combo))
names(OR_table)<-names(B_WHI)
for(I in names(AA)){
	OR_table[[I]]<-rbind(
	data.table(Dataset=c(rep("WHI_afr",5),rep("HRS_afr",3), rep('UKB_afr', 5), rep('HRS_eur',1), rep("JHS_afr",3)),
	EUR_ANC=c(na.omit(c(unlist(lapply(names(WOW[[I]][['WHI_afr']][['0.975perc']]), function(X) median(test[[I]][Dataset=='WHI_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='WHI_afr'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW[[I]][['HRS_afr']][['0.975perc']]), function(X) median(test[[I]][Dataset=='HRS_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='HRS_afr'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW[[I]][['UKB_afr']][['0.975perc']]), function(X) median(test[[I]][Dataset=='UKB_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='UKB_afr'][, EUR_ANC]))), 
	median(test[[I]][Dataset=='HRS_eur'][Quantile==1][, EUR_ANC]), 
	na.omit(c(unlist(lapply(names(WOW[[I]][['JHS_afr']][['0.975perc']]), function(X) median(test[[I]][Dataset=='JHS_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='JHS_afr'][, EUR_ANC])))), 
	Alpha=0.975, 
	OR=c(unlist(WOW[[I]][['WHI_afr']][['0.975perc']]), unlist(WOW[[I]][['HRS_afr']][['0.975perc']]), unlist(WOW[[I]][['UKB_afr']][['0.975perc']]),  unlist(WOW[[I]][['HRS_eur']][['0.975perc']]), unlist(WOW[[I]][['JHS_afr']][['0.975perc']])),
	Quantile=c(names(WOW[[I]][['WHI_afr']][['0.975perc']]), names(WOW[[I]][['HRS_afr']][['0.975perc']]),names(WOW[[I]][['UKB_afr']][['0.975perc']]),as.character(1), names(WOW[[I]][['JHS_afr']][['0.975perc']])), Prun=I),
	data.table(
	Dataset=c(rep("WHI_afr",5),rep("HRS_afr",3), rep('UKB_afr', 5), rep('HRS_eur',1),rep("JHS_afr",3)),
	EUR_ANC=c(na.omit(c(unlist(lapply(names(WOW[[I]][['WHI_afr']][['0.95perc']]), function(X) median(test[[I]][Dataset=='WHI_afr'][Quantile==X][, EUR_ANC]))),median(combo[[I]][Dataset=='WHI_afr'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW[[I]][['HRS_afr']][['0.95perc']]), function(X) median(test[[I]][Dataset=='HRS_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='HRS_afr'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW[[I]][['UKB_afr']][['0.95perc']]), function(X) median(test[[I]][Dataset=='UKB_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='UKB_afr'][, EUR_ANC]))),
	median(test[[I]][Dataset=='HRS_eur'][Quantile==1][, EUR_ANC]),
	na.omit(c(unlist(lapply(names(WOW[[I]][['JHS_afr']][['0.95perc']]), function(X) median(test[[I]][Dataset=='JHS_afr'][Quantile==X][, EUR_ANC]))),median(combo[[I]][Dataset=='JHS_afr'][, EUR_ANC])))),
	Alpha=0.95,
	OR=c(unlist(WOW[[I]][['WHI_afr']][['0.95perc']]), unlist(WOW[[I]][['HRS_afr']][['0.95perc']]), unlist(WOW[[I]][['UKB_afr']][['0.95perc']]), unlist(WOW[[I]][['HRS_eur']][['0.95perc']]), unlist(WOW[[I]][['JHS_afr']][['0.95perc']])),
	Quantile=c(names(WOW[[I]][['WHI_afr']][['0.95perc']]), names(WOW[[I]][['HRS_afr']][['0.95perc']]),names(WOW[[I]][['UKB_afr']][['0.95perc']]), as.character(1), names(WOW[[I]][['JHS_afr']][['0.95perc']])),Prun=I),
	data.table(
	Dataset=c(rep("WHI",5),rep("HRS_afr",3), rep('UKB_AFR', 5), rep('HRS_eur',1), rep("JHS_afr",3)),
	EUR_ANC=c(na.omit(c(unlist(lapply(names(WOW[[I]][['WHI_afr']][['0.9perc']]), function(X) median(test[[I]][Dataset=='WHI_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='WHI_afr'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW[[I]][['HRS_afr']][['0.9perc']]), function(X) median(test[[I]][Dataset=='HRS_afr'][Quantile==X][, EUR_ANC]))),median(combo[[I]][Dataset=='HRS_afr'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW[[I]][['UKB_afr']][['0.9perc']]), function(X) median(test[[I]][Dataset=='UKB_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='UKB_afr'][, EUR_ANC]))),
	median(test[[I]][Dataset=='HRS_eur'][Quantile==1][, EUR_ANC]),
	na.omit(c(unlist(lapply(names(WOW[[I]][['JHS_afr']][['0.9perc']]), function(X) median(test[[I]][Dataset=='JHS_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='JHS_afr'][, EUR_ANC])))),
	Alpha=0.9,
	OR=c(unlist(WOW[[I]][['WHI_afr']][['0.9perc']]), unlist(WOW[[I]][['HRS_afr']][['0.9perc']]), unlist(WOW[[I]][['UKB_afr']][['0.9perc']]), unlist(WOW[[I]][['HRS_eur']][['0.9perc']]),unlist(WOW[[I]][['JHS_afr']][['0.9perc']])),
	Quantile=c(names(WOW[[I]][['WHI_afr']][['0.9perc']]), names(WOW[[I]][['HRS_afr']][['0.9perc']]),names(WOW[[I]][['UKB_afr']][['0.9perc']]), as.character(1),names(WOW[[I]][['JHS_afr']][['0.9perc']])), Prun=I),
	data.table(
        Dataset=c(rep("WHI_afr",5),rep("HRS_afr",3), rep('UKB_afr', 5), rep('HRS_eur',1), rep("JHS_afr",3)),
        EUR_ANC=c(na.omit(c(unlist(lapply(names(WOW[[I]][['WHI_afr']][['0.85perc']]), function(X) median(test[[I]][Dataset=='WHI_afr'][Quantile==X][, EUR_ANC]))),  median(combo[[I]][Dataset=='WHI_afr'][, EUR_ANC]))),
        na.omit(c(unlist(lapply(names(WOW[[I]][['HRS_afr']][['0.85perc']]), function(X) median(test[[I]][Dataset=='HRS_afr'][Quantile==X][, EUR_ANC]))),median(combo[[I]][Dataset=='HRS_afr'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW[[I]][['UKB_afr']][['0.85perc']]), function(X) median(test[[I]][Dataset=='UKB_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='UKB_afr'][, EUR_ANC]))),
        median(test[[I]][Dataset=='HRS_eur'][Quantile==1][, EUR_ANC]),
	na.omit(c(unlist(lapply(names(WOW[[I]][['JHS_afr']][['0.85perc']]), function(X) median(test[[I]][Dataset=='JHS_afr'][Quantile==X][, EUR_ANC]))),  median(combo[[I]][Dataset=='JHS_afr'][, EUR_ANC])))),
        Alpha=0.85,
        OR=c(unlist(WOW[[I]][['WHI_afr']][['0.85perc']]), unlist(WOW[[I]][['HRS_afr']][['0.85perc']]), unlist(WOW[[I]][['UKB_afr']][['0.85perc']]), unlist(WOW[[I]][['HRS_eur']][['0.85perc']]), unlist(WOW[[I]][['JHS_afr']][['0.85perc']])),
        Quantile=c(names(WOW[[I]][['WHI_afr']][['0.85perc']]), names(WOW[[I]][['HRS_afr']][['0.85perc']]),names(WOW[[I]][['UKB_afr']][['0.85perc']]), as.character(1), names(WOW[[I]][['JHS_afr']][['0.85perc']])),Prun=I),
        data.table(
        Dataset=c(rep("WHI",5),rep("HRS_afr",3), rep('UKB_afr', 5),  rep('HRS_eur',1), rep("JHS_afr",3)),
        EUR_ANC=c(na.omit(c(unlist(lapply(names(WOW[[I]][['WHI_afr']][['0.8perc']]), function(X) median(test[[I]][Dataset=='WHI_afr'][Quantile==X][, EUR_ANC]))),median(combo[[I]][Dataset=='WHI_afr'][, EUR_ANC]))),
        na.omit(c(unlist(lapply(names(WOW[[I]][['HRS_afr']][['0.8perc']]), function(X) median(test[[I]][Dataset=='HRS_afr'][Quantile==X][, EUR_ANC]))),median(combo[[I]][Dataset=='HRS_afr'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW[[I]][['UKB_afr']][['0.8perc']]), function(X) median(test[[I]][Dataset=='UKB_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='UKB_afr'][, EUR_ANC]))),
        median(test[[I]][Dataset=='HRS_eur'][Quantile==1][, EUR_ANC]),
	na.omit(c(unlist(lapply(names(WOW[[I]][['JHS_afr']][['0.8perc']]), function(X) median(test[[I]][Dataset=='JHS_afr'][Quantile==X][, EUR_ANC]))),median(combo[[I]][Dataset=='JHS_afr'][, EUR_ANC])))),
        Alpha=0.8,
        OR=c(unlist(WOW[[I]][['WHI_afr']][['0.8perc']]), unlist(WOW[[I]][['HRS_afr']][['0.8perc']]), unlist(WOW[[I]][['UKB_afr']][['0.8perc']]), unlist(WOW[[I]][['HRS_eur']][['0.8perc']]), unlist(WOW[[I]][['JHS_afr']][['0.8perc']])),
        Quantile=c(names(WOW[[I]][['WHI_afr']][['0.8perc']]), names(WOW[[I]][['HRS_afr']][['0.8perc']]),names(WOW[[I]][['UKB_afr']][['0.8perc']]), as.character(1),names(WOW[[I]][['JHS_afr']][['0.8perc']])), Prun=I)
	)
	cat(I)
	cat('\n')
}

for(I in names(AA)){
OR_table[[I]]<-rbind(OR_table[[I]][Quantile!='all'],OR_table[[I]][Quantile=='all'][, Quantile:=EUR_ANC])
}
OR2_table<-vector('list', length(combo))
names(OR2_table)<-names(AA)
for(I in names(AA)){
        OR2_table[[I]]<-rbind(
        data.table(Dataset=c(rep("WHI_afr",5), rep('HRS_afr', 3), rep('UKB_afr', 5), rep('HRS_eur',1), rep("JHS_afr",3)),
        EUR_ANC=c(na.omit(c(unlist(lapply(names(WOW2[[I]][['WHI_afr']][['0.975perc']]), function(X) median(test2[[I]][Dataset=='WHI_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='WHI_afr'][, EUR_ANC]))),
        na.omit(c(unlist(lapply(names(WOW2[[I]][['HRS_afr']][['0.975perc']]), function(X) median(test2[[I]][Dataset=='HRS_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='HRS_afr'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW2[[I]][['UKB_afr']][['0.975perc']]), function(X) median(test2[[I]][Dataset=='UKB_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='UKB_afr'][, EUR_ANC]))),
        median(test2[[I]][Dataset=='HRS_eur'][Quantile==1][, EUR_ANC]),
	na.omit(c(unlist(lapply(names(WOW2[[I]][['JHS_afr']][['0.975perc']]), function(X) median(test2[[I]][Dataset=='JHS_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='JHS_afr'][, EUR_ANC])))),
        Alpha=0.975,
        OR=c(sapply(WOW2[[I]][['WHI_afr']][['0.975perc']], function(X) X$F),sapply(WOW2[[I]][['HRS_afr']][['0.975perc']], function(X) X$F), sapply(WOW2[[I]][['UKB_afr']][['0.975perc']], function(X) X$F), sapply(WOW2[[I]][['HRS_eur']][['0.975perc']], function(X) X$F), sapply(WOW2[[I]][['JHS_afr']][['0.975perc']], function(X) X$F)),
	Quantile=c(names(WOW2[[I]][['WHI_afr']][['0.975perc']]), names(WOW2[[I]][['HRS_afr']][['0.975perc']]),names(WOW2[[I]][['UKB_afr']][['0.975perc']]), as.character(1),names(WOW2[[I]][['JHS_afr']][['0.975perc']])), Prun=I),
	data.table(Dataset=c(rep("WHI",5),rep("HRS_afr",3), rep('UKB_afr', 5), rep('HRS_eur',1), rep("JHS_afr",3)),
        EUR_ANC=c(na.omit(c(unlist(lapply(names(WOW2[[I]][['WHI_afr']][['0.95perc']]), function(X) median(test2[[I]][Dataset=='WHI_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='WHI_afr'][, EUR_ANC]))),
        na.omit(c(unlist(lapply(names(WOW2[[I]][['HRS_afr']][['0.95perc']]), function(X) median(test2[[I]][Dataset=='HRS_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='HRS_afr'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW2[[I]][['UKB_afr']][['0.95perc']]), function(X) median(test2[[I]][Dataset=='UKB_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='UKB_afr'][, EUR_ANC]))),
        median(test2[[I]][Dataset=='HRS_eur'][Quantile==1][, EUR_ANC]),
	na.omit(c(unlist(lapply(names(WOW2[[I]][['JHS_afr']][['0.95perc']]), function(X) median(test2[[I]][Dataset=='JHS_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='JHS_afr'][, EUR_ANC])))),
        Alpha=0.95,	
        OR=c(sapply(WOW2[[I]][['WHI_afr']][['0.95perc']], function(X) X$F),sapply(WOW2[[I]][['HRS_afr']][['0.95perc']], function(X) X$F), sapply(WOW2[[I]][['UKB_afr']][['0.95perc']], function(X) X$F), sapply(WOW2[[I]][['HRS_eur']][['0.95perc']], function(X) X$F), sapply(WOW2[[I]][['JHS_afr']][['0.95perc']], function(X) X$F)),
        Quantile=c(names(WOW2[[I]][['WHI_afr']][['0.95perc']]), names(WOW2[[I]][['HRS_afr']][['0.95perc']]),names(WOW2[[I]][['UKB_afr']][['0.95perc']]), as.character(1),names(WOW2[[I]][['JHS_afr']][['0.95perc']])), Prun=I),
        data.table(Dataset=c(rep("WHI_afr",5),rep("HRS_afr",3), rep('UKB_afr', 5), rep('HRS_eur', 1), rep("JHS_afr",3)),
        EUR_ANC=c(na.omit(c(unlist(lapply(names(WOW2[[I]][['WHI_afr']][['0.9perc']]), function(X) median(test2[[I]][Dataset=='WHI_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='WHI_afr'][, EUR_ANC]))),
        na.omit(c(unlist(lapply(names(WOW2[[I]][['HRS_afr']][['0.9perc']]), function(X) median(test2[[I]][Dataset=='HRS_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='HRS_afr'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW2[[I]][['UKB_afr']][['0.9perc']]), function(X) median(test2[[I]][Dataset=='UKB_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='UKB_afr'][, EUR_ANC]))),
        median(test2[[I]][Dataset=='HRS_eur'][Quantile==1][, EUR_ANC]), 
	na.omit(c(unlist(lapply(names(WOW2[[I]][['JHS_afr']][['0.9perc']]), function(X) median(test2[[I]][Dataset=='JHS_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='JHS_afr'][, EUR_ANC])))), 
        Alpha=0.9,
        OR=c(sapply(WOW2[[I]][['WHI_afr']][['0.9perc']], function(X) X$F),sapply(WOW2[[I]][['HRS_afr']][['0.9perc']], function(X) X$F), sapply(WOW2[[I]][['UKB_afr']][['0.9perc']], function(X) X$F), sapply(WOW2[[I]][['HRS_eur']][['0.9perc']], function(X) X$F), sapply(WOW2[[I]][['JHS_afr']][['0.9perc']], function(X) X$F)),
        Quantile=c(names(WOW2[[I]][['WHI_afr']][['0.9perc']]), names(WOW2[[I]][['HRS_afr']][['0.9perc']]),names(WOW2[[I]][['UKB_afr']][['0.9perc']]), as.character(1),names(WOW2[[I]][['JHS_afr']][['0.9perc']])), Prun=I),
        data.table(Dataset=c(rep("WHI_afr",5),rep("HRS_afr",3), rep('UKB_afr', 5), rep('HRS_eur', 1), rep('JHS_afr',3)),
	EUR_ANC=c(na.omit(c(unlist(lapply(names(WOW2[[I]][['WHI_afr']][['0.85perc']]), function(X) median(test2[[I]][Dataset=='WHI_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='WHI_afr'][, EUR_ANC]))),
        na.omit(c(unlist(lapply(names(WOW2[[I]][['HRS_afr']][['0.85perc']]), function(X) median(test2[[I]][Dataset=='HRS_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='HRS_afr'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW2[[I]][['UKB_afr']][['0.85perc']]), function(X) median(test2[[I]][Dataset=='UKB_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='UKB_afr'][, EUR_ANC]))),
        median(test2[[I]][Dataset=='HRS_eur'][Quantile==1][, EUR_ANC]),
	na.omit(c(unlist(lapply(names(WOW2[[I]][['JHS_afr']][['0.85perc']]), function(X) median(test2[[I]][Dataset=='JHS_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='JHS_afr'][, EUR_ANC])))),
        Alpha=0.85,
        OR=c(sapply(WOW2[[I]][['WHI_afr']][['0.85perc']], function(X) X$F),sapply(WOW2[[I]][['HRS_afr']][['0.85perc']], function(X) X$F), sapply(WOW2[[I]][['UKB_afr']][['0.85perc']], function(X) X$F), sapply(WOW2[[I]][['HRS_eur']][['0.85perc']], function(X) X$F), sapply(WOW2[[I]][['JHS_afr']][['0.85perc']], function(X) X$F)),
        Quantile=c(names(WOW2[[I]][['WHI_afr']][['0.85perc']]), names(WOW2[[I]][['HRS_afr']][['0.85perc']]),names(WOW2[[I]][['UKB_afr']][['0.85perc']]), as.character(1),names(WOW2[[I]][['JHS_afr']][['0.85perc']])), Prun=I),
        data.table(Dataset=c(rep("WHI_afr",5),rep("HRS_afr",3), rep('UKB_afr', 5), rep('HRS_eur', 1), rep("JHS_afr",3)),
	EUR_ANC=c(na.omit(c(unlist(lapply(names(WOW2[[I]][['WHI_afr']][['0.8perc']]), function(X) median(test2[[I]][Dataset=='WHI_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='WHI_afr'][, EUR_ANC]))),
        na.omit(c(unlist(lapply(names(WOW2[[I]][['HRS_afr']][['0.8perc']]), function(X) median(test2[[I]][Dataset=='HRS_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='HRS_afr'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW2[[I]][['UKB_afr']][['0.8perc']]), function(X) median(test2[[I]][Dataset=='UKB_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='UKB_afr'][, EUR_ANC]))),
        median(test2[[I]][Dataset=='HRS_eur'][Quantile==1][, EUR_ANC]),
	na.omit(c(unlist(lapply(names(WOW2[[I]][['JHS_afr']][['0.8perc']]), function(X) median(test2[[I]][Dataset=='JHS_afr'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='JHS_afr'][, EUR_ANC])))),
        Alpha=0.8,
        OR=c(sapply(WOW2[[I]][['WHI_afr']][['0.8perc']], function(X) X$F),sapply(WOW2[[I]][['HRS_afr']][['0.8perc']], function(X) X$F), sapply(WOW2[[I]][['UKB_afr']][['0.8perc']], function(X) X$F), sapply(WOW2[[I]][['HRS_eur']][['0.8perc']], function(X) X$F), sapply(WOW2[[I]][['JHS_afr']][['0.8perc']], function(X) X$F)),
        Quantile=c(names(WOW2[[I]][['WHI_afr']][['0.8perc']]), names(WOW2[[I]][['HRS_afr']][['0.8perc']]),names(WOW2[[I]][['UKB_afr']][['0.8perc']]), as.character(1),names(WOW2[[I]][['JHS_afr']][['0.8perc']])), Prun=I)
	)
        cat(I)
        cat('\n')
}


for(I in names(AA)){
OR2_table[[I]]<-rbind(OR2_table[[I]][Quantile!='all'],OR2_table[[I]][Quantile=='all'][, Quantile:=EUR_ANC])
#OR2_table[[I]][EUR_ANC>=0.05]-> OR2_table[[I]] #revisit this later. 31/08
}

do.call(rbind, OR_table)-> OR_table2
do.call(rbind, OR2_table)-> OR_table3
fwrite(OR_table2, file='~/height_prediction/figs/Odds_Ratio_all_datasets.txt', quote=F, sep="\t")
fwrite(OR_table3, file='~/height_prediction/figs/Odds_Ratio_all_datasets_v2.txt', quote=F, sep="\t")
as.factor(OR_table[[63]]$Alpha)-> OR_table[[63]]$Alpha
OR_table[[63]][, logOR:=log(OR)]

as.factor(OR2_table[[63]]$Alpha)-> OR2_table[[63]]$Alpha
OR2_table[[63]][, logOR:=log(OR)]

colors<-brewer.pal(n = 5, name = 'RdBu')
factor(OR_table[[63]]$Alpha, levels=c(0.975,0.95,0.9, 0.85, 0.8))-> OR_table[[63]]$Alpha
png(paste0('~/height_prediction/figs/logOR_test_', names(AA)[63], '.png'))
one<-ggplot(OR_table[[63]], aes(x=EUR_ANC, y=logOR,  colour=Alpha)) +
geom_point(size=2) + labs(y="log(OR)", x="European Ancestry Proportion") + geom_hline(yintercept=0, linetype="dashed", color = "orange") + stat_smooth(method = "lm", se=F, aes(group=Alpha)) + theme(legend.text=element_text(size=12)) + scale_colour_manual(values=colors)
print(one)
dev.off()

png(paste0('~/height_prediction/figs/logOR_v2_test', names(AA)[63], '.png'))
one<-ggplot(OR2_table[[63]], aes(x=EUR_ANC, y=logOR,colour=Alpha)) +
geom_point(size=2) + labs( y="log(OR)", x="European Ancestry Proportion") + geom_hline(yintercept=0, linetype="dashed", color = "orange") + stat_smooth(method = "lm", se=F, aes(group=Alpha)) + theme(legend.text=element_text(size=12)) + scale_colour_manual(values=colors)
print(one)
dev.off()

png(paste0('~/height_prediction/figs_for_paper/figs/Fig_S8.png'), res=300, unit="in", height=8, width=7)
two<-ggplot(OR_table[[63]], aes(x=EUR_ANC, y=OR,colour=Alpha)) +
geom_point(size=2) + labs(y="OR", x="European Ancestry Proportion") + geom_hline(yintercept=1, linetype="dashed", color = "orange") + stat_smooth(method = "lm", se=F, aes(group=Alpha)) + theme(legend.key=element_blank(),legend.background=element_blank(),legend.title=element_blank(), legend.text=element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
scale_colour_manual(values=colors)
print(two)
dev.off()

png(paste0('~/height_prediction/figs/OR2_test_' , names(AA)[63],'_v2.png'))
two<-ggplot(OR2_table[[63]], aes(x=EUR_ANC, y=OR, shape=Dataset, colour=Alpha)) +
geom_point(size=3) + labs(title="Odds-ratio of P(>=Xth HEIGHT quantile|>=Xth PRS quantile)", y="OR") + geom_hline(yintercept=1, linetype="dashed", color = "orange") + stat_smooth(method = "lm", se=F, aes(group=Alpha)) + 
theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),  axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),legend.text=element_text(size=12)) + scale_colour_manual(values=colors)
print(two)
dev.off()

###

combo[[63]]-> tmp
melt(tmp[,.(PGS, Std.PRS, HEIGHTX, Res.Height, Dataset)])-> me
ggplot(me, aes(x=value, group=Dataset, color=Dataset))+ geom_density() + facet_wrap(~variable, nrow=4, scales='free') + 
theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),  axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text=element_text(size=12))
ggsave('~/height_prediction/figs/test_density.png')
#} #end of if(args[1]=='gwas')
txtStop()

