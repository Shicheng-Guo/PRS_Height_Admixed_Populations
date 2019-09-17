#!/usr/bin/env Rscript
############################
##preable
library("optparse")
library(data.table)
library(dplyr)
library(ggplot2);library(reshape2); library(wesanderson)
library(rlist)
library(asbio)
library(GGally)
library(tidyr)
library(hexbin)
library(psychometric)
library(boot)

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

for(I in names(B_HRS_afr)){ 
	ALL<-rbind(B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_HRS_eur[[I]])
	my_plot<-ggplot(ALL, aes(x=Med_Eur_Anc, y=R_sq, colour=Dataset, shape=Dataset)) +
	geom_point(size=1.5, fill="white") +
	geom_errorbar(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset,ymin=boots_perc_L, ymax=boots_perc_U), width=0.05, size=0.8) +
	geom_errorbarh(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset, xmin=HVB_L, xmax=HVB_U), width=0.05, size=0.5) +  scale_color_brewer(palette="Dark2") +
	ylab(expression(paste("Partial R"^"2"))) + xlab("European Ancestry Proportion") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),axis.text.x=element_text(size=9), axis.text.y=element_text(size=9),legend.key=element_blank(),legend.background=element_blank(),legend.title=element_blank())
	print(my_plot)
	ggsave(paste0('~/height_prediction/imputed/figs/error_bars_all_v2_', I, '.png'))
	}
#
ALL2<-vector('list', length(names(B_HRS_afr)))
names(ALL2)<-names(B_HRS_afr)

for(I in names(B_HRS_afr)){
	ALL2[[I]]<-rbind(B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_HRS_eur[[I]])
	tmp<-1/c(var(results.UKB_afr[[I]][[1]]$t),var(results.UKB_afr[[I]][[2]]$t), var(results.UKB_afr[[I]][[3]]$t), var(results.UKB_afr[[I]][[4]]$t), var(results.HRS_afr[[I]][[1]]$t), var(results.HRS_afr[[I]][[2]]$t), var(results.HRS_eur[[I]]$t))  #weighing lm by boostrap replicates.
	cbind(ALL2[[I]], W=tmp)-> ALL2[[I]]
	my_plot2<-ggplot(ALL2[[I]], aes(x=Med_Eur_Anc, y=R_sq)) +
		geom_point(size=1.5, shape=21, fill="white") + stat_smooth(method = "lm", mapping = aes(weight = W), col='black') +
		scale_color_brewer(palette="Dark2") +
		ylab(expression(paste("Partial R"^"2"))) + xlab("European Ancestry Proportion") +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),axis.text.x=element_text(size=9), axis.text.y=element_text(size=9), legend.key=element_blank(), legend.background=element_blank(),legend.title=element_blank())
	print(my_plot2)
	ggsave(paste0('~/height_prediction/imputed/figs/error_bars_all_v3_', I, '.png'))
	}
ALL2b<-vector('list', length(names(B_HRS_afr)))
names(ALL2b)<-names(B_HRS_afr)

for(I in names(B_HRS_afr)){
	ALL2b[[I]]<-rbind(B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_HRS_eur[[I]])
        tmp<-1/c(var(results.UKB_afr[[I]][[1]]$t),var(results.UKB_afr[[I]][[2]]$t), var(results.UKB_afr[[I]][[3]]$t), var(results.UKB_afr[[I]][[4]]$t), var(results.HRS_afr[[I]][[1]]$t), var(results.HRS_afr[[I]][[2]]$t),var(results.HRS_eur[[I]]$t))
	cbind(ALL2b[[I]], W=tmp)-> ALL2b[[I]]
        my_plot2<-ggplot(ALL2b[[I]], aes(x=Med_Eur_Anc, y=R_sq)) +
                geom_point(size=1.5, shape=21, fill="white") + stat_smooth(method = "lm", mapping = aes(weight = W), col='black') +
                scale_color_brewer(palette="Dark2") +
                ylab(expression(paste("Partial R"^"2"))) + xlab("European Ancestry Proportion") +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),axis.text.x=element_text(size=9), axis.text.y=element_text(size=9), legend.key=element_blank(), legend.background=element_blank(),legend.title=element_blank())
        print(my_plot2)
        ggsave(paste0('~/height_prediction/imputed/figs/error_bars_all_v3b_', I, '.png'))
	}

for(I in names(B_HRS_afr)){
        my_plot<-ggplot(ALL2b[[I]], aes(x=Med_Eur_Anc, y=R_sq,colour=Dataset)) +
        geom_point(aes(shape=Dataset), size=1.5, fill="white") + stat_smooth(data=ALL2[[I]],method = "lm", mapping = aes(weight = W), col='black') +
        geom_errorbar(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset,ymin=boots_perc_L, ymax=boots_perc_U), width=0.05, size=0.8) +
        #geom_line(color='lightgray')+
        geom_errorbarh(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset, xmin=HVB_L, xmax=HVB_U), width=0.05, size=0.5) +  scale_color_brewer(palette="Dark2")+
        ylab(expression(paste("Partial R"^"2"))) + xlab("European Ancestry Proportion") +
  	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),axis.text.x=element_text(size=9), axis.text.y=element_text(size=9), legend.key=element_blank(), legend.background=element_blank(), legend.title=element_blank()) 
	print(my_plot)
        ggsave(paste0('~/height_prediction/imputed/figs/error_bars_all_v4_', I, '.png'))
}
#

