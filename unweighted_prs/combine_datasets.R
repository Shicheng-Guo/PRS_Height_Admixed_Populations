#!/usr/bin/env Rscript
############################
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
library(psychometric)
library(boot)
library(RColorBrewer)
#args<-'gwas'
##############################################################
#combine all

##
readRDS('~/height_prediction/unweighted_prs/output/results.WHI.Rds')-> results.WHI
readRDS('~/height_prediction/unweighted_prs/output/results.JHS.Rds')-> results.JHS
readRDS('~/height_prediction/unweighted_prs/output/results.UKB_afr.Rds')-> results.UKB_afr
readRDS('~/height_prediction/unweighted_prs/output/results.HRS_eur.Rds')-> results.HRS_eur
readRDS('~/height_prediction/unweighted_prs/output/results.HRS_afr.Rds')-> results.HRS_afr
readRDS('~/height_prediction/unweighted_prs/output/B_WHI.Rds')->B_WHI
readRDS('~/height_prediction/unweighted_prs/output/B_JHS.Rds')->B_JHS
readRDS('~/height_prediction/unweighted_prs/output/B_UKB_afr.Rds')-> B_UKB_afr
readRDS('~/height_prediction/unweighted_prs/output/B_HRS_eur.Rds')-> B_HRS_eur
readRDS('~/height_prediction/unweighted_prs/output/B_HRS_afr.Rds')-> B_HRS_afr
readRDS('~/height_prediction/unweighted_prs/output/PGS3_UKB_afr.Rds')-> PGS3_UKB_afr
readRDS('~/height_prediction/unweighted_prs/output/PGS3_WHI.Rds')-> PGS3_WHI
readRDS('~/height_prediction/unweighted_prs/output/PGS3_JHS.Rds')-> PGS3_JHS
readRDS('~/height_prediction/unweighted_prs/output/PGS3_HRS_eur.Rds')-> PGS3_HRS_eur
readRDS('~/height_prediction/unweighted_prs/output/PGS3_HRS_afr.Rds')-> PGS3_HRS_afr

#if(args[1]=='sib_betas'){
#	readRDS(paste0('~/height_prediction/', args[1], '/ukb_eur/output/results.UKB_eur.Rds'))-> results.UKB_eur
#	readRDS(paste0('~/height_prediction/', args[1], '/ukb_eur/output/B_UKB_eur.Rds'))-> B_UKB_eur
#	readRDS(paste0('~/height_prediction/', args[1], '/ukb_eur/output/PGS3_UKB_eur.Rds'))-> PGS3_UKB_eur
#}
for(I in names(B_JHS)){ #JHS lacks the LD prunning methods
#	if(args[1]=='sib_betas'){
#		ALL<-rbind(B_JHS[[I]][1:2,][, Dataset:='JHS_afr'], B_WHI[[I]][1:4,][, Dataset:='WHI_afr'], B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_HRS_eur[[I]], B_UKB_eur[[I]])
#	} else{
	ALL<-rbind(B_JHS[[I]][1:2,][, Dataset:='JHS_afr'], B_WHI[[I]][1:4,][, Dataset:='WHI_afr'], B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_HRS_eur[[I]])
	ALL$Dataset<-factor(ALL$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "HRS_eur"))
#	}
	#rbind(ALL[!(Dataset %in% c('UKB_EUR', 'pennBB_EA'))], ALL[Dataset %in% c('UKB_EUR', 'pennBB_EA')][, Med_Eur_Anc:=1])-> ALL
	my_plot<-ggplot(ALL, aes(x=Med_Eur_Anc, y=R_sq, colour=Dataset, shape=Dataset)) +
	geom_point(size=1.5, fill="white", alpha=0.8) +
	geom_errorbar(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset,ymin=boots_perc_L, ymax=boots_perc_U), width=0.05, size=0.8) +
	#geom_line(color='lightgray')+
	geom_errorbarh(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset, xmin=HVB_L, xmax=HVB_U), width=0.05, size=0.5) +  scale_color_manual(values=c(brewer.pal(4, 'Set1'),"#101010")) +
	ylab(expression(paste("Partial R"^"2"))) + xlab("European Ancestry Proportion") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),legend.key=element_blank(),legend.background=element_blank(),legend.title=element_blank(), legend.text=element_text(size=12))
	print(my_plot)
	ggsave(paste0('~/height_prediction/unweighted_prs/figs/error_bars_all_v2_', I, '.png'))
}
#
ALL2<-vector('list', length(names(B_JHS)))
names(ALL2)<-names(B_JHS)

for(I in names(B_JHS)){
#	if(args[1]=='sib_betas'){
#		ALL2[[I]]<-rbind(B_JHS[[I]][1:2,][, Dataset:='JHS_afr'], B_WHI[[I]][1:4,][, Dataset:='WHI_afr'], B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_HRS_eur[[I]], B_UKB_eur[[I]])
#		tmp<-1/c(var(results.JHS[[I]][[1]]$t), var(results.JHS[[I]][[2]]$t), var(results.WHI[[I]][[1]]$t), var(results.WHI[[I]][[2]]$t), var(results.WHI[[I]][[3]]$t), var(results.WHI[[I]][[4]]$t), var(results.UKB_afr[[I]][[1]]$t),var(results.UKB_afr[[I]][[2]]$t), var(results.UKB_afr[[I]][[3]]$t), var(results.UKB_afr[[I]][[4]]$t), var(results.HRS_afr[[I]][[1]]$t), var(results.HRS_afr[[I]][[2]]$t), var(results.HRS_eur[[I]]$t), var(results.UKB_eur[[I]]$t))
#	} else{
		ALL2[[I]]<-rbind(B_JHS[[I]][1:2,][, Dataset:='JHS_afr'], B_WHI[[I]][1:4,][, Dataset:='WHI_afr'], B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_HRS_eur[[I]])
	#rbind(ALL2[[I]][!(Dataset %in% c('UKB_EUR', 'pennBB_EA'))], ALL2[[I]][Dataset %in% c('UKB_EUR', 'pennBB_EA')][, Med_Eur_Anc:=1])-> ALL2[[I]]
		tmp<-1/c(var(results.JHS[[I]][[1]]$t), var(results.JHS[[I]][[2]]$t), var(results.WHI[[I]][[1]]$t), var(results.WHI[[I]][[2]]$t), var(results.WHI[[I]][[3]]$t), var(results.WHI[[I]][[4]]$t), var(results.UKB_afr[[I]][[1]]$t),var(results.UKB_afr[[I]][[2]]$t), var(results.UKB_afr[[I]][[3]]$t), var(results.UKB_afr[[I]][[4]]$t), var(results.HRS_afr[[I]][[1]]$t), var(results.HRS_afr[[I]][[2]]$t), var(results.HRS_eur[[I]]$t))  #weighing lm by boostrap replicates.
#	}
	cbind(ALL2[[I]], W=tmp)-> ALL2[[I]]
	ALL2[[I]]$Dataset<-factor(ALL2[[I]]$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "HRS_eur"))
	my_plot2<-ggplot(ALL2[[I]], aes(x=Med_Eur_Anc, y=R_sq)) +
	geom_point(size=1.5, shape=21, fill="white", alpha=0.8) + stat_smooth(method = "lm", mapping = aes(weight = W), col='black') +
	scale_color_manual(values=c(brewer.pal(4, 'Set1'),"#101010")) +
	ylab(expression(paste("Partial R"^"2"))) + xlab("European Ancestry Proportion") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.key=element_blank(), legend.background=element_blank(),legend.title=element_blank(), legend.text=element_text(size=12))
	print(my_plot2)
	ggsave(paste0('~/height_prediction/unweighted_prs/figs/error_bars_all_v3_', I, '.png'))
}

a<-data.table(Name=names(B_JHS), Intercept=unlist(lapply(1:40, function(I) coef(lm(R_sq~Med_Eur_Anc, weights= W, data=ALL2[[I]]))[[1]])), Slope=unlist(lapply(1:40, function(I) coef(lm(R_sq~Med_Eur_Anc, weights= W, data=ALL2[[I]]))[[2]])), R_sq=unlist(lapply(1:40, function(I) summary(lm(R_sq~Med_Eur_Anc, weights= W, data=ALL2[[I]]))[9])), P=unlist(lapply(1:40, function(I) summary(lm(R_sq~Med_Eur_Anc, weights= W, data=ALL2[[I]]))$coefficients[8])))

fwrite(a, file="/~/height_prediction/unweighted_prs/figs/SM_Table1.txt", sep=",")
ALL2b<-vector('list', length(names(B_JHS)))
names(ALL2b)<-names(B_JHS)

for(I in names(B_JHS)){
#	if(args[1]=='sib_betas'){
 #       	ALL2b[[I]]<-rbind(B_JHS[[I]][1:2,][, Dataset:='JHS_afr'], B_WHI[[I]][1:4,][, Dataset:='WHI_afr'], B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_HRS_eur[[I]], B_UKB_eur[[I]])
 #              tmp<-1/c(var(results.JHS[[I]][[1]]$t), var(results.JHS[[I]][[2]]$t), var(results.WHI[[I]][[1]]$t), var(results.WHI[[I]][[2]]$t), var(results.WHI[[I]][[3]]$t), var(results.WHI[[I]][[4]]$t), var(results.UKB_afr[[I]][[1]]$t),var(results.UKB_afr[[I]][[2]]$t), var(results.UKB_afr[[I]][[3]]$t), var(results.UKB_afr[[I]][[4]]$t), var(results.HRS_afr[[I]][[1]]$t), var(results.HRS_afr[[I]][[2]]$t), var(results.HRS_eur[[I]]$t), var(results.UKB_eur[[I]]$t))
#	} else{
        	ALL2b[[I]]<-rbind(B_JHS[[I]][1:2,][, Dataset:='JHS_afr'], B_WHI[[I]][1:4,][, Dataset:='WHI_afr'], B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_HRS_eur[[I]])
                tmp<-1/c(var(results.JHS[[I]][[1]]$t), var(results.JHS[[I]][[2]]$t), var(results.WHI[[I]][[1]]$t), var(results.WHI[[I]][[2]]$t), var(results.WHI[[I]][[3]]$t), var(results.WHI[[I]][[4]]$t), var(results.UKB_afr[[I]][[1]]$t),var(results.UKB_afr[[I]][[2]]$t), var(results.UKB_afr[[I]][[3]]$t), var(results.UKB_afr[[I]][[4]]$t), var(results.HRS_afr[[I]][[1]]$t), var(results.HRS_afr[[I]][[2]]$t), var(results.HRS_eur[[I]]$t))  #weighing lm by boostrap replicates.
#        }
        cbind(ALL2b[[I]], W=tmp)-> ALL2b[[I]]
	ALL2b[[I]]$Dataset<-factor(ALL2b[[I]]$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "HRS_eur"))
	my_plot2<-ggplot(ALL2b[[I]], aes(x=Med_Eur_Anc, y=R_sq)) +
                geom_point(size=1.5, shape=21, fill="white", alpha=0.8) + stat_smooth(method = "lm", mapping = aes(weight = W), col='black') +
                scale_color_manual(values=c(brewer.pal(4, 'Set1'),"#101010")) +
		ylab(expression(paste("Partial R"^"2"))) + xlab("European Ancestry Proportion") +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.key=element_blank(), legend.background=element_blank(),legend.title=element_blank(), legend.text=element_text(size=12))
        print(my_plot2)
        ggsave(paste0('~/height_prediction/unweighted_prs/figs/error_bars_all_v3b_', I, '.png'))
}

for(I in names(B_JHS)){
	ALL2[[I]]$Dataset<-factor(ALL2[[I]]$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "HRS_eur"))
        my_plot<-ggplot(ALL2[[I]], aes(x=Med_Eur_Anc, y=R_sq,colour=Dataset)) +
        geom_point(aes(shape=Dataset), size=1.5, fill="white", alpha=0.8) + stat_smooth(data=ALL2[[I]],method = "lm", mapping = aes(weight = W), col='black') +
        geom_errorbar(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset,ymin=boots_perc_L, ymax=boots_perc_U), width=0.05, size=0.8) +
        #geom_line(color='lightgray')+
        geom_errorbarh(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset, xmin=HVB_L, xmax=HVB_U), width=0.05, size=0.5) + 
	scale_color_manual(values=c(brewer.pal(4, 'Set1'),"#101010")) + 
        ylab(expression(paste("Partial R"^"2"))) + xlab("European Ancestry Proportion") +
  	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.key=element_blank(), legend.background=element_blank(), legend.title=element_blank(), legend.text=element_text(size=12)) 
	print(my_plot)
        ggsave(paste0('~/height_prediction/unweighted_prs/figs/error_bars_all_v4_', I, '.png'))
}
#
#stop here 04/09/2019
ALL3<-vector('list', length(names(B_JHS)))
names(ALL3)<- names(B_JHS)

for (I in names(B_JHS)){
#	if(args[1]=='sib_betas'){
#		ALL3[[I]]<-rbind(B_JHS[[I]][1:2,][, Dataset:='JHS_afr'], B_WHI[[I]][1:4,][, Dataset:='WHI_afr'], B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_HRS_eur[[I]], B_UKB_eur[[I]])
#	 tmp<-lm(R_sq~Med_Eur_Anc,weights=1/
#        c(var(results.JHS[[I]][[1]]$t), var(results.JHS[[I]][[2]]$t), var(results.WHI[[I]][[1]]$t), var(results.WHI[[I]][[2]]$t), var(results.WHI[[I]][[3]]$t), var(results.WHI[[I]][[4]]$t), var(results.UKB_afr[[I]][[1]]$t),var(results.UKB_afr[[I]][[2]]$t), var(results.UKB_afr[[I]][[3]]$t), var(results.UKB_afr[[I]][[4]]$t), var(results.HRS_afr[[I]][[1]]$t), var(results.HRS_afr[[I]][[2]]$t),var(results.HRS_eur[[I]]$t), var(results.UKB_eur[[I]])), data=ALL3[[I]])
#	} else{
                ALL3[[I]]<-rbind(B_JHS[[I]][1:2,][, Dataset:='JHS_afr'], B_WHI[[I]][1:4,][, Dataset:='WHI_afr'], B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_HRS_eur[[I]])
	        tmp<-lm(R_sq~Med_Eur_Anc,weights=1/
        c(var(results.JHS[[I]][[1]]$t), var(results.JHS[[I]][[2]]$t), var(results.WHI[[I]][[1]]$t), var(results.WHI[[I]][[2]]$t), var(results.WHI[[I]][[3]]$t), var(results.WHI[[I]][[4]]$t), var(results.UKB_afr[[I]][[1]]$t),var(results.UKB_afr[[I]][[2]]$t), var(results.UKB_afr[[I]][[3]]$t), var(results.UKB_afr[[I]][[4]]$t), var(results.HRS_afr[[I]][[1]]$t), var(results.HRS_afr[[I]][[2]]$t),var(results.HRS_eur[[I]]$t)), data=ALL3[[I]])
#	}
	#rbind(ALL3[[I]][!(Dataset %in% c('UKB_EUR', 'pennBB_EA'))], ALL3[[I]][Dataset %in% c('UKB_EUR', 'pennBB_EA')][, Med_Eur_Anc:=1])-> ALL3[[I]]
	ALL3[[I]][,Set:=I]
	readRDS(paste0('~/height_prediction/gwas/WHI/output/Nr_SNPs_WHI.Rds'))[Name==I][, Nr]->a
	readRDS(paste0('~/height_prediction/gwas/ukb_afr/output/Nr_SNPs_UKB_afr.Rds'))[Name==I][, Nr]->b
	readRDS(paste0('~/height_prediction/gwas/JHS/output/Nr_SNPs_JHS.Rds'))[Name==I][, Nr]->d
	readRDS(paste0('~/height_prediction/gwas/HRS_eur/output/Nr_SNPs_HRS.Rds'))[Name==I][, Nr]->f
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

do.call(rbind,ALL3)[,.(Set,Intercept,Slope_Intercept, Slope, Nr_SNPs_WHI, Nr_SNPs_UKB, R_sq, Med_Eur_Anc)]->ALL4
#add nr of snps

combo<-vector('list', length(PGS3_JHS))
names(combo)<-names(PGS3_JHS)

for (I in names(PGS3_JHS)){
	rbind(PGS3_WHI[[I]][,.(SUBJID,AGE, age2, HEIGHTX,PGS, SEX,EUR_ANC)][,SUBJ_ID:=SUBJID][, age:=AGE][, sex:='FEMALE'][, SUBJID:=NULL][, AGE:=NULL][, SEX:=NULL][,Dataset:='WHI_afr'][, Res.Height:=resid(lm(HEIGHTX~age+age2))],
	PGS3_HRS_afr[[I]][,SUBJ_ID:=ID][, age:=AGE][, age2:=AGE2][, HEIGHTX:=HEIGHT][, sex:=SEX][,.(SUBJ_ID,age, age2, HEIGHTX,PGS, sex, EUR_ANC)][,Dataset:='HRS_afr'][, Res.Height:=resid(lm(HEIGHTX~sex+age+age2))],
	PGS3_HRS_eur[[I]][,SUBJ_ID:=ID][, age:=AGE][, age2:=AGE2][, HEIGHTX:=HEIGHT][, sex:=SEX][,.(SUBJ_ID,age, age2, HEIGHTX,PGS, sex, EUR_ANC)][,Dataset:='HRS_eur'][, Res.Height:=resid(lm(HEIGHTX~sex+age+age2))],
	PGS3_UKB_afr[[I]][, sex:=Sex][, HEIGHTX:=Height][,SUBJ_ID:=ID][,.(SUBJ_ID, Age, age2, HEIGHTX, PGS, sex, EUR_ANC)][, age:=Age][, Age:=NULL][,Dataset:='UKB_afr'][, Res.Height:=resid(lm(HEIGHTX~sex+age+age2))],
	PGS3_JHS[[I]][,.(SUBJID, age, age2, HEIGHTX, PGS, sex, EUR_ANC)][,SUBJ_ID:=SUBJID][, SUBJID:=NULL][,Dataset:='JHS_afr'][, Res.Height:=resid(lm(HEIGHTX~sex+age+age2))])[,Prun_Set:=I][,Nr_SNPs:=0]-> combo[[I]]
	readRDS(paste0('~/height_prediction/gwas/WHI/output/Nr_SNPs_WHI.Rds'))[Name==I][, Nr]->a
        readRDS(paste0('~/height_prediction/gwas/ukb_afr/output/Nr_SNPs_UKB_afr.Rds'))[Name==I][, Nr]->b
        readRDS(paste0('~/height_prediction/gwas/JHS/output/Nr_SNPs_JHS.Rds'))[Name==I][, Nr]->d
        readRDS(paste0('~/height_prediction/gwas/HRS_eur/output/Nr_SNPs_HRS.Rds'))[Name==I][, Nr]->f
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
fwrite(combo2, file="~/height_prediction/unweighted_prs/output/WHI_JHS_UKB_HRS_summ.txt", sep="\t")
