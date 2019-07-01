##preable
#combine all


#brary("optparse")
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
readRDS('WHI/results.WHI.Rds')-> results.WHI
readRDS('JHS/results.JHS.Rds')-> results.JHS
readRDS('pennBB_afr/results.pennBB_afr.Rds')-> results.pennBB_afr
readRDS('pennBB_eur/results.pennBB_eur.Rds')-> results.pennBB_eur
readRDS('ukb_afr/results.UKB_afr.Rds')-> results.UKB_afr
readRDS('ukb_eur/results.UKB_eur.Rds')-> results.UKB_eur
readRDS('WHI/B_WHI.Rds')->B_WHI
readRDS('JHS/B_JHS.Rds')->B_JHS
readRDS('pennBB_afr/B_pennBB_afr.Rds')-> B_pennBB_afr
readRDS('pennBB_eur/B_pennBB_eur_v2.Rds')-> B_pennBB_eur
readRDS('ukb_eur/B_UKB_eur_v2.Rds')->B_UKB_eur
readRDS('ukb_afr/B_UKB_afr.Rds')-> B_UKB_afr
readRDS('pennBB_afr/PGS3_pennBB_afr.Rds')-> PGS3_pennBB_afr
readRDS('pennBB_eur/PGS3_pennBB_eur.Rds')-> PGS3_pennBB_eur
readRDS('ukb_afr/PGS3_UKB_afr.Rds')-> PGS3_UKB_afr
readRDS('ukb_eur/PGS3_UKB_eur.Rds')-> PGS3_UKB_eur
readRDS('WHI/PGS3_WHI.Rds')-> PGS3_WHI
readRDS('JHS/PGS3_JHS.Rds')-> PGS3_JHS

for(I in names(B_JHS)){ #JHS lacks the LD prunning methods
	ALL<-rbind(B_JHS[[I]][1:2,], B_WHI[[I]][1:5,], B_UKB_afr[[I]][1:4,], B_pennBB_afr[[I]][1:2,], B_UKB_eur[[I]], B_pennBB_eur[[I]])
	rbind(ALL[!(Dataset %in% c('UKB_EUR', 'pennBB_EA'))], ALL[Dataset %in% c('UKB_EUR', 'pennBB_EA')][, Med_Eur_Anc:=1])-> ALL
	my_plot<-ggplot(ALL, aes(x=Med_Eur_Anc, y=R_sq,group=Dataset, colour=Dataset)) +
	geom_point(size=1.5, shape=21, fill="white") +
	geom_errorbar(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset,ymin=boots_perc_L, ymax=boots_perc_U), width=0.05, size=0.8) +
	#geom_line(color='lightgray')+
	geom_errorbarh(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset, xmin=HVB_L, xmax=HVB_U), width=0.05, size=0.5) +  scale_color_brewer(palette="Dark2")+
	labs(title = "All Datasets") + ylab("R-squared")+ xlab("European Ancestry Proportion") +
	theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=9), axis.text.y=element_text(size=9))
	print(my_plot)
	ggsave(paste0('figs/error_bars_all_v2_', I, '.png'))
}
#
ALL2<-vector('list', length(names(B_JHS)))
names(ALL2)<-names(B_JHS)

for(I in names(B_JHS)){
	ALL2[[I]]<-rbind(B_JHS[[I]][1:2,],B_WHI[[I]][1:5,], B_UKB_afr[[I]][1:4,], B_pennBB_afr[[I]][1:2,], B_UKB_eur[[I]], B_pennBB_eur[[I]])
	rbind(ALL2[[I]][!(Dataset %in% c('UKB_EUR', 'pennBB_EA'))], ALL2[[I]][Dataset %in% c('UKB_EUR', 'pennBB_EA')][, Med_Eur_Anc:=1])-> ALL2[[I]]
	tmp<-1/c(var(results.JHS[[I]][[1]]$t), var(results.JHS[[I]][[2]]$t), var(results.WHI[[I]][[1]]$t), var(results.WHI[[I]][[2]]$t), var(results.WHI[[I]][[3]]$t), var(results.WHI[[I]][[4]]$t), var(results.WHI[[I]][[5]]$t), var(results.UKB_afr[[I]][[1]]$t),var(results.UKB_afr[[I]][[2]]$t), var(results.UKB_afr[[I]][[3]]$t), var(results.UKB_afr[[I]][[4]]$t), var(results.pennBB_afr[[I]][[1]]$t), var(results.pennBB_afr[[I]][[2]]$t), var(results.UKB_eur[[I]]$total$t), var(results.pennBB_eur[[I]]$total$t))  #weighing lm by boostrap replicates.
	cbind(ALL2[[I]], W=tmp)-> ALL2[[I]]
	my_plot2<-ggplot(ALL2[[I]], aes(x=Med_Eur_Anc, y=R_sq)) +
		geom_point(size=1.5, shape=21, fill="white") + stat_smooth(method = "lm", mapping = aes(weight = W), col='black') +
		scale_color_brewer(palette="Dark2") +
		labs(title = "All Datasets") + ylab("R-squared") + xlab("European Ancestry Proportion")+
		theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),  axis.text.x=element_text(size=9), axis.text.y=element_text(size=9))
	print(my_plot2)
	ggsave(paste0('figs/error_bars_all_v3_', I, '.png'))
}

ALL2b<-vector('list', length(names(B_JHS)))
names(ALL2b)<-names(B_JHS)

for(I in names(B_JHS)){
        ALL2b[[I]]<-rbind(B_JHS[[I]][1:2,], B_WHI[[I]][1:5,], B_UKB_afr[[I]][1:4,], B_pennBB_afr[[I]][1:2,], B_pennBB_eur[[I]])
        rbind(ALL2b[[I]][!(Dataset %in% 'pennBB_EA')], ALL2b[[I]][Dataset %in% 'pennBB_EA'][, Med_Eur_Anc:=1])-> ALL2b[[I]]
        tmp<-1/c(var(results.JHS[[I]][[1]]$t), var(results.JHS[[I]][[2]]$t),var(results.WHI[[I]][[1]]$t), var(results.WHI[[I]][[2]]$t), var(results.WHI[[I]][[3]]$t), var(results.WHI[[I]][[4]]$t), var(results.WHI[[I]][[5]]$t), var(results.UKB_afr[[I]][[1]]$t),var(results.UKB_afr[[I]][[2]]$t), var(results.UKB_afr[[I]][[3]]$t), var(results.UKB_afr[[I]][[4]]$t),var(results.pennBB_afr[[I]][[1]]$t), var(results.pennBB_afr[[I]][[2]]$t),var(results.pennBB_eur[[I]]$total$t))  #weighing lm by boostrap replicates.
        cbind(ALL2b[[I]], W=tmp)-> ALL2b[[I]]
        my_plot2<-ggplot(ALL2b[[I]], aes(x=Med_Eur_Anc, y=R_sq)) +
                geom_point(size=1.5, shape=21, fill="white") + stat_smooth(method = "lm", mapping = aes(weight = W), col='black') +
                scale_color_brewer(palette="Dark2") +
                labs(title = "All Datasets") + ylab("R-squared") + xlab("European Ancestry Proportion")+
                theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),  axis.text.x=element_text(size=9), axis.text.y=element_text(size=9))
        print(my_plot2)
        ggsave(paste0('figs/error_bars_all_v3b_', I, '.png'))
}


for(I in names(B_JHS)){
        my_plot<-ggplot(ALL2[[I]], aes(x=Med_Eur_Anc, y=R_sq,colour=Dataset)) +
        geom_point(size=1.5, shape=21, fill="white") + stat_smooth(data=ALL2[[I]][Dataset != "UKB_EUR"],method = "lm", mapping = aes(weight = W), col='black') +
        geom_errorbar(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset,ymin=boots_perc_L, ymax=boots_perc_U), width=0.05, size=0.8) +
        #geom_line(color='lightgray')+
        geom_errorbarh(aes(x=Med_Eur_Anc, group=Dataset, colour=Dataset, xmin=HVB_L, xmax=HVB_U), width=0.05, size=0.5) +  scale_color_brewer(palette="Dark2")+
        labs(title = "All Datasets") + ylab("R-squared")+ xlab("European Ancestry Proportion") +
	theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=9), axis.text.y=element_text(size=9))
        print(my_plot)
        ggsave(paste0('figs/error_bars_all_v4_', I, '.png'))
}
#

ALL3<-vector('list', length(names(B_JHS)))
names(ALL3)<- names(B_JHS)

for (I in names(B_JHS)){
	ALL3[[I]]<-rbind(B_JHS[[I]][1:2,], B_WHI[[I]][1:5,],B_UKB_afr[[I]][1:4,], B_pennBB_afr[[I]][1:2,], B_UKB_eur[[I]], B_pennBB_eur[[I]])
	rbind(ALL3[[I]][!(Dataset %in% c('UKB_EUR', 'pennBB_EA'))], ALL3[[I]][Dataset %in% c('UKB_EUR', 'pennBB_EA')][, Med_Eur_Anc:=1])-> ALL3[[I]]
	ALL3[[I]][,Set:=I]
	tmp<-lm(R_sq~Med_Eur_Anc,weights=1/
	c(var(results.JHS[[I]][[1]]$t),var(results.JHS[[I]][[2]]$t), var(results.WHI[[I]][[1]]$t), var(results.WHI[[I]][[2]]$t), var(results.WHI[[I]][[3]]$t), var(results.WHI[[I]][[4]]$t), var(results.WHI[[I]][[5]]$t),var(results.UKB_afr[[I]][[1]]$t),var(results.UKB_afr[[I]][[2]]$t), var(results.UKB_afr[[I]][[3]]$t), var(results.UKB_afr[[I]][[4]]$t), var(results.pennBB_afr[[I]][[1]]$t),var(results.pennBB_afr[[I]][[2]]$t), var(results.UKB_eur[[I]]$total$t), var(results.pennBB_eur[[I]]$total$t)), data=ALL3[[I]])  #weight lm
	readRDS('Nr_SNPs_WHI.Rds')[Name==I][, Nr]->a
	readRDS('Nr_SNPs_UKB.Rds')[Name==I][, Nr]->b
	readRDS('Nr_SNPs_JHS.Rds')[Name==I][, Nr]->d
	readRDS('Nr_SNPs_pennBB_afr.Rds')[Name==I][, Nr]->f
	readRDS('Nr_SNPs_pennBB_eur.Rds')[Name==I][, Nr]->g
	ALL3[[I]][,Intercept:=coef(tmp)[[1]]][,Slope:=coef(tmp)[[2]]]
	ALL3[[I]][,Slope_Intercept:=sum(coef(tmp))]
	ALL3[[I]][, Nr_SNPs_WHI:=a]
	ALL3[[I]][, Nr_SNPs_UKB:=b]
	ALL3[[I]][, Nr_SNPs_JHS:=d]
	ALL3[[I]][, Nr_SNPs_pennBB_afr:=f]
	ALL3[[I]][, Nr_SNPs_pennBB_eur:=g]
	cat(I, ' \n')
}

do.call(rbind,ALL3)[,.(Set,Intercept,Slope_Intercept, Slope, Nr_SNPs_WHI, Nr_SNPs_UKB, R_sq, Med_Eur_Anc)]->ALL4
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

ggplot(dt_LD,aes(x=method, y=value, colour=window, shape=variable)) + geom_point(size=2.5, alpha=1) + theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=9), axis.text.y=element_text(size=9))
ggsave('figs/reg_rsq_eur_anc_LD.png')
#
#
melt(dt_genet)->dt_genet
dt_genet[, c("method", "window","p") := tstrsplit(Set, "_")]
as.factor(dt_genet$window)-> dt_genet$window
ggplot() + geom_point(data=dt_genet,aes(x=p, y=value, colour=window, shape=variable), size=2.5, alpha = 1) + theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=9), axis.text.y=element_text(size=9)) 
ggsave('figs/reg_rsq_eur_anc_genet.png')
#
#

melt(dt_phys)-> dt_phys
dt_phys[, c("method", "window","p") := tstrsplit(Set, "_")]
as.factor(dt_phys$window)-> dt_phys$window
factor(dt_phys$window, levels(dt_phys$window)[c(3,7,2,8,6,4,1,5)])-> dt_phys$window
ggplot() + geom_point(data=dt_phys,aes(x=p, y=value, colour=window, shape=variable), size=2.5, alpha = 1) + theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=9), axis.text.y=element_text(size=9))
ggsave('figs/reg_rsq_eur_anc_phys.png')

#plot as function of Nr_SNPs

data.table(Set=unique(ALL4$Set), Nr_SNPs=unique(ALL4$Nr_SNPs_UKB))->A1

ALL4[grep("phys",  ALL4$Set),][,.(Set,Intercept,Slope, Slope_Intercept, R_sq, Med_Eur_Anc)][Med_Eur_Anc==1][, Med_Eur_Anc:=NULL]->dt_phys
ALL4[grep("genet", ALL4$Set),][,.(Set,Intercept,Slope,Slope_Intercept, R_sq, Med_Eur_Anc)][Med_Eur_Anc==1][, Med_Eur_Anc:=NULL]->dt_genet
ALL4[grep("LD",    ALL4$Set),][,.(Set,Intercept,Slope, Slope_Intercept, R_sq, Med_Eur_Anc)][Med_Eur_Anc==1][, Med_Eur_Anc:=NULL]->dt_LD
dt_LD %>% group_by(Set) %>% dplyr::mutate(R_sq=min(R_sq)) %>% as.data.table-> dt_LD #remove ukb_eur
dt_phys %>% group_by(Set) %>% dplyr::mutate(R_sq=min(R_sq)) %>% as.data.table-> dt_phys
dt_genet %>% group_by(Set) %>% dplyr::mutate(R_sq=min(R_sq)) %>% as.data.table-> dt_genet

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
ggplot() + geom_point(data=dt_phys,aes(x=Nr_SNPs, y=value, colour=window, shape=variable), size=2.5, alpha = 1) + theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=9), axis.text.y=element_text(size=9))
ggsave('figs/reg_rsq_eur_anc_phys_v2.png')

ggplot() + geom_point(data=dt_genet,aes(x=Nr_SNPs, y=value, colour=window, shape=variable), size=2.5, alpha = 1) + theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=9), axis.text.y=element_text(size=9))
ggsave('figs/reg_rsq_eur_anc_genet_v2.png')

ggplot() + geom_point(data=dt_LD,aes(x=Nr_SNPs, y=value, colour=Set, shape=variable), size=2.5, alpha = 1) + theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=9), axis.text.y=element_text(size=9))
ggsave('figs/reg_rsq_eur_anc_LD_v2.png')

#try to use just one plot
dt_LD[, p:=0.01]
dt_LD[, .(Set, variable, value, method, window,p, Nr_SNPs)]-> dt_LD
rbind(dt_phys,dt_genet, dt_LD)-> dt

cat('STOP HERE STOP HERE\n')
ggplot() + geom_point(data=dt,aes(x=Nr_SNPs, y=value, colour=method, shape=p), size=1.2, alpha = 0.7) + theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=9), axis.text.y=element_text(size=9)) + facet_wrap(~variable, nrow=2, scales='free_y')
ggsave('figs/test.png')


ggplot() + geom_point(data=dt[variable=='R_sq'], ,aes(x=Nr_SNPs, y=value, colour=method, shape=p), size=1.2, alpha = 0.7) + theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=9), axis.text.y=element_text(size=9)) 
ggsave('figs/testb.png')

ggplot(dt, aes(x=Nr_SNPs, y=value, colour=method, shape=p)) + geom_point(size=1.2, alpha = 0.7) + scale_shape_manual(values=c(16,3,15,0,17,12)) + geom_line(alpha=0.4) + theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=9), axis.text.y=element_text(size=9)) +  facet_wrap(~variable, nrow=2, scales='free_y')
ggsave('figs/test2.png')


ggplot(dt[variable=='R_sq'], aes(x=Nr_SNPs, y=value, colour=method, shape=p)) + geom_point(size=1.2, alpha = 0.7) + scale_shape_manual(values=c(16,3,15,0,17,12)) + geom_line(alpha=0.4) + theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=9), axis.text.y=element_text(size=9))
ggsave('figs/test2b.png')

ggplot(dt, aes(x=Nr_SNPs, y=value, colour=window, shape=method)) + geom_point(size=1.2, alpha = 0.7) + scale_shape_manual(values=c(16,3,13,0,17,12)) + geom_line(alpha=0.4) + theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=9), axis.text.y=element_text(size=9)) +  facet_wrap(~variable, nrow=2, scales='free_y')
ggsave('figs/test3.png')


#ggplot(dt, aes(x=Nr_SNPs, y=value, shape=p, colour=int)) + geom_point(size=1.2, alpha = 0.7) + geom_line(alpha=0.4) +  facet_wrap(~variable, nrow=2, scales='free_y')
#ggsave('figs/test4.png')

###
combo<-vector('list', length(PGS3_JHS))
names(combo)<-names(PGS3_JHS)

for (I in names(PGS3_JHS)){
	rbind(PGS3_WHI[[I]][,.(SUBJID,AGE, age2, HEIGHTX,PGS, SEX,EUR_ANC)][,SUBJ_ID:=SUBJID][, age:=AGE][, sex:='FEMALE'][, SUBJID:=NULL][, AGE:=NULL][, SEX:=NULL][,Dataset:='WHI'][, Res.Height:=resid(lm(HEIGHTX~age+age2))],
	PGS3_pennBB_afr[[I]][,.(SUBJ_ID,age, age2, HEIGHTX,PGS, sex, EUR_ANC)][,Dataset:='pennBB_AA'][, Res.Height:=resid(lm(HEIGHTX~sex+age+age2))],
	PGS3_pennBB_eur[[I]][,.(SUBJ_ID,age, age2, HEIGHTX,PGS, sex, EUR_ANC)][,Dataset:='pennBB_EA'][, Res.Height:=resid(lm(HEIGHTX~sex+age+age2))],
	PGS3_UKB_afr[[I]][, sex:=Sex][, HEIGHTX:=Height][,SUBJ_ID:=ID][,.(SUBJ_ID, Age, age2, HEIGHTX, PGS, sex, EUR_ANC)][, age:=Age][, Age:=NULL][,Dataset:='UKB_AFR'][, Res.Height:=resid(lm(HEIGHTX~sex+age+age2))],
	PGS3_UKB_eur[[I]][, sex:=Sex][, HEIGHTX:=Height][,SUBJ_ID:=ID][,.(SUBJ_ID, Age, age2, HEIGHTX, PGS, sex, EUR_ANC)][, age:=Age][, Age:=NULL][,Dataset:='UKB_EUR'][, Res.Height:=resid(lm(HEIGHTX~sex+age+age2))],
	PGS3_JHS[[I]][,.(SUBJID, age, age2, HEIGHTX, PGS, sex, EUR_ANC)][,SUBJ_ID:=SUBJID][, SUBJID:=NULL][,Dataset:='JHS'][, Res.Height:=resid(lm(HEIGHTX~sex+age+age2))])[,Prun_Set:=I][,Nr_SNPs:=0]-> combo[[I]]
	readRDS('Nr_SNPs_WHI.Rds')[Name==I][, Nr]-> a
	readRDS('Nr_SNPs_UKB.Rds')[Name==I][, Nr]-> b
	readRDS('Nr_SNPs_JHS.Rds')[Name==I][, Nr]-> d
	readRDS('Nr_SNPs_pennBB_afr.Rds')[Name==I][, Nr]-> f
	combo[[I]][Dataset=='WHI']$Nr_SNPs<-a
	combo[[I]][Dataset=='JHS']$Nr_SNPs<-d
	combo[[I]][Dataset=='UKB_AFR']$Nr_SNPs<-b
	combo[[I]][Dataset=='UKB_EUR']$Nr_SNPs<-b
	combo[[I]][Dataset=='pennBB_EA']$Nr_SNPs<-f
	combo[[I]][Dataset=='pennBB_AA']$Nr_SNPs<-f
	combo[[I]][, Std.PRS:=scale(PGS), by=Dataset]-> combo[[I]]
	cat(I, ' done\n')
}
do.call(rbind,combo)-> combo2
combo2[Dataset=='UKB_EUR', EUR_ANC:=1]
#combo2[Dataset=='pennBB_EA', EUR_ANC:=1]
fwrite(combo2, file="WHI_pennBB_UKB_JHS_summ.txt", sep="\t")
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
z2<-quantile(data2[Dataset=='UKB_EUR'][,Std.PRS], probs=x)[[1]]
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

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F3_X=F3_x(x=X,data=combo[[Y]][Dataset=='WHI'],  data2=combo[[Y]]), Quantile=X))))-> AA
names(AA)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F2_X=F2_x(x=X,data=combo[[Y]][Dataset=='WHI']), Quantile=X))))-> AAA
names(AAA)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F3_X=F3_x(x=X,data=combo[[Y]][Dataset=='JHS'],  data2=combo[[Y]]), Quantile=X))))-> AJ
names(AJ)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F2_X=F2_x(x=X,data=combo[[Y]][Dataset=='JHS']), Quantile=X))))-> AAJ
names(AAJ)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F3_X=F3_x(x=X,data=combo[[Y]][Dataset=='UKB_AFR'],  data2=combo[[Y]]), Quantile=X))))-> AU
names(AU)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F2_X=F2_x(x=X,data=combo[[Y]][Dataset=='UKB_AFR']), Quantile=X))))-> AAU
names(AAU)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F3_X=F3_x(x=X,data=combo[[Y]][Dataset=='pennBB_AA'],  data2=combo[[Y]]), Quantile=X))))-> AP
names(AP)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F2_X=F2_x(x=X,data=combo[[Y]][Dataset=='pennBB_AA']), Quantile=X))))-> AAP
names(AAP)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F3_X=F3_x(x=X,data=combo[[Y]][Dataset=='pennBB_EA'],  data2=combo[[Y]]), Quantile=X))))-> APE
names(APE)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F2_X=F2_x(x=X,data=combo[[Y]][Dataset=='pennBB_EA']), Quantile=X))))-> AAPE
names(AAPE)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F2_X=F2_x(x=X,data=combo[[Y]][Dataset=='UKB_EUR']), Quantile=X))))-> AUE
names(AUE)<-names(B_JHS)


I<-names(AA)[63]

png(paste0('figs/OR_WHI_', I,  ".png"))
ggplot(AA[[I]], aes(x=Quantile, y=F3_X.F)) +
geom_point(size=2) + labs(title="Odds ratio of P(>=Xth HEIGHT quantile|>=Xth PRS quantile)", y="OR") + geom_hline(yintercept=1, linetype="dashed", color = "orange") +
theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
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
tmp11<-AUE[[I]]
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
tmp1[, F3_X.F:=NULL][, F3_X.G:=NULL][,F3_X.f:=NULL][,F3_X.g:=NULL][, Group:="WHI_Eur"]
tmp3[, F3_X.F:=NULL][, F3_X.G:=NULL][,F3_X.f:=NULL][,F3_X.g:=NULL][, Group:="JHS_Eur"]
tmp5[, F3_X.F:=NULL][, F3_X.G:=NULL][,F3_X.f:=NULL][,F3_X.g:=NULL][, Group:="UKB_AFR_Eur"]
tmp7[,F3_X.F:=NULL][, F3_X.G:=NULL][,F3_X.f:=NULL][,F3_X.g:=NULL][, Group:="pennBB_AA_Eur"]
tmp9[,F3_X.F:=NULL][, F3_X.G:=NULL][,F3_X.f:=NULL][,F3_X.g:=NULL][, Group:="pennBB_EA_Eur"]
#
setDT(tmp2)
tmp2[,OR:=F2_X][, F2_X:=NULL][, Group:='WHI_matched']
setDT(tmp4)
tmp4[,OR:=F2_X][, F2_X:=NULL][, Group:="JHS_matched"]
setDT(tmp6)
tmp6[,OR:=F2_X][, F2_X:=NULL][, Group:="UKB_AFR_matched"]
setDT(tmp8)
tmp8[,OR:=F2_X][, F2_X:=NULL][, Group:="pennBB_AA_matched"]
setDT(tmp10)
tmp10[,OR:=F2_X][, F2_X:=NULL][,Group:="pennBB_EA_matched"]
setDT(tmp11)
tmp11[,OR:=F2_X][, F2_X:=NULL][,Group:="UKB_EUR_matched"]


all<-rbind(tmp1,tmp2,tmp3,tmp4,tmp5, tmp6,tmp7, tmp8, tmp9, tmp10, tmp11)
all[, Group2:=c(rep("WHI", 100), rep("JHS", 100), rep("UKB_AFR", 100), rep("pennBB_AA",100), rep("pennBB_EA",100), rep("UKB_EUR",50))]
all[, Group3:=c(rep("Eur", 50), rep("Matched",50), rep("Eur",50), rep("Matched",50), rep("Eur",50), rep("Matched", 50), rep("Eur",50), rep("Matched", 50), rep("Eur", 50), rep("Matched", 50), rep("Matched", 50))]
all[, Prs.Quantile:=Quantile]

pdf(paste0('figs/OR_WHI_JHS_', I,  ".pdf"))
ggplot(all[Prs.Quantile<=0.975], aes(x=Prs.Quantile, y=OR, group=Group, color=Group2, linetype=Group3)) + geom_line() +
geom_point(size=0.7) + labs(y="OR", x="Quantile") +
theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
dev.off()


pdf(paste0('figs/OR_WHI_JHS_', I,  "_v2.pdf"))
ggplot(all[Prs.Quantile<=0.975], aes(x=Prs.Quantile, y=OR, group=Group, color=Group2, linetype=Group3)) + geom_smooth() +
geom_point(size=0.7) + labs(y="OR", x="Quantile") + 
theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
dev.off()
#lapply(1:length(names(A)),function(X) combo[[X]][,Quantile:= cut(EUR_ANC,breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.2), na.rm=TRUE),include.lowest=TRUE)])


WOW<-vector('list', length(combo));test<-vector('list', length(combo))
names(WOW)<-names(B_WHI);names(test)<-names(B_WHI)
for(I in names(AA)){
	a<-combo[[I]][Dataset=='WHI'][,Quantile:= cut(EUR_ANC,breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.2), na.rm=TRUE),include.lowest=TRUE)]
	b<-combo[[I]][Dataset=='pennBB_AA'][,Quantile:= cut(EUR_ANC,breaks=quantile(EUR_ANC, probs=seq(0,1, by=1/3), na.rm=TRUE),include.lowest=TRUE)]
	d<-combo[[I]][Dataset=='UKB_AFR'][,Quantile:= cut(EUR_ANC,breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.25), na.rm=TRUE),include.lowest=TRUE)]
	e<-combo[[I]][Dataset=='UKB_EUR'][, EUR_ANC:=1][,Quantile:= as.factor(as.character(1))];
	f<-combo[[I]][Dataset=='pennBB_EA'][,EUR_ANC:=1][,Quantile:=as.factor(as.character(1))]
	g<-combo[[I]][Dataset=='JHS'][,Quantile:= cut(EUR_ANC,breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.5), na.rm=TRUE),include.lowest=TRUE)]
	test[[I]]<-rbind(a,b,d,e,f,g);test[[I]]$Quantile<-as.character(test[[I]]$Quantile)
	WOW[[I]]<-vector('list', 6);names(WOW[[I]])<- c('WHI', 'pennBB_AA','UKB_AFR', 'UKB_EUR', 'pennBB_EA', 'JHS')
	lapply(levels(a$Quantile), function(X) F2_x(data=a[Quantile==X], x=0.975))-> WOW[[I]][['WHI']][[1]];list.append(WOW[[I]][['WHI']][[1]], F2_x(data=a, x=0.975))-> WOW[[I]][['WHI']][[1]]	
	lapply(levels(a$Quantile), function(X) F2_x(data=a[Quantile==X], x=0.95))-> WOW[[I]][['WHI']][[2]];list.append(WOW[[I]][['WHI']][[2]], F2_x(data=a, x=0.95))-> WOW[[I]][['WHI']][[2]]
	lapply(levels(a$Quantile), function(X) F2_x(data=a[Quantile==X], x=0.9))-> WOW[[I]][['WHI']][[3]];list.append(WOW[[I]][['WHI']][[3]], F2_x(data=a, x=0.9))-> WOW[[I]][['WHI']][[3]]
	lapply(levels(a$Quantile), function(X) F2_x(data=a[Quantile==X], x=0.85))-> WOW[[I]][['WHI']][[4]];list.append(WOW[[I]][['WHI']][[4]], F2_x(data=a, x=0.85))-> WOW[[I]][['WHI']][[4]]
	lapply(levels(a$Quantile), function(X) F2_x(data=a[Quantile==X], x=0.8))-> WOW[[I]][['WHI']][[5]];list.append(WOW[[I]][['WHI']][[5]], F2_x(data=a, x=0.8))-> WOW[[I]][['WHI']][[5]]
	names(WOW[[I]][['WHI']])<-c('0.975perc', '0.95perc', '0.9perc', '0.85perc', '0.8perc');names(WOW[[I]][['WHI']][['0.9perc']])<-c(levels(a$Quantile), 'all')
	names(WOW[[I]][['WHI']][['0.975perc']])<-c(levels(a$Quantile), 'all');names(WOW[[I]][['WHI']][['0.95perc']])<-c(levels(a$Quantile), 'all')
	names(WOW[[I]][['WHI']][['0.85perc']])<-c(levels(a$Quantile), 'all');names(WOW[[I]][['WHI']][['0.8perc']])<-c(levels(a$Quantile),'all') 
	lapply(levels(b$Quantile), function(X) F2_x(data=b[Quantile==X], x=0.975))-> WOW[[I]][['pennBB_AA']][[1]];list.append(WOW[[I]][['pennBB_AA']][[1]], F2_x(data=b, x=0.975))-> WOW[[I]][['pennBB_AA']][[1]]
	lapply(levels(b$Quantile), function(X) F2_x(data=b[Quantile==X], x=0.95))-> WOW[[I]][['pennBB_AA']][[2]];list.append(WOW[[I]][['pennBB_AA']][[2]], F2_x(data=b, x=0.95))-> WOW[[I]][['pennBB_AA']][[2]]
	lapply(levels(b$Quantile), function(X) F2_x(data=b[Quantile==X], x=0.9))-> WOW[[I]][['pennBB_AA']][[3]];list.append(WOW[[I]][['pennBB_AA']][[3]], F2_x(data=b, x=0.9))-> WOW[[I]][['pennBB_AA']][[3]]
	lapply(levels(b$Quantile), function(X) F2_x(data=b[Quantile==X], x=0.85))-> WOW[[I]][['pennBB_AA']][[4]];list.append(WOW[[I]][['pennBB_AA']][[4]], F2_x(data=b, x=0.85))-> WOW[[I]][['pennBB_AA']][[4]]
        lapply(levels(b$Quantile), function(X) F2_x(data=b[Quantile==X], x=0.8))-> WOW[[I]][['pennBB_AA']][[5]];list.append(WOW[[I]][['pennBB_AA']][[5]], F2_x(data=b, x=0.8))-> WOW[[I]][['pennBB_AA']][[5]]
	names(WOW[[I]][['pennBB_AA']])<-c('0.975perc', '0.95perc', '0.9perc', '0.85perc', '0.8perc');names(WOW[[I]][['pennBB_AA']][['0.9perc']])<-c(levels(b$Quantile), 'all')
	names(WOW[[I]][['pennBB_AA']][['0.95perc']])<-c(levels(b$Quantile), 'all');names(WOW[[I]][['pennBB_AA']][['0.975perc']])<-c(levels(b$Quantile), 'all')
	names(WOW[[I]][['pennBB_AA']][['0.85perc']])<-c(levels(b$Quantile), 'all');names(WOW[[I]][['pennBB_AA']][['0.8perc']])<-c(levels(b$Quantile), 'all')
	lapply(levels(d$Quantile), function(X) F2_x(data=d[Quantile==X], x=0.975))-> WOW[[I]][['UKB_AFR']][[1]];list.append(WOW[[I]][['UKB_AFR']][[1]], F2_x(data=d, x=0.975))-> WOW[[I]][['UKB_AFR']][[1]]
	lapply(levels(d$Quantile), function(X) F2_x(data=d[Quantile==X], x=0.95))-> WOW[[I]][['UKB_AFR']][[2]];list.append(WOW[[I]][['UKB_AFR']][[2]], F2_x(data=d, x=0.95))-> WOW[[I]][['UKB_AFR']][[2]]
	lapply(levels(d$Quantile), function(X) F2_x(data=d[Quantile==X], x=0.9))-> WOW[[I]][['UKB_AFR']][[3]];list.append(WOW[[I]][['UKB_AFR']][[3]], F2_x(data=d, x=0.9))-> WOW[[I]][['UKB_AFR']][[3]]
	lapply(levels(d$Quantile), function(X) F2_x(data=d[Quantile==X], x=0.85))-> WOW[[I]][['UKB_AFR']][[4]];list.append(WOW[[I]][['UKB_AFR']][[4]], F2_x(data=d, x=0.85))-> WOW[[I]][['UKB_AFR']][[4]]	
	lapply(levels(d$Quantile), function(X) F2_x(data=d[Quantile==X], x=0.8))-> WOW[[I]][['UKB_AFR']][[5]];list.append(WOW[[I]][['UKB_AFR']][[5]], F2_x(data=d, x=0.8))-> WOW[[I]][['UKB_AFR']][[5]]
	names(WOW[[I]][['UKB_AFR']])<-c('0.975perc', '0.95perc', '0.9perc','0.85perc', '0.8perc');names(WOW[[I]][['UKB_AFR']][['0.9perc']])<-c(levels(d$Quantile), 'all')
	names(WOW[[I]][['UKB_AFR']][['0.95perc']])<-c(levels(d$Quantile), 'all');names(WOW[[I]][['UKB_AFR']][['0.975perc']])<-c(levels(d$Quantile), 'all');names(WOW[[I]][['UKB_AFR']][['0.85perc']])<-c(levels(d$Quantile), 'all'); names(WOW[[I]][['UKB_AFR']][['0.8perc']])<-c(levels(d$Quantile), 'all')
	lapply(unique(e$Quantile), function(X) F2_x(data=e[Quantile==X], x=0.975))-> WOW[[I]][['UKB_EUR']][[1]]
	lapply(unique(e$Quantile), function(X) F2_x(data=e[Quantile==X], x=0.95))-> WOW[[I]][['UKB_EUR']][[2]];lapply(unique(e$Quantile), function(X) F2_x(data=e[Quantile==X], x=0.9))-> WOW[[I]][['UKB_EUR']][[3]]
	lapply(unique(e$Quantile), function(X) F2_x(data=e[Quantile==X], x=0.85))-> WOW[[I]][['UKB_EUR']][[4]];lapply(unique(e$Quantile), function(X) F2_x(data=e[Quantile==X], x=0.8))-> WOW[[I]][['UKB_EUR']][[5]]
	names(WOW[[I]][['UKB_EUR']])<-c('0.975perc', '0.95perc', '0.9perc', '0.85perc', '0.8perc');names(WOW[[I]][['UKB_EUR']][['0.9perc']])<-unique(e$Quantile)
	names(WOW[[I]][['UKB_EUR']][['0.95perc']])<-unique(e$Quantile);names(WOW[[I]][['UKB_EUR']][['0.975perc']])<-unique(e$Quantile);names(WOW[[I]][['UKB_EUR']][['0.85perc']])<-unique(e$Quantile);names(WOW[[I]][['UKB_EUR']][['0.8perc']])<-unique(e$Quantile)
	lapply(unique(f$Quantile), function(X) F2_x(data=f[Quantile==X], x=0.975))-> WOW[[I]][['pennBB_EA']][[1]]
	lapply(unique(f$Quantile), function(X) F2_x(data=f[Quantile==X], x=0.95))-> WOW[[I]][['pennBB_EA']][[2]];lapply(unique(f$Quantile), function(X) F2_x(data=f[Quantile==X], x=0.9))-> WOW[[I]][['pennBB_EA']][[3]]
	lapply(unique(f$Quantile), function(X) F2_x(data=f[Quantile==X], x=0.85))-> WOW[[I]][['pennBB_EA']][[4]];lapply(unique(f$Quantile), function(X) F2_x(data=f[Quantile==X], x=0.8))-> WOW[[I]][['pennBB_EA']][[5]]
	names(WOW[[I]][['pennBB_EA']])<-c('0.975perc', '0.95perc', '0.9perc','0.85perc', '0.8perc');names(WOW[[I]][['pennBB_EA']][['0.9perc']])<-unique(f$Quantile)
	names(WOW[[I]][['pennBB_EA']][['0.95perc']])<-unique(f$Quantile);names(WOW[[I]][['pennBB_EA']][['0.975perc']])<-unique(f$Quantile)
	names(WOW[[I]][['pennBB_EA']][['0.85perc']])<-unique(f$Quantile);names(WOW[[I]][['pennBB_EA']][['0.8perc']])<-unique(f$Quantile)
	lapply(levels(g$Quantile), function(X) F2_x(data=g[Quantile==X], x=0.975))-> WOW[[I]][['JHS']][[1]];list.append(WOW[[I]][['JHS']][[1]], F2_x(data=d, x=0.975))-> WOW[[I]][['JHS']][[1]]
        lapply(levels(g$Quantile), function(X) F2_x(data=g[Quantile==X], x=0.95))-> WOW[[I]][['JHS']][[2]];list.append(WOW[[I]][['JHS']][[2]], F2_x(data=d, x=0.95))-> WOW[[I]][['JHS']][[2]]
        lapply(levels(g$Quantile), function(X) F2_x(data=g[Quantile==X], x=0.9))-> WOW[[I]][['JHS']][[3]];list.append(WOW[[I]][['JHS']][[3]], F2_x(data=d, x=0.9))-> WOW[[I]][['JHS']][[3]]
        lapply(levels(g$Quantile), function(X) F2_x(data=g[Quantile==X], x=0.85))-> WOW[[I]][['JHS']][[4]];list.append(WOW[[I]][['JHS']][[4]], F2_x(data=d, x=0.85))-> WOW[[I]][['JHS']][[4]]
        lapply(levels(g$Quantile), function(X) F2_x(data=g[Quantile==X], x=0.8))-> WOW[[I]][['JHS']][[5]];list.append(WOW[[I]][['JHS']][[5]], F2_x(data=d, x=0.8))-> WOW[[I]][['JHS']][[5]]
        names(WOW[[I]][['JHS']])<-c('0.975perc', '0.95perc', '0.9perc','0.85perc', '0.8perc');names(WOW[[I]][['JHS']][['0.9perc']])<-c(levels(g$Quantile), 'all')
        names(WOW[[I]][['JHS']][['0.95perc']])<-c(levels(g$Quantile), 'all');names(WOW[[I]][['JHS']][['0.975perc']])<-c(levels(g$Quantile), 'all');names(WOW[[I]][['JHS']][['0.85perc']])<-c(levels(g$Quantile), 'all'); names(WOW[[I]][['JHS']][['0.8perc']])<-c(levels(g$Quantile), 'all')
	cat(I);cat('\n')
}
	

WOW2<-vector('list', length(B_WHI));test2<-vector('list', length(B_WHI))
names(WOW2)<-names(B_WHI);names(test2)<-names(B_WHI)

for(I in names(AA)){
        a<-combo[[I]][Dataset=='WHI'][,Quantile:= cut(EUR_ANC,breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.2), na.rm=TRUE),include.lowest=TRUE)]
        b<-combo[[I]][Dataset=='pennBB_AA'][,Quantile:= cut(EUR_ANC,breaks=quantile(EUR_ANC, probs=seq(0,1, by=1/3), na.rm=TRUE),include.lowest=TRUE)]
        d<-combo[[I]][Dataset=='UKB_AFR'][,Quantile:= cut(EUR_ANC,breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.25), na.rm=TRUE),include.lowest=TRUE)]
        e<-combo[[I]][Dataset=='UKB_EUR'][, EUR_ANC:=1][,Quantile:= as.factor(as.character(1))];
	f<-combo[[I]][Dataset=='pennBB_EA'][,EUR_ANC:=1][,Quantile:=as.factor(as.character(1))]
	g<-combo[[I]][Dataset=='JHS'][,Quantile:= cut(EUR_ANC,breaks=quantile(EUR_ANC, probs=seq(0,1, by=0.5), na.rm=TRUE),include.lowest=TRUE)]
        test2[[I]]<-rbind(a,b,d,e,f,g);test2[[I]]$Quantile<-as.character(test2[[I]]$Quantile)
        WOW2[[I]]<-vector('list', 6);names(WOW2[[I]])<- c('WHI', 'pennBB_AA','UKB_AFR', 'UKB_EUR', 'pennBB_EA', 'JHS')
        lapply(levels(a$Quantile), function(X) F3_x(data=a[Quantile==X], x=0.975,  data2=combo[[I]]))-> WOW2[[I]][['WHI']][[1]];list.append(WOW2[[I]][['WHI']][[1]], F3_x(data=a, x=0.975,data2=combo[[I]]))-> WOW2[[I]][['WHI']][[1]]
        lapply(levels(a$Quantile), function(X) F3_x(data=a[Quantile==X], x=0.95, data2=combo[[I]]))-> WOW2[[I]][['WHI']][[2]];list.append(WOW2[[I]][['WHI']][[2]], F3_x(data=a, x=0.95,data2=combo[[I]]))-> WOW2[[I]][['WHI']][[2]]
        lapply(levels(a$Quantile), function(X) F3_x(data=a[Quantile==X], x=0.9, data2=combo[[I]]))-> WOW2[[I]][['WHI']][[3]];list.append(WOW2[[I]][['WHI']][[3]], F3_x(data=a, x=0.9,data2=combo[[I]]))-> WOW2[[I]][['WHI']][[3]]
        lapply(levels(a$Quantile), function(X) F3_x(data=a[Quantile==X], x=0.85, data2=combo[[I]]))-> WOW2[[I]][['WHI']][[4]];list.append(WOW2[[I]][['WHI']][[4]], F3_x(data=a, x=0.85,data2=combo[[I]]))-> WOW2[[I]][['WHI']][[4]]
	lapply(levels(a$Quantile), function(X) F3_x(data=a[Quantile==X], x=0.8, data2=combo[[I]]))-> WOW2[[I]][['WHI']][[5]];list.append(WOW2[[I]][['WHI']][[5]], F3_x(data=a, x=0.8,data2=combo[[I]]))-> WOW2[[I]][['WHI']][[5]]
        names(WOW2[[I]][['WHI']])<-c('0.975perc', '0.95perc', '0.9perc', '0.85perc', '0.8perc');names(WOW2[[I]][['WHI']][['0.9perc']])<-c(levels(a$Quantile), 'all')
        names(WOW2[[I]][['WHI']][['0.95perc']])<-c(levels(a$Quantile), 'all');names(WOW2[[I]][['WHI']][['0.975perc']])<-c(levels(a$Quantile), 'all')
 	names(WOW2[[I]][['WHI']][['0.85perc']])<-c(levels(a$Quantile), 'all');names(WOW2[[I]][['WHI']][['0.8perc']])<-c(levels(a$Quantile), 'all')
        lapply(levels(b$Quantile), function(X) F3_x(data=b[Quantile==X], x=0.975,  data2=combo[[I]]))-> WOW2[[I]][['pennBB_AA']][[1]];list.append(WOW2[[I]][['pennBB_AA']][[1]], F3_x(data=b, x=0.975,data2=combo[[I]]))-> WOW2[[I]][['pennBB_AA']][[1]]
        lapply(levels(b$Quantile), function(X) F3_x(data=b[Quantile==X], x=0.95, data2=combo[[I]]))-> WOW2[[I]][['pennBB_AA']][[2]];list.append(WOW2[[I]][['pennBB_AA']][[2]], F3_x(data=b, x=0.95,data2=combo[[I]]))-> WOW2[[I]][['pennBB_AA']][[2]]
        lapply(levels(b$Quantile), function(X) F3_x(data=b[Quantile==X], x=0.9, data2=combo[[I]]))-> WOW2[[I]][['pennBB_AA']][[3]];list.append(WOW2[[I]][['pennBB_AA']][[3]], F3_x(data=b, x=0.9,data2=combo[[I]]))-> WOW2[[I]][['pennBB_AA']][[3]]
        lapply(levels(b$Quantile), function(X) F3_x(data=b[Quantile==X], x=0.85, data2=combo[[I]]))-> WOW2[[I]][['pennBB_AA']][[4]];list.append(WOW2[[I]][['pennBB_AA']][[4]], F3_x(data=b, x=0.85,data2=combo[[I]]))-> WOW2[[I]][['pennBB_AA']][[4]]
	lapply(levels(b$Quantile), function(X) F3_x(data=b[Quantile==X], x=0.8, data2=combo[[I]]))-> WOW2[[I]][['pennBB_AA']][[5]];list.append(WOW2[[I]][['pennBB_AA']][[5]], F3_x(data=b, x=0.8,data2=combo[[I]]))-> WOW2[[I]][['pennBB_AA']][[5]]
        names(WOW2[[I]][['pennBB_AA']])<-c('0.975perc', '0.95perc', '0.9perc', '0.85perc', '0.8perc')
        names(WOW2[[I]][['pennBB_AA']][['0.9perc']])<-c(levels(b$Quantile), 'all')
        names(WOW2[[I]][['pennBB_AA']][['0.95perc']])<-c(levels(b$Quantile), 'all')
        names(WOW2[[I]][['pennBB_AA']][['0.975perc']])<-c(levels(b$Quantile), 'all')
        names(WOW2[[I]][['pennBB_AA']][['0.85perc']])<-c(levels(b$Quantile), 'all')
	names(WOW2[[I]][['pennBB_AA']][['0.8perc']])<-c(levels(b$Quantile), 'all')
        lapply(levels(d$Quantile), function(X) F3_x(data=d[Quantile==X], x=0.975,  data2=combo[[I]]))-> WOW2[[I]][['UKB_AFR']][[1]];list.append(WOW2[[I]][['UKB_AFR']][[1]], F3_x(data=d, x=0.975,data2=combo[[I]]))-> WOW2[[I]][['UKB_AFR']][[1]]
        lapply(levels(d$Quantile), function(X) F3_x(data=d[Quantile==X], x=0.95, data2=combo[[I]]))-> WOW2[[I]][['UKB_AFR']][[2]];list.append(WOW2[[I]][['UKB_AFR']][[2]], F3_x(data=d, x=0.95,data2=combo[[I]]))-> WOW2[[I]][['UKB_AFR']][[2]]
        lapply(levels(d$Quantile), function(X) F3_x(data=d[Quantile==X], x=0.9, data2=combo[[I]]))-> WOW2[[I]][['UKB_AFR']][[3]];list.append(WOW2[[I]][['UKB_AFR']][[3]], F3_x(data=d, x=0.9,data2=combo[[I]]))-> WOW2[[I]][['UKB_AFR']][[3]]
        lapply(levels(d$Quantile), function(X) F3_x(data=d[Quantile==X], x=0.85, data2=combo[[I]]))-> WOW2[[I]][['UKB_AFR']][[4]];list.append(WOW2[[I]][['UKB_AFR']][[4]], F3_x(data=d, x=0.85,data2=combo[[I]]))-> WOW2[[I]][['UKB_AFR']][[4]]
	lapply(levels(d$Quantile), function(X) F3_x(data=d[Quantile==X], x=0.8, data2=combo[[I]]))-> WOW2[[I]][['UKB_AFR']][[5]];list.append(WOW2[[I]][['UKB_AFR']][[5]], F3_x(data=d, x=0.8,data2=combo[[I]]))-> WOW2[[I]][['UKB_AFR']][[5]]
        names(WOW2[[I]][['UKB_AFR']])<-c('0.975perc', '0.95perc', '0.9perc','0.85perc',  '0.8perc')
        names(WOW2[[I]][['UKB_AFR']][['0.9perc']])<-c(levels(d$Quantile), 'all')
        names(WOW2[[I]][['UKB_AFR']][['0.95perc']])<-c(levels(d$Quantile), 'all')
        names(WOW2[[I]][['UKB_AFR']][['0.975perc']])<-c(levels(d$Quantile), 'all')
        names(WOW2[[I]][['UKB_AFR']][['0.85perc']])<-c(levels(d$Quantile), 'all')
	names(WOW2[[I]][['UKB_AFR']][['0.8perc']])<-c(levels(d$Quantile), 'all')
        lapply(unique(e$Quantile), function(X) F3_x(data=e[Quantile==X], x=0.975,  data2=combo[[I]]))-> WOW2[[I]][['UKB_EUR']][[1]];
        lapply(unique(e$Quantile), function(X) F3_x(data=e[Quantile==X], x=0.95, data2=combo[[I]]))-> WOW2[[I]][['UKB_EUR']][[2]];
        lapply(unique(e$Quantile), function(X) F3_x(data=e[Quantile==X], x=0.9, data2=combo[[I]]))-> WOW2[[I]][['UKB_EUR']][[3]];
        lapply(unique(e$Quantile), function(X) F3_x(data=e[Quantile==X], x=0.85, data2=combo[[I]]))-> WOW2[[I]][['UKB_EUR']][[4]];
	lapply(unique(e$Quantile), function(X) F3_x(data=e[Quantile==X], x=0.8, data2=combo[[I]]))-> WOW2[[I]][['UKB_EUR']][[5]];
        names(WOW2[[I]][['UKB_EUR']])<-c('0.975perc', '0.95perc', '0.9perc', '0.85perc', '0.8perc')
        names(WOW2[[I]][['UKB_EUR']][['0.9perc']])<-unique(e$Quantile)
        names(WOW2[[I]][['UKB_EUR']][['0.95perc']])<-unique(e$Quantile);names(WOW2[[I]][['UKB_EUR']][['0.975perc']])<-unique(e$Quantile)
        names(WOW2[[I]][['UKB_EUR']][['0.85perc']])<-unique(e$Quantile);names(WOW2[[I]][['UKB_EUR']][['0.8perc']])<-unique(e$Quantile)
        lapply(unique(f$Quantile), function(X) F3_x(data=f[Quantile==X], x=0.975,  data2=combo[[I]]))-> WOW2[[I]][['pennBB_EA']][[1]];
        lapply(unique(f$Quantile), function(X) F3_x(data=f[Quantile==X], x=0.95, data2=combo[[I]]))-> WOW2[[I]][['pennBB_EA']][[2]];
        lapply(unique(f$Quantile), function(X) F3_x(data=f[Quantile==X], x=0.9, data2=combo[[I]]))-> WOW2[[I]][['pennBB_EA']][[3]];
        lapply(unique(f$Quantile), function(X) F3_x(data=f[Quantile==X], x=0.85, data2=combo[[I]]))-> WOW2[[I]][['pennBB_EA']][[4]];
	lapply(unique(f$Quantile), function(X) F3_x(data=f[Quantile==X], x=0.8, data2=combo[[I]]))-> WOW2[[I]][['pennBB_EA']][[5]];
        names(WOW2[[I]][['pennBB_EA']])<-c('0.975perc', '0.95perc', '0.9perc','0.85perc', '0.8perc');names(WOW2[[I]][['pennBB_EA']][['0.9perc']])<-unique(f$Quantile)
        names(WOW2[[I]][['pennBB_EA']][['0.95perc']])<-unique(f$Quantile);names(WOW2[[I]][['pennBB_EA']][['0.975perc']])<-unique(f$Quantile)
        names(WOW2[[I]][['pennBB_EA']][['0.85perc']])<-unique(f$Quantile);names(WOW2[[I]][['pennBB_EA']][['0.8perc']])<-unique(f$Quantile)
	lapply(levels(g$Quantile), function(X) F3_x(data=g[Quantile==X], x=0.975,data2=combo[[I]]))-> WOW2[[I]][['JHS']][[1]];list.append(WOW2[[I]][['JHS']][[1]], F3_x(data=d, x=0.975,data2=combo[[I]]))-> WOW2[[I]][['JHS']][[1]]
        lapply(levels(g$Quantile), function(X) F3_x(data=g[Quantile==X], x=0.95,data2=combo[[I]]))-> WOW2[[I]][['JHS']][[2]];list.append(WOW2[[I]][['JHS']][[2]], F3_x(data=d, x=0.95,data2=combo[[I]]))-> WOW2[[I]][['JHS']][[2]]
        lapply(levels(g$Quantile), function(X) F3_x(data=g[Quantile==X], x=0.9,data2=combo[[I]]))-> WOW2[[I]][['JHS']][[3]];list.append(WOW2[[I]][['JHS']][[3]], F3_x(data=d, x=0.9,data2=combo[[I]]))-> WOW2[[I]][['JHS']][[3]]
        lapply(levels(g$Quantile), function(X) F3_x(data=g[Quantile==X], x=0.85,data2=combo[[I]]))-> WOW2[[I]][['JHS']][[4]];list.append(WOW2[[I]][['JHS']][[4]], F3_x(data=d, x=0.85,data2=combo[[I]]))-> WOW2[[I]][['JHS']][[4]]
        lapply(levels(g$Quantile), function(X) F3_x(data=g[Quantile==X], x=0.8,data2=combo[[I]]))-> WOW2[[I]][['JHS']][[5]];list.append(WOW2[[I]][['JHS']][[5]], F3_x(data=d, x=0.8,data2=combo[[I]]))-> WOW2[[I]][['JHS']][[5]]
        names(WOW2[[I]][['JHS']])<-c('0.975perc', '0.95perc', '0.9perc','0.85perc', '0.8perc');names(WOW2[[I]][['JHS']][['0.9perc']])<-c(levels(g$Quantile), 'all')
        names(WOW2[[I]][['JHS']][['0.95perc']])<-c(levels(g$Quantile), 'all');names(WOW2[[I]][['JHS']][['0.975perc']])<-c(levels(g$Quantile), 'all');names(WOW2[[I]][['JHS']][['0.85perc']])<-c(levels(g$Quantile), 'all'); names(WOW2[[I]][['JHS']][['0.8perc']])<-c(levels(g$Quantile), 'all')
        cat(I)
        cat('\n')
}

OR_table<-vector('list', length(combo))
names(OR_table)<-names(B_WHI)
for(I in names(AA)){
	OR_table[[I]]<-rbind(
	data.table(Dataset=c(rep("WHI",6),rep("pennBB_AA",4), rep('UKB_AFR', 5), rep('UKB_EUR', 1), rep('pennBB_EA',1), rep("JHS",3)),
	EUR_ANC=c(na.omit(c(unlist(lapply(names(WOW[[I]][['WHI']][['0.975perc']]), function(X) median(test[[I]][Dataset=='WHI'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='WHI'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW[[I]][['pennBB_AA']][['0.975perc']]), function(X) median(test[[I]][Dataset=='pennBB_AA'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='pennBB_AA'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW[[I]][['UKB_AFR']][['0.975perc']]), function(X) median(test[[I]][Dataset=='UKB_AFR'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='UKB_AFR'][, EUR_ANC]))), 
	median(test[[I]][Dataset=='UKB_EUR'][Quantile==1][, EUR_ANC]), median(test[[I]][Dataset=='pennBB_EA'][Quantile==1][, EUR_ANC]), 
	na.omit(c(unlist(lapply(names(WOW[[I]][['JHS']][['0.975perc']]), function(X) median(test[[I]][Dataset=='JHS'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='JHS'][, EUR_ANC])))), 
	Alpha=0.975, 
	OR=c(unlist(WOW[[I]][['WHI']][['0.975perc']]), unlist(WOW[[I]][['pennBB_AA']][['0.975perc']]), unlist(WOW[[I]][['UKB_AFR']][['0.975perc']]), unlist(WOW[[I]][['UKB_EUR']][['0.975perc']]), unlist(WOW[[I]][['pennBB_EA']][['0.975perc']]), unlist(WOW[[I]][['JHS']][['0.975perc']])),
	Quantile=c(names(WOW[[I]][['WHI']][['0.975perc']]), names(WOW[[I]][['pennBB_AA']][['0.975perc']]),names(WOW[[I]][['UKB_AFR']][['0.975perc']]), as.character(1),as.character(1), names(WOW[[I]][['JHS']][['0.975perc']])), Prun=I),
	data.table(
	Dataset=c(rep("WHI",6,rep("pennBB_AA",4), rep('UKB_AFR', 5), rep('UKB_EUR', 1), rep('pennBB_EA',1)),rep("JHS",3)),
	EUR_ANC=c(na.omit(c(unlist(lapply(names(WOW[[I]][['WHI']][['0.95perc']]), function(X) median(test[[I]][Dataset=='WHI'][Quantile==X][, EUR_ANC]))),median(combo[[I]][Dataset=='WHI'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW[[I]][['pennBB_AA']][['0.95perc']]), function(X) median(test[[I]][Dataset=='pennBB_AA'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='pennBB_AA'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW[[I]][['UKB_AFR']][['0.95perc']]), function(X) median(test[[I]][Dataset=='UKB_AFR'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='UKB_AFR'][, EUR_ANC]))),
	median(test[[I]][Dataset=='UKB_EUR'][Quantile==1][, EUR_ANC]), median(test[[I]][Dataset=='pennBB_EA'][Quantile==1][, EUR_ANC]), 
	na.omit(c(unlist(lapply(names(WOW[[I]][['JHS']][['0.95perc']]), function(X) median(test[[I]][Dataset=='JHS'][Quantile==X][, EUR_ANC]))),median(combo[[I]][Dataset=='JHS'][, EUR_ANC])))),
	Alpha=0.95,
	OR=c(unlist(WOW[[I]][['WHI']][['0.95perc']]), unlist(WOW[[I]][['pennBB_AA']][['0.95perc']]), unlist(WOW[[I]][['UKB_AFR']][['0.95perc']]), unlist(WOW[[I]][['UKB_EUR']][['0.95perc']]), unlist(WOW[[I]][['pennBB_EA']][['0.95perc']]), unlist(WOW[[I]][['JHS']][['0.95perc']])),
	Quantile=c(names(WOW[[I]][['WHI']][['0.95perc']]), names(WOW[[I]][['pennBB_AA']][['0.95perc']]),names(WOW[[I]][['UKB_AFR']][['0.95perc']]), as.character(1),as.character(1), names(WOW[[I]][['JHS']][['0.97perc']])),Prun=I),
	data.table(
	Dataset=c(rep("WHI",6),rep("pennBB_AA",4), rep('UKB_AFR', 5), rep('UKB_EUR', 1), rep('pennBB_EA',1), rep("JHS",3)),
	EUR_ANC=c(na.omit(c(unlist(lapply(names(WOW[[I]][['WHI']][['0.9perc']]), function(X) median(test[[I]][Dataset=='WHI'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='WHI'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW[[I]][['pennBB_AA']][['0.9perc']]), function(X) median(test[[I]][Dataset=='pennBB_AA'][Quantile==X][, EUR_ANC]))),median(combo[[I]][Dataset=='pennBB_AA'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW[[I]][['UKB_AFR']][['0.9perc']]), function(X) median(test[[I]][Dataset=='UKB_AFR'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='UKB_AFR'][, EUR_ANC]))),
	median(test[[I]][Dataset=='UKB_EUR'][Quantile==1][, EUR_ANC]), median(test[[I]][Dataset=='pennBB_EA'][Quantile==1][, EUR_ANC]),
	na.omit(c(unlist(lapply(names(WOW[[I]][['JHS']][['0.9perc']]), function(X) median(test[[I]][Dataset=='JHS'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='JHS'][, EUR_ANC])))),
	Alpha=0.9,
	OR=c(unlist(WOW[[I]][['WHI']][['0.9perc']]), unlist(WOW[[I]][['pennBB_AA']][['0.9perc']]), unlist(WOW[[I]][['UKB_AFR']][['0.9perc']]), unlist(WOW[[I]][['UKB_EUR']][['0.9perc']]), unlist(WOW[[I]][['pennBB_EA']][['0.9perc']]),unlist(WOW[[I]][['JHS']][['0.9perc']])),
	Quantile=c(names(WOW[[I]][['WHI']][['0.9perc']]), names(WOW[[I]][['pennBB_AA']][['0.9perc']]),names(WOW[[I]][['UKB_AFR']][['0.9perc']]), as.character(1),as.character(1),names(WOW[[I]][['JHS']][['0.9perc']])), Prun=I),
	data.table(
        Dataset=c(rep("WHI",6),rep("pennBB_AA",4), rep('UKB_AFR', 5), rep('UKB_EUR', 1), rep('pennBB_EA',1), rep("JHS",3)),
        EUR_ANC=c(na.omit(c(unlist(lapply(names(WOW[[I]][['WHI']][['0.85perc']]), function(X) median(test[[I]][Dataset=='WHI'][Quantile==X][, EUR_ANC]))),  median(combo[[I]][Dataset=='WHI'][, EUR_ANC]))),
        na.omit(c(unlist(lapply(names(WOW[[I]][['pennBB_AA']][['0.85perc']]), function(X) median(test[[I]][Dataset=='pennBB_AA'][Quantile==X][, EUR_ANC]))),median(combo[[I]][Dataset=='pennBB_AA'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW[[I]][['UKB_AFR']][['0.85perc']]), function(X) median(test[[I]][Dataset=='UKB_AFR'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='UKB_AFR'][, EUR_ANC]))),
        median(test[[I]][Dataset=='UKB_EUR'][Quantile==1][, EUR_ANC]), median(test[[I]][Dataset=='pennBB_EA'][Quantile==1][, EUR_ANC]),
	na.omit(c(unlist(lapply(names(WOW[[I]][['JHS']][['0.85perc']]), function(X) median(test[[I]][Dataset=='JHS'][Quantile==X][, EUR_ANC]))),  median(combo[[I]][Dataset=='JHS'][, EUR_ANC])))),
        Alpha=0.85,
        OR=c(unlist(WOW[[I]][['WHI']][['0.85perc']]), unlist(WOW[[I]][['pennBB_AA']][['0.85perc']]), unlist(WOW[[I]][['UKB_AFR']][['0.85perc']]), unlist(WOW[[I]][['UKB_EUR']][['0.85perc']]), unlist(WOW[[I]][['pennBB_EA']][['0.85perc']]), unlist(WOW[[I]][['JHS']][['0.85perc']])),
        Quantile=c(names(WOW[[I]][['WHI']][['0.85perc']]), names(WOW[[I]][['pennBB_AA']][['0.85perc']]),names(WOW[[I]][['UKB_AFR']][['0.85perc']]), as.character(1),as.character(1), names(WOW[[I]][['JHS']][['0.85perc']])),Prun=I),
        data.table(
        Dataset=c(rep("WHI",6),rep("pennBB_AA",4), rep('UKB_AFR', 5), rep('UKB_EUR', 1), rep('pennBB_EA',1), rep("JHS",3)),
        EUR_ANC=c(na.omit(c(unlist(lapply(names(WOW[[I]][['WHI']][['0.8perc']]), function(X) median(test[[I]][Dataset=='WHI'][Quantile==X][, EUR_ANC]))),median(combo[[I]][Dataset=='WHI'][, EUR_ANC]))),
        na.omit(c(unlist(lapply(names(WOW[[I]][['pennBB_AA']][['0.8perc']]), function(X) median(test[[I]][Dataset=='pennBB_AA'][Quantile==X][, EUR_ANC]))),median(combo[[I]][Dataset=='pennBB_AA'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW[[I]][['UKB_AFR']][['0.8perc']]), function(X) median(test[[I]][Dataset=='UKB_AFR'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='UKB_AFR'][, EUR_ANC]))),
        median(test[[I]][Dataset=='UKB_EUR'][Quantile==1][, EUR_ANC]), median(test[[I]][Dataset=='pennBB_EA'][Quantile==1][, EUR_ANC]),
	na.omit(c(unlist(lapply(names(WOW[[I]][['JHS']][['0.8perc']]), function(X) median(test[[I]][Dataset=='JHS'][Quantile==X][, EUR_ANC]))),median(combo[[I]][Dataset=='JHS'][, EUR_ANC])))),
        Alpha=0.8,
        OR=c(unlist(WOW[[I]][['WHI']][['0.8perc']]), unlist(WOW[[I]][['pennBB_AA']][['0.8perc']]), unlist(WOW[[I]][['UKB_AFR']][['0.8perc']]), unlist(WOW[[I]][['UKB_EUR']][['0.8perc']]), unlist(WOW[[I]][['pennBB_EA']][['0.8perc']]), unlist(WOW[[I]][['JHS']][['0.8perc']])),
        Quantile=c(names(WOW[[I]][['WHI']][['0.8perc']]), names(WOW[[I]][['pennBB_AA']][['0.8perc']]),names(WOW[[I]][['UKB_AFR']][['0.8perc']]), as.character(1),as.character(1),names(WOW[[I]][['JHS']][['0.8perc']])), Prun=I)
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
        data.table(Dataset=c(rep("WHI",6),rep("pennBB_AA",4), rep('UKB_AFR', 5), rep('UKB_EUR', 1), rep('pennBB_EA',1), rep("JHS",3)),
        EUR_ANC=c(na.omit(c(unlist(lapply(names(WOW2[[I]][['WHI']][['0.975perc']]), function(X) median(test2[[I]][Dataset=='WHI'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='WHI'][, EUR_ANC]))),
        na.omit(c(unlist(lapply(names(WOW2[[I]][['pennBB_AA']][['0.975perc']]), function(X) median(test2[[I]][Dataset=='pennBB_AA'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='pennBB_AA'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW2[[I]][['UKB_AFR']][['0.975perc']]), function(X) median(test2[[I]][Dataset=='UKB_AFR'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='UKB_AFR'][, EUR_ANC]))),
        median(test2[[I]][Dataset=='UKB_EUR'][Quantile==1][, EUR_ANC]), median(test2[[I]][Dataset=='pennBB_EA'][Quantile==1][, EUR_ANC]),
	na.omit(c(unlist(lapply(names(WOW2[[I]][['JHS']][['0.975perc']]), function(X) median(test2[[I]][Dataset=='JHS'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='JHS'][, EUR_ANC])))),
        Alpha=0.975,
        OR=c(sapply(WOW2[[I]][['WHI']][['0.975perc']], function(X) X$F),sapply(WOW2[[I]][['pennBB_AA']][['0.975perc']], function(X) X$F), sapply(WOW2[[I]][['UKB_AFR']][['0.975perc']], function(X) X$F), sapply(WOW2[[I]][['UKB_EUR']][['0.975perc']], function(X) X$F), sapply(WOW2[[I]][['pennBB_EA']][['0.975perc']], function(X) X$F), sapply(WOW2[[I]][['JHS']][['0.975perc']], function(X) X$F)),
	Quantile=c(names(WOW2[[I]][['WHI']][['0.975perc']]), names(WOW2[[I]][['pennBB_AA']][['0.975perc']]),names(WOW2[[I]][['UKB_AFR']][['0.975perc']]), as.character(1),as.character(1), names(WOW2[[I]][['JHS']][['0.975perc']])), Prun=I),
	data.table(Dataset=c(rep("WHI",6),rep("pennBB_AA",4), rep('UKB_AFR', 5), rep('UKB_EUR', 1), rep('pennBB_EA',1), rep("JHS",3)),
        EUR_ANC=c(na.omit(c(unlist(lapply(names(WOW2[[I]][['WHI']][['0.95perc']]), function(X) median(test2[[I]][Dataset=='WHI'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='WHI'][, EUR_ANC]))),
        na.omit(c(unlist(lapply(names(WOW2[[I]][['pennBB_AA']][['0.95perc']]), function(X) median(test2[[I]][Dataset=='pennBB_AA'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='pennBB_AA'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW2[[I]][['UKB_AFR']][['0.95perc']]), function(X) median(test2[[I]][Dataset=='UKB_AFR'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='UKB_AFR'][, EUR_ANC]))),
        median(test2[[I]][Dataset=='UKB_EUR'][Quantile==1][, EUR_ANC]), median(test2[[I]][Dataset=='pennBB_EA'][Quantile==1][, EUR_ANC]),
	na.omit(c(unlist(lapply(names(WOW2[[I]][['JHS']][['0.95perc']]), function(X) median(test2[[I]][Dataset=='JHS'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='JHS'][, EUR_ANC])))),
        Alpha=0.95,	
        OR=c(sapply(WOW2[[I]][['WHI']][['0.95perc']], function(X) X$F),sapply(WOW2[[I]][['pennBB_AA']][['0.95perc']], function(X) X$F), sapply(WOW2[[I]][['UKB_AFR']][['0.95perc']], function(X) X$F), sapply(WOW2[[I]][['UKB_EUR']][['0.95perc']], function(X) X$F), sapply(WOW2[[I]][['pennBB_EA']][['0.95perc']], function(X) X$F), sapply(WOW2[[I]][['JHS']][['0.95perc']], function(X) X$F)),
        Quantile=c(names(WOW2[[I]][['WHI']][['0.95perc']]), names(WOW2[[I]][['pennBB_AA']][['0.95perc']]),names(WOW2[[I]][['UKB_AFR']][['0.95perc']]), as.character(1),as.character(1), names(WOW2[[I]][['JHS']][['0.95perc']])), Prun=I),
        data.table(Dataset=c(rep("WHI",6),rep("pennBB_AA",4), rep('UKB_AFR', 5), rep('UKB_EUR', 1), rep('pennBB_EA',1), rep("JHS",3)),
        EUR_ANC=c(na.omit(c(unlist(lapply(names(WOW2[[I]][['WHI']][['0.9perc']]), function(X) median(test2[[I]][Dataset=='WHI'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='WHI'][, EUR_ANC]))),
        na.omit(c(unlist(lapply(names(WOW2[[I]][['pennBB_AA']][['0.9perc']]), function(X) median(test2[[I]][Dataset=='pennBB_AA'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='pennBB_AA'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW2[[I]][['UKB_AFR']][['0.9perc']]), function(X) median(test2[[I]][Dataset=='UKB_AFR'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='UKB_AFR'][, EUR_ANC]))),
        median(test2[[I]][Dataset=='UKB_EUR'][Quantile==1][, EUR_ANC]), median(test2[[I]][Dataset=='pennBB_EA'][Quantile==1][, EUR_ANC]),
	na.omit(c(unlist(lapply(names(WOW2[[I]][['JHS']][['0.9perc']]), function(X) median(test2[[I]][Dataset=='JHS'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='JHS'][, EUR_ANC])))), 
        Alpha=0.9,
        OR=c(sapply(WOW2[[I]][['WHI']][['0.9perc']], function(X) X$F),sapply(WOW2[[I]][['pennBB_AA']][['0.9perc']], function(X) X$F), sapply(WOW2[[I]][['UKB_AFR']][['0.9perc']], function(X) X$F), sapply(WOW2[[I]][['UKB_EUR']][['0.9perc']], function(X) X$F), sapply(WOW2[[I]][['pennBB_EA']][['0.9perc']], function(X) X$F), sapply(WOW2[[I]][['JHS']][['0.9perc']], function(X) X$F)),
        Quantile=c(names(WOW2[[I]][['WHI']][['0.9perc']]), names(WOW2[[I]][['pennBB_AA']][['0.9perc']]),names(WOW2[[I]][['UKB_AFR']][['0.9perc']]), as.character(1),as.character(1), names(WOW2[[I]][['JHS']][['0.9perc']])), Prun=I),
        data.table(Dataset=c(rep("WHI",6),rep("pennBB_AA",4), rep('UKB_AFR', 5), rep('UKB_EUR', 1), rep('pennBB_EA',1), rep("JHS",3)),
	EUR_ANC=c(na.omit(c(unlist(lapply(names(WOW2[[I]][['WHI']][['0.85perc']]), function(X) median(test2[[I]][Dataset=='WHI'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='WHI'][, EUR_ANC]))),
        na.omit(c(unlist(lapply(names(WOW2[[I]][['pennBB_AA']][['0.85perc']]), function(X) median(test2[[I]][Dataset=='pennBB_AA'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='pennBB_AA'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW2[[I]][['UKB_AFR']][['0.85perc']]), function(X) median(test2[[I]][Dataset=='UKB_AFR'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='UKB_AFR'][, EUR_ANC]))),
        median(test2[[I]][Dataset=='UKB_EUR'][Quantile==1][, EUR_ANC]), median(test2[[I]][Dataset=='pennBB_EA'][Quantile==1][, EUR_ANC]),
	na.omit(c(unlist(lapply(names(WOW2[[I]][['JHS']][['0.85perc']]), function(X) median(test2[[I]][Dataset=='JHS'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='JHS'][, EUR_ANC])))),
        Alpha=0.85,
        OR=c(sapply(WOW2[[I]][['WHI']][['0.85perc']], function(X) X$F),sapply(WOW2[[I]][['pennBB_AA']][['0.85perc']], function(X) X$F), sapply(WOW2[[I]][['UKB_AFR']][['0.85perc']], function(X) X$F), sapply(WOW2[[I]][['UKB_EUR']][['0.85perc']], function(X) X$F), sapply(WOW2[[I]][['pennBB_EA']][['0.85perc']], function(X) X$F), sapply(WOW2[[I]][['JHS']][['0.85perc']], function(X) X$F)),
        Quantile=c(names(WOW2[[I]][['WHI']][['0.85perc']]), names(WOW2[[I]][['pennBB_AA']][['0.85perc']]),names(WOW2[[I]][['UKB_AFR']][['0.85perc']]), as.character(1),as.character(1), names(WOW2[[I]][['JHS']][['0.85perc']])), Prun=I),
        data.table(Dataset=c(rep("WHI",6),rep("pennBB_AA",4), rep('UKB_AFR', 5), rep('UKB_EUR', 1), rep('pennBB_EA',1), rep("JHS",3)),
	EUR_ANC=c(na.omit(c(unlist(lapply(names(WOW2[[I]][['WHI']][['0.8perc']]), function(X) median(test2[[I]][Dataset=='WHI'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='WHI'][, EUR_ANC]))),
        na.omit(c(unlist(lapply(names(WOW2[[I]][['pennBB_AA']][['0.8perc']]), function(X) median(test2[[I]][Dataset=='pennBB_AA'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='pennBB_AA'][, EUR_ANC]))),
	na.omit(c(unlist(lapply(names(WOW2[[I]][['UKB_AFR']][['0.8perc']]), function(X) median(test2[[I]][Dataset=='UKB_AFR'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='UKB_AFR'][, EUR_ANC]))),
        median(test2[[I]][Dataset=='UKB_EUR'][Quantile==1][, EUR_ANC]), median(test2[[I]][Dataset=='pennBB_EA'][Quantile==1][, EUR_ANC]),
	na.omit(c(unlist(lapply(names(WOW2[[I]][['JHS']][['0.8perc']]), function(X) median(test2[[I]][Dataset=='JHS'][Quantile==X][, EUR_ANC]))), median(combo[[I]][Dataset=='JHS'][, EUR_ANC])))),
        Alpha=0.8,
        OR=c(sapply(WOW2[[I]][['WHI']][['0.8perc']], function(X) X$F),sapply(WOW2[[I]][['pennBB_AA']][['0.8perc']], function(X) X$F), sapply(WOW2[[I]][['UKB_AFR']][['0.8perc']], function(X) X$F), sapply(WOW2[[I]][['UKB_EUR']][['0.8perc']], function(X) X$F), sapply(WOW2[[I]][['pennBB_EA']][['0.8perc']], function(X) X$F), sapply(WOW2[[I]][['JHS']][['0.8perc']], function(X) X$F)),
        Quantile=c(names(WOW2[[I]][['WHI']][['0.8perc']]), names(WOW2[[I]][['pennBB_AA']][['0.8perc']]),names(WOW2[[I]][['UKB_AFR']][['0.8perc']]), as.character(1),as.character(1), names(WOW2[[I]][['JHS']][['0.8perc']])), Prun=I)
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
fwrite(OR_table2, file='Odds_Ratio_all_datasets.txt', quote=F, sep="\t")
fwrite(OR_table3, file='Odds_Ratio_all_datasets_v2.txt', quote=F, sep="\t")
as.factor(OR_table[[63]]$Alpha)-> OR_table[[63]]$Alpha
OR_table[[63]][, logOR:=log(OR)]

as.factor(OR2_table[[63]]$Alpha)-> OR2_table[[63]]$Alpha
OR2_table[[63]][, logOR:=log(OR)]

png(paste0('figs/logOR_test_', names(AA)[63], '.png'))
one<-ggplot(OR_table[[63]], aes(x=EUR_ANC, y=logOR, shape=Dataset, colour=Alpha)) +
geom_point(size=2) + labs(title="Log odds-ratio of P(>=Xth HEIGHT quantile|>=Xth PRS quantile)", y="log(OR)") + geom_hline(yintercept=0, linetype="dashed", color = "orange") + stat_smooth(method = "lm", se=F, aes(group=Alpha))
print(one)
dev.off()

png(paste0('figs/logOR_v2_test', names(AA)[63], '.png'))
one<-ggplot(OR2_table[[63]], aes(x=EUR_ANC, y=logOR, shape=Dataset, colour=Alpha)) +
geom_point(size=2) + labs(title="Log odds-ratio of P(>=Xth HEIGHT quantile|>=Xth PRS quantile)", y="log(OR)") + geom_hline(yintercept=0, linetype="dashed", color = "orange") + stat_smooth(method = "lm", se=F, aes(group=Alpha))
print(one)
dev.off()

png(paste0('figs/OR_test_', names(AA)[63], '.png'))
two<-ggplot(OR_table[[63]], aes(x=EUR_ANC, y=OR, shape=Dataset, colour=Alpha)) +
geom_point(size=2) + labs(title="Odds-ratio of P(>=Xth HEIGHT quantile|>=Xth PRS quantile)", y="OR") + geom_hline(yintercept=1, linetype="dashed", color = "orange") + stat_smooth(method = "lm", se=F, aes(group=Alpha))
print(two)
dev.off()

png(paste0('figs/OR2_test_' , names(AA)[63],'_v2.png'))
two<-ggplot(OR2_table[[63]], aes(x=EUR_ANC, y=OR, shape=Dataset, colour=Alpha)) +
geom_point(size=3) + labs(title="Odds-ratio of P(>=Xth HEIGHT quantile|>=Xth PRS quantile)", y="OR") + geom_hline(yintercept=1, linetype="dashed", color = "orange") + stat_smooth(method = "lm", se=F, aes(group=Alpha)) + 
theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),  axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
print(two)
dev.off()

###area under the curve
I<-names(AA)[63]
lapply(seq(0.2, 0.98, by=0.01), function(X) F3_x(data=combo[[I]][Dataset=='WHI'],data2=combo[[I]], x=X))-> my.auc.WHI
lapply(seq(0.2, 0.98, by=0.01), function(X) F3_x(data=combo[[I]][Dataset=='JHS'],data2=combo[[I]], x=X))-> my.auc.JHS
lapply(seq(0.2, 0.98, by=0.01), function(X) F3_x(data=combo[[I]][Dataset=='UKB_EUR'],data2=combo[[I]], x=X))-> my.auc.UKB_EUR
lapply(seq(0.2, 0.98, by=0.01), function(X) F3_x(data=combo[[I]][Dataset=='UKB_AFR'],data2=combo[[I]], x=X))-> my.auc.UKB_AFR
lapply(seq(0.2, 0.98, by=0.01), function(X) F3_x(data=combo[[I]][Dataset=='pennBB_AA'],data2=combo[[I]], x=X))-> my.auc.pennBB_AFR
lapply(seq(0.2, 0.98, by=0.01), function(X) F3_x(data=combo[[I]][Dataset=='pennBB_EA'],data2=combo[[I]], x=X))-> my.auc.pennBB_EUR
my.dt.WHI<-data.table(F=unlist(lapply(my.auc.WHI, function(x) x$F)), G=unlist(lapply(my.auc.WHI, function(x) x$G)), TP=unlist(lapply(my.auc.WHI, function(x) x$f)), FP=unlist(lapply(my.auc.WHI, function(x) x$g)), Alpha=seq(0.2, 0.98, by=0.01), Dataset='WHI', Ref="UKB_EUR")
my.dt.JHS<-data.table(F=unlist(lapply(my.auc.JHS, function(x) x$F)), G=unlist(lapply(my.auc.JHS, function(x) x$G)), TP=unlist(lapply(my.auc.JHS, function(x) x$f)), FP=unlist(lapply(my.auc.JHS, function(x) x$g)), Alpha=seq(0.2, 0.98, by=0.01), Dataset='JHS', Ref="UKB_EUR")
my.dt.UKB_EUR<-data.table(F=unlist(lapply(my.auc.UKB_EUR, function(x) x$F)), G=unlist(lapply(my.auc.UKB_EUR, function(x) x$G)), TP=unlist(lapply(my.auc.UKB_EUR, function(x) x$f)), FP=unlist(lapply(my.auc.UKB_EUR, function(x) x$g)), Alpha=seq(0.2, 0.98, by=0.01), Dataset='UKB_EUR', Ref="UKB_EUR")
my.dt.UKB_AFR<-data.table(F=unlist(lapply(my.auc.UKB_AFR, function(x) x$F)), G=unlist(lapply(my.auc.UKB_AFR, function(x) x$G)), TP=unlist(lapply(my.auc.UKB_AFR, function(x) x$f)), FP=unlist(lapply(my.auc.UKB_AFR, function(x) x$g)), Alpha=seq(0.2, 0.98, by=0.01), Dataset='UKB_AFR', Ref="UKB_EUR")
my.dt.pennBB_EUR<-data.table(F=unlist(lapply(my.auc.pennBB_EUR, function(x) x$F)), G=unlist(lapply(my.auc.pennBB_EUR, function(x) x$G)), TP=unlist(lapply(my.auc.pennBB_EUR, function(x) x$f)), FP=unlist(lapply(my.auc.pennBB_EUR, function(x) x$g)), Alpha=seq(0.2, 0.98, by=0.01), Dataset='pennBB_EA', Ref="UKB_EUR")
my.dt.pennBB_AFR<-data.table(F=unlist(lapply(my.auc.pennBB_AFR, function(x) x$F)), G=unlist(lapply(my.auc.pennBB_AFR, function(x) x$G)), TP=unlist(lapply(my.auc.pennBB_AFR, function(x) x$f)), FP=unlist(lapply(my.auc.pennBB_AFR, function(x) x$g)), Alpha=seq(0.2, 0.98, by=0.01), Dataset='pennBB_AA', Ref="UKB_EUR")
rbind(my.dt.WHI, my.dt.JHS,my.dt.UKB_EUR, my.dt.UKB_AFR, my.dt.pennBB_AFR, my.dt.pennBB_EUR)-> my.dt
#rbind(my.dt.WHI, my.dt.UKB_EUR, my.dt.UKB_AFR)-> my.dt
#

lapply(seq(0.2, 0.98, by=0.01), function(X) F4_x(data=combo[[I]][Dataset=='WHI'], x=X))-> my.auc.WHI.f4
lapply(seq(0.2, 0.98, by=0.01), function(X) F4_x(data=combo[[I]][Dataset=='JHS'], x=X))-> my.auc.JHS.f4
lapply(seq(0.2, 0.98, by=0.01), function(X) F4_x(data=combo[[I]][Dataset=='UKB_EUR'], x=X))-> my.auc.UKB_EUR.f4
lapply(seq(0.2, 0.98, by=0.01), function(X) F4_x(data=combo[[I]][Dataset=='UKB_AFR'], x=X))-> my.auc.UKB_AFR.f4
lapply(seq(0.2, 0.98, by=0.01), function(X) F4_x(data=combo[[I]][Dataset=='pennBB_AA'], x=X))-> my.auc.pennBB_AFR.f4
lapply(seq(0.2, 0.98, by=0.01), function(X) F4_x(data=combo[[I]][Dataset=='pennBB_EA'], x=X))-> my.auc.pennBB_EUR.f4
my.dt.WHI.f4<-data.table(F=unlist(lapply(my.auc.WHI.f4, function(x) x$F)), G=unlist(lapply(my.auc.WHI.f4, function(x) x$G)), TP=unlist(lapply(my.auc.WHI.f4, function(x) x$f)), FP=unlist(lapply(my.auc.WHI.f4, function(x) x$g)), Alpha=seq(0.2, 0.98, by=0.01), Dataset='WHI', Ref="Matched")
my.dt.JHS.f4<-data.table(F=unlist(lapply(my.auc.JHS.f4, function(x) x$F)), G=unlist(lapply(my.auc.JHS.f4, function(x) x$G)), TP=unlist(lapply(my.auc.JHS.f4, function(x) x$f)), FP=unlist(lapply(my.auc.JHS.f4, function(x) x$g)), Alpha=seq(0.2, 0.98, by=0.01), Dataset='JHS', Ref="Matched")
my.dt.UKB_EUR.f4<-data.table(F=unlist(lapply(my.auc.UKB_EUR.f4, function(x) x$F)), G=unlist(lapply(my.auc.UKB_EUR.f4, function(x) x$G)), TP=unlist(lapply(my.auc.UKB_EUR.f4, function(x) x$f)), FP=unlist(lapply(my.auc.UKB_EUR.f4, function(x) x$g)), Alpha=seq(0.2, 0.98, by=0.01), Dataset='UKB_EUR', Ref="Matched")
my.dt.UKB_AFR.f4<-data.table(F=unlist(lapply(my.auc.UKB_AFR.f4, function(x) x$F)), G=unlist(lapply(my.auc.UKB_AFR.f4, function(x) x$G)), TP=unlist(lapply(my.auc.UKB_AFR.f4, function(x) x$f)), FP=unlist(lapply(my.auc.UKB_AFR.f4, function(x) x$g)), Alpha=seq(0.2, 0.98, by=0.01), Dataset='UKB_AFR', Ref="Matched")
my.dt.pennBB_EUR.f4<-data.table(F=unlist(lapply(my.auc.pennBB_EUR.f4, function(x) x$F)), G=unlist(lapply(my.auc.pennBB_EUR.f4, function(x) x$G)), TP=unlist(lapply(my.auc.pennBB_EUR.f4, function(x) x$f)), FP=unlist(lapply(my.auc.pennBB_EUR.f4, function(x) x$g)), Alpha=seq(0.2, 0.98, by=0.01), Dataset='pennBB_EA', Ref="Matched")
my.dt.pennBB_AFR.f4<-data.table(F=unlist(lapply(my.auc.pennBB_AFR.f4, function(x) x$F)), G=unlist(lapply(my.auc.pennBB_AFR.f4, function(x) x$G)), TP=unlist(lapply(my.auc.pennBB_AFR.f4, function(x) x$f)), FP=unlist(lapply(my.auc.pennBB_AFR.f4, function(x) x$g)), Alpha=seq(0.2, 0.98, by=0.01), Dataset='pennBB_AA', Ref="Matched")
rbind(my.dt.WHI.f4, my.dt.JHS.f4,  my.dt.UKB_EUR.f4, my.dt.UKB_AFR.f4, my.dt.pennBB_AFR.f4, my.dt.pennBB_EUR.f4)-> my.dt.f4
#rbind(my.dt.WHI.f4, my.dt.UKB_EUR.f4, my.dt.UKB_AFR.f4)-> my.dt.f4
rbind(my.dt, my.dt.f4)-> my.dt


ggplot(my.dt[Alpha>=0.50], aes(x=FP, y=TP, colour=Dataset, group=Alpha)) + geom_point(size=0.3)  + facet_wrap(~Ref, nrow=2, scales='free_y')

ggsave('figs/auc_test.pdf')

#this didn't work well. let's try something else
#I<-names(A)[67]
#my.alpha<-0.9
#data<-combo[[I]][Dataset=='WHI']
#data2<-combo[[I]]
#z<-quantile(data[,HEIGHTX], probs=my.alpha)[[1]]
#z2<-quantile(data2[Dataset=='UKB_EUR'][,PGS], probs=my.alpha)[[1]]
#tp<-combo[[I]][Dataset=='WHI'][,Outcome:=ifelse(HEIGHTX>=z, "tall", "not-tall")]

#pdf('figs/roc_test.pdf')
#plot.roc(Outcome ~ HEIGHTX, tp)
#dev.off() #this doesn't make any sense either. need to continue tomorrow.



###

combo[[63]]-> tmp
melt(tmp[,.(PGS, Std.PRS, HEIGHTX, Res.Height, Dataset)])-> me
ggplot(me, aes(x=value, group=Dataset, color=Dataset))+ geom_density() + facet_wrap(~variable, nrow=4, scales='free') + 
theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15),  axis.text.x=element_text(size=9), axis.text.y=element_text(size=9))
ggsave('figs/test_density.png')

