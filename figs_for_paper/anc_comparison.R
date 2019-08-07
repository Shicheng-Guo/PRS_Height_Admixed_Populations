#compare local ancestry and admixture estimates for all individuals

library(data.table)
library(reshape)
library(ggplot2)
library(dplyr)

#read local ancestry

#HRS_afr

dtsets<-c('HRS_afr', 'WHI', 'JHS', 'ukb_afr')

rfmix_anc<-vector('list', 4)
names(rfmix_anc)<-dtsets

for(I in dtsets){
	ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/', I, '/rfmix_anc_chr', X, '.txt'))))
	rfmix_anc[[I]]<-ancestry %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual
	rfmix_anc[[I]][,Method:='RFmix']
	rfmix_anc[[I]][,Dataset:=I]
	cat(I , ' done\n')
}

rfmix_anc[['WHI']]<-rfmix_anc[['WHI']][, Dataset:='WHI_afr']
rfmix_anc[['JHS']]<-rfmix_anc[['JHS']][, Dataset:='JHS_afr']
rfmix_anc[['ukb_afr']]<-rfmix_anc[['ukb_afr']][, Dataset:='UKB_afr']
#read admixture

admix_anc<-vector('list', 4)
names(admix_anc)<-dtsets

tp<-fread('/project/mathilab/data/HRS/admixture/HRS_AFR_b37_strand_prune_include.2.Q')
colnames(tp)<-c('EUR_ANC', 'AFR_ANC')
tp[,SUBJID:=fread('/project/mathilab/data/HRS/data/HRS_AFR_b37_strand_prune_include.fam')[,V2]]
tp[,Method:='Admixture']
tp[,Dataset:='HRS_afr']
admix_anc[['HRS_afr']]<-tp

tp<-fread('/project/mathilab/data/WHI/admxiture/WHI_b37_strand_prune_include.2.Q')
colnames(tp)<-c('AFR_ANC', 'EUR_ANC')
tp[,SUBJID:=fread('/project/mathilab/data/WHI/data/WHI_b37_strand_prune_include.fam')[,V2]]
tp[,SUBJID:=paste0("0_", SUBJID)]
tp[,Method:='Admixture']
tp[,Dataset:='WHI_afr']
admix_anc[['WHI']]<-tp


tp<-fread('/project/mathilab/data/JHS/admixture/JHS_b37_strand_prune.2.Q')
colnames(tp)<-c('EUR_ANC', 'AFR_ANC')
tp[,SUBJID:=paste0("0_", fread('/project/mathilab/data/JHS/data/JHS_b37_strand_prune.fam')[,V1], ":", fread('/project/mathilab/data/JHS/data/JHS_b37_strand_prune.fam')[,V2])]
tp[,Method:='Admixture']
tp[,Dataset:='JHS_afr']
admix_anc[['JHS']]<-tp

tp<-fread('/project/mathilab/data/UKB/admixture/UKB_AFR_prune.2.Q')
colnames(tp)<-c('AFR_ANC', 'EUR_ANC')
tp[,SUBJID:=rfmix_anc[['ukb_afr']]$SUBJID]
tp[,Method:='Admixture']
tp[,Dataset:='UKB_afr']
admix_anc[['ukb_afr']]<-tp

##combine

do.call(rbind, rfmix_anc)-> rfmix_anc
do.call(rbind, admix_anc)-> admix_anc
setkey(rfmix_anc, SUBJID)
setkey(admix_anc, SUBJID)

rfmix_anc[admix_anc][,AFR_ANC:=NULL][,i.AFR_ANC:=NULL][,i.Method:=NULL][,i.Dataset:=NULL]-> dt
#data.table(SUBJID=rfmix_anc$SUBJID, RFMix=rfmix_anc[, EUR_ANC], Admixture=admix_anc[,EUR_ANC], Dataset=rfmix_anc[,Dataset])
colnames(dt)[5]<-'Admixture'
colnames(dt)[2]<-'RFMix'
na.omit(dt)-> dt
ggplot(dt, aes(x=Admixture, y=RFMix, color=Dataset)) + geom_point(cex=0.2, alpha=0.5) + geom_smooth(method='lm')
ggsave('~/height_prediction/figs_for_paper/anc_comp.pdf')

