#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
fread('output/HRS_AFR.ldpred.bim')-> hrs_afr
colnames(hrs_afr)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')

#fwrite(final2, file="output/HRS_AFR.ldpred.bim", sep="\t", col.names=F)
final_res<-vector('list',22)
for(chr in 1:22){
	#afr<-final2[CHR==chr]
	afr<-hrs_afr[CHR==chr]
	fread(paste0('/project/mathilab/data/UKB/imputed/ukb_imp_chr', chr, '_eur.bim'))-> eur
	colnames(eur)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')
	final_res[[chr]]<-select(merge(afr, eur, by=c('CHR', 'POS')), CHR, SNP.y, V3.x, POS) %>% as.data.table
	cat(chr, '\n')
	gc()
}

do.call(rbind, final_res)-> final_res
colnames(final_res)<-c('CHR', 'SNP_ID', 'V3', 'POS')
final_res[, SNP:=paste0(CHR, ":", POS)]

fread('zcat output/Height.QC.gz')-> sum_stats
merge(sum_stats, final_res, by=c('CHR', 'POS', 'SNP')) %>% as.data.table-> final
select(final, CHR, SNP, V3.x, POS, A2, A1) %>% as.data.table-> final2


fwrite(final2, file="output/HRS_AFR_UKB_EUR.ldpred.bim", sep="\t", col.names=F)

system('cp output/UKB_EUR.ldpred.fam output/HRS_AFR_UKB_EUR.ldpred.fam')
system('cp output/UKB_EUR.ldpred.bed output/HRS_AFR_UKB_EUR.ldpred.bed')

