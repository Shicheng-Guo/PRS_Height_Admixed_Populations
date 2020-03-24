#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
fread('~/height_prediction/ldpred/output/UKB_EUR.ldpred.bim')-> ukb_eur
colnames(ukb_eur)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')

#fwrite(final2, file="output/HRS_AFR.ldpred.bim", sep="\t", col.names=F)
final_res<-vector('list',22)
for(chr in 22:1){
        #afr<-final2[CHR==chr]
        eur2<-ukb_eur[CHR==chr]
        fread(paste0('/project/mathilab/data/UKB/imputed/ukb_imp_chr', chr, '_eur.bim'))-> eur
        colnames(eur)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')
        final_res[[chr]]<-select(merge(eur2, eur, by=c('CHR', 'POS')), CHR, SNP.y, V3.x, POS) %>% as.data.table
	ukb_eur<-ukb_eur[CHR!=chr]
	remove(eur2, eur)
        cat(chr, '\n')
        gc()
}

do.call(rbind, final_res)-> final_res #744077 SNPs
colnames(final_res)<-c('CHR', 'SNP_ID', 'V3', 'POS')
final_res[, SNP:=paste0(CHR, ":", POS)]

fread('zcat output/Height.QC.gz')-> sum_stats
merge(sum_stats, final_res, by=c('CHR', 'POS', 'SNP')) %>% as.data.table-> final #688076 SNPs
select(final, CHR, SNP, V3, POS, A2, A1) %>% as.data.table-> final2

cat('checkpoint\n')
fwrite(final2, file="~/height_prediction/ldpred/output/UKB_EUR_UKB_EUR.ldpred.bim", sep="\t", col.names=F)

system('cp ~/height_prediction/ldpred/output/UKB_EUR.ldpred.fam output/UKB_EUR_UKB_EUR.ldpred.fam')
system('cp ~/height_prediction/ldpred/output/UKB_EUR.ldpred.bed output/UKB_EUR_UKB_EUR.ldpred.bed')
######################
######################
remove(list=ls())
fread('~/height_prediction/ldpred/output/UKB_AFR.ldpred.bim')-> ukb_afr
colnames(ukb_afr)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')

#fwrite(final2, file="output/HRS_AFR.ldpred.bim", sep="\t", col.names=F)
final_res<-vector('list',22)
for(chr in 22:1){
        #afr<-final2[CHR==chr]
        afr2<-ukb_afr[CHR==chr]
        fread(paste0('/project/mathilab/data/UKB/imputed/ukb_imp_chr', chr, '_eur.bim'))-> eur
        colnames(eur)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')
        final_res[[chr]]<-select(merge(afr2, eur, by=c('CHR', 'POS')), CHR, SNP.y, V3.x, POS) %>% as.data.table
        ukb_afr<-ukb_afr[CHR!=chr]
        remove(afr2, eur)
        cat(chr, '\n')
        gc()
}

do.call(rbind, final_res)-> final_res #744077 SNPs
colnames(final_res)<-c('CHR', 'SNP_ID', 'V3', 'POS')
final_res[, SNP:=paste0(CHR, ":", POS)]

fread('zcat output/Height.QC.gz')-> sum_stats
merge(sum_stats, final_res, by=c('CHR', 'POS', 'SNP')) %>% as.data.table-> final #688076 SNPs
select(final, CHR, SNP, V3, POS, A2, A1) %>% as.data.table-> final2

cat('checkpoint\n')
fwrite(final2, file="~/height_prediction/ldpred/output/UKB_AFR_UKB_EUR.ldpred.bim", sep="\t", col.names=F)

system('cp ~/height_prediction/ldpred/output/UKB_AFR.ldpred.fam output/UKB_AFR_UKB_EUR.ldpred.fam')
system('cp ~/height_prediction/ldpred/output/UKB_AFR.ldpred.bed output/UKB_AFR_UKB_EUR.ldpred.bed')
#######################
#######################
remove(list=ls())
fread('~/height_prediction/ldpred/output/HRS_AFR.ldpred.bim')-> hrs_afr
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

do.call(rbind, final_res)-> final_res #2100121 SNPs
colnames(final_res)<-c('CHR', 'SNP_ID', 'V3', 'POS')
final_res[, SNP:=paste0(CHR, ":", POS)]

fread('zcat output/Height.QC.gz')-> sum_stats
merge(sum_stats, final_res, by=c('CHR', 'POS', 'SNP')) %>% as.data.table-> final # 1523862 SNPs
select(final, CHR, SNP, V3, POS, A2, A1) %>% as.data.table-> final2

cat('checkpoint\n')
fwrite(final2, file="~/height_prediction/ldpred/output/HRS_AFR_UKB_EUR.ldpred.bim", sep="\t", col.names=F)

system('cp ~/height_prediction/ldpred/output/HRS_AFR.ldpred.fam output/HRS_AFR_UKB_EUR.ldpred.fam')
system('cp ~/height_prediction/ldpred/output/HRS_AFR.ldpred.bed output/HRS_AFR_UKB_EUR.ldpred.bed')
##############################################
##############################################
remove(list=ls())
fread('~/height_prediction/ldpred/output/HRS_EUR.ldpred.bim')-> hrs_eur
colnames(hrs_eur)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')
#
cat('another check\n')
final_res<-vector('list',22)
for(chr in 22:1){
        #afr<-final2[CHR==chr]
        eur2<-hrs_eur[CHR==chr]
        fread(paste0('/project/mathilab/data/UKB/imputed/ukb_imp_chr', chr, '_eur.bim'))-> eur
        colnames(eur)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')
        final_res[[chr]]<-select(merge(eur2, eur, by=c('CHR', 'POS')), CHR, SNP.y, V3.x, POS) %>% as.data.table
        cat(chr, '\n')
	remove(eur2, eur)
	hrs_eur<-hrs_eur[CHR!=chr]
        gc()
}
cat('checkpoint EUR\n')
do.call(rbind, final_res)-> final_res #2048211 SNPs
colnames(final_res)<-c('CHR', 'SNP_ID', 'V3', 'POS')
final_res[, SNP:=paste0(CHR, ":", POS)]

fread('zcat ~/height_prediction/ldpred/output/Height.QC.gz')-> sum_stats
merge(sum_stats, final_res, by=c('CHR', 'POS', 'SNP')) %>% as.data.table-> final #1527648 SNPs
select(final, CHR, SNP, V3, POS, A2, A1) %>% as.data.table-> final2 


fwrite(final2, file="~/height_prediction/ldpred/output/HRS_EUR_UKB_EUR.ldpred.bim", sep="\t", col.names=F)

system('cp ~/height_prediction/ldpred/output/HRS_EUR.ldpred.fam output/HRS_EUR_UKB_EUR.ldpred.fam')
system('cp ~/height_prediction/ldpred/output/HRS_EUR.ldpred.bed output/HRS_EUR_UKB_EUR.ldpred.bed')
#
#WHI
remove(list=ls())
fread('~/height_prediction/ldpred/output/WHI.ldpred.bim')-> whi
colnames(whi)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')
whi<-whi[CHR %in% 1:22]
whi$CHR<-as.integer(whi$CHR)
#
cat('another check\n')
final_res<-vector('list',22)
for(chr in 22:1){
        #afr<-final2[CHR==chr]
        afr2<-whi[CHR==chr]
        fread(paste0('/project/mathilab/data/UKB/imputed/ukb_imp_chr', chr, '_eur.bim'))-> eur
        colnames(eur)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')
        final_res[[chr]]<-select(merge(afr2, eur, by=c('CHR', 'POS')), CHR, SNP.y, V3.x, POS) %>% as.data.table
        cat(chr, '\n')
        remove(afr2, eur)
        whi<-whi[CHR!=chr]
        gc()
}

cat('checkpoint EUR\n')
do.call(rbind, final_res)-> final_res #864655 SNPs
colnames(final_res)<-c('CHR', 'SNP_ID', 'V3', 'POS')
final_res[, SNP:=paste0(CHR, ":", POS)]

fread('zcat ~/height_prediction/ldpred/output/Height.QC.gz')-> sum_stats # 748723 SNPs
merge(sum_stats, final_res, by=c('CHR', 'POS', 'SNP')) %>% as.data.table-> final
select(final, CHR, SNP, V3, POS, A2, A1) %>% as.data.table-> final2

fwrite(final2, file="~/height_prediction/ldpred/output/WHI_UKB_EUR.ldpred.bim", sep="\t", col.names=F)

system('cp ~/height_prediction/ldpred/output/WHI.ldpred.fam output/WHI_UKB_EUR.ldpred.fam')
system('cp ~/height_prediction/ldpred/output/WHI.ldpred.bed output/WHI_UKB_EUR.ldpred.bed')
##
#JHS
#WHI
remove(list=ls())
fread('~/height_prediction/ldpred/output/JHS.ldpred.bim')-> jhs
colnames(jhs)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')
jhs<-jhs[CHR %in% 1:22]
jhs$CHR<-as.integer(jhs$CHR)
#
cat('another check\n')
final_res<-vector('list',22)
for(chr in 1:22){
        #afr<-final2[CHR==chr]
        afr2<-jhs[CHR==chr]
        fread(paste0('/project/mathilab/data/UKB/imputed/ukb_imp_chr', chr, '_eur.bim'))-> eur
        colnames(eur)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')
        final_res[[chr]]<-select(merge(afr2, eur, by=c('CHR', 'POS')), CHR, SNP.y, V3.x, POS) %>% as.data.table
        cat(chr, '\n')
        remove(afr2, eur)
        jhs<-jhs[CHR!=chr]
        gc()
}
cat('checkpoint EUR\n')
do.call(rbind, final_res)-> final_res #821551 SNPs
colnames(final_res)<-c('CHR', 'SNP_ID', 'V3', 'POS')
final_res[, SNP:=paste0(CHR, ":", POS)]

fread('zcat ~/height_prediction/ldpred/output/Height.QC.gz')-> sum_stats
merge(sum_stats, final_res, by=c('CHR', 'POS', 'SNP')) %>% as.data.table-> final # 706341 SNPs
select(final, CHR, SNP, V3, POS, A2, A1) %>% as.data.table-> final2


fwrite(final2, file="~/height_prediction/ldpred/output/JHS_UKB_EUR.ldpred.bim", sep="\t", col.names=F)

system('cp ~/height_prediction/ldpred/output/JHS.ldpred.fam output/JHS_UKB_EUR.ldpred.fam')
system('cp ~/height_prediction/ldpred/output/JHS.ldpred.bed output/JHS_UKB_EUR.ldpred.bed')
