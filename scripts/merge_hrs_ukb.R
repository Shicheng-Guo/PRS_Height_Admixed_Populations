#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
fread('~/height_prediction/ldpred/output/1000g_all.bim')-> g1000
colnames(g1000)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')
gc()
fread('zcat output/Height.QC.gz')-> sum_stats
sum_stats[, SNP:=paste0(CHR, "_", POS)]
gc()
fread('~/height_prediction/ldpred/output/UKB_EUR.ldpred.bim')-> ukb_eur
colnames(ukb_eur)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')

setkey(g1000, CHR, POS)
setkey(ukb_eur, CHR, POS)
g1000[ukb_eur, nomatch=0]-> out
out[,.(CHR, SNP, V3, POS, REF, ALT)]-> out
merge(sum_stats, out, by=c('CHR', 'POS', 'SNP')) %>% as.data.table-> final #638811
select(final, CHR, SNP, V3, POS, A2, A1) %>% as.data.table-> final2

cat('checkpoint\n')
fwrite(final2, file="~/height_prediction/ldpred/output/UKB_EUR_g1000.ldpred.bim", sep="\t", col.names=F)

system('cp ~/height_prediction/ldpred/output/UKB_EUR.ldpred.fam output/UKB_EUR_g1000.ldpred.fam')
system('cp ~/height_prediction/ldpred/output/UKB_EUR.ldpred.bed  output/UKB_EUR_g1000.ldpred.bed')
######################
######################
remove(ukb_eur, final, out, final2)
fread('~/height_prediction/ldpred/output/UKB_AFR.ldpred.bim')-> ukb_afr
colnames(ukb_afr)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')

setkey(g1000, CHR, POS)
setkey(ukb_afr, CHR, POS)
g1000[ukb_afr, nomatch=0]-> out
out[,.(CHR, SNP, V3, POS, REF, ALT)]-> out
merge(sum_stats, out, by=c('CHR', 'POS', 'SNP')) %>% as.data.table-> final #638811
select(final, CHR, SNP, V3, POS, A2, A1) %>% as.data.table-> final2

cat('checkpoint\n')
fwrite(final2, file="~/height_prediction/ldpred/output/UKB_AFR_g1000.ldpred.bim", sep="\t", col.names=F)

system('cp ~/height_prediction/ldpred/output/UKB_AFR.ldpred.fam output/UKB_AFR_g1000.ldpred.fam')
system('cp ~/height_prediction/ldpred/output/UKB_AFR.ldpred.bed output/UKB_AFR_g1000.ldpred.bed')
#######################
#######################
remove(ukb_afr, final, final2, out)
fread('~/height_prediction/ldpred/output/HRS_AFR.ldpred.bim')-> hrs_afr
colnames(hrs_afr)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')

setkey(g1000, CHR, POS)
setkey(hrs_afr, CHR, POS)
g1000[hrs_afr, nomatch=0]-> out
out[,.(CHR, SNP, V3, POS, REF, ALT)]-> out
merge(sum_stats, out, by=c('CHR', 'POS', 'SNP')) %>% as.data.table-> final #1421790
select(final, CHR, SNP, V3, POS, A2, A1) %>% as.data.table-> final2

fwrite(final2, file="~/height_prediction/ldpred/output/HRS_AFR_g1000.ldpred.bim", sep="\t", col.names=F)

system('cp ~/height_prediction/ldpred/output/HRS_AFR.ldpred.fam output/HRS_AFR_g1000.ldpred.fam')
system('cp ~/height_prediction/ldpred/output/HRS_AFR.ldpred.bed output/HRS_AFR_g1000.ldpred.bed')
##############################################
##############################################
remove(out, final, final2, hrs_afr)
fread('~/height_prediction/ldpred/output/HRS_EUR.ldpred.bim')-> hrs_eur
colnames(hrs_eur)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')
#
setkey(g1000, CHR, POS)
setkey(hrs_eur, CHR, POS)
g1000[hrs_eur, nomatch=0]-> out
out[,.(CHR, SNP, V3, POS, REF, ALT)]-> out
merge(sum_stats, out, by=c('CHR', 'POS', 'SNP')) %>% as.data.table-> final #1421790
select(final, CHR, SNP, V3, POS, A2, A1) %>% as.data.table-> final2


fwrite(final2, file="~/height_prediction/ldpred/output/HRS_EUR_g1000.ldpred.bim", sep="\t", col.names=F)

system('cp ~/height_prediction/ldpred/output/HRS_EUR.ldpred.fam output/HRS_EUR_g1000.ldpred.fam')
system('cp ~/height_prediction/ldpred/output/HRS_EUR.ldpred.bed output/HRS_EUR_g1000.ldpred.bed')
#
#WHI
remove(out, final, final2, hrs_eur)
fread('~/height_prediction/ldpred/output/WHI.ldpred.bim')-> whi
colnames(whi)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')
whi<-whi[CHR %in% 1:22]
whi$CHR<-as.integer(whi$CHR)
#

setkey(g1000, CHR, POS)
setkey(whi, CHR, POS)
g1000[whi, nomatch=0]-> out
out[,.(CHR, SNP, V3, POS, REF, ALT)]-> out
merge(sum_stats, out, by=c('CHR', 'POS', 'SNP')) %>% as.data.table-> final #708356
select(final, CHR, SNP, V3, POS, A2, A1) %>% as.data.table-> final2

fwrite(final2, file="~/height_prediction/ldpred/output/WHI_g1000.ldpred.bim", sep="\t", col.names=F)

system('cp ~/height_prediction/ldpred/output/WHI.ldpred.fam output/WHI_g1000.ldpred.fam')
system('cp ~/height_prediction/ldpred/output/WHI.ldpred.bed output/WHI_g1000.ldpred.bed')
##
#JHS
#WHI
remove(out, final, final2, whi)
fread('~/height_prediction/ldpred/output/JHS.ldpred.bim')-> jhs
colnames(jhs)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')
jhs<-jhs[CHR %in% 1:22]
jhs$CHR<-as.integer(jhs$CHR)
#
setkey(g1000, CHR, POS)
setkey(jhs, CHR, POS)
g1000[jhs, nomatch=0]-> out
out[,.(CHR, SNP, V3, POS, REF, ALT)]-> out
merge(sum_stats, out, by=c('CHR', 'POS', 'SNP')) %>% as.data.table-> final #668260
select(final, CHR, SNP, V3, POS, A2, A1) %>% as.data.table-> final2
fwrite(final2, file="~/height_prediction/ldpred/output/JHS_g1000.ldpred.bim", sep="\t", col.names=F)

system('cp ~/height_prediction/ldpred/output/JHS.ldpred.fam output/JHS_g1000.ldpred.fam')
system('cp ~/height_prediction/ldpred/output/JHS.ldpred.bed output/JHS_g1000.ldpred.bed')
