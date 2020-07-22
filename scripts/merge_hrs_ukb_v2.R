#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
fread('~/height_prediction/ldpred/output/UKB_EUR.ldpred.bim')-> ukb_eur
colnames(ukb_eur)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')

#fwrite(final2, file="output/HRS_AFR.ldpred.bim", sep="\t", col.names=F)
ukb_eur[,SNP:=paste0(CHR,"_", POS)]

fwrite(ukb_eur, file="~/height_prediction/ldpred/output/UKB_EUR.ldpred.bim", sep="\t", col.names=F)

######################
######################
fread('~/height_prediction/ldpred/output/UKB_AFR.ldpred.bim')-> ukb_afr
colnames(ukb_afr)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')

ukb_afr[, SNP:=paste0(CHR, "_", POS)]

fwrite(ukb_afr, file="~/height_prediction/ldpred/output/UKB_AFR.ldpred.bim", sep="\t", col.names=F)
#######################
#######################
fread('~/height_prediction/ldpred/output/HRS_AFR.ldpred.bim')-> hrs_afr
colnames(hrs_afr)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')

hrs_afr[, SNP:=paste0(CHR, "_", POS)]

fwrite(hrs_afr, file="~/height_prediction/ldpred/output/HRS_AFR.ldpred.bim", sep="\t", col.names=F)

##############################################
##############################################
fread('~/height_prediction/ldpred/output/HRS_EUR.ldpred.bim')-> hrs_eur
colnames(hrs_eur)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')
#
hrs_eur[, SNP:=paste0(CHR, "_", POS)]

fwrite(hrs_eur, file="~/height_prediction/ldpred/output/HRS_EUR.ldpred.bim", sep="\t", col.names=F)

#
#WHI
fread('~/height_prediction/ldpred/output/WHI.ldpred.bim')-> whi
colnames(whi)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')
#whi<-whi[CHR %in% 1:22]
#whi$CHR<-as.integer(whi$CHR)
#
whi[, SNP:=paste0(CHR, "_", POS)]

fwrite(whi, file="~/height_prediction/ldpred/output/WHI.ldpred.bim", sep="\t", col.names=F)

##
#JHS
fread('~/height_prediction/ldpred/output/JHS.ldpred.bim')-> jhs
colnames(jhs)<-c('CHR', 'SNP', 'V3', 'POS', 'REF', 'ALT')
#jhs<-jhs[CHR %in% 1:22]
#jhs$CHR<-as.integer(jhs$CHR)
#
jhs[, SNP:=paste0(CHR, "_", POS)]

fwrite(jhs, file="~/height_prediction/ldpred/output/JHS.ldpred.bim", sep="\t", col.names=F)
#####


