#!/usr/bin/env Rscript

old <- Sys.time()
#args<-c('phys_100000_0.0005', 21, 0.10)
## Load libraries
library(optparse)
library(data.table)
library(dplyr)
library(readr)
library(tidyr)
library(parallel)
options(scipen=999)
source('~/height_prediction/scripts/mclapply2.R')


cat('checkpoint number 1\n')
plink<-fread('~/height_prediction/runSmartpCA-master/UKB_AFR_imputed/association_v3.Res.Height.glm.linear.adjusted', header=T, fill=T)
colnames(plink)[2]<-'MarkerName'
colnames(plink)[1]<-'CHR'

plink2<-fread('~/height_prediction/runSmartpCA-master/UKB_AFR_imputed/test3.txt', fill=T)
colnames(plink2)<-c("CHR","POS", "MarkerName","REF","ALT","Effect_Allele_plink","TEST","OBS_CT","PLINK", "SE_plink","T_STAT", "UNADJ") #plink is BETA
plink2<-select(plink2, -c("REF", "ALT"))
gc()
select(merge(plink, plink2, by=c("CHR", "MarkerName")), CHR, MarkerName, POS, Effect_Allele_plink, PLINK, SE_plink)-> final_plink #3388119
remove(plink, plink2)
gc()
final_plink$POS<-as.numeric(final_plink$POS)
gc()
final_plink$CHR<-as.numeric(final_plink$CHR)
gc()
final_plink<-na.omit(final_plink) #3382050
arrange(final_plink, CHR, POS) %>% as.data.table-> final_plink
gc()


saveRDS(final_plink, file='~/height_prediction/loc_anc_analysis/output/final_plink.Rds')
final_plink[, PLINK:=scale(PLINK)]
saveRDS(final_plink, file='~/height_prediction/loc_anc_analysis/output/final_plink_v2.Rds')
new <- Sys.time() - old # calculate difference
print(new) # print in nice format
