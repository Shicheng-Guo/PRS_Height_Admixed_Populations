#!/usr/bin/env Rscript
#######################################################################################################################################################
#Script Name    :make_all_tables.R
#Description    :Generate all supplementary tables present in the paper.
#Args           :NA
#Author         :Bárbara D Bitarello
#Email          :barbarabitarello@gmail.com
#Date created   :22-07-2020
#Requirements   :all analyses have been performed and output saved
#USage          :Rscript --vanills scripts/make_all_tables.R
########################################################################################################################################################
#save start time
ptm <- proc.time()
#load packages
suppressPackageStartupMessages({
        library(data.table)
        library(dplyr)
        library(optparse)
	library(xlsx)
})
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
}
############################
options(scipen=999)

readRDS('/home/bbita/height_prediction/gwas/WHI/output/Nr_SNPs_WHI.Rds')[, Dataset:='WHI_afr']-> WHI
readRDS('/home/bbita/height_prediction/gwas/JHS/output/Nr_SNPs_JHS.Rds')[, Dataset:='JHS_afr']-> JHS
readRDS('/home/bbita/height_prediction/gwas/ukb_afr/output/Nr_SNPs_UKB_afr.Rds')[, Dataset:='UKB_afr']-> UKB_AFR
readRDS('/project/mathilab/bbita/gwas_admix/new_height/Nr_SNPs_UKB.Rds')[, Dataset:='UKB_eur']-> UKB_EUR
readRDS('/home/bbita/height_prediction/gwas/HRS_afr/output/Nr_SNPs_HRS_afr.Rds')[, Dataset:='HRS_afr']-> HRS_afr
readRDS('/home/bbita/height_prediction/gwas/HRS_eur/output/Nr_SNPs_HRS.Rds')[, Dataset:='HRS_eur']-> HRS_eur

readRDS('/home/bbita/height_prediction/sib_betas/WHI/output/Nr_SNPs_WHI.Rds')[, Dataset:='WHI_afr(sib)']-> WHI_sib
readRDS('/home/bbita/height_prediction/sib_betas/JHS/output/Nr_SNPs_JHS.Rds')[, Dataset:='JHS_afr(sib)']-> JHS_sib
readRDS('/home/bbita/height_prediction/sib_betas/ukb_afr/output/Nr_SNPs_UKB_afr.Rds')[, Dataset:='UKB_afr(sib)']-> UKB_AFR_sib
readRDS('/home/bbita/height_prediction/sib_betas/ukb_eur/output/Nr_SNPs_UKB_eur.Rds')[, Dataset:='UKB_eur(sib)']-> UKB_EUR_sib
readRDS('/home/bbita/height_prediction/sib_betas/HRS_afr/output/Nr_SNPs_HRS_afr.Rds')[, Dataset:='HRS_afr(sib)']-> HRS_afr_sib
readRDS('/home/bbita/height_prediction/sib_betas/HRS_eur/output/Nr_SNPs_HRS.Rds')[, Dataset:='HRS_eur(sib)']-> HRS_eur_sib

#Table S1
write.xlsx(x = WHI, file = "/home/bbita/height_prediction/tables/Table_S1.xlsx",
        sheetName = "1", row.names = FALSE)
write.xlsx(x = JHS, file = "/home/bbita/height_prediction/tables/Table_S1.xlsx",
        sheetName = "2", row.names = FALSE, append=T)
write.xlsx(x =UKB_AFR, file = "/home/bbita/height_prediction/tables/Table_S1.xlsx",
        sheetName = "3", row.names = FALSE, append=T)
write.xlsx(x = UKB_EUR, file = "/home/bbita/height_prediction/tables/Table_S1.xlsx",
        sheetName = "4", row.names = FALSE, append=T)
write.xlsx(x = HRS_afr, file = "/home/bbita/height_prediction/tables/Table_S1.xlsx",
        sheetName = "5", row.names = FALSE, append=T)
write.xlsx(x = HRS_eur, file = "/home/bbita/height_prediction/tables/Table_S1.xlsx",
        sheetName = "6", row.names = FALSE, append=T)

#Table S2
write.xlsx(x = WHI_sib, file = "/home/bbita/height_prediction/tables/Table_S1.xlsx",
        sheetName = "7", row.names = FALSE, append=T)
write.xlsx(x = JHS_sib, file = "/home/bbita/height_prediction/tables/Table_S1.xlsx",
        sheetName = "8", row.names = FALSE, append=T)
write.xlsx(x =UKB_AFR_sib, file = "/home/bbita/height_prediction/tables/Table_S1.xlsx",
        sheetName = "9", row.names = FALSE, append=T)
write.xlsx(x = UKB_EUR_sib, file = "/home/bbita/height_prediction/tables/Table_S1.xlsx",
        sheetName = "10", row.names = FALSE, append=T)
write.xlsx(x = HRS_afr_sib, file = "/home/bbita/height_prediction/tables/Table_S1.xlsx",
        sheetName = "11", row.names = FALSE, append=T)
write.xlsx(x = HRS_eur_sib, file = "/home/bbita/height_prediction/tables/Table_S1.xlsx",
        sheetName = "12", row.names = FALSE, append=T)

###Imputed

names_imp<-paste(rep(c(5000,10000,25000,50000,75000,100000, 500000), 5), rep(c(0.0005, 0.00005, 0.000005, 0.0000005, 0.00000005), each=7), sep="_")
imp<-vector('list', 2)
names(imp)<-c('UKB', 'HRS')
imp[['UKB']]<-vector('list', 35)
names(imp[['UKB']])<-names_imp
imp[['HRS']]<-vector('list', 35)
names(imp[['HRS']])<-names_imp
for(I in names_imp){
	imp[['UKB']][[I]]<-nrow(do.call(rbind, readRDS(paste0('/home/bbita/height_prediction/imputed/output/UKB_vec_all_', I,'.Rds'))))
	imp[['HRS']][[I]]<-nrow(do.call(rbind, readRDS(paste0('/home/bbita/height_prediction/imputed/output/vec_all_', I,'.Rds'))))
}

imputed<-data.table(Set=paste0("phys_",names_imp), UKB=unlist(imp[['UKB']]), HRS=unlist(imp[['HRS']]))
imputed[, Name:=Set][,Nr_SNPs_UKB:=UKB][,Nr_SNPs_HRS:=HRS][, UKB:=NULL][,HRS:=NULL][,Set:=NULL]
merge(imputed,select(readRDS('~/height_prediction/imputed/output/comparison.Rds'), HRS_afr_imp, HRS_eur_imp, Name), by="Name")-> imputed
imputed<-merge(imputed, select(readRDS('~/height_prediction/imputed/output/comparison_ukb.Rds'), UKB_eur, UKB_afr, Name))
write.xlsx(x = imputed, file = "/home/bbita/height_prediction/tables/Table_S1.xlsx",
        sheetName = "13", row.names = FALSE, append=T)

###
select(readRDS('/home/bbita/height_prediction/unweighted_prs/output/Nr_SNPs_WHI.Rds')[, Dataset:='WHI_afr'],Name, Nr, Part_R2, Dataset)-> WHI_unw
select(readRDS('/home/bbita/height_prediction/unweighted_prs/output/Nr_SNPs_JHS.Rds')[, Dataset:='JHS_afr'],Name, Nr, Part_R2, Dataset)-> JHS_unw
select(readRDS('/home/bbita/height_prediction/unweighted_prs/output/Nr_SNPs_UKB_afr.Rds')[, Dataset:='UKB_afr'],Name, Nr, Part_R2, Dataset)-> UKB_AFR_unw
select(readRDS('/home/bbita/height_prediction/unweighted_prs/output/Nr_SNPs_UKB_eur.Rds')[, Dataset:='UKB_eur'],Name, Nr, Part_R2, Dataset)-> UKB_EUR_unw
select(readRDS('/home/bbita/height_prediction/unweighted_prs/output/Nr_SNPs_HRS_afr.Rds')[, Dataset:='HRS_afr'],Name, Nr, Part_R2, Dataset)-> HRS_afr_unw
select(readRDS('/home/bbita/height_prediction/unweighted_prs/output/Nr_SNPs_HRS.Rds')[, Dataset:='HRS_eur'], Name, Nr, Part_R2, Dataset)-> HRS_eur_unw


write.xlsx(x = WHI_unw, file = "/home/bbita/height_prediction/tables/Table_S1.xlsx",
        sheetName = "14", row.names = FALSE, append=T)
write.xlsx(x = JHS_unw, file = "/home/bbita/height_prediction/tables/Table_S1.xlsx",
        sheetName = "15", row.names = FALSE, append=T)
write.xlsx(x = UKB_AFR_unw, file = "/home/bbita/height_prediction/tables/Table_S1.xlsx",
        sheetName = "16", row.names = FALSE,append=T)
write.xlsx(x = UKB_EUR_unw, file = "/home/bbita/height_prediction/tables/able_S1.xlsx",
        sheetName = "17", row.names = FALSE,append=T)
write.xlsx(x = HRS_afr_unw, file = "/home/bbita/height_prediction/tables/Table_S1.xlsx",
        sheetName = "18", row.names = FALSE,append=T)
write.xlsx(x = HRS_eur_unw, file = "/home/bbita/height_prediction/tables/Table_S1.xlsx",
        sheetName = "19", row.names = FALSE,append=T)


hrs_afr<-levels(readRDS('~/height_prediction/strat_prs/output/rec_quant_HRS_afr_gwas_20000_AA_v2.Rds')$Quantile)
hrs_eur<-levels(readRDS('~/height_prediction/strat_prs/output/rec_quant_HRS_eur_gwas_20000_AA_v2.Rds')$Quantile)
jhs<-levels(readRDS('~/height_prediction/strat_prs/output/rec_quant_JHS_gwas_20000_AA_v2.Rds')$Quantile)
whi<-levels(readRDS('~/height_prediction/strat_prs/output/rec_quant_WHI_gwas_20000_AA_v2.Rds')$Quantile)
ukb_afr<-levels(readRDS('~/height_prediction/strat_prs/output/rec_quant_ukb_afr_gwas_20000_AA_v2.Rds')$Quantile)

A<-data.table(Dataset=c('JHS_afr', 'WHI_afr', 'UKB_afr', 'HRS_afr', 'HRS_eur'), 
Q1=c(jhs[1],whi[1], ukb_afr[1],hrs_afr[1],hrs_eur[1]),
Q2=c(jhs[2],whi[2], ukb_afr[2],hrs_afr[2],hrs_eur[2]),
Q3=c(jhs[3],whi[3], ukb_afr[3],hrs_afr[3],hrs_eur[3]),
Q4=c(jhs[4],whi[4], ukb_afr[4],hrs_afr[4],hrs_eur[4]))

write.xlsx(x = A, file = "/home/bbita/height_prediction/tables/Table_S5.xlsx",
        sheetName = "Rec_Quantiles_AA_Map", row.names = FALSE)


hrs_afr<-levels(readRDS('~/height_prediction/strat_prs/output/rec_quant_HRS_afr_gwas_20000_AA_v2.Rds')$Quantile)

B<-rbind(unlist(readRDS('~/height_prediction/strat_prs/output/part_R2_JHS_gwas_20000_AA_v2.Rds')),
unlist(readRDS('~/height_prediction/strat_prs/output/part_R2_WHI_gwas_20000_AA_v2.Rds')),
unlist(readRDS('~/height_prediction/strat_prs/output/part_R2_ukb_afr_gwas_20000_AA_v2.Rds')),
unlist(readRDS('~/height_prediction/strat_prs/output/part_R2_HRS_afr_gwas_20000_AA_v2.Rds')),
unlist(readRDS('~/height_prediction/strat_prs/output/part_R2_HRS_eur_gwas_20000_AA_v2.Rds')))
B<-as.data.frame(B)
cbind(data.table(Dataset=c('JHS_afr', 'WHI_afr', 'UKB_afr', 'HRS_afr', 'HRS_eur')), B)-> B
write.xlsx(x = B, file = "/home/bbita/height_prediction/tables/Table_S5.xlsx",
        sheetName = "Partial_R2_quantile", row.names = FALSE, append=T)

#
tim<-proc.time()[[3]] - ptm[[3]]
cat('Total elapsed time is ', tim, ' seconds\n')
#End

