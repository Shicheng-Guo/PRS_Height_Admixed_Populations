#!/usr/bin/env Rscript
############################
library(data.table)
library(dplyr)
library(xlsx)
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
write.xlsx(x = WHI, file = "/home/bbita/height_prediction/figs_for_paper/tables/Table_S1.xlsx",
        sheetName = "WHI_afr", row.names = FALSE)
write.xlsx(x = JHS, file = "/home/bbita/height_prediction/figs_for_paper/tables/Table_S1.xlsx",
        sheetName = "JHS_afr", row.names = FALSE, append=T)
write.xlsx(x =UKB_AFR, file = "/home/bbita/height_prediction/figs_for_paper/tables/Table_S1.xlsx",
        sheetName = "UKB_afr", row.names = FALSE, append=T)
write.xlsx(x = UKB_EUR, file = "/home/bbita/height_prediction/figs_for_paper/tables/Table_S1.xlsx",
        sheetName = "UKB_eur", row.names = FALSE, append=T)
write.xlsx(x = HRS_afr, file = "/home/bbita/height_prediction/figs_for_paper/tables/Table_S1.xlsx",
        sheetName = "HRS_afr", row.names = FALSE, append=T)
write.xlsx(x = HRS_eur, file = "/home/bbita/height_prediction/figs_for_paper/tables/Table_S1.xlsx",
        sheetName = "HRS_eur", row.names = FALSE, append=T)

#Table S2
write.xlsx(x = WHI_sib, file = "/home/bbita/height_prediction/figs_for_paper/tables/Table_S2.xlsx",
        sheetName = "WHI_afr", row.names = FALSE)
write.xlsx(x = JHS_sib, file = "/home/bbita/height_prediction/figs_for_paper/tables/Table_S2.xlsx",
        sheetName = "JHS_afr", row.names = FALSE, append=T)
write.xlsx(x =UKB_AFR_sib, file = "/home/bbita/height_prediction/figs_for_paper/tables/Table_S2.xlsx",
        sheetName = "UKB_afr", row.names = FALSE, append=T)
write.xlsx(x = UKB_EUR_sib, file = "/home/bbita/height_prediction/figs_for_paper/tables/Table_S2.xlsx",
        sheetName = "UKB_eur", row.names = FALSE, append=T)
write.xlsx(x = HRS_afr_sib, file = "/home/bbita/height_prediction/figs_for_paper/tables/Table_S2.xlsx",
        sheetName = "HRS_afr", row.names = FALSE, append=T)
write.xlsx(x = HRS_eur_sib, file = "/home/bbita/height_prediction/figs_for_paper/tables/Table_S2.xlsx",
        sheetName = "HRS_eur", row.names = FALSE, append=T)

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
write.xlsx(x = imputed, file = "/home/bbita/height_prediction/figs_for_paper/tables/Table_S3.xlsx",
        sheetName = "HRS_and_UKB", row.names = FALSE)


#fwrite(data.table(Set=names_imp, UKB=unlist(imp[['UKB']]), HRS=unlist(imp[['HRS']])), "~/height_prediction/figs_for_paper/figs/imputed.txt")
#fwrite(readRDS('~/height_prediction/imputed/output/comparison.Rds'), "~/height_prediction/figs_for_paper/figs/imp_comparison.txt")

##unweighted PRS:




###
select(readRDS('/home/bbita/height_prediction/unweighted_prs/output/Nr_SNPs_WHI.Rds')[, Dataset:='WHI_afr'],Name, Nr, Part_R2)-> WHI_unw
select(readRDS('/home/bbita/height_prediction/unweighted_prs/output/Nr_SNPs_JHS.Rds')[, Dataset:='JHS_afr'],Name, Nr, Part_R2)-> JHS_unw
select(readRDS('/home/bbita/height_prediction/unweighted_prs/output/Nr_SNPs_UKB_afr.Rds')[, Dataset:='UKB_afr'],Name, Nr, Part_R2)-> UKB_AFR_unw
select(readRDS('/home/bbita/height_prediction/unweighted_prs/output/Nr_SNPs_UKB_eur.Rds')[, Dataset:='UKB_eur'],Name, Nr, Part_R2)-> UKB_EUR_unw
select(readRDS('/home/bbita/height_prediction/unweighted_prs/output/Nr_SNPs_HRS_afr.Rds')[, Dataset:='HRS_afr'],Name, Nr, Part_R2)-> HRS_afr_unw
select(readRDS('/home/bbita/height_prediction/unweighted_prs/output/Nr_SNPs_HRS.Rds')[, Dataset:='HRS_eur'], Name, Nr, Part_R2)-> HRS_eur_unw


write.xlsx(x = WHI_unw, file = "/home/bbita/height_prediction/figs_for_paper/tables/Table_S4.xlsx",
        sheetName = "WHI_afr", row.names = FALSE)
write.xlsx(x = JHS_unw, file = "/home/bbita/height_prediction/figs_for_paper/tables/Table_S4.xlsx",
        sheetName = "JHS_afr", row.names = FALSE, append=T)
write.xlsx(x = UKB_AFR_unw, file = "/home/bbita/height_prediction/figs_for_paper/tables/Table_S4.xlsx",
        sheetName = "UKB_afr", row.names = FALSE,append=T)
write.xlsx(x = HRS_afr_unw, file = "/home/bbita/height_prediction/figs_for_paper/tables/Table_S4.xlsx",
        sheetName = "HRS_afr", row.names = FALSE,append=T)
write.xlsx(x = HRS_eur_unw, file = "/home/bbita/height_prediction/figs_for_paper/tables/Table_S4.xlsx",
        sheetName = "HRS_eur", row.names = FALSE,append=T)
write.xlsx(x = UKB_EUR_unw, file = "/home/bbita/height_prediction/figs_for_paper/tables/Table_S4.xlsx",
        sheetName = "UKB_eur", row.names = FALSE,append=T)



