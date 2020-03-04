library(data.table)
library(dplyr)
library(ggplot2)
library(xlsx)
options(scipen=999)

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

write.xlsx(x = A, file = "/home/bbita/height_prediction/figs_for_paper/tables/Table_S5.xlsx",
        sheetName = "Rec_Quantiles_AA_Map", row.names = FALSE)


hrs_afr<-levels(readRDS('~/height_prediction/strat_prs/output/rec_quant_HRS_afr_gwas_20000_AA_v2.Rds')$Quantile)

B<-rbind(unlist(readRDS('~/height_prediction/strat_prs/output/part_R2_JHS_gwas_20000_AA_v2.Rds')),
unlist(readRDS('~/height_prediction/strat_prs/output/part_R2_WHI_gwas_20000_AA_v2.Rds')),
unlist(readRDS('~/height_prediction/strat_prs/output/part_R2_ukb_afr_gwas_20000_AA_v2.Rds')),
unlist(readRDS('~/height_prediction/strat_prs/output/part_R2_HRS_afr_gwas_20000_AA_v2.Rds')),
unlist(readRDS('~/height_prediction/strat_prs/output/part_R2_HRS_eur_gwas_20000_AA_v2.Rds')))
B<-as.data.frame(B)
cbind(data.table(Dataset=c('JHS_afr', 'WHI_afr', 'UKB_afr', 'HRS_afr', 'HRS_eur')), B)-> B
write.xlsx(x = B, file = "/home/bbita/height_prediction/figs_for_paper/tables/Table_S5.xlsx",
        sheetName = "Partial_R2_quantile", row.names = FALSE, append=T)

