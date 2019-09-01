library(data.table)

readRDS('/home/bbita/height_prediction/gwas/WHI/output/Nr_SNPs_WHI.Rds')[, Dataset:='WHI_afr']-> WHI
readRDS('/home/bbita/height_prediction/gwas/JHS/output/Nr_SNPs_JHS.Rds')[, Dataset:='JHS_afr']-> JHS
readRDS('/home/bbita/height_prediction/gwas/ukb_afr/output/Nr_SNPs_UKB_afr.Rds')[, Dataset:='UKB_afr']-> UKB_AFR
readRDS('/project/mathilab/bbita/gwas_admix/new_height/Nr_SNPs_UKB.Rds')[, Dataset:='UKB_eur']-> UKB_EUR
readRDS('/home/bbita/height_prediction/gwas/HRS_afr/output/Nr_SNPs_HRS_afr.Rds')[, Dataset:='HRS_afr']-> HRS_afr
readRDS('/home/bbita/height_prediction/gwas/HRS_eur/output/Nr_SNPs_HRS.Rds')[, Dataset:='HRS_eur']-> HRS_eur


fwrite(WHI, '~/height_prediction/figs_for_paper/figs/WHI.txt')
fwrite(JHS, '~/height_prediction/figs_for_paper/figs/JHS.txt')
fwrite(UKB_AFR, '~/height_prediction/figs_for_paper/figs/UKB_AFR.txt')
fwrite(UKB_EUR, '~/height_prediction/figs_for_paper/figs/UKB_EUR.txt')
fwrite(HRS_afr, '~/height_prediction/figs_for_paper/figs/HRS_afr.txt')
fwrite(HRS_eur, '~/height_prediction/figs_for_paper/figs/HRS_eur.txt')


readRDS('/home/bbita/height_prediction/sib_betas/WHI/output/Nr_SNPs_WHI.Rds')[, Dataset:='WHI_afr(sib)']-> WHI_sib
readRDS('/home/bbita/height_prediction/sib_betas/JHS/output/Nr_SNPs_JHS.Rds')[, Dataset:='JHS_afr(sib)']-> JHS_sib
readRDS('/home/bbita/height_prediction/sib_betas/ukb_afr/output/Nr_SNPs_UKB_afr.Rds')[, Dataset:='UKB_afr(sib)']-> UKB_AFR_sib
readRDS('/home/bbita/height_prediction/sib_betas/HRS_afr/output/Nr_SNPs_HRS_afr.Rds')[, Dataset:='HRS_afr(sib)']-> HRS_afr_sib
readRDS('/home/bbita/height_prediction/sib_betas/HRS_eur/output/Nr_SNPs_HRS.Rds')[, Dataset:='HRS_eur(sib)']-> HRS_eur_sib
#TO DO UKB EUR

fwrite(WHI_sib, '~/height_prediction/figs_for_paper/figs/WHI_sib.txt')
fwrite(JHS_sib, '~/height_prediction/figs_for_paper/figs/JHS_sib.txt')
fwrite(UKB_AFR_sib, '~/height_prediction/figs_for_paper/figs/UKB_AFR_sib.txt')
fwrite(HRS_afr_sib, '~/height_prediction/figs_for_paper/figs/HRS_afr_sib.txt')
fwrite(HRS_eur_sib, '~/height_prediction/figs_for_paper/figs/HRS_eur_sib.txt')


