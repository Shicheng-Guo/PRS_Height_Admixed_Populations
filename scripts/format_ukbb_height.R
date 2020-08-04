library(data.table)
library(biomaRt)

fread('zcat ../input/50_raw_filtered.txt.gz')-> ukb_height #read in GWAS summary statistics for height from the UK Biobank
ukb_height[,c("CHR", "POS","Allele2","Allele1") := tstrsplit(variant, ":", fixed=TRUE)][,variant:=NULL]  #fix columns. In the UKB, variants are listed as CHR:POS: REF:ALT, where ALT is the effect allele. So, in order to be compatilble with my scripts (where allele 1 is effect all
ukb_height[, N:=n_complete_samples][, AC:=NULL][, BETA:=beta][,P:=pval]
ukb_height[,n_complete_samples:=NULL][, beta:=NULL][, pval:=NULL][, SE:=se]
ukb_height[,.(Allele1,Allele2, BETA, SE, P, N, CHR, POS)]-> ukb_height
ukb_height[CHR %in% 1:22]-> ukb_height
na.omit(ukb_height)-> ukb_height
ukb_height[,CHR:=as.integer(CHR)]
ukb_height[,POS:=as.integer(POS)]
setkey(ukb_height, CHR, POS)
gc()

JHS<-fread('../input/JHS/JHS_b37_strand.bim')[,c(1,2,4)]
colnames(JHS)<-c('CHR', 'SNP', 'POS')
setkey(JHS, CHR, POS)
ukb_height[JHS,nomatch=0][,.(SNP,Allele1,Allele2, BETA, SE, P, N, CHR, POS)]-> JHS1
fwrite(JHS1,'../gwas/JHS/plink2/output/betas_for_plink.txt', sep = "\t")
remove(JHS, JHS1)
#WHI
WHI<-fread('../input/WHI/WHI_b37_strand_include.bim')[,c(1,2,4)]
colnames(WHI)<-c('CHR', 'SNP', 'POS')
setkey(WHI, CHR, POS)
ukb_height[WHI,nomatch=0][,.(SNP,Allele1,Allele2, BETA, SE, P, N, CHR, POS)]-> WHI1
fwrite(WHI1,'../gwas/WHI/plink2/output/betas_for_plink.txt', sep = "\t")
remove(WHI, WHI1)
#HRS_afr
HRS_afr<-fread('../input/HRS_afr/HRS_AFR_b37_strand_include.bim')[,c(1,2,4)]
colnames(HRS_afr)<-c('CHR', 'SNP', 'POS')
setkey(HRS_afr, CHR, POS)
ukb_height[HRS_afr,nomatch=0][,.(SNP,Allele1,Allele2, BETA, SE, P, N, CHR, POS)]-> HRS_afr1
fwrite(HRS_afr1,'../gwas/HRS_afr/plink2/output/betas_for_plink.txt', sep = "\t")
remove(HRS_afr, HRS_afr1)
#
HRS_eur<-fread('../input/HRS_eur/HRS_EUR_b37_strand_include.bim')[,c(1,2,4)]
colnames(HRS_eur)<-c('CHR', 'SNP', 'POS')
setkey(HRS_eur, CHR, POS)
ukb_height[HRS_eur,nomatch=0][,.(SNP,Allele1,Allele2, BETA, SE, P, N, CHR, POS)]-> HRS_eur1
fwrite(HRS_eur1,'../gwas/HRS_eur/plink2/output/betas_for_plink.txt', sep = "\t")
remove(HRS_eur, HRS_eur1)
#
UKB_eur<-fread('../input/ukb_eur/UKB_EUR.bim')[,c(1,2,4)]
colnames(UKB_eur)<-c('CHR', 'SNP', 'POS')
setkey(UKB_eur, CHR, POS)
ukb_height[UKB_eur,nomatch=0][,.(SNP,Allele1,Allele2, BETA, SE, P, N, CHR, POS)]-> UKB_eur1
fwrite(UKB_eur1,'../gwas/ukb_eur/plink2/output/betas_for_plink.txt', sep = "\t")
remove(UKB_eur, UKB_eur1)
#
UKB_afr<-fread('../input/ukb_afr/UKB_AFR.bim')[,c(1,2,4)]
colnames(UKB_afr)<-c('CHR', 'SNP', 'POS')
setkey(UKB_afr, CHR, POS)
ukb_height[UKB_afr,nomatch=0][,.(SNP,Allele1,Allele2, BETA, SE, P, N, CHR, POS)]-> UKB_afr1
fwrite(UKB_afr1,'../gwas/ukb_afr/plink2/output/betas_for_plink.txt', sep = "\t")
remove(UKB_afr, UKB_afr1)

