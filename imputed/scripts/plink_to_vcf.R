#!/usr/bin/env Rscript
############################
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
}
options(scipen=10)
options(digits=10)
W<-args[1]
P<-args[2]

if(args[3]=='HRS'){
	for (chr in 22:1){
		system(paste0('plink2 --bfile /project/mathilab/data/HRS/data/imputed/HRS_AFR_imputed_chr', chr, ' --extract ~/height_prediction/imputed/output/chr_', chr,'_', W, '_', P, '.txt --recode vcf --out ~/height_prediction/imputed/output/HRS_afr_', W, '_', P, '_chr_', chr))
	#	system(paste0("grep -v '##' ~/height_prediction/imputed/output/HRS_afr_",W, '_', P, '_chr_', chr, '.vcf > ~/height_prediction/imputed/output/temp_', W, '_', P, '_chr_', chr, '.txt '))
		system(paste0("mv  ~/height_prediction/imputed/output/temp_", W, '_', P, '_chr_', chr, '.txt  ~/height_prediction/imputed/output/HRS_afr_',W, '_', P, '_chr_', chr, '.vcf'))
		system(paste0('plink2 --bfile /project/mathilab/data/HRS/data/imputed/HRS_EUR_imputed_chr', chr, ' --extract ~/height_prediction/imputed/output/chr_', chr,'_', W, '_', P, '.txt --recode vcf --out ~/height_prediction/imputed/output/HRS_eur_', W, '_', P, '_chr_', chr))
	#	system(paste0("grep -v '##' ~/height_prediction/imputed/output/HRS_eur_",W, '_', P, '_chr_', chr, '.vcf > ~/height_prediction/imputed/output/temp_', W, '_', P, '_chr_', chr, '.txt '))
       		system(paste0("mv  ~/height_prediction/imputed/output/temp_", W, '_', P, '_chr_', chr, '.txt  ~/height_prediction/imputed/output/HRS_eur_',W, '_', P, '_chr_', chr, '.vcf'))
        	system(paste0('bgzip ~/height_prediction/imputed/output/HRS_afr_', W, '_', P, '_chr_', chr, '.vcf'))
		system(paste0('bgzip ~/height_prediction/imputed/output/HRS_eur_', W, '_', P, '_chr_', chr, '.vcf'))
       		cat(chr, 'r')
	}
}

#UKB
if(args[3]=='UKB'){
	for (chr in 22:1){
        	system(paste0('plink2 --bfile /project/mathilab/data/UKB/imputed/ukb_imp_chr', chr, '_afr --extract ~/height_prediction/imputed/output/UKB_chr_', chr,'_', W, '_', P, '.txt --recode vcf --out ~/height_prediction/imputed/output/UKB_afr_', W, '_', P, '_chr_', chr))
        #	system(paste0("grep -v '##' ~/height_prediction/imputed/output/UKB_afr_",W, '_', P, '_chr_', chr, '.vcf > ~/height_prediction/imputed/output/temp_', W, '_', P, '_chr_', chr, '.txt '))
        	system(paste0("mv  ~/height_prediction/imputed/output/temp_", W, '_', P, '_chr_', chr, '.txt  ~/height_prediction/imputed/output/UKB_afr_',W, '_', P, '_chr_', chr, '.vcf'))
        	system(paste0('plink2 --bfile /project/mathilab/data/UKB/imputed/ukb_imp_chr', chr, '_eur --extract ~/height_prediction/imputed/output/UKB_chr_', chr,'_', W, '_', P, '.txt --recode vcf --out ~/height_prediction/imputed/output/UKB_eur_', W, '_', P, '_chr_', chr))
        #	system(paste0("grep -v '##' ~/height_prediction/imputed/output/UKB_eur_",W, '_', P, '_chr_', chr, '.vcf > ~/height_prediction/imputed/output/temp_', W, '_', P, '_chr_', chr, '.txt '))
        	system(paste0("mv  ~/height_prediction/imputed/output/temp_", W, '_', P, '_chr_', chr, '.txt  ~/height_prediction/imputed/output/UKB_eur_',W, '_', P, '_chr_', chr, '.vcf'))
        	system(paste0('bgzip ~/height_prediction/imputed/output/UKB_afr_', W, '_', P, '_chr_', chr, '.vcf'))
        	system(paste0('bgzip ~/height_prediction/imputed/output/UKB_eur_', W, '_', P, '_chr_', chr, '.vcf'))
       		cat(chr, 'r')
	}
}

