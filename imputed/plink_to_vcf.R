#!/usr/bin/env Rscript
############################
options(scipen=10)
options(digits=10)

for(W in c(5000,10000,25000,50000,75000,100000,500000)){
	for(P in c(0.0005,0.00005, 0.000005,0.0000005,0.00000005)){
		for (chr in 22:1){
        		system(paste0('plink2 --bfile /project/mathilab/data/HRS/data/imputed/HRS_AFR_imputed_chr', chr, ' --extract ~/height_prediction/imputed/output/chr_', chr,'_', W, '_', P, '.txt --recode vcf --out ~/height_prediction/imputed/output/HRS_afr_', W, '_', P, '_chr_', chr))
	        	system(paste0('plink2 --bfile /project/mathilab/data/HRS/data/imputed/HRS_EUR_imputed_chr', chr, ' --extract ~/height_prediction/imputed/output/chr_', chr,'_', W, '_', P, '.txt --recode vcf --out ~/height_prediction/imputed/output/HRS_eur_', W, '_', P, '_chr_', chr))
        		system(paste0('bgzip ~/height_prediction/imputed/output/HRS_afr_', W, '_', P, '_chr_', chr, '.vcf'))
			system(paste0('bgzip ~/height_prediction/imputed/output/HRS_eur_', W, '_', P, '_chr_', chr, '.vcf'))
       			cat(chr, 'r')
		}
		cat(P, ' done\n')	
	}
	cat(W, ' done\n')
}

#UKB has problems with chr 1 and 2, so far now I won't run this:

#for(chr in 22:1){
#        system(paste0('plink2 --bgen /project/mathilab/data/UKB/imputed/ukb_imp_chr',chr,'_afr.bgen --extract ~/height_prediction/imputed/output/chr_', chr,'.txt --recode vcf --out ~/height_prediction/imputed/output/UKB_afr_chr_', chr))
#	system(paste0('plink2 --bgen /project/mathilab/data/UKB/imputed/ukb_imp_chr',chr,'_eur.bgen --extract ~/height_prediction/imputed/output/chr_', chr,'.txt --recode vcf --out ~/height_prediction/imputed/output/UKB_eur_chr_', chr))
#        system(paste0('bgzip ~/height_prediction/imputed/output/UKB_eur_chr_', chr, '.vcf'))
#	system(paste0('bgzip ~/height_prediction/imputed/output/UKB_afr_chr_', chr, '.vcf'))
#        cat(chr, 'done\n')
#}


