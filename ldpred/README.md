##Try to replicate R2xEUR_ANC plots using LDpred.
##FOr now, using this tutorial as reference: https://choishingwan.github.io/PRS-Tutorial/ldpred/
##And also the ldpred tutorial here: 
source ~/ENV/bin/activate

1. Preprocessing the base data file and coordinating the data
The first step is a data synchronization step, where two or three data sets, genotypes and summary statistics are synchronized. This generates a HDF5 file which contains the synchronized genotypes. This step can be done by running

ldpred coord

use --help for detailed options. This step requires at least one genotype file (the LD reference genotypes), where we recommend at least 1000 unrelated individuals with the same ancestry make-up as the individuals for which summary statistics datasets are obtained from. Another genotype file can also be given if the user intends to validate the predictions using a separate set of genotypes.

```
awk '{print $2,$1,$3,$4}' /project/mathilab/data/UKB/UKB_EUR_pheno.txt > My_Pheno.txt
Rscript --vanilla R_script.R
Rscript --vanilla format_sumstat.R #format summary statistics file

#/project/mathilab/data/HRS/data/imputed/

plink2 \
    --bfile /project/mathilab/data/UKB/UKB_EUR \ 
    --pheno Pheno.txt \
    --keep /project/mathilab/data/UKB/UKB_EUR_IDS \
    --make-bed \
    --out ~/height_prediction/ldpred/output/EUR.ldpred 
#UKB_AFR
awk '{print $2,$1,$3,$4}' /project/mathilab/data/UKB/UKB_AFR_pheno.txt > My_Pheno.txt
Rscript --vanilla R_script.R

plink2 \
    --bfile /project/mathilab/data/UKB/UKB_AFR \
    --pheno Pheno.txt \
    --keep /project/mathilab/data/UKB/UKB_AFR_IDS \
    --make-bed \
    --out ~/height_prediction/ldpred/output/UKB_AFR.ldpred

#HRS_AFR

awk '{print $2,$1,$3,$4}' /project/mathilab/data/HRS/data/HRS_AFR_phenotypes.txt |sed s/SEX/Sex/|sed  s/HEIGHT/Height/ |sed s/AGE/Age/ > My_Pheno.txt
Rscript --vanilla R_script.R

plink2 \
    --bfile /project/mathilab/data/HRS/data/HRS_AFR_b37_strand_include \
    --keep /project/mathilab/data/HRS/data/HRS_AFR_IDS.fam \
    --mac 1 \
    --make-bed \
    --out ~/height_prediction/ldpred/output/HRS_AFR.ldpred

Rscript --vanilla merge_hrs_ukb.R
```

```
# There are 336,474 samples in the Height GWAS
ldpred coord \
    --rs SNP_ID \
    --A1 A1 \
    --A2 A2 \
    --pos POS \
    --chr CHR \
    --pval PVAL \
    --eff BETA \
    --ssf-format CUSTOM \
    --eff_type LINREG\
    --N 336474 \
    --ssf output/Height.QC.gz \ 
    --out output/EUR.coord \
    --gf output/EUR.ldpred

#HRS_AFR

ldpred coord \
    --rs SNP_ID \
    --A1 A1 \
    --A2 A2 \
    --pos POS \
    --chr CHR \
    --pval PVAL \
    --eff BETA \
    --ssf-format CUSTOM \
    --eff_type LINREG\
    --N 336474 \
    --ssf output/Height.QC.gz \
    --out output/HRS_AFR.coord \
    --gf output/HRS_AFR_UKB_EUR.ldpred > logs/logcoordhrsafr.log
```

2. Adjust the effect size estimates

```
# LDpred recommend radius to be Total number of SNPs in target / 3000= 589696/3000=196
#Regarding choice of the LD panel, its LD structure should ideally be similar to the training data for which the summary statistics are calculated.(UKB_EUR)
bsub -M 90000 -o ~/height_prediction/ldpred/logs/loggibbs  -e ~/height_prediction/ldpred/logs/loggibbs ~/.local/bin/ldpred gibbs  --cf ~/height_prediction/ldpred/output/EUR.coord  --ldr 196 --ldf ~/height_prediction/ldpred/output/EUR.ld --out ~/height_prediction/ldpred/output/EUR.weight --N 336474

```


3. Calculate the PRS
```
#For just PRS (no validation), run below code (UKB_EUR)
ldpred score \
    --gf output/EUR.ldpred \
    --rf output/EUR.weight \
    --out output/EUR.score \
    --only-score \
    --pf-format LSTANDARD > logs/logscore

#ukb afr
ldpred score \
    --gf output/UKB_AFR.ldpred \
    --rf output/EUR.weight \
    --out output/UKB_AFR.score \
    --only-score \
    --pf-format LSTANDARD > logs/logscore_ukb_afr

#WHI



#JHS


#HRS_afr
ldpred score \
    --gf output/HRS_AFR.ldpred \
    --rf output/EUR.weight \
    --out output/HRS_AFR.score \
    --only-score \
    --pf-format LSTANDARD > logs/logscore_hrs_afr

#HRS_eur
ldpred score \
    --gf output/HRS_EUR.ldpred \
    --rf output/EUR.weight \
    --out output/HRS_EUR.score \
    --only-score \
    --pf-format LSTANDARD > logs/logscore_hrs_eur

```
