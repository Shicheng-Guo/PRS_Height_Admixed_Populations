##Try to replicate R2xEUR_ANC plots using LDpred.
##FOr now, using this tutorial as reference: https://choishingwan.github.io/PRS-Tutorial/ldpred/
##And also the ldpred tutorial here: 
source ~/ENV/bin/activate

1. Preprocessing the base data file and coordinating the data
The first step is a data synchronization step, where two or three data sets, genotypes and summary statistics are synchronized. This generates a HDF5 file which contains the synchronized genotypes. This step can be done by running

```
ldpred coord
```

use --help for detailed options. This step requires at least one genotype file (the LD reference genotypes), where we recommend at least 1000 unrelated individuals with the same ancestry make-up as the individuals for which summary statistics datasets are obtained from. Another genotype file can also be given if the user intends to validate the predictions using a separate set of genotypes.

```
1000G EUR (except FIN) as LD panel

grep EUR /project/mathilab/data/1kg/20130502_phase3_final/integrated_call_samples_v3.20130502.ALL.panel |grep -v FIN|awk 'OFS="\t"{print $1, $1}' > EUR_samples.txt

for chr in {1..22};
do
plink2 --vcf /project/mathilab/data/1kg/20130502_phase3_final/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --keep EUR_samples.txt --make-bed --out output/1000g_chr${chr}
done
#merge all vcfs into one

#UKB_EUR
awk '{print $2,$1,$3,$4}' /project/mathilab/data/UKB/UKB_EUR_pheno.txt > My_Pheno_UKB_eur.txt
Rscript --vanilla format_pheno.R My_Pheno_UKB_eur.txt
Rscript --vanilla format_sumstat.R #format summary statistics file

#LDpred does not support filtering of samples and SNPs, so therefore we must generate a new QC'ed genotype file using plink:
plink2 \
    --bfile /project/mathilab/data/UKB/UKB_EUR \ 
    --pheno Pheno_UKB_eur.txt \
    --keep /project/mathilab/data/UKB/UKB_EUR_IDS \
    --make-bed \
    --out ~/height_prediction/ldpred/output/UKB_EUR.ldpred 
#
#UKB_AFR
awk '{print $2,$1,$3,$4}' /project/mathilab/data/UKB/UKB_AFR_pheno.txt > My_Pheno_UKB_afr.txt
Rscript --vanilla format_pheno.R My_Pheno_UKB_afr.txt

plink2 \
    --bfile /project/mathilab/data/UKB/UKB_AFR \
    --pheno Pheno_UKB_afr.txt \
    --keep /project/mathilab/data/UKB/UKB_AFR_IDS \
    --make-bed \
    --out ~/height_prediction/ldpred/output/UKB_AFR.ldpred
#
#HRS_AFR
awk '{print $2,$1,$3,$4}' /project/mathilab/data/HRS/data/HRS_AFR_phenotypes.txt |sed s/SEX/Sex/|sed  s/HEIGHT/Height/ |sed s/AGE/Age/ > My_Pheno_HRS_afr.txt
Rscript --vanilla format_pheno.R My_Pheno_HRS_afr.txt

plink2 \
    --bfile /project/mathilab/data/HRS/data/HRS_AFR_b37_strand_include \
    --keep /project/mathilab/data/HRS/data/HRS_AFR_IDS.fam \
    --pheno Pheno_HRS_afr.txt \
    --mac 1 \
    --make-bed \
    --out ~/height_prediction/ldpred/output/HRS_AFR.ldpred
#
#HRS_EUR
awk '{print $2,$1,$3,$4}' /project/mathilab/data/HRS/data/HRS_EUR_phenotypes.txt |sed s/SEX/Sex/|sed  s/HEIGHT/Height/ |sed s/AGE/Age/ > My_Pheno_HRS_eur.txt
Rscript --vanilla format_pheno.R My_Pheno_HRS_eur.txt

plink2 \
    --bfile /project/mathilab/data/HRS/data/HRS_EUR_b37_strand_include \
    --keep /project/mathilab/data/HRS/data/HRS_EUR_IDS.fam \
    --pheno Pheno_HRS_eur.txt \
    --mac 1 \
    --make-bed \
    --out ~/height_prediction/ldpred/output/HRS_EUR.ldpred
#WHI
awk '{print $5,$1,$87, $9}' /project/mathilab/data/WHI/data/WHI_phenotypes.txt |sed  s/HEIGHTX/Height/ |sed s/AGE/Age/|sed s/SEX/Sex/|sed s/SUBJID/ID/ > My_Pheno_WHI.txt
Rscript --vanilla format_pheno.R My_Pheno_HRS_eur.txt

plink2 \
    --bfile /project/mathilab/data/WHI/data/WHI_b37_strand_include \
    --pheno Pheno_WHI.txt \
    --mac 1 \
    --make-bed \
    --out ~/height_prediction/ldpred/output/WHI.ldpred

#JHS
cp /project/mathilab/data/JHS/data/JHS_phenotypes.txt .
sed -i 's/African American/African_American/' JHS_phenotypes.txt
awk '{print $5, $1, $2,$15, $10}' JHS_phenotypes.txt|sed 's/SUBJID/ID/'|sed 's/height_baseline/Height/'|sed 's/age_baseline/Age/'|sed 's/SEX/Sex/' > My_Pheno_JHS.txt
cp /project/mathilab/data/JHS/data/JHS_b37_strand* .
awk '{print $2, $1, $3, $4, $5, $6}' JHS_b37_strand.fam > tmp && mv tmp JHS_b37_strand.fam

plink2 \
    --bfile JHS_b37_strand \
    --pheno Pheno_JHS.txt \
    --mac 1 \
    --make-bed \
    --out ~/height_prediction/ldpred/output/JHS.ldpred

```
**Merge**
```
Rscript --vanilla merge_hrs_ukb.R
```
**COORD**
```
#Preprocessing the base data file:

# There are 361,194 samples in the Height GWAS
#UKB_EUR
ldpred coord \
    --rs SNP \
    --A1 A1 \
    --A2 A2 \
    --pos POS \
    --chr CHR \
    --pval PVAL \
    --eff BETA \
    --ssf-format CUSTOM \
    --eff_type LINREG\
    --N 361194  \
    --ssf output/Height.QC.gz \ 
    --out output/UKB_EUR.coord \
    --gf output/UKB_EUR_UKB_EUR.ldpred

#HRS_AFR
ldpred coord \
    --rs SNP \
    --A1 A1 \
    --A2 A2 \
    --pos POS \
    --chr CHR \
    --pval PVAL \
    --eff BETA \
    --ssf-format CUSTOM \
    --eff_type LINREG\
    --N 361194 \
    --vbim output/HRS_AFR_UKB_EUR.ldpred.bim \
    --ssf output/Height.QC.gz \
    --out output/UKB_EUR_HRS_AFR.coord \
    --gf output/UKB_EUR_UKB_EUR.ldpred

#HRS_EUR
ldpred coord \
    --rs SNP \
    --A1 A1 \
    --A2 A2 \
    --pos POS \
    --chr CHR \
    --pval PVAL \
    --eff BETA \
    --ssf-format CUSTOM \
    --eff_type LINREG\
    --N 361194 \
    --vbim output/HRS_EUR_UKB_EUR.ldpred.bim \
    --ssf output/Height.QC.gz \
    --out output/UKB_EUR_HRS_EUR.coord \
    --gf output/UKB_EUR_UKB_EUR.ldpred
#UKB_AFR
ldpred coord \
    --rs SNP \
    --A1 A1 \
    --A2 A2 \
    --pos POS \
    --chr CHR \
    --pval PVAL \
    --eff BETA \
    --ssf-format CUSTOM \
    --eff_type LINREG\
    --N 361194 \
    --vbim output/UKB_AFR_UKB_EUR.ldpred.bim \
    --ssf output/Height.QC.gz \
    --out output/UKB_EUR_UKB_AFR.coord \
    --gf output/UKB_EUR_UKB_EUR.ldpred 
#WHI
ldpred coord \
    --rs SNP \
    --A1 A1 \
    --A2 A2 \
    --pos POS \
    --chr CHR \
    --pval PVAL \
    --eff BETA \
    --ssf-format CUSTOM \
    --eff_type LINREG\
    --N 361194 \
    --vbim output/WHI_UKB_EUR.ldpred.bim \
    --ssf output/Height.QC.gz \
    --out output/UKB_EUR_WHI.coord \
    --gf output/UKB_EUR_UKB_EUR.ldpred
#JHS
ldpred coord \
    --rs SNP \
    --A1 A1 \
    --A2 A2 \
    --pos POS \
    --chr CHR \
    --pval PVAL \
    --eff BETA \
    --ssf-format CUSTOM \
    --eff_type LINREG\
    --N 361194 \
    --vbim output/JHS_UKB_EUR.ldpred.bim \
    --ssf output/Height.QC.gz \
    --out output/UKB_EUR_JHS.coord \
    --gf output/UKB_EUR_UKB_EUR.ldpred
```

2. Adjust the effect size estimates

```
# LDpred recommend radius to be Total number of SNPs in target / 3000. E.g. 784256/3000=261
#Regarding choice of the LD panel, its LD structure should ideally be similar to the training data for which the summary statistics are calculated.(UKB_EUR)

#UKB_EUR 558565/3000~186
bsub -q denovo "ldpred gibbs  --cf ~/height_prediction/ldpred/output/UKB_EUR.coord  --ldr 250 --ldf ~/height_prediction/ldpred/output/UKB_EUR.ld --out ~/height_prediction/ldpred/output/UKB_EUR.weight --no-ld-compression --n-iter 1000 --n-burn-in 50 > ~/height_prediction/ldpred/logs/gibbs_ukb_eur.log"
#UKB_AFR  same as UKB_Eur
#HRS_eur 215580/3000~72
bsub -q denovo "ldpred gibbs  --cf ~/height_prediction/ldpred/output/UKB_EUR_HRS_EUR.coord  --ldr 250 --ldf ~/height_prediction/ldpred/output/HRS_EUR.ld --out ~/height_prediction/ldpred/output/HRS_EUR.weight --no-ld-compression --n-iter 1000 --n-burn-in 50 >  ~/height_prediction/ldpred/logs/gibbs_hrs_eur.log"
#HRS_afr 215504/3000~72
bsub -q denovo "ldpred gibbs  --cf ~/height_prediction/ldpred/output/UKB_EUR_HRS_AFR.coord  --ldr 250 --ldf ~/height_prediction/ldpred/output/HRS_AFR.ld --out ~/height_prediction/ldpred/output/HRS_AFR.weight --no-ld-compression --n-iter 1000 --n-burn-in 50 > ~/height_prediction/ldpred/logs/gibbs_hrs_afr.log"
#WHI 72425/3000 ~ 24
bsub -q denovo "ldpred gibbs  --cf ~/height_prediction/ldpred/output/UKB_EUR_WHI.coord  --ldr 250 --ldf ~/height_prediction/ldpred/output/WHI.ld --out ~/height_prediction/ldpred/output/WHI.weight --no-ld-compression --n-iter 1000 --n-burn-in 50 > ~/height_prediction/ldpred/logs/gibbs_whi.log"
#JHS 68988/3000 23
bsub -q denovo "ldpred gibbs  --cf ~/height_prediction/ldpred/output/UKB_EUR_JHS.coord --ldr 250 --ldf ~/height_prediction/ldpred/output/JHS.ld --out ~/height_prediction/ldpred/output/JHS.weight --no-ld-compression --n-iter 1000 --n-burn-in 50 > ~/height_prediction/ldpred/logs/gibbs_jhs.log"

#P+T
ldpred p+t --cf ~/height_prediction/ldpred/output/UKB_EUR.coord  --ldr 250 --out ~/height_prediction/ldpred/output/UKB_EUR_pt.weight > logs/pt_ukb_eur.log
ldpred p+t --cf ~/height_prediction/ldpred/output/UKB_EUR_HRS_EUR.coord  --ldr 250 --out ~/height_prediction/ldpred/output/HRS_EUR_pt.weight  > logs/pt_hrs_eur.log
ldpred p+t --cf ~/height_prediction/ldpred/output/UKB_EUR_HRS_AFR.coord  --ldr 250 --out ~/height_prediction/ldpred/output/HRS_AFR_pt.weight > logs/pt_hrs_afr.log
ldpred p+t --cf ~/height_prediction/ldpred/output/UKB_EUR_WHI.coord  --ldr 250 --out ~/height_prediction/ldpred/output/WHI_pt.weight > logs/pt_whi.log
ldpred p+t --cf ~/height_prediction/ldpred/output/UKB_EUR_JHS.coord --ldr 250 --out ~/height_prediction/ldpred/output/JHS_pt.weight > logs/pt_JHS.log
ldpred p+t --cf ~/height_prediction/ldpred/output/UKB_EUR_UKB_AFR.coord --ldr 250 --out ~/height_prediction/ldpred/output/UKB_AFR_pt.weight  > logs/pt_ukb_afr.log
```

3. Calculate the PRS
```
#For just PRS (no validation), run below code (UKB_EUR)
ldpred score \
    --gf output/UKB_EUR_UKB_EUR.ldpred \
    --rf output/UKB_EUR.weight \
    --out output/UKB_EUR.score \
    --only-score \
     --summary-file 
    --pf-format LSTANDARD 

ldpred score \
    --gf output/UKB_EUR_UKB_EUR.ldpred \
    --rf output/UKB_EUR_pt.weight \
    --out output/UKB_EUR_pt.score \
    --only-score \
    --pf-format LSTANDARD > logs/score_pt_ukb_eur.log

#ukb afr
ldpred score \
    --gf output/UKB_AFR_UKB_EUR.ldpred \
    --rf output/UKB_AFR.weight \
    --out output/UKB_AFR.score \
    --only-score \
    --pf-format LSTANDARD > logs/logscore_ukb_afr
ldpred score \
    --gf output/UKB_AFR_UKB_EUR.ldpred \
    --rf output/UKB_AFR_pt.weight \
    --out output/UKB_AFR_pt.score \
    --only-score \
    --pf-format LSTANDARD > logs/score_pt_ukb_afr.log

#HRS_EUR
ldpred score \
    --gf output/HRS_EUR_UKB_EUR.ldpred \
    --rf output/HRS_EUR.weight \
    --out output/HRS_EUR.score \
    --only-score \
    --pf-format LSTANDARD > logs/score_hrs_eur.log
ldpred score \
    --gf output/HRS_EUR_UKB_EUR.ldpred \
    --rf output/HRS_EUR_pt.weight \
    --out output/HRS_EUR_pt.score \
    --only-score \
    --pf-format LSTANDARD > logs/score_pt_hrs_eur.log

#HRS_afr
ldpred score \
    --gf output/HRS_AFR_UKB_EUR.ldpred \
    --rf output/HRS_AFR.weight \
    --out output/HRS_AFR.score \
    --only-score \
    --pf-format LSTANDARD > logs/score_hrs_afr.log
ldpred score \
    --gf output/HRS_AFR_UKB_EUR.ldpred \
    --rf output/HRS_AFR_pt.weight \
    --out output/HRS_AFR_pt.score \
    --only-score \
    --pf-format LSTANDARD > logs/score_pt_hrs_afr.log
#WHI
ldpred score \
    --gf output/WHI_UKB_EUR.ldpred \
    --rf output/WHI.weight \
    --out output/WHI.score \
    --only-score \
    --pf-format LSTANDARD 
ldpred score \
    --gf output/WHI_UKB_EUR.ldpred \
    --rf output/WHI_pt.weight \
    --out output/WHI_pt.score \
    --only-score \
    --pf-format LSTANDARD > logs/score_pt_whi.log
#JHS
ldpred score \
    --gf output/JHS_UKB_EUR.ldpred \
    --rf output/JHS.weight \
    --out output/JHS.score \
    --only-score \
    --pf-format LSTANDARD
ldpred score \
    --gf output/JHS_UKB_EUR.ldpred \
    --rf output/JHS_pt.weight \
    --out output/JHS_pt.score \
    --only-score \
    --pf-format LSTANDARD > logs/score_pt_jhs.log
```
