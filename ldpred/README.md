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
cp /home/bbita/height_prediction/runSmartpCA-master/UKB_AFR/R_script.R .
Rscript --vanilla R_script.r
Rscript --vanilla format_sumstat.R #format summary statistics file

plink2 \
    --bfile /project/mathilab/data/UKB/UKB_EUR \ 
    --pheno Pheno.txt \
    --keep /project/mathilab/data/UKB/UKB_EUR_IDS \
    --make-bed \
    --out ~/height_prediction/ldpred/output/EUR.ldpred 
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
```

2. Adjust the effect size estimates

```
# LDpred recommend radius to be Total number of SNPs in target / 3000= 589696/3000=196

bsub -M 90000 -o ~/height_prediction/ldpred/logs/loggibbs  -e ~/height_prediction/ldpred/logs/loggibbs ~/.local/bin/ldpred gibbs  --cf ~/height_prediction/ldpred/output/EUR.coord  --ldr 196 --ldf ~/height_prediction/ldpred/output/EUR.ld --out ~/height_prediction/ldpred/output/EUR.weight --N 336474
#bsub -M 90000 -o ~/height_prediction/ldpred/logs/loggibbs_1E-04.weight -e ~/height_prediction/ldpred/logs/loggibbs_1E-04.weight ldpred gibbs --cf ~/height_prediction/ldpred/output/EUR.coord  --ldr 196 --ldf ~/height_prediction/ldpred/output/EUR.ld --out ~/height_prediction/ldpred/output/EUR_1E-04.weight --N 336474 --f 0.0001 --h2 0.8
#bsub -M 90000 -o ~/height_prediction/ldpred/logs/loggibbs_1E-02.weight -e ~/height_prediction/ldpred/logs/loggibbs_1E-02.weight ldpred gibbs --cf ~/height_prediction/ldpred/output/EUR.coord  --ldr 196 --ldf ~/height_prediction/ldpred/output/EUR.ld --out ~/height_prediction/ldpred/output/EUR_1E-02.weight --N 336474 --f 0.01 --h2 0.8
#bsub -M 90000 -o ~/height_prediction/ldpred/logs/loggibbs_1E-03.weight -e ~/height_prediction/ldpred/logs/loggibbs_1E-03.weight ldpred gibbs --cf ~/height_prediction/ldpred/output/EUR.coord  --ldr 196 --ldf ~/height_prediction/ldpred/output/EUR.ld --out ~/height_prediction/ldpred/output/EUR_1E-03.weight --N 336474 --f 0.001 --h2 0.8

```


3. Calculate the PRS
```
ldpred score \
    --gf output/EUR.ldpred \
    --rf output/EUR_1E-02.weights_LDpred \
    --out output/EUR_1E-02.score \
    --pf ~/height_prediction/input/ukb_eur/UKB_EUR_pheno.txt \
    --pf-format LSTANDARD 
#the above would be for validation. For just PRS, run below code.
ldpred score \
    --gf output/EUR.ldpred \
    --rf output/EUR_1E-02.weight \
    --out output/EUR_1E-02.score \
    --only-score \
    --pf-format LSTANDARD

ldpred score \
    --gf output/EUR.ldpred \
    --rf output/EUR_1E-03.weight \
    --out output/EUR_1E-03.score \
    --only-score \
    --pf-format LSTANDARD


```
