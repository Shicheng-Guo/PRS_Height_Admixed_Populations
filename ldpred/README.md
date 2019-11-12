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
    --bfile /project/mathilab/data/UKB/UKB_EUR \ #QC'd EUROPEAN bfiles
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
    --N 336474 \
    --ssf output/Height.QC.gz \ 
    --out output/EUR.coord \
    --gf output/EUR.ldpred
```

2. Adjust the effect size estimates

```
# LDpred recommend radius to be Total number of SNPs in target / 3000= 589696/3000=196

bsub -M 90000 -o ~/height_prediction/ldpred/loggibbs  -e ~/height_prediction/ldpred/loggibbs ~/.local/bin/ldpred gibbs  --cf ~/height_prediction/ldpred/output/EUR.coord  --ldr 196 --ldf ~/height_prediction/ldpred/output/EUR.ld --out ~/height_prediction/ldpred/output/EUR.weight --N 336474
```


3. Calculate the PRS
```
ldpred score \
    --gf output/EUR.ldpred \
    --rf output/EUR.weight \
    --out output/EUR.score \
    --pf output/EUR.height \
    --pf-format LSTANDARD 
```


