#Use imputed HRS  data

#1. might require additional filterning, say for genotyping quality and MAF
#2. will probably use rbgen R package

```
for i in 5000 10000 25000 50000 75000 100000 500000 1000000;
do
for j in 0.0005 0.00005 0.000005 0.0000005 0.00000005;
do
bsub -M 100000 Rscript --vanilla ~/height_prediction/imputed/scripts/prune_imputed.R ${i} ${j} HRS
bsub -M 100000 Rscript --vanilla ~/height_prediction/imputed/scripts/prune_imputed.R ${i} ${j} UKB
done
done


for i in 5000 10000 25000 50000 75000 100000 500000 1000000;
do
for j in 0.0005 0.00005 0.000005 0.0000005 0.00000005;
do
bsub -M 100000 -e ~/height_prediction/imputed/logs/plink_${i}_${j} -o ~/height_prediction/imputed/logs/plink_${i}_${j} Rscript --vanilla ~/height_prediction/imputed/scripts/plink_to_vcf.R ${i} ${j} HRS
bsub -M 100000 -e ~/height_prediction/imputed/logs/plink_${i}_${j} -o ~/height_prediction/imputed/logs/plink_${i}_${j} Rscript --vanilla ~/height_prediction/imputed/scripts/plink_to_vcf.R ${i} ${j} UKB
done
done


for i in 5000 10000 25000 50000 75000 100000 500000 1000000;
do
for j in 0.0005 0.00005 0.000005 0.0000005 0.00000005;
do
bsub -M 100000 -e ~/height_prediction/imputed/logs/prs_${i}_${j} -o ~/height_prediction/imputed/logs/prs_${i}_${j} Rscript --vanilla ~/height_prediction/imputed/scripts/run_PRS.R ${i} ${j}
bsub -M 100000 -e ~/height_prediction/imputed/logs/UKB_prs_${i}_${j} -o ~/height_prediction/imputed/logs/UKB_prs_${i}_${j} Rscript --vanilla ~/height_prediction/imputed/scripts/run_PRS_UKB.R ${i} ${j}
done
done


Rscript --vanilla scripts/pheno_pred.R
Rscript --vanilla scripts/pheno_pred_UKB.R


#headers

gunzip -c ~/height_prediction/imputed/output/HRS_afr_100000_0.0005_chr_22.vcf.gz|head -7 > ~/height_prediction/imputed/header_HRS_afr.txt
gunzip -c ~/height_prediction/imputed/output/HRS_eur_100000_0.0005_chr_22.vcf.gz|head -7 > ~/height_prediction/imputed/header_HRS_eur.txt
gunzip -c ~/height_prediction/imputed/output/UKB_afr_100000_0.0005_chr_22.vcf.gz|head -7 > ~/height_prediction/imputed/header_UKB_afr.txt
gunzip -c ~/height_prediction/imputed/output/UKB_eur_100000_0.0005_chr_22.vcf.gz|head -7 > ~/height_prediction/imputed/header_UKB_eur.txt
```
