*Prepare vcfs*

```
bsub -M 90000 Rscript --vanilla ~/height_prediction/gwas/HRS_eur/plink2/R_script.R
for chr in {1..22};do
plink --noweb --vcf output/chr${chr}.vcf.gz --double-id --make-bed --out output/chr${chr}
echo ${chr}
done
```

