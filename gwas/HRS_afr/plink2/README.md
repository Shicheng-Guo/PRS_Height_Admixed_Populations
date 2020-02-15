*Prepare vcfs*
```
Rscript --vanilla R_script.R
for chr in {1..22};do
plink --noweb --vcf output/chr${chr}.vcf.gz --double-id --make-bed --out output/chr${chr}
echo ${chr}
done
```
