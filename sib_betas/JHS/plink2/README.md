*Prepare vcfs*
```
Rscript --vanilla R_script.R
for chr in {1..22};do
plink --noweb --vcf chr${chr}.vcf.gz --double-id --make-bed --out chr${chr}
echo ${chr}
done
```
*Prepare p values file*
cp ../../input/sibestimates_50.tsv output/
sed -i 's/b/BETA/' output/sibestimates_50.tsv
sed -i 's/p/P/' output/sibestimates_50.tsv
#sed -i 's/se/SE/' 50.assoc.tsv
sed -i 's/MarkerName/SNP/' output/sibestimates_50.tsv
sed -i 's/BETAeta/BETA/' output/sibestimates_50.tsv
