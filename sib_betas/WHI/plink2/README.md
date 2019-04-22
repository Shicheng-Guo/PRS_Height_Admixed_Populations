*Prepare vcfs*
```
Rscript --vanilla R_script.R
for chr in {1..22};do
plink --noweb --vcf chr${chr}.vcf.gz --double-id --make-bed --out chr${chr}
echo ${chr}
done
```
*Prepare p values file*
```
cp ../../input/sib_beta_gwas_p.txt output/
sed -i 's/b/BETA/' sib_beta_gwas_p.txt
sed -i 's/p/P/' sib_beta_gwas_p.txt
#sed -i 's/se/SE/' 50.assoc.tsv
sed -i 's/MarkerName/SNP/' output/sib_beta_gwas_p.txt
```
