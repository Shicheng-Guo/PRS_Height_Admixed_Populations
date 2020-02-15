```
Rscript --vanilla R_script.R

for chr in {1..22};do 
plink --noweb --vcf output/chr${chr}.vcf.gz --double-id --make-bed --out output/chr${chr} 
echo ${chr} 
done

cd output/

zcat ~/height_prediction/gwas/input/50_raw_filtered.txt.gz > 50.assoc.tsv

sed -i 's/beta/BETA/' 50.assoc.tsv 
sed -i 's/pval/P/' 50.assoc.tsv 
sed -i 's/se/SE/' 50.assoc.tsv 
sed -i 's/rsid/SNP/' 50.assoc.tsv 
sed -i 's/nCompleteSamples/N/' 50.assoc.tsv 
awk 'OFS="\t"{print $2, $1, $6, $7, $9, $3}' 50.assoc.tsv > tmp && mv tmp 50.assoc.tsv

awk -F":" '$1=$1' OFS="\t" 50.assoc.tsv|awk 'OFS="\t"{print $1,$4,$5,$6,$7, $8, $9}'| sed "1 s/^.*$/SNP\tAllele1\tAllele2\tBETA\tSE\tP\tN/" > tmp mv tmp 50.assoc.tsv
```

