
Run R_script.R

for chr in {1..22};do
plink --noweb --vcf chr${chr}.vcf.gz --double-id --make-bed --out chr${chr}
echo ${chr}
done


cp /project/mathilab/bbita/gwas_admix/sib_gwas/sib_beta_gwas_p.txt .
sed -i 's/b/BETA/' sib_beta_gwas_p.txt
sed -i 's/p/P/' sib_beta_gwas_p.txt
#sed -i 's/se/SE/' 50.assoc.tsv
sed -i 's/MarkerName/SNP/' sib_beta_gwas_p.txt

