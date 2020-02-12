*Get SNP CHR and POS for each dataset*

```
cat ~/height_prediction/input/JHS/JHS_b37_strand.bas.vcf |grep -v "^#"|awk 'OFS="\t"{print $1,$2}' > JHS.txt
cat ~/height_prediction/input/WHI/WHI_b37_strand_include.bas.vcf |grep -v "^#"|awk 'OFS="\t"{print $1,$2}' > WHI.txt
cat ~/height_prediction/input/HRS_afr/HRS_AFR_b37_strand_include.bas.vcf |grep -v "^#"|awk 'OFS="\t"{print $1,$2}' > HRS_afr.txt
cat ~/height_prediction/input/ukb_afr/UKB_AFR.bas.vcf |grep -v "^#"|awk 'OFS="\t"{print $1,$2}' > UKB_afr.txt
```

*Combine them into one file*

```
cat WHI.txt >> all.txt
cat JHS.txt >> all.txt
cat HRS_afr.txt >> all.txt
cat UKB_afr.txt >> all.txt
```

*Clean*
```
awk -F"\t" '!seen[$1, $2]++' all.txt > sorted.txt #list of SNPs to run the gwas for #3290580  SNPs
rm all.txt
```

*Create UKB_AFR samples list for UKB AFR Imputed*

```
echo 'ID_1 ID_2 missing sex'> my_samples.txt
echo '0 0 0 D' >> my_samples.txt
grep -F -f /project/mathilab/data/UKB/imputed/UKB_AFR_IDS /project/mathilab/data/UKB/imputed/ukb33923_imp_v3_s487320.sample|grep -v '6007195' >> ~/height_prediction/runSmartpCA-master/UKB_AFR_imputed/my_samples.txt #for some reason this only has 8,809 instead of 8,813. Note that we remove individual 6007195 which has been removed from UKB
```


*Make VCF files from Imputed UKB_AFR data [keep only bi-allelic SNPs and preserve REF/ALT info*

Only SNPs present in the sorted.txt file are kept.
```
for chr in {1..22};
do
awk '$1=='${chr}'{print $2}' sorted.txt > tmp${chr}
grep -F -f tmp${chr} /project/mathilab/data/UKB/imputed/ukb_imp_chr${chr}_afr.bim  > ~/height_prediction/runSmartpCA-master/UKB_AFR_imputed/out${chr}.txt
awk '{print $2}' out${chr}.txt |sort|uniq -u > tmp_${chr}.txt #keep only bi-allelic SNPs
grep -F -f tmp_${chr}.txt out${chr}.txt| awk 'length($5)<2'|awk 'OFS="\t"{print $2,$5}' > ref_chr${chr}.txt
plink2 --bgen /project/mathilab/data/UKB/imputed/ukb_imp_chr${chr}_afr.bgen --sample ~/height_prediction/runSmartpCA-master/UKB_AFR_imputed/my_samples.sample --extract tmp_${chr}.txt  --ref-allele force ~/height_prediction/runSmartpCA-master/UKB_AFR_imputed/ref_chr${chr}.txt --recode vcf --out ~/height_prediction/runSmartpCA-master/UKB_AFR_imputed/chr${chr}
bgzip ~/height_prediction/runSmartpCA-master/UKB_AFR_imputed/chr${chr}.vcf
echo ${chr}
echo 'done'
done
```
*Combine individual chromosomes VCFs into one VCF*

```
bcftools concat chr1.vcf.gz chr2.vcf.gz -Oz -o chr1_2.vcf.gz
bcftools concat chr3.vcf.gz chr4.vcf.gz -Oz -o chr3_4.vcf.gz
bcftools concat chr5.vcf.gz chr6.vcf.gz -Oz -o chr5_6.vcf.gz
bcftools concat chr7.vcf.gz chr8.vcf.gz -Oz -o chr7_8.vcf.gz
bcftools concat chr9.vcf.gz chr10.vcf.gz -Oz -o chr9_10.vcf.gz
bcftools concat chr11.vcf.gz chr12.vcf.gz -Oz -o chr11_12.vcf.gz
bcftools concat chr13.vcf.gz chr14.vcf.gz chr15.vcf.gz chr16.vcf.gz chr17.vcf.gz chr18.vcf.gz chr19.vcf.gz chr20.vcf.gz chr21.vcf.gz chr22.vcf.gz -Oz -o chr13_22.vcf.gz
bcftools concat chr1_2.vcf.gz chr3_4.vcf.gz chr5_6.vcf.gz chr7_8.vcf.gz chr9_10.vcf.gz chr11_12.vcf.gz chr13_22.vcf.gz  -Oz -o chr1_22.vcf.gz
tabix -p vcf chr1_22.vcf.gz
```

*Remove intermediate vcf files*
```
for chr in {1..22};
do
rm chr${chr}.vcf.gz
done
rm chr1_2.vcf.gz chr3_4.vcf.gz chr5_6.vcf.gz chr7_8.vcf.gz chr9_10.vcf.gz chr11_12.vcf.gz chr13_122.vcf.gz
```

*Get SNP IDs and clean up temp files*
```
awk '{print $2}' tmp_*.txt |sort|uniq > SNP_ids.txt
rm out*.txt
rm tmp*
```
*Set up phenotype file and copy header and PCA files*

```
awk '{print $2,$1,$3,$4}' ~/height_prediction/input/ukb_afr/UKB_AFR_pheno.txt > My_Pheno.txt
cp ~/height_prediction/runSmartpCA-master/UKB_AFR/header.txt .
cp ~/height_prediction/runSmartpCA-master/UKB_AFR/UKB_AFR.bas.evec .
sort -nms My_Pheno.txt > MyPheno.txt
echo "SUBJID PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 POP" > tmp
cat UKB_AFR.bas.evec >> tmp
Rscript --vanilla ~/height_prediction/runSmartpCA-master/UKB_AFR_imputed/R_script.R #further configure phenotype file
#and then
awk '{print $1,$2,$18}' Pheno.txt > PHENO.txt #use residual height
```

*Convert back to plink format and keep only required columns in .fam file*
```
plink --vcf ~/height_prediction/runSmartpCA-master/UKB_AFR_imputed/chr1_22.vcf.gz --out ~/height_prediction/runSmartpCA-master/UKB_AFR_imputed/chr1_22
awk 'OFS="";{print $1, "_", $1, "\t",$2, "_", $2, "\t", $3, "\t", $4, "\t", $5, "\t", $6}' chr1_22.fam > temp
mv temp chr1_22.fam
```

*Run GWAS and clean output file*

```
plink2 --bfile ~/height_prediction/runSmartpCA-master/UKB_AFR_imputed/chr1_22  --pheno PHENO.txt --allow-no-sex --covar Pheno.txt --covar-name Age2,Sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 -out ~/height_prediction/runSmartpCA-master/UKB_AFR_imputed/association_v3 --glm --adjust
grep ADD association_v3.Res.Height.glm.linear > plink_ukb_afr_height_glm_linear.txt
paste <(awk 'OFS="\t"{print $2,$1,$3,$4}' association_v3.Res.Height.glm.linear.adjusted|sort -k 1) <(sort -k 3 <(grep ADD association_v3.Res.Height.glm.linear)|awk 'OFS="\t"{print $3,$1,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12}')|awk '$1==$5;OFS="\t"{print $1,$2,$3,$4,$7,$8,$9,$10,$12,$13,$14,$15,$16}' > tmp1
cat <(echo -e "ID\tCHR\tUNADJ\tGC\tPOS\tREF\tALT\tA1\tOBS_CT\tBETA\tSE\tT_STAT\tP") tmp1 > ukb_afr_gwas_AA_v3_imputed.txt
rm tmp1
```
