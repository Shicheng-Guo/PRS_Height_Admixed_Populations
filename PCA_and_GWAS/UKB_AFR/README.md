*Prepare input files for UKB_AFR PCA*

```
awk '{b=$1"_"$1;print b,"UKA"}' /project/mathilab/data/UKB/UKB_AFR.fam  > panel_file2.txt
touch "pop2.txt"
echo "UKA" > pop2.txt
cp ~/height_prediction/input/ukb_afr/UKB_AFR* .
#edit arguments file
#remove individual '6007195': TO DO TO DO TO DO
awk 'NR==28{print $0}' UKB_AFR.bas.vcf|tr "\t" "\n"|grep -n '6007195' #line 8801
awk '{$8801="";print $0; OFS="\t"}' UKB_AFR.bas.vcf > test.vcf
bgzip test.vcf

bsub -n 10 -M 40000 -R "span [hosts=1] rusage [mem=40480]" -o logPCA -e logPCA bash ~/height_prediction/runSmartpCA-master/generatePCA.sh ~/height_prediction/runSmartpCA-master/UKB_AFR/generatePCA_arg.txt
```

cp ~/height_prediction/input/ukb_afr/UKB_AFR_pheno.txt .
awk '{print $2,$1,$3,$4}' UKB_AFR_pheno.txt > My_Pheno.txt
head -1 UKB_AFR.bas.evec  > header.txt

sed -i "1d" UKB_AFR.bas.evec
sort UKB_AFR.bas.evec > tp
mv tp UKB_AFR.bas.evec
#sed -i 's/0_//g' UKB_AFR.bas.evec
sort -nms My_Pheno.txt > MyPheno.txt

echo "SUBJID PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 POP" > tmp
cat UKB_AFR.bas.evec >> tmp

#make Pheno.txt
R_script --vanilla R_script.R

#and then
awk '{print $1,$2,$18}' Pheno.txt > PHENO.txt #use residual height


awk 'a=$1"_"$1{print a,a,$3,$4,$5,$6}' UKB_AFR.fam > ww
mv ww UKB_AFR.fam

plink2 --bfile UKB_AFR --pheno PHENO.txt --allow-no-sex --covar Pheno.txt --covar-name Age2,Sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --out association --linear --adjust

grep ADD association.Res.Height.glm.linear  > test.txt

paste <(awk 'OFS="\t"{print $2,$1,$3,$4}' association.Height.glm.linear.adjusted|sort -k 1) <(sort -k 3 <(grep ADD association.Height.glm.linear)|awk 'OFS="\t"{print $3,$1,$2,$4,$5,$6,$7,$8,$9,$10,$11,$
12}')|awk '$1==$5;OFS="\t"{print $1,$2,$3,$4,$7,$8,$9,$10,$12,$13,$14,$15,$16}' > tmp1
cat <(echo -e "ID\tCHR\tUNADJ\tGC\tPOS\tREF\tALT\tA1\tOBS_CT\tBETA\tSE\tT_STAT\tP") tmp1 > ukb_afr_gwas_AA.txt
rm tmp1
rm tmp
#p value in GWAS



##a test

plink2 --bfile UKB_AFR --pheno PHENO.txt --allow-no-sex --covar Pheno.txt --covar-name Age2,Sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 -out association_v3 --glm --adjust

grep ADD association_v3.Res.Height.glm.linear > test3.txt

paste <(awk 'OFS="\t"{print $2,$1,$3,$4}' association_v3.Res.Height.glm.linear.adjusted|sort -k 1) <(sort -k 3 <(grep ADD association_v3.Res.Height.glm.linear)|awk 'OFS="\t"{print $3,$1,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12}')|awk '$1==$5;OFS="\t"{print $1,$2,$3,$4,$7,$8,$9,$10,$12,$13,$14,$15,$16}' > tmp1
cat <(echo -e "ID\tCHR\tUNADJ\tGC\tPOS\tREF\tALT\tA1\tOBS_CT\tBETA\tSE\tT_STAT\tP") tmp1 > ukb_afr_gwas_AA_v3.txt
rm tmp1
