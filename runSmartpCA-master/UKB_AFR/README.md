#UKB_AFR

awk '{b=$1"_"$1;print b,"UKA"}' ../../ukb_afr/UKB_AFR.fam > panel_file2.txt

touch "pop2.txt"
echo "UKA" > pop2.txt

cp ../../ukb_afr/UKB_AFR* .
bgzip UKB_AFR.bas.vcf
#edit arguments file
bsub -n 10 -M 40000 -R "span [hosts=1] rusage [mem=40480]" -o logPCA -e logPCA bash /project/mathilab/bbita/gwas_admix/new_height/runSmartpCA-master/generatePCA.sh /project/mathilab/bbita/gwas_admix/new_height/runSmartpCA-master/UKB_AFR/generatePCA_arg.txt


cp /project/mathilab/data/UKB/UKB_AFR_pheno.txt .
awk '{print $2,$1,$3,$4}' UKB_AFR_pheno.txt > My_Pheno.txt
head -1 UKB_AFR.bas.evec  > header.txt

sed -i "1d" UKB_AFR.bas.evec
sort UKB_AFR.bas.evec > tp
mv tp UKB_AFR.bas.evec
#sed -i 's/0_//g' UKB_AFR.bas.evec
sort -nms My_Pheno.txt > MyPheno.txt

echo "SUBJID PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 POP" > tmp
cat UKB_AFR.bas.evec >> tmp

#open R and combine the two and make
Pheno.txt





#and then
awk '{print $1,$2,$3,$4,$5}' Pheno.txt > PHENO.txt


awk 'a=$1"_"$1{print a,a,$3,$4,$5,$6}' UKB_AFR.fam > ww
mv ww UKB_AFR.fam

plink2 --bfile UKB_AFR --pheno PHENO.txt --allow-no-sex --covar Pheno.txt --covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,Age --out association --linear --adjust

grep ADD association.Height.glm.linear > test.txt

paste <(awk 'OFS="\t"{print $2,$1,$3,$4}' association.Height.glm.linear.adjusted|sort -k 1) <(sort -k 3 <(grep ADD association.Height.glm.linear)|awk 'OFS="\t"{print $3,$1,$2,$4,$5,$6,$7,$8,$9,$10,$11,$
12}')|awk '$1==$5;OFS="\t"{print $1,$2,$3,$4,$7,$8,$9,$10,$12,$13,$14,$15,$16}' > tmp1
cat <(echo -e "ID\tCHR\tUNADJ\tGC\tPOS\tREF\tALT\tA1\tOBS_CT\tBETA\tSE\tT_STAT\tP") tmp1 > ukb_afr_gwas_AA.txt
rm tmp1
rm tmp
#p value in GWAS



##a test

plink2 --bfile UKB_AFR --pheno PHENO.txt --allow-no-sex --covar Pheno.txt --covar-name Sex, Age, PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --out association_v2 --linear --adjust
plink2 --bfile UKB_AFR --pheno PHENO.txt --allow-no-sex --covar Pheno.txt --covar-name Sex,Age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --out association_v3 --glm --adjust



grep ADD association_v3.Height.glm.linear > test.txt

paste <(awk 'OFS="\t"{print $2,$1,$3,$4}' association_v3.Height.glm.linear.adjusted|sort -k 1) <(sort -k 3 <(grep ADD association_v3.Height.glm.linear)|awk 'OFS="\t"{print $3,$1,$2,$4,$5,$6,$7,$8,$9,$10,$11,$
12}')|awk '$1==$5;OFS="\t"{print $1,$2,$3,$4,$7,$8,$9,$10,$12,$13,$14,$15,$16}' > tmp1
cat <(echo -e "ID\tCHR\tUNADJ\tGC\tPOS\tREF\tALT\tA1\tOBS_CT\tBETA\tSE\tT_STAT\tP") tmp1 > ukb_afr_gwas_AA_v3.txt
rm tmp1
rm tmp

