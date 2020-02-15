RUN A GWAS ON INDIVIDUALS WITH SOME LEVEL OF AFRICAN ANCESTRY

*WHI, ukb_afr datasets
*use plink to run gwas
*goal it to compare estimated 'African' betas with those we estimated with the local ancestry analyses


#WHI
awk '{b="0_"$2;print b,"AA"}' ../WHI/WHI_b37_strand_include.fam > panel_file.txt
touch "pop.txt"
echo "AA" > pop.txt

cp ../WHI/WHI_b37_strand_include.bas.vcf.gz .
bash /project/mathilab/bbita/gwas_admix/new_height/PCA_and_GWAS/generatePCA.sh /project/mathilab/bbita/gwas_admix/new_height/runSmartpCA-master/generatePCA_arg.txt

*copy files to my local directory*
```
cp /project/mathilab/data/WHI/data/WHI_b37_strand_include* .
cp /project/mathilab/data/WHI/data/WHI_phenotypes.txt .
awk '{print $2,$1,$3,$4, $5,$87, $9}' WHI_phenotypes.txt > My_Pheno.txt
head -1 WHI_b37_strand_include.bas.evec > header.txt

sed -i "1d" WHI_b37_strand_include.bas.evec
sort WHI_b37_strand_include.bas.evec > tp
mv tp WHI_b37_strand_include.bas.evec
sed -i 's/0_//g' WHI_b37_strand_include.bas.evec
sort -nms My_Pheno.txt > MyPheno.txt

echo "SUBJID PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 POP" > tmp
cat WHI_b37_strand_include.bas.evec >> tmp

#open R and combine the two and make 
Pheno.txt

and then
awk '{print $1,$2,$3,$4,$5,$6}' Pheno.txt > PHENO.txt

plink2 --bfile WHI_b37_strand_include --pheno PHENO.txt --allow-no-sex --covar Pheno.txt --covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,AGE --out association --linear --adjust

paste <(awk 'OFS="\t"{print $2,$1,$3,$4}' association.HEIGHTX.glm.linear.adjusted|sort -k 1) <(sort -k 3 <(grep ADD association.HEIGHTX.glm.linear)|awk 'OFS="\t"{print $3,$1,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12}')|awk '$1==$5;OFS="\t"{print $1,$2,$3,$4,$7,$8,$9,$10,$12,$13,$14,$15,$16}' > tmp1
cat <(echo -e "ID\tCHR\tUNADJ\tGC\tPOS\tREF\tALT\tA1\tOBS_CT\tBETA\tSE\tT_STAT\tP") tmp1 > whi_gwas_AA.txt
rm tmp1
rm tmp
#p value in GWAS




https://www.researchgate.net/post/How_do_I_adjust_the_p_value_of_SNPs-phenotype_association_study_for_multiple_testing
