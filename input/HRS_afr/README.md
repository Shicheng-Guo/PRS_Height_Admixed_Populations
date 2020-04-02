
*copy files to my local directory*
```
cp /project/mathilab/data/HRS/data/HRS_AFR_b37_strand_include* .
#cp /project/mathilab/data/HRS/admixture/HRS_AFR_b37_strand_prune_include.2.Q . ##order comes from file below
#cp /project/mathilab/data/HRS/data/HRS_EUR_b37_strand_include.fam .  #the order of samples
cp /project/mathilab/data/HRS/data/HRS_AFR_phenotypes.txt .
```
*convert to vcf format*

```
plink --allow-extra-chr --allow-no-sex --autosome \
--bfile /project/mathilab/data/HRS/data/HRS_AFR_b37_strand_include \
--recode vcf \
--keep-allele-order \
--out ~/height_prediction/input/HRS_afr/HRS_AFR_b37_strand_include.bas
```

```
awk '/^ *#/ { print; }' HRS_AFR_b37_strand_include.bas.vcf > header_HRS_afr.txt

for chr in {1..22};
do
touch chr${chr}_bas.vcf
cat header_HRS_afr.txt >  chr${chr}_bas.vcf
awk -v i=$chr '$1==i' HRS_AFR_b37_strand_include.bas.vcf >> chr${chr}_bas.vcf
echo 'chr'
echo $chr
echo 'done'
done
```

*compress/index*
```
for chr in {1..22};
do
bgzip -c chr${chr}_bas.vcf > chr${chr}_bas.vcf.gz &&
tabix -p vcf chr${chr}_bas.vcf.gz
done

bgzip -c HRS_AFR_b37_strand_include.bas.vcf > HRS_AFR_b37_strand_include.bas.vcf.gz
